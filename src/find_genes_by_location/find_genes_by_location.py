import argparse
from collections import defaultdict
import csv
from dataclasses import dataclass, field
from enum import Enum, unique, auto
import os
import sys
import tempfile
import yaml
import zipfile

import gffutils
from google.protobuf import json_format
from ncbi.datasets.v1alpha1 import dataset_catalog_pb2
from ncbi.datasets.v1alpha1.reports import assembly_pb2
from ncbi.datasets.reports.report_reader import DatasetsReportReader


def retrieve_assembly_report(zip_in, catalog, assm_acc: str) -> assembly_pb2.AssemblyDataReport:
    report_files = get_catalog_files_for_assembly(catalog, dataset_catalog_pb2.File.FileType.DATA_REPORT, assm_acc)
    for path in report_files:
        yaml = zip_in.read(path)
        rpt_rdr = DatasetsReportReader()
        return rpt_rdr.assembly_report(yaml)


def retrieve_data_catalog(zip_in) -> dataset_catalog_pb2.Catalog:
    catalog_json = zip_in.read('ncbi_dataset/data/dataset_catalog.json')
    return json_format.Parse(catalog_json, dataset_catalog_pb2.Catalog())


def get_catalog_files_for_assembly(catalog: dataset_catalog_pb2.Catalog, desired_filetype: dataset_catalog_pb2.File.FileType, assm_acc: str):
    report_files = get_catalog_files(catalog, desired_filetype, assm_acc)
    filepaths = []
    for assm_acc, paths in report_files.items():
        filepaths.extend(paths)
    return filepaths


def get_catalog_files(catalog: dataset_catalog_pb2.Catalog, desired_filetype: dataset_catalog_pb2.File.FileType, assm_acc: str = None):
    files = defaultdict(list)
    for assm in catalog.assemblies:
        acc = assm.accession
        if assm_acc and assm_acc != acc:
            continue
        for f in assm.files:
            filepath = os.path.join('ncbi_dataset', 'data', f.file_path)
            if f.file_type == desired_filetype:
                files[acc].append(filepath)
    return files


def get_zip_file_for_acc(acc, path):
    fname = os.path.join(path, f'{acc}.zip')
    if os.path.isfile(fname):
        return fname
    return None


@dataclass
class Gene:
    id: str
    feat_type: str
    name: str
    chrom: str
    strand: str
    range_start: int
    range_stop: int
    protein_accession: str = ""

    def get_fields(self):
        return [self.feat_type, self.name, self.range_start, self.range_stop, self.protein_accession]

    def name_val(self):
        return self.protein_accession if self.protein_accession else self.name


def find_genes_by_loc(gff3_db, csvout, assm_acc, seq_acc, start, stop, extra_fields):
    found_genes = []
    feat_types = ('gene', 'pseudogene')
    for gene in gff3_db.region(seqid=seq_acc, start=start, end=stop, featuretype=feat_types, completely_within=False):
        gene_name = gene.attributes.get('Name', None)[0]
        prot_acc = ""
        if gene.attributes['gene_biotype'][0] == 'protein_coding':
            cds = list(gff3_db.children(gene, featuretype='CDS'))
            prot_acc = cds[0].attributes.get('protein_id', None)[0]
        geneobj = Gene(
            gene.id,
            gene.featuretype,
            gene_name,
            gene.chrom,
            gene.strand,
            gene.start,
            gene.stop,
            prot_acc,
        )
        csvout.writerow([assm_acc, seq_acc, start, stop, *extra_fields, *geneobj.get_fields()])
        found_genes.append(geneobj)
    return found_genes


class FindGenesByLoc:
    default_packages_dir = os.path.join('var', 'data', 'packages')

    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--packages-dir', type=str, default=self.default_packages_dir,
                            help=f'root of input data directory [{self.default_packages_dir}]')
        parser.add_argument('--locs', type=str, help='file containing genomic locations')
        self.args = parser.parse_args()
        self.writer = csv.writer(sys.stdout, dialect='excel-tab')

    def read_data(self):
        for row in csv.reader(iter(sys.stdin.readline, ''), dialect='excel-tab'):
            yield row

    def run(self):
        for assm_acc, seq_acc, start, stop, *extra in self.read_data():
            self.find_all_for_location(assm_acc, seq_acc, start, stop, extra)

    def process_loc_for_gff(self, zin, gff_fname, assm_acc, seq_acc, start, stop, extra_fields):
        with tempfile.NamedTemporaryFile() as tmpfile:
            tmpfile.write(zin.read(gff_fname))
            db = gffutils.create_db(
                tmpfile.name,
                dbfn=':memory:',
                force=True,
                keep_order=True,
                merge_strategy='merge',
                sort_attribute_values=True
            )
            find_genes_by_loc(db, self.writer, assm_acc, seq_acc, start, stop, extra_fields)

    def find_all_for_location(self, assm_acc, seq_acc, start, stop, extra_fields):
        zip_file = get_zip_file_for_acc(assm_acc, self.args.packages_dir)
        try:
            with zipfile.ZipFile(zip_file, 'r') as zin:
                catalog = retrieve_data_catalog(zin)
                gff_files = get_catalog_files(catalog, dataset_catalog_pb2.File.FileType.GFF3)
                for assm_acc, gff_files in gff_files.items():
                    report = retrieve_assembly_report(zin, catalog, assm_acc)
                    for gff_fname in gff_files:
                        self.process_loc_for_gff(zin, gff_fname, assm_acc, seq_acc, start, stop, extra_fields)
        except zipfile.BadZipFile:
            print(f'{zip_file} is not a zip file')


if __name__ == '__main__':
    FindGenesByLoc().run()

