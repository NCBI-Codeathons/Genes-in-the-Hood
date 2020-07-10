import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from enum import Enum, unique, auto
from pprint import PrettyPrinter
import json
import os
import re
import tempfile
from typing import List
import zipfile

import gffutils
from google.protobuf import json_format
import ncbi.datasets
from ncbi.datasets.v1alpha1 import dataset_catalog_pb2
from ncbi.datasets.v1alpha1.reports import assembly_pb2
from ncbi.datasets.reports.report_reader import DatasetsReportReader


def retrieve_assembly_report(zip_in, catalog) -> assembly_pb2.AssemblyDataReport:
    report_files = get_files(catalog, dataset_catalog_pb2.File.FileType.DATA_REPORT)
    if not report_files:
        return None
    for assm_acc, files in report_files.items():
        for fname in files:
            yaml = zip_in.read(fname)
            rpt_rdr = DatasetsReportReader()
            return rpt_rdr.assembly_report(yaml)


def retrieve_data_catalog(zip_in) -> dataset_catalog_pb2.Catalog:
    catalog_json = zip_in.read('ncbi_dataset/data/dataset_catalog.json')
    return json_format.Parse(catalog_json, dataset_catalog_pb2.Catalog())


def get_files(data_catalog: dataset_catalog_pb2.Catalog, desired_filetype: dataset_catalog_pb2.File.FileType):
    files = defaultdict(list)
    for assm in data_catalog.assemblies:
        acc = assm.accession
        for f in assm.files:
            filepath = os.path.join('ncbi_dataset', 'data', f.file_path)
            if f.file_type == desired_filetype:
                files[acc].append(filepath)
    return files


def get_zip_files_for_accs(accs, path):
    zip_files = []
    for acc in accs:
        fname = os.path.join(path, f'{acc}.zip')
        if os.path.isfile(fname):
                zip_files.append(fname)
    return zip_files


def get_zip_files(accfile_path, path):
    with open(accfile_path, 'r') as accfile:
        return get_zip_files_for_accs([acc.rstrip('\n') for acc in accfile], path)
    return None


def get_all_files(path):
    for root, dirs, files in os.walk(path):
        zip_files.extend([os.path.join(root, x) for x in files if x.endswith('.zip')])
    return zip_files


@dataclass
class Gene:
    id: str
    feat_type: str
    name: str
    chrom: str
    strand: str
    range_start: int
    range_stop: int
    protein_accession: str = None

    def __str__(self):
        return f'{self.chrom}\t{self.feat_type}\t{self.name}\t{self.range_start}\t{self.range_stop}\t{self.protein_accession}\n'

    def name_val(self):
        return self.protein_accession if self.protein_accession else self.name


@unique
class Justification(Enum):
    LEFT = auto()
    RIGHT = auto()


def justify_vals(vals, size, default_val, justify: Justification = Justification.LEFT):
    add = size - len(vals)
    if add < 1:
        # return a subset of vals
        if justify == Justification.LEFT:
            return vals[0:size]
        elif justify == Justification.RIGHT:
            return vals[-size:]
    # pad vals with additional elements
    if justify == Justification.LEFT:
        vals.extend([default_val] * add)
    elif justify == Justification.RIGHT:
        vals = ([default_val] * add) + vals
    return vals


@dataclass
class Assembly:
    assm_acc: str = ""
    tax_id: int = 0


@dataclass
class Neighbors:
    gene: Gene
    upstream: List[Gene] = field(default_factory=list)
    downstream: List[Gene] = field(default_factory=list)
    assm: Assembly = field(default_factory=Assembly)

    def _get_table(self, upstream_count=10, downstream_count=10):
        u = [g.name_val() for g in self.upstream]
        d = [g.name_val() for g in self.downstream]
        # justify
        u = justify_vals(u, upstream_count, None, Justification.RIGHT)
        d = justify_vals(d, downstream_count, None, Justification.LEFT)
        return [*u, self.gene.name, *d]

    def get_neighborhood(self, upstream_count=10, downstream_count=10):
        return self._get_table(upstream_count, downstream_count)

    def to_text(self):
        return f'{self.assm} -- {self.get_neighborhood()}\n'


class NeighborhoodFreqs:
    def __init__(self, upstream_count=10, downstream_count=10):
        self.all_freqs = Freqs()
        size = upstream_count + downstream_count + 1
        self.freqs = [Freqs() for n in range(size)]

    def add_terms(self, terms):
        for t, f in zip(terms, self.freqs):
            self.all_freqs.add_term(t)
            f.add_term(t)

    def get_freq_table(self, topn=10):
        return [f.get_relative_freqs(topn) for f in self.freqs]

    def get_global_freqs(self, topn=10):
        return self.all_freqs.get_relative_freqs(topn)


class Freqs(Counter):
    def add_term(self, term):
        if term:
            self[term] += 1

    def get_relative_freqs(self, topn=10):
        total_count = sum(self.values())
        return [(t, count / total_count) for t, count in self.most_common(topn)]


def extract_genes(gff3_db, desired_gene):
    neighbors = []
    gene_neighborhood = None
    found_gene = None
    genes_by_chrom = defaultdict(list)

    for feat_type in ['gene', 'pseudogene']:
        for gene in gff3_db.features_of_type(feat_type):
            gene_name = gene.attributes.get('Name', None)[0]
            prot_acc = None
            if gene.attributes['gene_biotype'][0] == 'protein_coding':
                cds = list(gff3_db.children(gene, featuretype='CDS'))
                prot_acc = cds[0].attributes.get('protein_id', None)[0]

            geneobj = Gene(
                gene.id,
                feat_type,
                gene_name,
                gene.chrom,
                gene.strand,
                gene.start,
                gene.stop,
                prot_acc,
            )
            genes_by_chrom[gene.chrom].append(geneobj)
            if gene_name == desired_gene:
                gene_neighborhood = Neighbors(geneobj)
    if gene_neighborhood:
        get_neighborhood_by_count(gene_neighborhood, genes_by_chrom[gene_neighborhood.gene.chrom], 10)
    return gene_neighborhood


def get_neighborhood_by_count(neighborhood, genes, count):
    sorted_genes = sorted(genes, key=lambda g: g.range_start)
    idx = sorted_genes.index(neighborhood.gene)
    offset = count + 1
    if idx > 0:
        idx_start = idx - offset if idx >= offset else 0
        neighborhood.upstream = sorted_genes[idx_start:idx - 1]
    if idx < len(sorted_genes):
        idx_stop = idx + offset
        if idx_stop >= len(sorted_genes):
            idx_stop = len(sorted_genes) - 1
        neighborhood.downstream = sorted_genes[idx + 1:idx_stop]


def get_neighborhood_by_pos(genes, seq, lo, hi):
    # TODO: here we want to feed gff_lines into gffutils.create_db()
    # not sure if it needs to be on disk or if we can use io.StringIO

    # Doing something quick-and-dirty for now
    genes = []
    for line in gff_lines:
        col = line.split('\t')
        if col[0] == seq:
            pos1, pos2, gene_type = int(col[3]), int(col[4]), col[2]
            if gene_type in ('gene', 'pseudogene') and (pos1 in range(lo, hi) or pos2 in range(lo, hi)):
                attr = col9_attributes(col[8])
                gene = attr.get('gene') or attr.get('product') or attr.get('Name') or attr.get('locus_tag')
                if gene_type == 'pseudogene':
                    gene += '(ps)'
                genes.append(gene)
    return '-'.join(genes)


def col9_attributes(col_9):
    parts = [x.split('=') for x in col_9.split(';')]
    return {x: y for (x, y) in parts}



class ThisApp:
    default_input_accfile_path = '/usr/local/data/cas9_gene.gcf'
    default_input_path = os.path.join('var', 'data', 'packages')
    default_output_path = os.path.join('var', 'data', 'neighborhood')
    default_gene = 'cas9'
    defult_window_bp = 10000

    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-I', '--input-path', type=str, default=self.default_input_path,
                            help=f'root of input data directory [{self.default_input_path}]')
        parser.add_argument('-O', '--output-path', type=str, default=self.default_output_path,
                            help=f'root of output data directory [{self.default_output_path}]')
        parser.add_argument('-g', '--gene', type=str, default=self.default_gene,
                            help=f'gene symbol [{self.default_gene}]')
        parser.add_argument('-w', '--window', type=int, default=self.defult_window_bp,
                            help=f'gene symbol [{self.defult_window_bp}]')
        parser.add_argument('-a', '--accession', type=str, default=None,
                            help='accession - limit the analysis to this accession [none]')
        parser.add_argument('-A', '--accession-file', type=str, default=self.default_input_accfile_path,
                            help='file of accessions - limit the analysis to these accessions')
        self.args = parser.parse_args()
        self.freqs = NeighborhoodFreqs()

    def run(self):
        zip_files = get_zip_files(self.args.accession_file, self.args.input_path)
        self.process_zip_files(zip_files)

    def process_zip_file(self, zip_file, gene, fout1, fout2):
        try:
            with zipfile.ZipFile(zip_file, 'r') as zin:
                catalog = retrieve_data_catalog(zin)
                report = retrieve_assembly_report(zin, catalog)
                gff_files = get_files(catalog, dataset_catalog_pb2.File.FileType.GFF3)
                for assm_acc, gff_files in gff_files.items():
                    print(f'processing assembly {assm_acc} for {report.species_name} {report.strain}')
                    for fname in gff_files:
                        with tempfile.NamedTemporaryFile() as tmpfile:
                            tmpfile.write(zin.read(fname))
                            db = gffutils.create_db(
                                tmpfile.name,
                                dbfn = ':memory:',
                                force=True,
                                keep_order=True,
                                merge_strategy='merge',
                                sort_attribute_values=True
                            )
                            found_gene = extract_genes(db, gene)
                            found_gene.assm.assm_acc = assm_acc
                            found_gene.assm.tax_id = report.tax_id
                            self.freqs.add_terms(found_gene.get_neighborhood())
                            fout1.write(str(found_gene.gene))
                            fout2.write(found_gene.to_text())
                            # seq, pos1, pos2, gene_type = find_gene(gene, gff_lines)
        except zipfile.BadZipFile:
            print(f'{zip_file} is not a zip file')
            fout1.write('--ERROR--')
            fout2.write('--ERROR--')

    def process_zip_files(self, zip_files, accessions=None):
        print(f'Processing {len(zip_files)} assemblies ...')
        gene = self.args.gene
        window = self.args.window
        report_file_1 = os.path.join(self.args.output_path, f'assembly_{gene}_report.txt')
        report_file_2 = os.path.join(self.args.output_path, f'neighborhood_{gene}_report.txt')
        with open(report_file_1, 'w') as fout1, open(report_file_2, 'w') as fout2:
            for zip_file in zip_files[:4]:
                self.process_zip_file(zip_file, gene, fout1, fout2)
        pp = PrettyPrinter()
        pp.pprint(self.freqs.get_freq_table())
        pp.pprint(self.freqs.get_global_freqs())


if __name__ == '__main__':
    ThisApp().run()
