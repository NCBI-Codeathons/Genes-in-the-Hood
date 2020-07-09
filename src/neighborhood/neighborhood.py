import argparse
import gffutils
import os
import zipfile


class ThisApp:
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
        parser.add_argument('-A', '--accession-file', type=str, default=None,
                            help='file of accessions - limit the analysis to these accessions')
        self.args = parser.parse_args()

    def run(self):
        zip_files = self.get_zip_files(self.args.input_path)
        accessions = self.get_desired_accessions()
        self.process_zip_files(zip_files, accessions)

    def get_zip_files(self, path):
        all_zip_files = []
        for root, dirs, files in os.walk(path):
            all_zip_files += [os.path.join(root, x) for x in files if x.endswith('.zip')]
        return all_zip_files

    def get_desired_accessions(self):
        desired = []
        if self.args.accession:
            desired = [self.args.accession]
        elif self.args.accession_file:
            with open(self.args.accession_file, 'r') as fin:
                desired = fin.read().split('\n')
        return desired

    def process_zip_files(self, zip_files, accessions=None):
        count = len(accessions) if accessions else len(zip_files)
        if count == 0:
            print('There is nothing to do!')
            return
        print(f'Processing {count} assemblies ...')
        gene = self.args.gene
        window = self.args.window
        report_file_1 = os.path.join(self.args.output_path, f'assembly_{gene}_report.txt')
        report_file_2 = os.path.join(self.args.output_path, f'neighborhood_{gene}_report.txt')
        with open(report_file_1, 'w') as fout1, open(report_file_2, 'w') as fout2:
            for zip_file in zip_files:
                acc = self.get_accession(zip_file)
                if accessions and acc not in accessions:
                    continue
                try:
                    with zipfile.ZipFile(zip_file, 'r') as zin:
                        gff_lines = self.read_file_lines(zin, f'ncbi_dataset/data/{acc}/genomic.gff')
                        seq, pos1, pos2, gene_type = self.find_gene(gene, gff_lines)
                        fout1.write(f'{acc}\t{gene_type}\n')
                        if gene_type is not None:
                            print(f'{acc}: found {gene} {gene_type}')
                            lo, hi = pos1 - window, pos2 + window
                            hood = self.get_neighborhood(gff_lines, seq, lo, hi)
                            fout2.write(f'{acc}\t{hood}\n')
                except PermissionError as err:
                    print(f'{acc} : {err}')
                    fout1.write(f'{acc}\tError\n')
                except zipfile.BadZipFile:
                    print(f'{zip_file} is not a zip file')
                    fout1.write(f'{acc}\tError\n')
                except KeyError as err:
                    print(f'{acc} : {err}')
                    fout1.write(f'{acc}\tError\n')

    def get_accession(self, file_name):
        # For now, extract from filename,  Ideally, it should come from a file within
        # the package, like ncbi_dataset/data/dataset_catalog.json
        parts = os.path.split(file_name)
        x = parts[-1].split('.')
        return f'{x[0]}.{x[1]}'

    def read_file_lines(self, zin, file_name):
        # print(f'Extracting {file_name}')
        with zin.open(file_name, 'r') as fin:
            lines = [x.decode() for x in fin.readlines()]
            return [x for x in lines if not x[0] == '#']

    def find_gene(self, gene, gff_lines):
        gene_attr = f';gene={gene};'
        for line in gff_lines:
            col = line.split('\t')
            if col[2] in ('gene', 'pseudogene') and gene_attr in line:
                return col[0], int(col[3]), int(col[4]), col[2]
        return None, None, None, None

    def get_neighborhood(self, gff_lines, seq, lo, hi):
        # TODO: here we want to feed gff_lines into gffutils.create_db()
        # not sure if it needs to be on disk or if we can use io.StringIO

        # Doing something quick-and-dirty for now
        genes = []
        for line in gff_lines:
            col = line.split('\t')
            if col[0] == seq:
                pos1, pos2, gene_type = int(col[3]), int(col[4]), col[2]
                if gene_type in ('gene', 'pseudogene') and (pos1 in range(lo, hi) or pos2 in range(lo, hi)):
                    attr = self.attributes(col[8])
                    gene = attr.get('gene') or attr.get('product') or attr.get('Name') or attr.get('locus_tag')
                    if gene_type == 'pseudogene':
                        gene += '(ps)'
                    genes.append(gene)
        return '-'.join(genes)

    def attributes(self, col_9):
        parts = [x.split('=') for x in col_9.split(';')]
        return {x: y for (x, y) in parts}


if __name__ == '__main__':
    ThisApp().run()
