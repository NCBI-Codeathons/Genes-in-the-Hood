import argparse
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
        self.args = parser.parse_args()

    def run(self):
        zip_files = self.get_zip_files(self.args.input_path)
        self.process_zip_files(zip_files)

    def get_zip_files(self, path):
        zip_files = []
        for root, dirs, files in os.walk(path):
            zip_files += [os.path.join(root, x) for x in files if x.endswith('.zip')]
        return zip_files

    def process_zip_files(self, zip_files):
        print(f'Processing {len(zip_files)} assemblies ...')
        gene = self.args.gene
        report_file_1 = os.path.join(self.args.output_path, f'assembly_{gene}_report.txt')
        with open(report_file_1, 'w') as fout:
            for zip_file in zip_files:
                # print('ZIP', zip_file) ...
                acc = self.get_accession(zip_file)
                try:
                    with zipfile.ZipFile(zip_file, 'r') as zin:
                        gff_lines = self.read_file_lines(zin, f'ncbi_dataset/data/{acc}/genomic.gff')
                        seq, pos1, pos2, gene_type = self.find_gene(gene, gff_lines)
                        fout.write(f'{acc}\t{gene_type}\t{seq}\t{pos1}\t{pos2}\n')
                except zipfile.BadZipFile:
                    print(f'{zip_file} is not a zip file')
                    fout.write(f'{acc}\tError\n')
                except KeyError as err:
                    print(f'{acc} : {err}')
                    fout.write(f'{acc}\tError\n')

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


if __name__ == '__main__':
    ThisApp().run()
