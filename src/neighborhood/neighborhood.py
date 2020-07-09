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
        for zip_file in zip_files:
            # print(zip_file)
            acc = self.get_accession(zip_file)
            try:
                with zipfile.ZipFile(zip_file, 'r') as zin:
                    gff_lines = self.read_file_lines(zin, f'ncbi_dataset/data/{acc}/genomic.gff')
            except zipfile.BadZipFile:
                print(f'{zip_file} is not a zip file')

    def get_accession(self, file_name):
        # For now, extract from filename,  Ideally, it should come from a file within
        # the package, like ncbi_dataset/data/dataset_catalog.json
        parts = os.path.split(file_name)
        x = parts[-1].split('.')
        return f'{x[0]}.{x[1]}'

    def read_file_lines(self, zin, file_name):
        print(f'Extracting {file_name}')
        with zin.open(file_name) as fin:
            return fin.read()


if __name__ == '__main__':
    ThisApp().run()
