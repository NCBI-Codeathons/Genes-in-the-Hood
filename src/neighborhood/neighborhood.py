import argparse
import os
import zipfile


class ThisApp:
    input_path = os.path.join('var', 'data')
    output_path = os.path.join('var', 'data', 'neighborhood')
    gene = 'cas9'

    def run(self):
        zip_files = self.get_zip_files(self.input_path)
        for zip_file in zip_files:
            try:
                with zipfile.ZipFile(zip_file, 'r') as zin:
                    print(f'Reading {zip_file} ...')
            except zipfile.BadZipFile:
                print(f'{zip_file} is not a zip file')

    def get_zip_files(self, path):
        """
        Get a list of ZIP packages, each expected to contain ncbi_dataset for one NCBI assembly

        :param path: str - root of input directory
        :return: list of str - file paths for discovered ZIP files
        """
        zip_files = []
        for root, dirs, files in os.walk(path):
            zip_files += [os.path.join(root, x) for x in files if x.endswith('.zip')]
        return zip_files


if __name__ == '__main__':
    ThisApp().run()
