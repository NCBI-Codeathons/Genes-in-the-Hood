import os
import yaml


class ThisApp:
    infile = os.path.join('data', 'cas9_neighborhood.yaml')
    outfile = 'cas9-global-freqs.md'

    def run(self):
        data = self.read_data(self.infile)
        self.write_markdown(self.outfile, data)

    def read_data(self, file_name):
        with open(file_name, 'r') as fin:
            print(f'Reading file {file_name} ...')
            data = yaml.load(fin)
            return data

    def global_levels(self, data):
        freqs = data['freqs']['global']
        top = freqs[0]['freq']
        return {x['term']: int(10 * (x['freq']/top)) for x in freqs}

    def write_markdown(self, file_name, data):
        with open(file_name, 'w') as fout:
            print(f'Writing file {file_name} ...')

            fout.write('# Genes in the (cas9) Hood\n\n')
            fout.write('The shading level indicates the global frequency in the cas9 gene neighborhood.\n\n')

            fout.write('|   |   |   |   |   |   |   |   |   |   | 0 |   |   |   |   |   |   |   |   |   |   |\n')
            fout.write('|---|---|---|---|---|---|---|---|---|---|----|---|---|---|---|---|---|---|---|---|---|\n')

            hoods = data['neighbors']
            levels = self.global_levels(data)
            for hood in hoods:
                for i in range(10 - len(hood['upstream'])):
                    fout.write('| ')
                for gene in hood['upstream']:
                    fout.write(self.gene_md(gene, levels))
                fout.write(self.gene_md(hood['gene'], levels))
                for gene in hood['downstream']:
                    fout.write(self.gene_md(gene, levels))
                fout.write('|\n')

    def gene_md(self, gene, levels):
        # GitHub won't show an md file if it is too big.  Adding the link puts it over the limit.
        # link = f'https://www.ncbi.nlm.nih.gov/protein/'
        name = gene['name']
        color = levels.get(name, 0) or levels.get(gene['protein_accessio'], 0)
        image = f'c{color}.png'
        return f'| ![alt text](img/{image} "{name}") '


if __name__ == '__main__':
    ThisApp().run()
