import csv
import os
import re
import sys
import pandas as pd

def tsv2html(tsv_file_name, html_file_name):
    df = pd.read_csv(tsv_file_name,sep='\t', header=0)
    df.rename(columns={"#assembly": "assembly"}, inplace=True)
    df["assembly"] = df["assembly"].apply( # insert links
        lambda x: "<a href='https://www.ncbi.nlm.nih.gov/assembly/{}'>{}</a>".format(x, x)
    )
    df["accession"] = df["accession"].apply( # insert links
        lambda x: "<a href='https://www.ncbi.nlm.nih.gov/nucleotide/{}'>{}</a>".format(x, x)
    )
    pd.set_option('colheader_justify', 'center') 
    with open(html_file_name,'w') as html_file:
        html_file.write(
            df.to_html(
                index=False,
                render_links=True,
                escape=False,
            )
        )

if __name__ == "__main__":
   tsv2html(sys.argv[1], sys.argv[2])
