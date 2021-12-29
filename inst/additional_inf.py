import argparse
import pandas as pd
from Bio import Entrez
from Bio import Medline

def main_foo(pmids_lst):
    handle = Entrez.efetch(db='pubmed', id=pmids_lst, rettype='medline', retmode='text')
    records = Medline.parse(handle)

    df = pd.DataFrame(records)
    return(df)


if __name__ == '__main__':
    #print(help("modules"))
    Entrez.email = "Your.Name.Here@example.org"

