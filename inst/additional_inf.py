import argparse
import pandas as pd
from Bio import Entrez
from Bio import Medline


def main_foo(pmids_lst, output_filename):
    with open(pmids_lst) as f:
        pmid_list = f.read().splitlines()

    handle = Entrez.efetch(db='pubmed', id=pmid_list, rettype='medline', retmode='text')
    records = Medline.parse(handle)

    df = pd.DataFrame(records)
    df.to_csv(output_filename, sep="\t", index=False)



if __name__ == '__main__':
    Entrez.email = "Your.Name.Here@example.org"
    
    ap = argparse.ArgumentParser(description="Creates .tsv file containing additional info about PMID list.")
    ap.add_argument("-f", help="Path to PMID list file", type=str, required=True)
    ap.add_argument("-o", help="Path to output .tsv file", type=str, required=True)

    opts = ap.parse_args()
    PMID_file = opts.f
    output_file = opts.o

    main_foo(PMID_file, output_file)
    print("Parse Medline info completed")
