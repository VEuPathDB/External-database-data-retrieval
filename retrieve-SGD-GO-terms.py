import pandas as pd
from itertools import takewhile
import gzip
import argparse
from urllib.request import urlretrieve
from datetime import datetime
from tqdm import tqdm
tqdm.pandas()

def download_gaf(url):
    date_time = datetime.now().strftime("%Y-%m-%d") 
    filename='sgd'+date_time+'.gaf'
    urlretrieve(url, filename)
    return filename

def clean_sgd_gaf(file, outname):
    df = pd.read_csv(file,sep='\t', comment="!", header=None)

    identifiers = df[1]
    synonyms = list(df[10].str.split('|'))

    # save the comments 
    with open(file, 'r') as fobj:
        # takewhile returns an iterator over all the lines 
        # that start with the comment string
        headiter = takewhile(lambda s: s.startswith('!'), fobj) 
        # convert it to a list
        header = list(headiter)

    # Extract the veupath identifiers from the 11th column (index 10).
    # The current pattern is that it is 7 characters long and begins with Y. 
    # It should also not end with lowercase 'a' to remove synonyms ending with 'beta'
    # Check that this is still right.
    veu_id = [[ele for ele in sub if ((len(ele)==7) & (ele[0]=='Y') & (ele[-1]!='a'))] for sub in synonyms]

    sub_list=[]
    for sub in synonyms:
        ele_list=[]
        for ele in sub:
            if (len(ele)==7) & (ele[0]=='Y') & (ele[-1]!='a'):
                ele_list.append(ele)
            elif (len(ele)==9) & (ele[0]=='Y'):
                if (ele[7]=='-'):
                    ele_list.append(ele)
            elif (len(sub)==1) and (ele not in ele_list):
                ele_list.append('')
        sub_list.append(ele_list)

    # remove longer identifiers where two or more exist
    # and remove additional identifiers 
    # (these are usually aliased so should be accessible that way)
    def pop_off(sub_list):
        for sub in sub_list:
            if len(sub)>1:
                if (len(sub[0])==7) and (len(sub[0])!=7):
                    sub.pop(1)
                elif (len(sub[0])!=7) and (len(sub[0])==7):
                    sub.pop(0)
                else:
                    sub.pop(0)
        length = len([ele for ele in sub_list if len(ele)>1])
        return length

    popped = pop_off(sub_list)
    while popped>0:
        pop_off(sub_list)

    sub_list = [''.join(ele) for ele in sub_list]

    # map SGD and veupath identifiers
    # this might take a minute if it's a big gaf file
    mapping = dict(zip(identifiers,sub_list))
    mapping = {k:v for k,v in mapping.items() if v != ''}
    for key, value in tqdm(mapping.items()):
        df.loc[df[1]==key,10] = df.loc[df[1]==key,1].str.replace(key,value)
        df.loc[df[1]==key,1] = value


    # write the comments and dataframe to a tsv
    def write_gaf(df, header, outname):
        f = open(outname+"sgd.gaf", 'a', newline='',encoding="utf-8")
        for i in header:
            f.write(str(i))
        df.to_csv(f, sep='\t', index=None, header=None)
        f.close()
        return
    name="sgd-edited"
    write_gaf(df, header, outname)
    with open(outname+"sgd.gaf", 'rb') as src, gzip.open(outname+".gaf.gz", 'wb') as dst:
        dst.writelines(src)
    return

def main():
    parser = argparse.ArgumentParser(
    description="Reformat the SGDs GO annotations",
    epilog="written by Helen Davison")
    parser.add_argument('--gaf', \
                        help="The gaf file of SGD GO annotations",
                        required=True)
    parser.add_argument('-o','--outname', \
                        help="The name prefix to give your outputs",
                        required=True)
    args = parser.parse_args()

    file = args.gaf
    outname = args.outname
    print("\033[33m {}\033[0;0m".format("***REMINDER***"))
    print('Always check that the identifiers still match the pattern:')
    print('1. begins with capital Y')
    print('2. is seven characters long (unless it has a "-A" or "-B" suffix')
    clean_sgd_gaf(file, outname)
    return