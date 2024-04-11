## retrieve data from the databases

from bs4 import BeautifulSoup
import urllib.request as urllib2
import re
import gzip
import shutil
import csv
import  pandas as pd
from datetime import datetime
from pathlib import Path


'''
Uses the 'database-retrieval-organism-list.xlsx' 
Make sure its all up to date!
CGD GO term retrieval specifically uses the 'veupath_name' and 'taxon_id' columns
'''

## FUNCTIONS 

# this function ca retrieve gff information automatically from CGD but is no 
# longer needed as Veupath only takes structural information from INSDC

def retrieve_cgd_gffs(df):
    external_db_names = df[df["retrieve_annotation"]==True]["external_db_name"]

    for external_db_name in external_db_names:
        url = "http://www.candidagenome.org/download/gff/"+external_db_name+"/"
        html = urllib2.urlopen(url).read()
        soup = BeautifulSoup(html)
        files = soup.find_all('a', href=True)

        if "_current_features.gff" not in str(files): 
            last_link = soup.find_all('a', href=True)[-2] # this avoids acrchive files
            latest_content = urllib2.urlopen(url+last_link['href']).read()
            latest_soup = BeautifulSoup(latest_content)
            files = latest_soup.find_all('a', href=True)

        gff_file = [s for s in files if "_current_features.gff" in str(s)][0]['href']

        html = url+last_link['href']+gff_file

        urllib2.urlretrieve(html, gff_file)

        # Candida albicans SC5314
        # product descriptions also have ORF identifiers in the “Note=” attribute. 
        # Remove the orf IDs that are between “%28” and “%29".
        if external_db_name =="C_albicans_SC5314":
            with open(gff_file) as fin, open("prepped-"+gff_file, "w+") as fout:
                for line in fin:
                    line = re.sub(r'%28.*%29', r'', line)
                    fout.write(line)
    return

# GO term retrieval 

def retrieve_cgd_GO(df):
    url = "http://current.geneontology.org/annotations/cgd.gaf.gz"
    names = df[(df["retrieve_GO"]==True)&(df["database"]=='CGD')][["veupath_name","taxon_id"]]
    names['taxon_id'] = names['taxon_id'].str.replace(' ','')
    urllib2.urlretrieve(url, "cgd.gaf.gz")
    with gzip.open("cgd.gaf.gz", 'rb') as f_in:
        with open("cgd.gaf", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    # load gaf as a dataframe
    gaf_df = pd.read_csv("cgd.gaf", sep="\t", comment="!", dtype='str', header=None).dropna(axis=1, how='all')
    gaf_df = gaf_df[gaf_df[0]=='CGD']
    # Move FungiDB gene IDs from column 11 to 2, and move the CGD IDs to column 11
    fungiDB_ids = []
    for index, row in gaf_df.iterrows():
        if str(row[10]) != "nan" : 
            fungiDB_id = row[10].split("|")[0]
            row[10] = row[10].replace(fungiDB_id, str(row[1]))
            fungiDB_ids.append(fungiDB_id)
        elif str(row[10]) == "nan":
            fungiDB_id = row[1]
            fungiDB_ids.append(fungiDB_id)
    gaf_df[1] = fungiDB_ids
    # remove taxon IDs that do not exist in veupath
    gaf_df = gaf_df.loc[gaf_df[12].isin(list(names['taxon_id']))]
    # split gaf file by organism using taxID (column 12)
    dfs = [y for x,y in gaf_df.groupby(12)]
    # make a directory for the files
    dirname = 'CGD-GO-terms_' + datetime.today().strftime('%Y-%b-%d')
    directory_path = Path.cwd()/dirname
    directory_path.mkdir()
    print(f"Successfully made the '{directory_path}' directory.")
    # write df to a tab delimited file with writer.writerow(['!gaf-version: 2.2']) at the start
    for index, row in names.iterrows():
        name=row['veupath_name']
        taxa=row['taxon_id']
        for i in dfs:
            if i[12].unique()==taxa:
                file= name+'_GO-terms_'+datetime.today().strftime('%Y-%b-%d')+'.gaf'
                filepath = directory_path/file
                gzfile = file+'.gz'
                filepathgz = directory_path/gzfile
                with open(filepath, 'w') as f:
                    f.write('!gaf-version: 2.2\n')
                i.to_csv(filepath, sep='\t', index=None, header=False, mode='a')
                with open(filepath, 'rb') as f_in:
                    with gzip.open(filepathgz, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
    return

def main():
    # list of organisms to reference when retrieving annotation information
    print("This script uses the 'database-retrieval-organism-list.xlsx'")
    print("Make sure its all up to date and is in your current directory!")
    print("CGD GO term retrieval specifically uses the 'veupath_name' and 'taxon_id' columns")
    print()
    print("If you have checked, enter 'y':")
    answer = input()
    if answer == 'y':
        df = pd.read_excel('database-retrieval-organism-list.xlsx', sheet_name='database-retrieval-organism-lis')
        retrieve_cgd_GO(df)
    else:
        print("\nYou did not answer yes :( ")
        print("See you later!")
    return

main()