# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 16:01:40 2024

@author: telma
"""
import argparse
import pandas as pd
import selenium
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
import pyperclip as pc
from datetime import datetime
from pathlib import Path
from tabulate import tabulate


# --------------------------------------------------------------

def get_shared_ids(phi_genes,veu):
    # get lists of matching species taxon IDs
    shared_taxonID = set.intersection(set(veu['Species NCBI taxon ID']),set(phi_genes['Pathogen_NCBI_species_Taxonomy ID']))

    # Get records with shared species taxon IDs
    shared_IDs = phi_genes.loc[phi_genes['Pathogen_NCBI_species_Taxonomy ID'].isin(shared_taxonID),['Pathogen_NCBI_species_Taxonomy ID','Pathogen_species','ProteinID','GeneLocusID', 'PHI_MolConn_ID']]

    # Add veupath project
    project_IDs = veu[['VEuPathDB Project','Species NCBI taxon ID']]
    project_IDs = project_IDs.rename(columns={'Species NCBI taxon ID':'Pathogen_NCBI_species_Taxonomy ID'})
    shared_IDs = shared_IDs.merge(project_IDs, how='left').drop_duplicates() 

    # make an index for each record
    index = []
    for i in range(1,len(shared_IDs)+1):
        index.append("PHI-Base_v4-16_"+str(i))
    
    shared_IDs['index']=index

    # Save the gene records
    shared_IDs.set_index('index').to_csv('PHIbase-VeupathDB-gene-records-to-check_'+datetime.today().strftime('%Y-%b-%d')+'.csv')
    print("Shared IDs have been saved to:")
    print('PHIbase-VeupathDB-gene-records-to-check_'+datetime.today().strftime('%Y-%b-%d')+'.csv\n')

    return shared_IDs
# --------------------------------------------------------------

# OPTIONAL use selenium to retrieve Gene ID matches from VeupathDB 
# otherise just copy-paste the protein and gene IDs from shared_IDs into Veupath manually

def automate_id_retrieval(shared_IDs):
    def get_veupath_IDs(keys):  
        # Using Chrome to access web
        url = "https://veupathdb.org/veupathdb/app/search/transcript/GeneByLocusTag"
        driver = webdriver.Firefox()
        
        # Open the website
        driver.get(url)
        driver.implicitly_wait(3)
        
        # Select the id box
        # If needed, use  “ctrl + shift + i” on the webpage to find and copy the CSS address for the textarea element
        id_box = driver.find_element(by=By.CSS_SELECTOR, value="html body.vpdb-Body div.main-stack div#wdk-container div.vpdb-RootContainer.vpdb-RootContainer__header-expanded.vpdb-RootContainer__news-collapsed.EuPathDB main.ebrc-Main.vpdb-Main div.wdk-QuestionForm div.wdk-TabsContainer div.wdk-TabContent form div div.wdk-QuestionFormParameterList div.wdk-QuestionFormParameterControl div.wdk-DatasetParam ul.wdk-DatasetParamSectionList li.wdk-DatasetParamSection div.wdk-DatasetParamControl div.wdk-DatasetParamIdList textarea")
        
        # put a list of IDs in the search box
        id_box.clear()
        pc.copy(keys)
        id_box.send_keys(Keys.CONTROL, 'v')
        
        
        # submit the search
        submit_button = driver.find_element(by=By.CSS_SELECTOR,value="#wdk-container > div > main > div > div.wdk-TabsContainer > div.wdk-TabContent > form > div.wdk-QuestionFormSubmitSection > button")
        submit_button.click()
        driver.implicitly_wait(5)
        
        # click download
        download_button = driver.find_element(by=By.XPATH,value='//*[@id="wdk-container"]/div/main/div/div[3]/div[2]/div[2]/div/div[2]/div/div[2]/div[1]/div/span[1]/div/a')
        download_button.click()
        
        # click tab/csv report 
        driver.implicitly_wait(3)
        choose_report = driver.find_element(by=By.CSS_SELECTOR, value="#wdk-container > div > main > div > div:nth-child(2) > div:nth-child(2) > ul > li:nth-child(1)")
        choose_report.click()
        
        # select download type (csv)
        driver.implicitly_wait(5)
        download_type = driver.find_element(by=By.CSS_SELECTOR, value="#wdk-container > div > main > div > div:nth-child(3) > div > div.eupathdb-ReporterFormWrapper > div.eupathdb-ReporterForm > div.eupathdb-ReporterFormGroup.eupathdb-ReporterFormGroup__otherOptions > div:nth-child(2) > div > ul > li:nth-child(2) > label > input[type=radio]")
        download_type.click()
        
        # download the data
        driver.implicitly_wait(3)
        get_genes = driver.find_element(by=By.CSS_SELECTOR, value="#wdk-container > div > main > div > div:nth-child(3) > div > div.eupathdb-ReporterFormWrapper > div.eupathdb-ReporterFormSubmit > div > button")
        driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        get_genes.click()
        driver.implicitly_wait(5)
        
        # End the driver
        # driver.quit()
        return
    # then get a list of IDs and input them into the function
    ProteinIDs = '\n'.join(str(x) for x in list(shared_IDs['ProteinID'].dropna()))
    GeneIDs = '\n'.join(str(x) for x in list(shared_IDs['GeneLocusID'].dropna()))
    get_veupath_IDs(ProteinIDs)
    get_veupath_IDs(GeneIDs)
    return

# --------------------------------------------------------------
 
# CROSS REFERENCE BACK TO VEUPATH
# matching taxonID records between PHIbase and Veupath that need checking:

def load_search_results():
    # load the veupath gene id search results
    print("Provide path for Protein ID search results:")
    print("Example: GeneByLocusTag_Summary.csv")
    protein_path = input()
    protein_path = Path(protein_path)
    print("Provide path for GeneLocus ID search results:")
    print("Example: GeneByLocusTag_Summary(1).csv")
    genelocus_path = input()
    genelocus_path = Path(genelocus_path)

    ProteinID_matches = pd.read_csv(protein_path)
    GeneID_matches = pd.read_csv(genelocus_path)
    return ProteinID_matches, GeneID_matches

def clean_data(ProteinID_matches, GeneID_matches, shared_IDs):

    # some results include two IDs in the input IDs. Separate them and strip blankspace
    
    ProteinID_matches['Input ID'] = ProteinID_matches['Input ID'].str.split(',')
    ProteinID_matches = ProteinID_matches.explode('Input ID')
    ProteinID_matches['Input ID'] = ProteinID_matches['Input ID'].str.strip()
    
    GeneID_matches['Input ID'] = GeneID_matches['Input ID'].str.split(',')
    GeneID_matches = GeneID_matches.explode('Input ID')
    GeneID_matches['Input ID'] = GeneID_matches['Input ID'].str.strip()


    # merge data sets and drop any without results
    #first match to protein results
    protein_merge = shared_IDs.merge(ProteinID_matches[['Gene ID', 'source_id', 'Organism','Input ID']], left_on='ProteinID', right_on='Input ID', how='outer')
    
    # then fill in the gaps with the gene ID results
    gene_merge = shared_IDs.merge(GeneID_matches[['Gene ID', 'source_id', 'Organism','Input ID']], left_on='GeneLocusID', right_on='Input ID', how='outer')
    
    # finally merge all the data and drop empty records     
    all_merge = protein_merge.merge(gene_merge, how='outer').dropna(subset=['Gene ID','index']).sort_values(by='Organism')
    all_merge.set_index('index').to_csv('PHIbase-veupath-IDs-with-matches_'+datetime.today().strftime('%Y-%b-%d')+'.csv')
    
    return all_merge

def make_linkout_list(all_merge):

    # MAKE A LIST OF LINKOUTS
    # Produce an updated data file (.txt file) with mapped VEuPathDB gene ID. \
    # This file can be used for extraction of data for loading of both phenotypes and link outs.
    Veu_PHI = all_merge[['Gene ID', 'PHI_MolConn_ID', 'VEuPathDB Project']]
    Veu_PHI = Veu_PHI.rename(columns={'Gene ID':'Locus ID', 'PHI_MolConn_ID':'PHI MolConn ID'})
    Veu_PHI.drop_duplicates(inplace=True)

    ## ADD MISSING LINKOUTS
    # cross reference this list with the old PHI linkouts and check discrepancies
    
    # These are the known discrepancies as of 15th Feb 2024:
    # - Rhizopus delemar does not exist on PHI-base, but Rhizopus arrhizus does. The opposite is true for Veupath. Protein identifers match between the Gene IDs and Uniprot lists them as synonyms.
    # - There are conflicts with Fusarium asiaticum and Fusarium graminearum. The phibase IDs for F. asiaticum are based on sequence dat a that technically comes from F. graminearum PH-1.
    
    print("\nProvide a list of the last set of linkouts as a tab delimited text file.")
    print("Example: PHI-base_linkouts_v4.12.txt")
    previous_linkouts_file = input()
    previous_linkouts = pd.read_csv(previous_linkouts_file, sep='\t')
    Unshared_identifiers = Veu_PHI.merge(previous_linkouts, left_on=['Locus ID', 'PHI MolConn ID'], right_on=list(previous_linkouts.columns), how='right')
    Unshared_identifiers = Unshared_identifiers.loc[Unshared_identifiers['Locus ID'].isna()]
    Unshared_identifiers.dropna(how='all').to_csv('Unshared_identifiers_'+datetime.today().strftime('%Y-%b-%d')+'.csv', index=None)
    
    
    print("\nProvide a list of known missing linkouts as a tab delimited text file.")
    print("Example: missing_linkouts_2024-Apr-11.tsv")
    missing_linkouts_file = input()
    missing_linkouts = pd.read_csv(missing_linkouts_file, sep='\t')
    
    # missing_linkouts = [('FGRAMPH1_01G01797','PHI:10551','FungiDB'),
    #('FGRAMPH1_01G13535','PHI:5353', 'FungiDB'),
    #('FGRAMPH1_01G17601','PHI:6603', 'FungiDB'),
    #('FGRAMPH1_01G24679','PHI:6604', 'FungiDB'),
    #('FGRAMPH1_01G26035','PHI:10550', 'FungiDB'),
    #('RO3G_03471','PHI:4844', 'FungiDB'),
    #('RO3G_05087','PHI:4842', 'FungiDB'),
    #('RO3G_08029','PHI:4142', 'FungiDB'),
    #('RO3G_11000','PHI:4843', 'FungiDB'),
    #('RO3G_11882','PHI:4143', 'FungiDB')]
    
    print('''These are the known discrepancies as of 15th Feb 2024:
    - Rhizopus delemar does not exist on PHI-base, but Rhizopus arrhizus does. The opposite is true for Veupath. Protein identifers match between the Gene IDs and Uniprot lists them as synonyms.
    - There are conflicts with Fusarium asiaticum and Fusarium graminearum. The phibase IDs for F. asiaticum are based on sequence dat a that technically comes from F. graminearum PH-1.
    ''')
    print()
    print("The following are additional linkouts that need to be added.")
    print("Please add to this list or remove as needed\n")
    print(tabulate(missing_linkouts, headers="keys", showindex=False))
    print("\nMake sure to check the gene IDs in 'Unshared_identifierd_"+datetime.today().strftime('%Y-%b-%d')+".csv'")
    print("Do you need to add anything to this list? (y/n)")
    answer = input()
    if answer =='y':
        print('''Provide a tab delimited text file of additional linkouts 
        with the same column headers as the table above:''')
        answer = input()
        missing_linkouts = missing_linkouts.append(pd.read_csv(answer, sep='\t'))
    elif answer == 'n':
        print("Ok, continuing...")

    Veu_PHI = Veu_PHI.merge(pd.DataFrame(missing_linkouts, columns=['Locus ID', 'PHI MolConn ID', 'VEuPathDB Project']), how = 'outer')
    Veu_PHI.drop_duplicates(inplace=True)
    Veu_PHI.to_csv('PHI-base_linkouts_v4-16_'+datetime.today().strftime('%Y-%b-%d')+'.tsv', sep='\t', index=None)

    # SPLIT UP BY PROJECT
    projects = Veu_PHI['VEuPathDB Project'].unique()
    filenames = []
    for i in projects:
        Veu_PHI.loc[Veu_PHI['VEuPathDB Project']==i,['Locus ID', 'PHI MolConn ID']].to_csv(i+'-phibase-v4.16-linkouts_'+datetime.today().strftime('%Y-%b-%d')+'.tsv',sep='\t', index=None)
        filenames.append([i+'-phibase-v4.16-linkouts_'+datetime.today().strftime('%Y-%b-%d')+'.tsv'])
    print("Linkout lists are separated by database and are saved as:")
    print(tabulate(filenames,headers=['']))
    return

def main():
    parser = argparse.ArgumentParser(
    description="Walks you through processing PHIbase data to making linkouts for Veupath sites",
    epilog="written by Helen Davison")

    print('''Download the current release of phibase data from: 
          https://github.com/PHI-base/data/tree/master/releases''')
    print("Provide the file path:")
    phi_path = Path(input()) # phi-base_v4-16_2023-11-30.csv
    phi_genes = pd.read_csv(phi_path) 
    
    print('''\nRetrieve a list of VEuPathDB species, strain and their taxon IDs accross all of VEuPathDB from: 
          https://veupathdb.org/veupathdb/app/search/organism/GenomeDataTypes/result''')
    print("Provide the file path:")
    veu_path = Path(input()) # GenomeDataTypes_Summary.csv
    veu = pd.read_csv(veu_path)
    
    shared_IDs = get_shared_ids(phi_genes,veu)
    ProteinID_matches, GeneID_matches = load_search_results()
    all_merge = clean_data(ProteinID_matches, GeneID_matches, shared_IDs)
    make_linkout_list(all_merge)
    
    return

main()