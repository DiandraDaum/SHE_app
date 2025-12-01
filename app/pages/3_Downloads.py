import streamlit as st
import pandas as pd
import numpy as np
import time
import csv
from ete3 import Tree, TreeStyle, NodeStyle
#from phytreeviz import TreeViz
import sys
import pathlib
import base64
import textwrap
import time
import os
import glob
import shutil
from pathlib import Path

import foldcomp as fc

# Make download button for all comparisons
# Traffic analysis code
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
#@st.cache_data

def zip_download(matching_paths, zip_file_name, label):
    target_folder_path = matching_paths[0] 
    #zip_file_name = f"{proteome}_foldseek"
    output_base_dir = "temp_archives" 
    os.makedirs(output_base_dir, exist_ok=True) 
    shutil.make_archive(
        base_name=os.path.join(output_base_dir, zip_file_name),
        format='zip',
        root_dir=os.path.dirname(target_folder_path), # Start zipping from the parent directory
        base_dir=os.path.basename(target_folder_path)  # Only include the target folder itself
    )
    # Get the full path to the newly created ZIP file and make download availabel
    zip_file_path = os.path.join(output_base_dir, f"{zip_file_name}.zip")
    with open(zip_file_path, "rb") as file:
        st.download_button(
            label=label,
            data=file, # The actual byte content of the zip file
            file_name=f"{zip_file_name}.zip", # The name the user will see
            mime="application/zip", 
            icon=":material/download:",
        )
    # Clean up the generated zip file after use
    os.remove(zip_file_path)




##########page-----------------------------------------------
st.set_page_config(
    page_title="SHE: Downloads",
    page_icon= "data/tree/logo_tree.png", 
    layout="centered"
)

st.title('SHE: Structural History of Eukarya')
#----sidebar-------
st.sidebar.write("Do you like it?")
sentiment_mapping = [":material/thumb_down:", ":material/thumb_up:"]
selected = st.sidebar.feedback("thumbs")
if selected is not None:
    st.sidebar.markdown(f"You selected: {sentiment_mapping[selected]}")

st.subheader('Downloads')
st.write('Here it is possible to download the data (structures and trees) from our analysis.')
st.write('Please, select the type of file you want to download from the sidebar.')

st.sidebar.subheader('Downloads', divider="yellow")


if st.sidebar.checkbox('Complete Foldcomp results by Proteome'):
    st.subheader('Foldcomp results for all Proteome')
    st.write('Here it is possible to download the Foldseek results per Proteome. To read the files, see https://github.com/steineggerlab/foldcomp and our code on https://github.com/DiandraDaum/SHE_app/')
    proteome = st.sidebar.text_input(
    'Which proteome: UniProt ID?', key="proteome")
    st.write('Proteome selected for download: ', proteome)
    st.page_link("https://www.uniprot.org/proteomes/" + proteome, label="UniProt Proteomes", icon="🌎")
    if proteome == '':
        st.warning('Please input the UniProt proteome ID (i.e. UP000005640)')
    else:
        st.warning('This function is going to be released soon')
        subset_dir='data/subset/subset'
        proteome_dir = f"{subset_dir}/{proteome}_*"  # temporarely unziped proteome
        # temp zip proteome folder for 
        matching_paths = glob.glob(proteome_dir)
        if matching_paths:
            st.info("Zipping ... Please wait some seconds for the download button to appear :)")
            st.snow()
            zip_file_name = f"{proteome}_foldseek"
            label1="Download Foldseek results for " + proteome + " (.zip)"
            zip_download(matching_paths, zip_file_name, label1)
        else:
            st.error(f"Proteome {proteome} not found.")

if st.sidebar.checkbox('Alignment scores by Organism'):
    st.subheader('Protein-protein alignment bitscores by Organism')
    st.write('Here it is possible to download the Foldseek proteome-proteome comparison bitscores (.csv) only by organism')
    proteome = st.sidebar.text_input(
    'Which proteome: UniProt ID?', key="proteome")
    st.write('Proteome selected for download: ', proteome)
    st.page_link("https://www.uniprot.org/proteomes/" + proteome, label="UniProt Proteomes", icon="🌎")
    if proteome == '':
        st.warning('Please input the UniProt proteome ID (i.e. UP000005640)')
    else: 
        root_dir_bitscores='data/decompressed' #350 gb
        file_pattern= f"{root_dir_bitscores}/{proteome}_*"
        matching_paths = glob.glob(file_pattern)
        if matching_paths:
            st.info("Zipping all the comparison bitscore csvs for your proteome ... Please wait some seconds for the download button to appear :)")
            st.snow()
            label2="Download all bitscores for " + proteome + " (.zip)"
            zip_file_name = f"{proteome}_bitscores"
            zip_download(matching_paths, zip_file_name, label2)
        else:
            st.error("Proteome not found.")




if st.sidebar.checkbox('Protein structure predictions'):
    st.write('Here it is possible to download the predicted protein structures (.pdb) by protein ID, originally retrieved from AFDB v4')
    proteome = st.sidebar.text_input(
    'Which proteome: UniProt ID?', key="proteome")
    st.write('Proteome selected for download: ', proteome)
    st.page_link("https://www.uniprot.org/proteomes/" + proteome, label="UniProt Proteomes", icon="🌎")
            
    if st.sidebar.toggle('By protein'):
        st.subheader('Protein structure predictions by protein')
        protein = st.sidebar.text_input(
        'Which protein: UniProt ID?', key="protein")
        st.write('Protein structure selected for download: ', protein)
        st.page_link("https://www.uniprot.org/uniprotkb/" + protein + "/entry", label="UniProtKB", icon="🌎")
        if protein == '':
            st.warning('Please input the UniProt proteome ID (i.e. UP000000226) and the UniProt protein ID (i.e. V7ATE2) in the respective spaces')
        else:
            #find protein index in proteome ---------------------------------------------------------------
            index_dir='data/subset/index' #140 gb zipped
            proteome_indexes= f"{index_dir}/{proteome}_*_fc_id_mapping.csv"
            matching_proteome = glob.glob(proteome_indexes)
            if not matching_proteome:
                st.error(f"Error: No matching proteome file found for {proteome}")
            index_file_path = matching_proteome[0]
            proteome_df = pd.read_csv(index_file_path) # columns: 'fc_db_index' and 'afdb_id'
            af_protein=f"AF-{protein}-F1-model_v4"
            matching_rows = proteome_df[proteome_df['afdb_id'] == af_protein]
            if matching_rows.empty:
                st.error(f"Structure prediction for protein {af_protein} not found in {proteome}.")
            
            retrieved_index = matching_rows['fc_db_index'].iloc[0]     #index in 'fc_db_index' column of the first matching row
            
            #from  index get protein pdb --------------------------------------------------------------------
            subset_dir='data/subset/subset'
            wild_proteome_dir = f"{subset_dir}/{proteome}_*"
            proteome_dirs = glob.glob(wild_proteome_dir)
            if proteome_dirs:
                proteome_dir = proteome_dirs[0]
            else:
                st.error(f"Error: Could not find any directory matching the pattern: {wild_proteome_dir}")
                sys.exit(1)
            db_dir = f"{proteome_dir}/extracted_pdbs" 
            db_dir_path = Path(db_dir)

            # Locate the full database file
            db_name_wild = f"{proteome_dir}/{proteome}_*"
            db_name = glob.glob(db_name_wild)[0]
            
            os.makedirs(db_dir, exist_ok=True)
            target_name = None
            target_pdb_filepath = None
            protein_found = False
            #st.info(f"Looking for protein {af_protein} (database index {retrieved_index}). Please wait...")
            protein_index = 0
            st.info("Please wait for the download button to appear :)")
            st.snow()
            with fc.open(db_name) as db: 
                # Iterate through the database entries
                for (name_full, pdb_dc) in db:
                    if protein_index == retrieved_index:
                        target_name = name_full.split('.')[0] 
                        target_pdb_filepath = os.path.join(db_dir, target_name + ".pdb")
                        
                        # B. Check if the file already exists on disk
                        if Path(target_pdb_filepath).exists():
                            st.info("Protein file already exists. Preparing download...")
                        else:
                            # C. If it doesn't exist, extract and write it to the db_dir folder
                            st.info("Extracting and writing file to disk...")
                            # Note: The 'pdb_dc' content is assumed to be the PDB/structure text
                            with open(target_pdb_filepath, "w") as pdb_file:
                                pdb_file.write(pdb_dc)
                        
                        protein_found = True
                        break # Stop iterating immediately after the target is found
                    protein_index += 1
            
            if protein_found:
                target_pdb_filepath = os.path.join(db_dir, af_protein + ".pdb")
                label4="Download protein structure prediction for " + protein
                # Download the file using the dynamically determined path and name
                with open(target_pdb_filepath, "rb") as file:
                    st.download_button(
                        label=label4,
                        data=file,
                        file_name=target_name + ".pdb", # Use the dynamically found name
                        mime="application/octet-stream",
                        icon=":material/download:",
                    )
            shutil.rmtree(db_dir)


if st.sidebar.checkbox('Tree by Kingdom'):
    st.subheader('Tree (newick file) by Kingdom')
    st.write('Here it is possible to download the tree file (.newick) for each kingdom. To open it, try uploading it on iTOL (https://itol.embl.de/upload.cgi)!')
    df = pd.DataFrame({
        'kingdom': ["All", "Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"],
        'number_species': [10, 20, 30, 40, 50, 70, 80] # change to actual number of species
        })

    for kingdom in df['kingdom']:
        file = "data/tree/" + kingdom + "/" + "kingdoms_" + kingdom + "_tree.newick"
        file_name= kingdom + "_tree.newick"
        label1="Download Tree " + kingdom
        label2="Download Tree Eukarya"
        if kingdom == "All":
            with open(file, "rb") as file:
                st.download_button(
                label=label2,
                data=file,
                file_name=file_name,
                mime="image/png",
                icon=":material/download:",
            )
        else:
            #print(f"Downloading file: {file}")
            with open(file, "rb") as file:
                st.download_button(
                    label=label1,
                    data=file,
                    file_name=file_name,
                    mime="image/png",
                    icon=":material/download:",
                )




