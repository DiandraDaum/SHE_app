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

# Make download button for all comparisons
# Traffic analysis code
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
#@st.cache_data

st.set_page_config(
    page_title="SHE: Downloads",
    page_icon= "data/tree/logo_tree.png", 
    layout="centered"
)

st.title('SHE: Structural History of Eukarya')
#-----------
st.sidebar.write("Do you like it?")
sentiment_mapping = [":material/thumb_down:", ":material/thumb_up:"]
selected = st.sidebar.feedback("thumbs")
if selected is not None:
    st.sidebar.markdown(f"You selected: {sentiment_mapping[selected]}")

st.subheader('Downloads')
st.write('Here it is possible to download the data (structures and trees) from our analysis.')
st.write('Please, select the type of file you want to download from the sidebar.')

st.sidebar.subheader('Downloads', divider="yellow")


if st.sidebar.checkbox('Complete Foldcomb results by Proteome'):
    st.subheader('Foldcomb results for all Proteome')
    st.write('Here it is possible to download the Foldseek results per Proteome')
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
            target_folder_path = matching_paths[0] 
            zip_file_name = f"{proteome}_foldseek"
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
                    label="Download all bitscores for " + proteome + " (.zip)",
                    data=file, # The actual byte content of the zip file
                    file_name=f"{zip_file_name}.zip", # The name the user will see
                    mime="application/zip", 
                    icon=":material/download:",
                )
            # Clean up the generated zip file after use
            os.remove(zip_file_path)
        else:
            st.error(f"Proteome {proteome} not found.")

if st.sidebar.checkbox('Alignment scores by Organism'):
    st.subheader('Protein-protein alignments by Organism')
    st.write('Here it is possible to download the Foldseek proteome-proteome comparison bitscores (.csv) only by organism')
    proteome = st.sidebar.text_input(
    'Which proteome: UniProt ID?', key="proteome")
    st.write('Proteome selected for download: ', proteome)
    st.page_link("https://www.uniprot.org/proteomes/" + proteome, label="UniProt Proteomes", icon="🌎")
    if proteome == '':
        st.warning('Please input the UniProt proteome ID (i.e. UP000005640)')
    else: 
        if st.sidebar.toggle('Bitscore only'):
            root_dir_bitscores='data/decompressed' #350 gb
            file_pattern= f"{root_dir_bitscores}/{proteome}_*"
            matching_paths = glob.glob(file_pattern)
            if matching_paths:
                st.info("Zipping all the comparison bitscore csvs for your proteome ... Please wait some seconds for the download button to appear :)")
                target_folder_path = matching_paths[0] 
                zip_file_name = f"{proteome}_bitscores"
                # output path: './temp_archives/UP000005640_bitscores.zip'
                output_base_dir = "temp_archives" 
                os.makedirs(output_base_dir, exist_ok=True) 
                # results in 'temp_archives/UP000005640_bitscores.zip'
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
                        label="Download all bitscores for " + proteome + " (.zip)",
                        data=file, # The actual byte content of the zip file
                        file_name=f"{zip_file_name}.zip", # The name the user will see
                        mime="application/zip", # Critical: Set the correct MIME type for zip files
                        icon=":material/download:",
                    )
                # Clean up the generated zip file after use
                os.remove(zip_file_path)
            else:
                st.error("Proteome not found.")




if st.sidebar.checkbox('Protein structure predictions'):
    st.write('Here it is possible to download the predicted protein structures (.pdb) per Proteome')
    proteome = st.sidebar.text_input(
    'Which proteome: UniProt ID?', key="proteome")
    st.write('Proteome selected for download: ', proteome)
    st.page_link("https://www.uniprot.org/proteomes/" + proteome, label="UniProt Proteomes", icon="🌎")
    if st.sidebar.toggle('By proteome'):
        st.subheader('Protein structure predictions for all Proteome')
        if proteome == '':
            st.warning('Please input the UniProt proteome ID (i.e. UP000005640)')
        else:
            st.info("Zipping the proteome predictions folder ... Please wait some seconds for the download button to appear :)")
            subset_dir='data/subset/subset'
            proteome_dir = f"{subset_dir}/{proteome}_*/extracted_pdbs" 
            matching_paths = glob.glob(proteome_dir)
            target_folder_path = matching_paths[0]
            # temp zip pdbs folder for download
            zip_file_name = f"{proteome}_preds"
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
                label2="Download protein structure prediction for " + proteome
                st.download_button(
                    label="Download predictions for " + proteome + " (.zip)",
                    data=file, # The actual byte content of the zip file
                    file_name=f"{zip_file_name}.zip", # The name the user will see
                    mime="application/zip", # Critical: Set the correct MIME type for zip files
                    icon=":material/download:",
                )
                # Clean up the generated zip file after use
                os.remove(zip_file_path)
            else:
                st.error(f"Predictions for Proteome {proteome} not found.")

            
    elif st.sidebar.toggle('By protein'):
        st.subheader('Protein structure predictions by protein')
        st.write('Here it is possible to download a single protein structure prediction (.pdb), obtained from AFDB v4')
        protein = st.sidebar.text_input(
        'Which protein: UniProt ID?', key="protein")
        st.write('Protein structure selected for download: ', protein)
        st.page_link("https://www.uniprot.org/uniprotkb/" + protein + "/entry", label="UniProtKB", icon="🌎")
        if protein == '':
            st.warning('Please input the UniProt protein ID (i.e. A6NI61)')
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
            if not matching_rows.empty:
                retrieved_index = matching_rows['fc_db_index'].iloc[0]     #index in 'fc_db_index' column of the first matching row
            else:
                st.error(f"Structure prediction for protein {af_protein} not found.")
            
            #from  index get protein pdb --------------------------------------------------------------------
            subset_dir='data/subset/subset'
            pdb_dir = f"{subset_dir}/{proteome}_*/extracted_pdbs" 
            matching_paths = glob.glob(proteome_dir)
            # extract pdbs
            protein_pdb=pdb_dir[0].iloc[retrieved_index] #.pdb file
            label3="Download protein structure prediction for " + protein
            with open(protein_pdb, "rb") as file:
                st.download_button(
                    label=label3,
                    data=file,
                    file_name=af_protein + ".pdb",
                    icon=":material/download:",
                )


if st.sidebar.checkbox('Tree by Kingdom'):
    st.subheader('Tree (newick file) by Kingdom')
    st.write('Here it is possible to download the tree file (.newick) for each kigdom. To open it, try uploading it on iTOL (https://itol.embl.de/upload.cgi)!')
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




