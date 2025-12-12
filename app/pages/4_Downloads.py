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



@st.cache_data
def set_background(image_path):
    """
    Adds a background image to the main Streamlit app container.
    """
    if not os.path.exists(image_path):
        st.error(f"Background image not found at: {image_path}")
        return

    with open(image_path, "rb") as f:
        data = f.read()
    
    # Convert image to Base64
    encoded_image = base64.b64encode(data).decode()

    # Inject CSS with the Base64 image
    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url("data:image/png;base64,{encoded_image}");
            background-size: cover;        /* Scale image to cover entire screen */
            background-position: center;   /* Center the image */
            background-repeat: no-repeat;  /* Do not repeat */
            background-attachment: fixed;  /* Keep background fixed when scrolling */
        }}
        </style>
        """,
        unsafe_allow_html=True
    )
##########page-----------------------------------------------
st.set_page_config(
    page_title="SHE: Downloads",
    page_icon= "data/tree/logo_tree.png", 
    layout="wide"
)
set_background('data/logos/SHE.png')
#st.title(':green[SHE: Structural History of Eukarya]')
st.subheader('Downloads')
st.write('Here it is possible to download the data (structures and trees) from our analysis.')
st.write('Please, below select the type of file you want to download.')
#----sidebar-------
st.sidebar.subheader('Downloads', divider="yellow")
st.sidebar.write("Do you like it?")
sentiment_mapping = [":material/thumb_down:", ":material/thumb_up:"]
selected = st.sidebar.feedback("thumbs")
if selected is not None:
    st.sidebar.markdown(f"You selected: {sentiment_mapping[selected]}")


#expanders
st.markdown("""
<style>
/* 1. Target the Expander Container (The 'details' HTML tag) */
[data-testid="stExpander"] {
    border: 1px solid rgba(250, 250, 250, 0.2); /* Subtle default border */
    border-radius: 12px;                        /* Rounded corners */
    transition: all 0.3s ease-in-out;           /* Smooth animation */
}

/* 2. HOVER STATE: The "Green Aura" */
[data-testid="stExpander"]:hover {
    border-color: #22c55e;   /* Brighter green border */
    
    /* The Aura Effect: Green Glow */
    box-shadow: 0 0 20px rgba(74, 222, 128, 0.4), 
                0 0 10px rgba(74, 222, 128, 0.2);
}

/* 3. OPTIONAL: Make the Title Green on Hover too */
[data-testid="stExpander"]:hover summary {
    color: #4ade80 !important;
}

/* 4. OPTIONAL: Customize the internal icon color */
[data-testid="stExpander"] summary svg {
    fill: #9ca3af !important; /* Default gray */
    transition: fill 0.3s ease-in-out;
}
[data-testid="stExpander"]:hover summary svg {
    fill: #4ade80 !important; /* Green on hover */
}
</style>
""", unsafe_allow_html=True)

with st.expander('Complete Foldcomp results by Proteome'):
    st.subheader('Foldcomp results for all Proteome')
    st.write('Here it is possible to download the Foldcomp results per Proteome. To read the files, see https://github.com/steineggerlab/foldcomp and our code on https://github.com/DiandraDaum/SHE_app/')
    proteome = st.text_input(
    'Which proteome: UniProt ID?', key="proteome1")
    st.write('Proteome selected for download: ', proteome)
    st.page_link("https://www.uniprot.org/proteomes/" + proteome, label="UniProt Proteomes", icon="🌎")
    if proteome == '':
        st.warning('Please input the UniProt proteome ID (i.e. UP000005640)')
    else:
        #st.warning('This function is going to be released soon')
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
            st.error(f"Proteome {proteome} not found. The **UniProt Proteomes** link above might help you find the ID you are looking for. Try again :)")

with st.expander('Alignment scores by Organism'):
    st.subheader('Protein-protein alignment bitscores by Organism')
    st.write('Here it is possible to download the Foldseek proteome-proteome comparison bitscores (.csv) only by organism')
    proteome = st.text_input(
    'Which proteome: UniProt ID?', key="proteome2")
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
            st.error(f"Proteome {proteome} not found. The **UniProt Proteomes** link above might help you find the ID you are looking for. Try again :)")




with st.expander('Protein Structure Predictions'):
    st.write('Here it is possible to download the predicted protein structures (.pdb) by protein ID, originally retrieved from AFDB v4')
    proteome = st.text_input(
    'Which proteome: UniProt ID?', key="proteome3")
    st.write('Proteome selected for download: ', proteome)
    st.page_link("https://www.uniprot.org/proteomes/" + proteome, label="UniProt Proteomes", icon="🌎")
    st.subheader('Protein structure predictions by protein')
    protein = st.text_input(
    'Which protein: UniProt ID?', key="protein")
    st.write('Protein structure selected for download: ', protein)
    st.page_link("https://www.uniprot.org/uniprotkb/" + protein + "/entry", label="UniProtKB", icon="🌎")
    if protein == '':
        st.warning('Please input the UniProt proteome ID (i.e. UP000000226) and the UniProt protein ID (i.e. V7ATE2) in the respective spaces')
    else:
        #find protein index in proteome ---------------------------------------------------------------
        index_dir='data/subset/index' 
        subset_dir='data/subset/subset'
        # --- 1. Find the Index File and Retrieved Index ---
        proteome_indexes = f"{index_dir}/{proteome}_*_fc_id_mapping.csv"
        matching_proteome = glob.glob(proteome_indexes)
        if not matching_proteome:
            st.error(f"No match found for proteome {proteome}. The **UniProt Proteomes** link above might help you find the ID you are looking for. Try again :)")
            st.stop()
        index_file_path = matching_proteome[0]
        
        try:
            proteome_df = pd.read_csv(index_file_path) # columns: 'fc_db_index' and 'afdb_id'
        except Exception as e:
            print(f"Error reading index file {index_file_path}: {e}")
            
        af_protein = f"AF-{protein}-F1-model_v4"
        matching_rows = proteome_df[proteome_df['afdb_id'] == af_protein]
        
        if matching_rows.empty:
            st.error(f"Structure prediction for protein {af_protein} not found in {proteome}. The **UniProtKB** link above might help you find the ID you are looking for. Try again :)")
            st.stop()
        
        # This is the correct database index
        retrieved_index = matching_rows['fc_db_index'].iloc[0] 
        
        # --- 2. Locate the Database Path ---
        wild_proteome_dir = f"{subset_dir}/{proteome}_*" #no .file
        proteome_dirs = glob.glob(wild_proteome_dir)
        if not proteome_dirs:
            print(f"Error: Could not find any directory matching the pattern: {wild_proteome_dir}")
            
        proteome_dir = proteome_dirs[0]
        db_dir = f"{proteome_dir}/extracted_pdbs" 
        os.makedirs(db_dir, exist_ok=True)
        full_proteome_id = Path(proteome_dir).name

        # Locate the full database file (e.g., UP000005640_9606_fc_db)
        db_name_wild = f"{proteome_dir}/{full_proteome_id}"
        db_name_list = glob.glob(db_name_wild)
        if not db_name_list:
            print(f"Error: Could not find any database file matching: {db_name_wild}")
        db_name = db_name_list[0]

        # --- 3. Directly Retrieve PDB using the Index ---
        st.info(f"c. Loading protein {af_protein} (database index {retrieved_index}). Please wait...")
        protein_found = False
        with fc.open(db_name) as db: 
            # Direct access using the index!
            target_name, pdb_dc = db[retrieved_index] 
        
        target_name = target_name.split('.')[0] 
        target_pdb_filepath = os.path.join(db_dir, target_name + ".pdb")

        #st.info("Extracting and writing file to disk...")
        with open(target_pdb_filepath, "w") as pdb_file:
            pdb_file.write(pdb_dc)
                    
        protein_found = True  # Stop iterating immediately after the target is found
        
        if protein_found:
            st.snow()
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


#if st.checkbox('Tree by Kingdom'):
with st.expander('Tree by Kingdom'):
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

# --- CITATION
st.subheader('Cite', divider="rainbow")
st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')

colm1, colm2, colm3, colm4 = st.columns(4)
with colm2:
    st.image("data/logos/scilifelab.png", use_container_width=False)
with colm3:
    st.image("data/logos/su.png", use_container_width=False)

