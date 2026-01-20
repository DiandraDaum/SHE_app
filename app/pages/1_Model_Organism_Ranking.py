import os
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
import streamlit as st
import pandas as pd
import numpy as np
import time
import csv
import sys
import pathlib
import base64
import textwrap
import time
from collections import defaultdict
from tqdm import tqdm
from datetime import datetime
import glob
from urllib.parse import urlparse
import py3Dmol
from stmol import showmol
from stmol import *
from IPython.display import HTML
from scipy.spatial.distance import pdist, squareform
from ete3 import Tree, TreeStyle, NodeStyle
import plotly.graph_objects as go
import plotly.express as px
#from ete3.webplugin import WebTreeApplication
#from ete3.treeview.web import WebTreeApplication
from Bio import Phylo
from Bio.Phylo import Consensus
from biotite.sequence.phylo import neighbor_joining
from io import StringIO
import copy
import shutil
import re
from pathlib import Path
import foldcomp as fc
import random
import requests
import subprocess
import gzip
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Structure, Model, PDBExceptions
from Bio.SVDSuperimposer import SVDSuperimposer

#Code to analyse:
# model organism search
#gene set search --> tree build
# Define navigation items for the custom header
NAV_ITEMS = [
    {"label": "SHE", "icon": "🏠", "path": "app.py"},
    {"label": "Model Organisms Ranking", "icon": "🐁", "path": "pages/1_Model_Organism_Ranking.py"},
    {"label": "Explore SHE tree", "icon": "🌳", "path": "pages/2_Explore_Tree.py"},
    {"label": "Interactive SHE Tree", "icon": "🔎", "path": "pages/3_Interactive_Tree.py"},
    {"label": "Downloads", "icon": ":material/download:", "path": "pages/4_Downloads.py"},
]
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
# save user email
EMAIL_FILE = 'data/email_results.csv'
def save_email(email):
    """Saves the email and current timestamp to the CSV file."""
    try:
        os.makedirs(os.path.dirname(EMAIL_FILE), exist_ok=True)
        file_exists = os.path.exists(EMAIL_FILE)
        
        if file_exists:
            with open(EMAIL_FILE, 'r', encoding='utf-8') as f:
                if email in f.read():
                    return "Exists" # Return a special flag for duplicates

        with open(EMAIL_FILE, 'a', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            if not file_exists:
                writer.writerow(['email', 'timestamp'])
            
            timestamp = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
            writer.writerow([email, timestamp])
        return True
    except Exception as e:
        st.error(f"An error occurred while saving: {e}")
        return False

def is_valid_email(email):
    """
    Checks if the email string matches a standard email pattern.
    """
    # This regex checks for: character sequence + @ + domain name + . + extension
    pattern = r"^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$"
    return re.match(pattern, email) is not None

@st.cache_data
def get_base64_of_bin_file(bin_file):
    """
    Reads a binary file and returns its Base64 encoded string.
    """
    try:
        with open(bin_file, 'rb') as f:
            data = f.read()
        return base64.b64encode(data).decode()
    except FileNotFoundError:
        st.error(f"Error: Video file not found at path: {bin_file}. Please check the path.")
        return None
    
@st.cache_data
def embed_video_background(video_path):
    """
    Embeds a local MP4 video as a fixed background using Base64 encoding.
    """
    video_base64 = get_base64_of_bin_file(video_path)
    
    if video_base64 is None:
        return # Stop execution if file wasn't found

    # 1. Construct the Base64 URL for the video source
    video_url = f"data:video/mp4;base64,{video_base64}"

    # 2. Inject the HTML/CSS for the background video
    st.markdown(
        f"""
        <style>
        /* This hides all Streamlit content (which we will put back later) */
        .stApp > header, .stApp > div:not(.video-container) {{
            z-index: 10; /* Bring content above the video */
            position: relative; 
            background: none; /* Make content background transparent */
        }}

        /* Video element to cover the entire screen */
        .video-background {{
            position: fixed;
            top: 0;
            left: 0;
            min-width: 100%;
            min-height: 100%;
            width: auto;
            height: auto;
            z-index: -100; /* Puts the video behind all other content */
            overflow: hidden;
        }}
        </style>
        
        <div class="video-container">
            <video class="video-background" autoplay muted loop>
                <source src="{video_url}" type="video/mp4">
                Your browser does not support the video tag.
            </video>
        </div>
        """,
        unsafe_allow_html=True
    )
video_path_bg = "data/tree/bg2.mp4" 
# Call the function to set the background
#embed_video_background(video_path_bg)

# find model organisms --------------------------------------------------------------------------------------------------------------------------
@st.cache_data
def load_mapping(mapping_file):
    try:
        mapping_df = pd.read_csv(mapping_file)
        return dict(zip(mapping_df['species_name'], mapping_df['id']))
    except Exception as e:
        print(f"Error reading mapping file: {e}")
        st.stop()
        sys.exit(1)
@st.cache_data
def load_mapping_proteome(mapping_file):
    try:
        mapping_df = pd.read_csv(mapping_file)
        return dict(zip(mapping_df['species_name'], mapping_df['proteome_id']))
    except Exception as e:
        print(f"Error reading mapping file: {e}")
        st.stop()
        sys.exit(1)
@st.cache_data
def find_comparison_csv(target_dir, compare_dir_id):
    """Finds the Human_vs_Species or Species_vs_Human CSV file."""
    search_pattern = os.path.join(target_dir, '**', f'*vs{compare_dir_id}_top_hits.csv')
    csv_files = glob.glob(search_pattern, recursive=True)
    
    if not csv_files:
        search_pattern_alt = os.path.join(target_dir, '**', f'*vs{compare_dir_id}*.csv')
        csv_files = glob.glob(search_pattern_alt, recursive=True)

    return csv_files[0] if csv_files else None

#@st.cache_data
def process_orthologs(csv_file, uniprot_ids, search_col='proteinA'):
    """
    Calculates scores for a comparison species.
    Returns: (Total Score, Presence Map, Score Vector)
    """
    total_score = 0.0
    presence_map = {uid: -1 for uid in uniprot_ids}
    score_map = {uid: 0.0 for uid in uniprot_ids} # Default 0.0 if missing

    best_hit_data=[]
    try:
        df = pd.read_csv(csv_file, header=None, usecols=[0, 1, 4], on_bad_lines='warn')
        df.columns = ['proteinA', 'proteinB', 'bit_score']
        
        for uniprot_id in uniprot_ids:
            # Search logic (your original feature)
            hits = df[df[search_col].str.contains(f"-{uniprot_id}-", na=False)]
            if hits.empty:
                hits = df[df[search_col].str.contains(uniprot_id, na=False)]
            
            if not hits.empty:
                best_score = hits['bit_score'].max()
                total_score += best_score

                # 2. Find the row(s) with the best score and take the first one
                best_hit_row = hits[hits['bit_score'] == best_score].iloc[0]
                # 3. Append the best hit row (which is a Series with 'proteinA', 'proteinB', 'bit_score' indices)
                best_hit_data.append(best_hit_row)
                
                # For iTOL
                if best_score > 0:
                    presence_map[uniprot_id] = 1
                
                # For Tree (Store actual score)
                score_map[uniprot_id] = best_score


    except Exception as e:
        print("Your input protein(s) cannot be found in the proteome of the input target organism... Try again with the right target protein(s): the UniProt ID has to be found in the proteome of the target species selected above.")
        st.stop()
        
    # Convert map to list (Vector) ensures order is identical for all species
    score_vector = [score_map[uid] for uid in uniprot_ids]
    best_hits_df = pd.DataFrame(best_hit_data)
    
    return total_score, presence_map, score_vector, best_hits_df

@st.cache_data
def process_orthologs2(csv_file, uniprot_ids, search_col='proteinA'):
    """
    Optimized: Calculates scores using vectorized grouping instead of iterative search.
    """
    # 1. Convert IDs to a set for instant O(1) lookups
    target_ids_set = set(uniprot_ids)
    
    # Initialize default results
    # We use a dictionary for fast updates, then convert to list at the end
    # Default score is 0.0, Default presence is -1 (absent)
    final_scores = {uid: 0.0 for uid in uniprot_ids}
    
    total_score = 0.0

    try:
        # Read the file
        # We only need the search column (protein name) and the score
        # Note: We read the whole file. If files are massive (>1GB), keep the chunking logic.
        # For standard blast outputs, reading the whole file is usually faster than chunking overhead.
        df = pd.read_csv(csv_file, header=None, usecols=[0, 1, 4], on_bad_lines='warn')
        df.columns = ['proteinA', 'proteinB', 'bit_score']
        
        # 2. Extract the ID from the column
        # Assuming format "AF-ID-..." -> Split by '-' and take index 1
        # This runs in C-speed (very fast)
        df['extracted_id'] = df[search_col].astype(str).str.split('-').str[1]
        
        # 3. Filter: Keep only rows where the ID is in our target list
        # This drastically reduces the data we need to process
        relevant_rows = df[df['extracted_id'].isin(target_ids_set)]
        
        if not relevant_rows.empty:
            # 4. Group by ID and find the max bit_score for each protein
            # This replaces the manual "if score > best_score" check
            max_scores = relevant_rows.groupby('extracted_id')['bit_score'].max()
            
            # 5. Update our results
            # Convert the pandas series to a dictionary and merge it
            found_scores = max_scores.to_dict()
            final_scores.update(found_scores)
            
            # Calculate total score sum
            total_score = max_scores.sum()

    except Exception as e:
        print(f"Error processing {csv_file}: {e}")
        # Return empty/default values on error
        return 0.0, {uid: -1 for uid in uniprot_ids}, [0.0] * len(uniprot_ids)
        
    # 6. Formatting for Output
    # Create the strictly ordered list (Vector) for the tree
    score_vector = [final_scores[uid] for uid in uniprot_ids]
    
    # Create the presence map for iTOL
    # Logic: 1 if score > 0, else -1
    presence_map = {uid: (1 if final_scores[uid] > 0 else -1) for uid in uniprot_ids}
    
    return total_score, presence_map, score_vector

def get_human_self_vector(target_dir, uniprot_ids):
    """
    [NEW FEATURE]
    Calculates the 'Perfect' vector for Human by finding the max self-score
    across all files. This allows Human to be plotted on the tree.
    """
    #st.info("Calculating Human self-scores for phylogenetic tree...")
    
    # Initialize with 0.0
    self_scores_map = {uid: 0.0 for uid in uniprot_ids}
    
    search_pattern = os.path.join(target_dir, '**', f'*vs*_top_hits.csv')
    all_csv_files = glob.glob(search_pattern, recursive=True)
    
    if not all_csv_files:
        return 0.0, [0.0]*len(uniprot_ids)

    # Optimization: We just need to fill the map with the max score found anywhere
    for csv_file in all_csv_files:
        try:
            # Read in chunks to save memory
            df_chunk = pd.read_csv(csv_file, header=None, usecols=[0, 4], on_bad_lines='warn',
                                 chunksize=100000, dtype={0: 'str', 4: 'float'})
            for chunk in df_chunk:
                chunk.columns = ['proteinA', 'bit_score']
                
                for uniprot_id in uniprot_ids:
                    # Check if we already have a high score (heuristic optimization)
                    # but to be safe and match original logic, we scan.
                    hits = chunk[chunk['proteinA'].str.contains(f"-{uniprot_id}-", na=False)]
                    if hits.empty:
                        hits = chunk[chunk['proteinA'].str.contains(uniprot_id, na=False)]
                    
                    if not hits.empty:
                        local_max = hits['bit_score'].max()
                        if local_max > self_scores_map[uniprot_id]:
                            self_scores_map[uniprot_id] = local_max
        except Exception:
            continue
            
    # Create the vector and total sum
    human_vector = [self_scores_map[uid] for uid in uniprot_ids]
    total_self_score = sum(human_vector)
    
    return total_self_score, human_vector

@st.cache_data
def sum_scores_for_orthologs(csv_file, uniprot_ids, search_col='proteinA'):
    """
    Loads a single CSV file and sums the *best hit* bit scores
    for each target UniProt ID.
    
    [FIX] Added 'search_col' to specify which column to search
    for the UniProt ID (for bi-directional search).
    """
    total_score = 0.0
    try:
        df = pd.read_csv(csv_file, header=None, usecols=[0, 1, 4], on_bad_lines='warn')
        df.columns = ['proteinA', 'proteinB', 'bit_score']
        
        for uniprot_id in uniprot_ids:
            # [FIX] Use the dynamic 'search_col' argument
            # This is 'proteinA' for a forward search
            # and 'proteinB' for a reverse search.
            hits = df[df[search_col].str.contains(f"-{uniprot_id}-", na=False)]
            
            if hits.empty:
                hits = df[df[search_col].str.contains(uniprot_id, na=False)]
            
            if not hits.empty:
                best_score = hits['bit_score'].max()
                total_score += best_score
            #else:
                #st.info(f"  > No hit found for {uniprot_id} in {csv_file}")

    except pd.errors.EmptyDataError:
        print(f"CSV file {csv_file} is empty. Skipping.")
    except Exception as e:
        print(f"Error processing CSV {csv_file}: {e}")
        
    return total_score

@st.cache_data
def get_self_scores(target_dir, target_dir_id, uniprot_ids):
    """
    Calculates the "self-score" for the target species, which will be
    the 1.0 baseline.
    
    [FIX] This function no longer looks for a "self-CSV".
    Instead, it scans all available CSVs to find the
    absolute max score for each target protein, which serves
    as the robust self-score baseline.
    """
    #st.info(f"Calculating self-scores for target {target_dir_id}...")
    
    # Find ALL comparison CSVs in the target directory
    search_pattern = os.path.join(target_dir, '**', f'*vs*_top_hits.csv')
    all_csv_files = glob.glob(search_pattern, recursive=True)
    
    if not all_csv_files:
        st.error(f"No comparison CSVs found at all in {target_dir}. Cannot calculate self-scores.")
        return 0.0

    total_self_score = 0.0
    
    for uniprot_id in uniprot_ids:
        max_score_for_this_protein = 0.0
        
        # Search all CSVs for the highest score for this one protein
        for csv_file in all_csv_files:
            try:
                # Load only the proteinA and bit_score columns
                df_chunk = pd.read_csv(csv_file, header=None, usecols=[0, 4], on_bad_lines='warn',
                                       chunksize=100000, dtype={0: 'str', 4: 'float'})
                
                for chunk in df_chunk:
                    chunk.columns = ['proteinA', 'bit_score']
                    
                    # Find hits for the UniProt ID
                    hits = chunk[chunk['proteinA'].str.contains(f"-{uniprot_id}-", na=False)]
                    if hits.empty:
                        hits = chunk[chunk['proteinA'].str.contains(uniprot_id, na=False)]
                    
                    if not hits.empty:
                        local_max = hits['bit_score'].max()
                        if local_max > max_score_for_this_protein:
                            max_score_for_this_protein = local_max
                            
            except pd.errors.EmptyDataError:
                # This file is empty, skip it
                continue
            except Exception as e:
                # Log other errors (e.g., file corrupted) but continue
                print(f"Could not fully read {csv_file} for self-score: {e}")
                continue
        
        if max_score_for_this_protein == 0.0:
            print(f"No score found for {uniprot_id} in any comparison file.")
        
        #st.info(f"  > Max score (baseline) for {uniprot_id} is: {max_score_for_this_protein}")
        total_self_score += max_score_for_this_protein

    return total_self_score


# reteive protein pdbs
@st.cache_data(show_spinner=False)
def get_pdb(protein, proteome):
    index_dir='data/subset/index' 
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
        print(f"Error: Could not find any directory matching the pattern: {wild_proteome_dir}")
        st.stop()
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
    st.info(f"b. Loading protein {af_protein} (database index {retrieved_index}). Please wait...")
    protein_index = 0
    with fc.open(db_name) as db: 
        # Iterate through the database entries
        for (name_full, pdb_dc) in db:
            if protein_index == retrieved_index:
                target_name = name_full.split('.')[0] 
                target_pdb_filepath = os.path.join(db_dir, target_name + ".pdb")
                # C. If it doesn't exist, extract and write it to the db_dir folder
                #st.info("Extracting and writing file to disk...")
                # Note: The 'pdb_dc' content is assumed to be the PDB/structure text
                with open(target_pdb_filepath, "w") as pdb_file:
                    pdb_file.write(pdb_dc)
                    st.success(f"Protein structure saved to {target_pdb_filepath}.")
                protein_found = True
                break # Stop iterating immediately after the target is found
            protein_index += 1
    if protein_found:
        #target_pdb_filepath = os.path.join(db_dir, af_protein + ".pdb")
        #st.success(f"Protein {af_protein} found and saved to {target_pdb_filepath}.")
        return target_pdb_filepath
    else:
        st.error(f"Protein with index {retrieved_index} not found in the database {target_pdb_filepath}.")
        return None


@st.cache_data(show_spinner=False)
def get_pdb2(protein, proteome):
    index_dir='data/subset/index' 
    subset_dir='data/subset/subset'

    # --- 1. Find the Index File and Retrieved Index ---
    proteome_indexes = f"{index_dir}/{proteome}_*_fc_id_mapping.csv"
    matching_proteome = glob.glob(proteome_indexes)
    if not matching_proteome:
        st.error(f"Error: No matching proteome index file found for {proteome}")
        return None
    index_file_path = matching_proteome[0]
    
    try:
        proteome_df = pd.read_csv(index_file_path) # columns: 'fc_db_index' and 'afdb_id'
    except Exception as e:
        st.error(f"Error reading index file {index_file_path}: {e}")
        return None
        
    af_protein = f"AF-{protein}-F1-model_v4"
    matching_rows = proteome_df[proteome_df['afdb_id'] == af_protein]
    
    if matching_rows.empty:
        st.error(f"Structure prediction for protein {af_protein} not found in {proteome}.")
        return None
    
    # This is the correct database index
    retrieved_index = matching_rows['fc_db_index'].iloc[0] 
    
    # --- 2. Locate the Database Path ---
    wild_proteome_dir = f"{subset_dir}/{proteome}_*" #no .file
    proteome_dirs = glob.glob(wild_proteome_dir)
    if not proteome_dirs:
        st.error(f"Error: Could not find any directory matching the pattern: {wild_proteome_dir}")
        return None
        
    proteome_dir = proteome_dirs[0]
    db_dir = f"{proteome_dir}/extracted_pdbs" 
    os.makedirs(db_dir, exist_ok=True)
    full_proteome_id = Path(proteome_dir).name

    # Locate the full database file (e.g., UP000005640_9606_fc_db)
    db_name_wild = f"{proteome_dir}/{full_proteome_id}"
    db_name_list = glob.glob(db_name_wild)
    if not db_name_list:
         st.error(f"Error: Could not find any database file matching: {db_name_wild}")
         return None
    db_name = db_name_list[0]

    # --- 3. Directly Retrieve PDB using the Index ---
    #st.info(f"b. Loading protein {af_protein} (database index {retrieved_index}). Please wait...")
    # Create a placeholder for the temporary message
    status_placeholder = st.empty()
    # Display the loading message in the placeholder
    status_placeholder.info(f"b. Loading protein {af_protein} (database index {retrieved_index})...")
    try:
        with fc.open(db_name) as db: 
            # Direct access using the index!
            target_name, pdb_dc = db[retrieved_index] 
        
        target_name = target_name.split('.')[0] 
        target_pdb_filepath = os.path.join(db_dir, target_name + ".pdb")

        #st.info("Extracting and writing file to disk...")
        with open(target_pdb_filepath, "w") as pdb_file:
            pdb_file.write(pdb_dc)
        status_placeholder.empty()
        return target_pdb_filepath
        
    except IndexError:
        status_placeholder.error(f"Error: Database index {retrieved_index} is out of bounds for the database file {db_name}.")
        return None
    except Exception as e:
        status_placeholder.error(f"An unexpected error occurred during FoldComp retrieval: {e}")
        return None

@st.cache_data(show_spinner=False)
def render_3dmol_pdb(color, structure_path_or_url: str):
    """
    Loads and renders a protein structure, supporting both local file paths 
    (PDB/CIF) and remote AlphaFold CIF URLs.
    """
    #st.info("Preparing protein structure for 3D rendering...")
    structure_data = None
    structure_format = None
    
    # 1. Determine source (URL vs. Local File) and load data
    if structure_path_or_url.startswith(('http://', 'https://')):
        # --- Handle Remote URL (Original Logic) ---
        st.info("Fetching structure from remote URL...")
        try:
            cif_response = requests.get(structure_path_or_url, timeout=30)
            cif_response.raise_for_status() 
            structure_data = cif_response.text
            structure_format = 'cif'
        except requests.exceptions.RequestException as e:
            st.error(f"Network error during structure download: {e}")
            return
            
    else:
        # --- Handle Local File Path ---
        local_filepath = structure_path_or_url
        if not os.path.exists(local_filepath):
            st.error(f"Local file not found at path: {local_filepath}")
            return

        # Determine format based on extension
        if local_filepath.lower().endswith(('.cif', '.mmcif')):
            structure_format = 'cif'
        elif local_filepath.lower().endswith('.pdb'):
            structure_format = 'pdb' # Use 'pdb' format for py3Dmol
        else:
            st.error("Unsupported local file format. Must be .pdb or .cif.")
            return

        try:
            # Read the content of the local file
            with open(local_filepath, 'r') as f:
                structure_data = f.read()
        except Exception as e:
            st.error(f"Error reading local structure file: {e}")
            return
    
    if not structure_data:
        st.error("Structure data is empty. Cannot render.")
        return

    # 2. Initialize the py3Dmol viewer object
    view = py3Dmol.view(width=800, height=500, js='https://3dmol.org/build/3Dmol.js')
    view.setBackgroundColor('black')
    # 3. Add the model using the determined format
    view.addModel(structure_data, structure_format) # Use the determined format
    #view.zoomTo()

    # 4. Set the coloring scheme (lDDT for AlphaFold)
    if color == "plDDT":
        view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
        #plot_plddt_legend().show()
    elif color == "chain":
        view.setStyle({'chain': 'A'}, {'cartoon': {'color': 'cyan'}})
        view.setStyle({'chain': 'B'}, {'cartoon': {'color': 'yellow'}})
        view.setStyle({'chain': 'C'}, {'cartoon': {'color': 'magenta'}})
    
    #view.replicateModel() 
    # --- The Fix for Centering ---
    #view.resetView()
    view.zoomTo()

    # 5. Render the view using the stmol component
    try:
        #showmol(view, height=510, width=800)
        showmol(view, height=600, width=600)
    except Exception as e:
         st.error(f"Error rendering 3Dmol viewer with stmol: {e}")


# -----protein alignment viewer -----
#FOLDSEEK_EXEC_PATH = "./foldseek_local"
FOLDSEEK_EXEC_PATH = os.path.abspath("./foldseek_local")

@st.cache_data(show_spinner=False)
def run_foldseek_alignment_only(query_path, target_path):
    """
    Runs the standalone Foldseek binary.
    """
    
    # 1. Setup Directories
    temp_OUTPUT_DIR = "temp_archives/temp_protein_alignments"
    if not os.path.exists(temp_OUTPUT_DIR):
        os.makedirs(temp_OUTPUT_DIR)
    temp_m8 = os.path.join(temp_OUTPUT_DIR, "temp_aln.m8")
    
    abs_query_path = os.path.abspath(query_path)
    abs_target_path = os.path.abspath(target_path)

    # 2. Verify Binary Exists (Sanity Check)
    if not os.path.exists(FOLDSEEK_EXEC_PATH):
        st.error("CRITICAL ERROR: 'foldseek_local' not found.")
        return None

    # 3. Define Command
    cmd = [
        FOLDSEEK_EXEC_PATH, "easy-search", 
        abs_query_path, abs_target_path, temp_m8, "tmp_fs",
        "-s", "7.5", 
        "--format-output", "qstart,tstart,qaln,taln,alntmscore",
        "--threads", "1", "-v", "0"
    ]
    
    #st.info(f"Running Foldseek alignment...")
    status_placeholder = st.empty()
    # Display the loading message in the placeholder
    status_placeholder.info(f"Running Foldseek alignment...")
    
    try:
        # 4. Run Subprocess (No 'env' parameter needed for standalone binary)
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300) 
        
        # 5. Check Results
        if not os.path.exists(temp_m8) or os.path.getsize(temp_m8) == 0:
            #status_placeholder.warning("Foldseek ran but found **NO STRUCTURAL HITS**.")
            return None
            
        df = pd.read_csv(temp_m8, sep="\t", names=["qstart", "tstart", "qaln", "taln", "tm"])
        
        if df.empty:
            #status_placeholder.warning("Foldseek returned an empty alignment table.")
            return None
            
        #st.success("Foldseek alignment completed successfully!")
        status_placeholder.empty()
        return df.iloc[0]
        
    except subprocess.CalledProcessError as e:
        #status_placeholder.error(f"FOLDSEEK FAILED.\n\nStderr: {e.stderr}")
        return None
    except Exception as e:
        #status_placeholder.error(f"An unexpected error occurred: {e}")
        return None

@st.cache_data
def get_format_from_path(filepath):
    """Simple helper to determine PDB/CIF format from file extension."""
    if filepath.lower().endswith(('.cif', '.mmcif')):
        return "cif"
    elif filepath.lower().endswith('.pdb'):
        return "pdb"
    return None

@st.cache_data(show_spinner=False)
def calculate_svd_and_merge(target_file, target_fmt, query_file, query_fmt, aln_data, out_name):
    """
    Uses Bio.SVDSuperimposer to align the Query (moving) structure to the Target (fixed) structure 
    based on Foldseek's alignment mapping.
    """
    # 1. Load Structures
    p_t = MMCIFParser(QUIET=True) if target_fmt == "cif" else PDBParser(QUIET=True)
    p_q = MMCIFParser(QUIET=True) if query_fmt == "cif" else PDBParser(QUIET=True)
    
    try:
        s_t = p_t.get_structure("T", target_file) # Target (fixed)
        s_q = p_q.get_structure("Q", query_file) # Query (moving)
        
        atoms_t = [r['CA'] for r in s_t.get_residues() if 'CA' in r]
        atoms_q = [r['CA'] for r in s_q.get_residues() if 'CA' in r]

        # 2. Map Atoms using Alignment Strings
        qaln = aln_data['qaln'] # Query Sequence (moving)
        taln = aln_data['taln'] # Target Sequence (fixed)
        
        # Start indices (Foldseek is 1-based, Python list is 0-based)
        q_idx = int(aln_data['qstart']) - 1
        t_idx = int(aln_data['tstart']) - 1
        
        ref_coords = [] # Fixed (Target)
        mov_coords = [] # Moving (Query)
        
        for q_char, t_char in zip(qaln, taln):
            current_q_atom = None
            current_t_atom = None
            
            # Advance pointers for Query (moving)
            if q_char != '-':
                if q_idx < len(atoms_q): current_q_atom = atoms_q[q_idx]
                q_idx += 1
            
            # Advance pointers for Target (fixed)
            if t_char != '-':
                if t_idx < len(atoms_t): current_t_atom = atoms_t[t_idx]
                t_idx += 1
            
            # If both exist and are not gaps, they are an aligned pair!
            if q_char != '-' and t_char != '-' and current_q_atom and current_t_atom:
                ref_coords.append(current_t_atom.get_coord())
                mov_coords.append(current_q_atom.get_coord())

        # 3. Run SVD Superimposer
        if len(ref_coords) < 3:
            # print("  Error: Not enough aligned residues for SVD.")
            return False, 0 # Return failure and 0 RMSD

        sup = SVDSuperimposer()
        sup.set(np.array(ref_coords), np.array(mov_coords))
        sup.run()
        rot, tra = sup.get_rotran()
        rms = sup.get_rms()

        # 4. Apply to Query (Whole structure)
        # Note: AlphaFold PDB/CIFs often have only one chain
        c_query = list(s_q.get_chains())[0]
        for atom in c_query.get_atoms():
            # SVDSuperimposer gives rot/tran suitable for: new = dot(old, rot) + tran
            atom.set_coord(np.dot(atom.get_coord(), rot) + tra)

        # 5. Save Merged (Target Chain A, Query Chain B)
        merged = Structure.Structure("Merged")
        model = Model.Model(0)
        merged.add(model)
        
        c_target = list(s_t.get_chains())[0]; c_target.id = "A" # Fixed Chain
        c_query.id = "B" # Moving Chain (now superimposed)
        
        model.add(c_target)
        model.add(c_query)
        
        io = PDBIO()
        io.set_structure(merged)
        io.save(out_name)
        return True, rms # Return success and RMSD

    except Exception as e:
        # st.error(f"SVD Error: {e}")
        return False, 0

@st.cache_data
def render_3dmol_alignment(color, structure_path, target_name, target_species, query_name, model_species, tm_score, rmsd, spin_on=True):
    """
    Renders the merged PDB file, applying different colors to the two chains 
    (Target: light blue, Query/Model: yellow).
    """
    if not os.path.exists(structure_path):
        st.error(f"Alignment file not found at path: {structure_path}")
        return

    # Load data from local file (format is always 'pdb' here after merge)
    with open(structure_path, 'r') as f:
        structure_data = f.read()
    
    # 1. Initialize the py3Dmol viewer object
    view = py3Dmol.view(width=800, height=500, js='https://3dmol.org/build/3Dmol.js')
    view.setBackgroundColor('black')
    # 2. Add the model
    view.addModel(structure_data, 'pdb')
    # 3. Set styles and colors for Chain A (Target) and Chain B (Query)
    #view.setStyle({'chain': 'A'}, {'cartoon': {'color': 'lightskyblue'}})
    #view.setStyle({'chain': 'B'}, {'cartoon': {'color': 'yellow'}})
    if color == "plDDT":
        view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
    #plot_plddt_legend().show()
    elif color == "chain":
        view.setStyle({'chain': 'A'}, {'cartoon': {'color': 'cyan'}})
        view.setStyle({'chain': 'B'}, {'cartoon': {'color': 'yellow'}})
        view.setStyle({'chain': 'C'}, {'cartoon': {'color': 'magenta'}})
    view.zoomTo() 
    if spin_on:
        # 'axis' can be 'x', 'y', or 'z'. 'speed' is the rotation speed.
        view.spin(True) 
    else:
        view.spin(False)
    # 4. Add labels/info
    st.markdown(f"""
        **Alignment Details:**
        * **TM-Score:** `{tm_score:.4f}` 
        * **RMSD:** `{rmsd:.4f}` Å
        * **Chain legend:** **{target_species}** **{target_name}** (**cyan**) vs. **{model_species}** **{query_name}** (**yellow**)
        * **plDDT legend: >90 (dark blue), 81-90 (cyan), 71-80 (green), 61-70 (yellow), 51-60 (orange), <51 (red)**
    """)

    # 5. Render the view
    try:
        showmol(view, height=510, width=800)
    except Exception as e:
         st.error(f"Error rendering 3Dmol viewer with stmol: {e}")


# fast change in display
@st.fragment
def render_interactive_structure(key_suffix, out_file, target_name, target_species, query_name, model_species, tm_score, rmsd):
    """
    This function handles ONLY the UI for the structure viewer.
    Interacting with widgets here will NOT reload the whole page.
    """
    c1, c2 = st.columns([1, 2])
    
    with c1:
        alignment_color = st.radio(
            "Color by:",
            ('chain', 'plDDT'),
            index=0, 
            horizontal=True,
            key=f"color_radio_{key_suffix}"
        )
        
        enable_spin = st.checkbox(
            "Spin Structure", 
            value=True, 
            key=f"spin_check_{key_suffix}"
        )

    # Call your existing 3Dmol render function
    # Note: Ensure render_3dmol_alignment handles the path string correctly
    render_3dmol_alignment(
        alignment_color,
        structure_path=out_file, 
        target_name=target_name, 
        target_species=target_species,
        query_name=query_name, 
        model_species=model_species,
        tm_score=tm_score, 
        rmsd=rmsd,
        spin_on=enable_spin
    )




@st.cache_data(show_spinner=False)
def prepare_ranking_data(target_species, model_organisms, protein_list_str):
    """
    Performs all data parsing, ranking, and merging.
    Returns the final DataFrames or None if an error occurs.
    """
    # 1. Setup
    uniprot_ids = protein_list_str.split()
    root_dir = 'data/decompressed/'
    mapping_file = 'data/core/model/extracted_and_renamed.csv'
    
    # Load mapping
    species_map = load_mapping(mapping_file)
    
    # Check Target
    target_dir_id = species_map.get(target_species)
    if not target_dir_id:
        return None, None, f"Target species '{target_species}' not found in mapping."

    target_dir = os.path.join(root_dir, target_dir_id)
    
    # Storage
    ranking_scores = {}
    itol_data = {}
    tree_data = {} 
    all_best_hits_dfs = [] 

    # 2. Calculate Human (Target) Data
    human_score, human_vector = get_human_self_vector(target_dir, uniprot_ids)
    ranking_scores[target_species] = human_score
    tree_data[target_species] = human_vector
    itol_data[target_species] = {uid: 1 for uid in uniprot_ids}

    # 3. Calculate Comparison Species
    for species in model_organisms:
        if species not in species_map: 
            continue
        
        compare_dir_id = species_map[species]
        csv_file = find_comparison_csv(target_dir, compare_dir_id)
        search_col = 'proteinA'
        
        # Try reverse direction if needed
        if not csv_file:
            reverse_dir = os.path.join(root_dir, compare_dir_id)
            if os.path.isdir(reverse_dir):
                csv_file = find_comparison_csv(reverse_dir, target_dir_id)
                search_col = 'proteinB'

        if csv_file:
            try:
                score, presence, vector, best_hits_df = process_orthologs2(csv_file, uniprot_ids, search_col)
            except Exception:
                best_hits_df = pd.DataFrame()
                st.warning(f"The protein ID **{', '.join(uniprot_ids)}** was not found in the reference proteome of **{target_species}**. This can happen if the protein ID is not reviewed in the reference proteome used for AFDB v4 predictions.")
                cols = st.columns(len(uniprot_ids))
                for i, uid in enumerate(uniprot_ids):
                    with cols[i]: # Prevent creating too many columns if list is long
                        st.page_link(f"https://www.uniprot.org/uniprotkb/{uid}/entry", label=f"Check {uid} existence in {target_species} proteome and potential updated reviewed ID on UniProt", icon="🌎")
                continue
            
            if search_col == 'proteinB':
                # This suggests the Input ID exists, but isn't in the provided AFDB subset.
                try:
                    temp_col = best_hits_df['proteinA'].copy()
                    best_hits_df['proteinA'] = best_hits_df['proteinB']
                    best_hits_df['proteinB'] = temp_col
                except Exception:
                    continue
            
            best_hits_df['Compared_Species'] = species 
            all_best_hits_dfs.append(best_hits_df)
            ranking_scores[species] = score
            itol_data[species] = presence
            tree_data[species] = vector
        else:
            ranking_scores[species] = 0.0
            itol_data[species] = {uid: -1 for uid in uniprot_ids}
            tree_data[species] = [0.0] * len(uniprot_ids)
            st.warning(f"No match found in {species}")

    # 4. Merge Results
    if all_best_hits_dfs:
        comparison_df = pd.concat(all_best_hits_dfs, ignore_index=True)
    #else: comparison_df = pd.DataFrame()

    # Build Rank DataFrame
    df_rank = pd.DataFrame.from_dict(ranking_scores, orient='index', columns=['Score'])
    df_target = df_rank.loc[[target_species]]
    df_other = df_rank.drop(target_species).sort_values(by='Score', ascending=False)
    
    df_rank_sorted = pd.concat([df_target, df_other])
    df_rank = df_rank_sorted.reset_index()
    df_rank.columns = ['Organism', 'Score']

    # Build Presence DataFrame
    presence_list = []
    for sp, pres in itol_data.items():
        row_data = {'Organism': sp}
        for uid in uniprot_ids:
            row_data[uid] = 1 if pres.get(uid, 0) == 1 else 0
        presence_list.append(row_data)

    df_presence = pd.DataFrame(presence_list)

    # Final Merge
    df_rank_filtered = df_rank[df_rank['Organism'] != target_species].copy()
    df_rank_filtered.reset_index(drop=True, inplace=True)
    df_rank_filtered.insert(0, 'Rank', df_rank_filtered.index + 1)
    
    df_final = pd.merge(df_rank_filtered, df_presence, on='Organism', how='left')
    
    return df_final, comparison_df, None  # None means no error

# hadles input protein mismatches
@st.cache_data(show_spinner=False)
def get_uniprot_suggestion(invalid_id, target_species_fmt):
    """
    Queries UniProt to find out why an ID is invalid and suggests a fix.
    Handles: Outdated IDs, Species Mismatches, and Unreviewed (TrEMBL) entries.
    """
    # Clean the target species name (e.g., "Homo_sapiens_Human" -> "Homo sapiens")
    clean_target = ' '.join(target_species_fmt.split('_')[:2])
    
    base_url = "https://rest.uniprot.org/uniprotkb"
    
    try:
        # 1. Check if the ID exists at all (Look up the ID)
        response = requests.get(f"{base_url}/{invalid_id}")
        
        if response.status_code == 404:
            return None, "ID does not exist in UniProt: please check the correctness of your UniProt ID."
        
        if response.status_code != 200:
            return None, "UniProt API error: please check the correctness of your UniProt ID."
            
        data = response.json()
        
        # 2. Check if it is a Secondary/Outdated Accession
        primary_accession = data.get('primaryAccession')
        if primary_accession and primary_accession != invalid_id:
            return primary_accession, f"The ID **{invalid_id}** is outdated. The current primary ID is **{primary_accession}**."

        # Extract Gene and Organism info for further checks
        input_gene = None
        if 'genes' in data and len(data['genes']) > 0:
            input_gene = data['genes'][0].get('geneName', {}).get('value', '')
            
        org_sci_name = data.get('organism', {}).get('scientificName', '')
        input_tax_id = data.get('organism', {}).get('taxonId')

        # 3. Check Species Mismatch (e.g., User put Mouse ID, Target is Human)
        if clean_target.lower() not in org_sci_name.lower():
            # The ID is valid, but it belongs to the wrong species.
            if input_gene:
                # Search for this gene in the intended Target Species (Reviewed only)
                query = f"gene_exact:{input_gene} AND organism_name:\"{clean_target}\" AND reviewed:true"
                search_url = f"{base_url}/search?query={query}&fields=accession&size=1"
                
                search_resp = requests.get(search_url)
                if search_resp.status_code == 200:
                    results = search_resp.json().get('results', [])
                    if results:
                        best_match = results[0]['primaryAccession']
                        return best_match, f"The ID **{invalid_id}** belongs to *{org_sci_name}*. Did you mean the **{clean_target}** ID **{best_match}** (Gene: {input_gene})?"
            
            return None, f"The ID **{invalid_id}** belongs to *{org_sci_name}*, not {clean_target}. Please change your target species or provide the correspondent ID for {clean_target}, accordingly."

        # 4. Check if it is Unreviewed (TrEMBL) -> Suggest Reviewed (Swiss-Prot)
        # If we reached here, the species is correct, but the ID is missing from AFDB.
        # It is likely a TrEMBL entry. Let's find the Swiss-Prot version.
        if input_gene and input_tax_id:
             query = f"gene_exact:{input_gene} AND organism_id:{input_tax_id} AND reviewed:true"
             search_url = f"{base_url}/search?query={query}&fields=accession&size=1"
             
             s_resp = requests.get(search_url)
             if s_resp.status_code == 200:
                 results = s_resp.json().get('results', [])
                 if results:
                     reviewed_id = results[0]['primaryAccession']
                     # Only suggest if it is actually a different ID
                     if reviewed_id != invalid_id:
                         return reviewed_id, f"The ID **{invalid_id}** is an unreviewed (TrEMBL) entry. The reviewed (Swiss-Prot) entry for gene **{input_gene}** is **{reviewed_id}**."

        # 5. Default Fallback
        return None, "ID exists but was probably not reviewed in the reference proteome used for AFDB v4 predictions. All its comparisons will be skipped."

    except Exception:
        return None, "Could not connect to UniProt."


#rank proteins
@st.cache_data(show_spinner=False)
def prepare_ranking_data2(target_species, model_organisms, protein_list_str):
    """
    Performs all data parsing, ranking, and merging.
    Handled partial mismatches by warning user but continuing with valid data.
    """
    # 1. Setup
    uniprot_ids = protein_list_str.split()
    root_dir = 'data/decompressed/'
    mapping_file = 'data/core/model/extracted_and_renamed.csv'
    
    # Load mapping
    try:
        species_map = load_mapping(mapping_file)
    except Exception:
        return None, None, "Could not load species mapping file."
    
    # Check Target
    target_dir_id = species_map.get(target_species) # returns ta
    #st.info(f"species_map: {target_dir_id}")
    if not target_dir_id:
        return None, None, f"Target species '{target_species}' not found in mapping."
    target_dir = os.path.join(root_dir, target_dir_id)
    
    # from target species to target proteome
    proteome_map = load_mapping_proteome(mapping_file)
    proteome_target_id = proteome_map.get(target_species)
    index_file_path = f"data/subset/index/{proteome_target_id}_*_fc_id_mapping.csv"
    matching_proteome = glob.glob(index_file_path)
    if not matching_proteome:
        print(f"Error: No matching proteome file found for {proteome_target_id}")
    else:
        print(f"Found proteome file for {proteome_target_id}")
    index_file_path = matching_proteome[0]
    proteome_df = pd.read_csv(index_file_path) # columns: 'fc_db_index' and 'afdb_id'
    
    # handle protein mismatch in target proteome
    uid_missing = []
    final_analysis_ids = []
    for uid in uniprot_ids:
        pid="AF-"+uid+"-F1-model_v4"
        #st.info(f"UID check: {uid}, PID check: {pid}")
        if pid in proteome_df['afdb_id'].values:
            final_analysis_ids.append(uid)
        else:
            uid_missing.append(uid)
            #st.error(f"The protein ID **{uid}** was not found in the reference proteome of **{target_species}**. This can happen if the protein ID is not reviewed in the reference proteome used for AFDB v4 predictions. All its comparisons will be skipped.")
            suggested_id, reason_msg = get_uniprot_suggestion(uid, target_species)
            if suggested_id:
                st.warning(f"**Mismatch:** Protein **{uid}** was not found in the reference proteome of your target species: **{target_species}**.")
                new_pid = "AF-" + suggested_id + "-F1-model_v4"
                if new_pid in proteome_df['afdb_id'].values:
                    # SUCCESS: Replace it
                    final_analysis_ids.append(suggested_id)
                    st.info(f"Replacing **{uid}** with **{suggested_id}** for this analysis...")
                    with st.expander(f"See details for {uid} replacement"):
                        st.info(f"{reason_msg}")
                        st.page_link(f"https://www.uniprot.org/uniprotkb/{uid}/entry", label=f"Check {uid} existence in {target_species} proteome and potential updated reviewed ID on UniProt", icon="🌎")
                else:
                    st.error(f"**Mismatch:** Protein **{uid}** was not found in the reference proteome of your target species: **{target_species}**.")
                    # Suggestion exists in UniProt but NOT in local AFDB (e.g. >2700aa)
                    with st.expander(f"See {uid} error details"):
                        st.warning(f"{reason_msg}. UniProt suggests **{suggested_id}**, but it is also missing from the AFDB v4 (it might be too large or new).")
                        st.page_link(f"https://www.uniprot.org/uniprotkb/{uid}/entry", label=f"Check {uid} existence in {target_species} proteome and potential updated reviewed ID on UniProt", icon="🌎")
            else:
                # If no suggestion
                if reason_msg:
                    with st.expander(f"See {uid} error details"):
                        st.warning(f"{reason_msg}")
                        st.page_link(f"https://www.uniprot.org/uniprotkb/{uid}/entry", label=f"Check {uid} existence in {target_species} proteome and potential updated reviewed ID on UniProt", icon="🌎")

    # Update the main list to use ONLY the validated/replaced IDs
    uniprot_ids = final_analysis_ids
    # Stop if we have no proteins left to analyze
    if not uniprot_ids:
        return None, None, "No valid proteins found in the target proteome to analyze."
    st.success(f"Proceeding with the following valid protein IDs: **{', '.join(uniprot_ids)}**")
    # Storage
    ranking_scores = {}
    itol_data = {}
    tree_data = {} 
    all_best_hits_dfs = [] 
    # 2. Calculate Human (Target) Data
    human_score, human_vector = get_human_self_vector(target_dir, uniprot_ids)
    ranking_scores[target_species] = human_score
    tree_data[target_species] = human_vector
    itol_data[target_species] = {uid: 1 for uid in uniprot_ids}
    # 3. Calculate Comparison Species
    for species in model_organisms:
        if species not in species_map: 
            continue
        
        compare_dir_id = species_map[species]
        csv_file = find_comparison_csv(target_dir, compare_dir_id)
        search_col = 'proteinA'
        
        # Try reverse direction if needed
        if not csv_file:
            reverse_dir = os.path.join(root_dir, compare_dir_id)
            if os.path.isdir(reverse_dir):
                csv_file = find_comparison_csv(reverse_dir, target_dir_id)
                search_col = 'proteinB'

        if csv_file:
            score, presence, vector, best_hits_df = process_orthologs(csv_file, uniprot_ids, search_col)

            # --- FIX: Handle Reverse Lookup (Column Swap) Robustly ---
            if search_col == 'proteinB':
                try:
                    temp_col = best_hits_df['proteinA'].copy()
                    best_hits_df['proteinA'] = best_hits_df['proteinB']
                    best_hits_df['proteinB'] = temp_col
                except Exception:
                    continue
            best_hits_df['Compared_Species'] = species 
            all_best_hits_dfs.append(best_hits_df)
            ranking_scores[species] = score
            itol_data[species] = presence
            tree_data[species] = vector
        else:
            ranking_scores[species] = 0.0
            itol_data[species] = {uid: -1 for uid in uniprot_ids}
            tree_data[species] = [0.0] * len(uniprot_ids)
            print(f"No comparison file found for {species}")

    # 4. Merge Results
    if all_best_hits_dfs:
        comparison_df = pd.concat(all_best_hits_dfs, ignore_index=True)
    else: 
        comparison_df = pd.DataFrame()
        st.error(f"No structural ortholog found in your model organisms.")
        for uid in uniprot_ids:
            st.page_link(f"https://www.uniprot.org/uniprotkb/{uid}/entry", label=f"Check {uid} existence in {target_species} proteome and potential updated reviewed ID on UniProt", icon="🌎")

    # Build Rank DataFrame
    df_rank = pd.DataFrame.from_dict(ranking_scores, orient='index', columns=['Score'])
    df_target = df_rank.loc[[target_species]]
    df_other = df_rank.drop(target_species).sort_values(by='Score', ascending=False)
    df_rank_sorted = pd.concat([df_target, df_other])
    df_rank = df_rank_sorted.reset_index()
    df_rank.columns = ['Organism', 'Score']

    # Build Presence DataFrame
    presence_list = []
    for sp, pres in itol_data.items():
        row_data = {'Organism': sp}
        for puid in uniprot_ids:
            row_data[puid] = 1 if pres.get(puid, 0) == 1 else 0
            if pres.get(puid, 0) != 1:
                if puid not in uid_missing:
                    st.warning(f"No structural match for {target_species} {puid} found in {sp}. Skipping.")
        presence_list.append(row_data)
    df_presence = pd.DataFrame(presence_list)

    # Final Merge
    df_rank_filtered = df_rank[df_rank['Organism'] != target_species].copy()
    df_rank_filtered.reset_index(drop=True, inplace=True)
    df_rank_filtered.insert(0, 'Rank', df_rank_filtered.index + 1)
    df_final = pd.merge(df_rank_filtered, df_presence, on='Organism', how='left')

    return df_final, comparison_df, None #(error_msg is None)


@st.cache_data
def get_base64_image(image_path):
    """Converts a local file to a base64 string for HTML embedding."""
    if not os.path.exists(image_path):
        return "" # Return empty if file not found to prevent errors
    with open(image_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode()
#---------------------------------------
st.markdown("""
    <style>
    /* Hide the link button (anchor) next to headers */
    [data-testid="stHeaderActionElements"] {
        display: none !important;
    }
    
    /* Optional: If the title still looks shifted, remove the padding usually reserved for the icon */
    h1, h2, h3, h4, h5, h6 {
        padding-top: 0rem !important;
    }
    </style>
""", unsafe_allow_html=True)
##MAIN##
st.set_page_config(
    page_title="SHE: Model",
    page_icon= "data/tree/logo_tree.png", 
    layout="wide"
)
#sidebar navigation
with st.sidebar:
    st.title("Browse")
    for item in NAV_ITEMS:
        st.page_link(item["path"], label=item["label"], icon=item["icon"])
st.sidebar.divider()

#st.title('Model Organism Ranking')
st.markdown("""
    <h2 style='
        text-align: center; 
        font-size: 3rem; 
        font-weight: bold; 
        /* Gradient Background */
        background: linear-gradient(to right, #00FF7F, #40E0D0); 
        /* Clip background to text */
        -webkit-background-clip: text; 
        -webkit-text-fill-color: transparent; 
        background-clip: text; 
        color: transparent;'>
        Model Organism Ranking
    </h2>
""", unsafe_allow_html=True)
#set_background('data/logos/model.png')
image_path = 'data/logos/model.png'
# 2. Check if file exists and encode
if os.path.exists(image_path):
    img_b64 = get_base64_image(image_path)

#st.write("Find the best model organism match for your biological pathway by combined structural matching for all pathway proteins.")

st.markdown(
    f"""
    <style>
    .banner-container {{
        /* Use the Base64 image as a background */
        background-image: url("data:image/png;base64,{img_b64}");
        background-size: cover; /* Ensures image fills the box */
        background-position: center;
        border-radius: 15px;
        padding: 100px 20px; /* Top/Bottom padding creates the height */
        text-align: center;
        position: relative;
        margin-bottom: 20px;
    }}
    
    /* Adds a dark overlay so text pops */
    .banner-overlay {{
        background-color: rgba(0, 0, 0, 0.4); /* 40% opaque black */
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        border-radius: 15px;
        z-index: 1;
    }}

    /* The text styling */
    .banner-text {{
        position: relative;
        z-index: 2; /* Ensures text is above the dark overlay */
        color: white;
        transition: all 0.3s ease-in-out;
        text-shadow: 2px 2px 5px rgba(0,0,0,0.8);
    }}

    /* The Glow Effect on Hover */
    .banner-text:hover {{
        text-shadow: 0 0 10px #fff, 0 0 20px #fff, 0 0 30px #00e6b8ff;
        transform: scale(1.02); /* Slight zoom for effect */
    }}
    </style>

    <div class="banner-container">
        <div class="banner-overlay"></div>
        <div class="banner-text">
            <h2>Find the best model organism to study your biological pathway with</h2>
        </div>
    </div>
    """,
    unsafe_allow_html=True
)



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
#st.title(':green[SHE: Structural History of Eukarya]')
#st.image('data/tree/model.png',ccaption='Protein-protein structure comparison to rank the best model organism for your target organism')
with st.expander(":violet[See Example Usage]"):
    st.subheader(":violet[Follow the following steps to discover the top model organism to study part of the human Fanconi anemia (FA) DNA repair pathway]")
    st.subheader(':green[Step 1. Choose Target Organim]')
    st.selectbox('Target organism:', 'Homo sapiens', key = 'ex1')
    st.write('Target organism: Homo sapiens')
    st.subheader(':green[Step 2. Select Model Organims]')
    st.selectbox('Model organisms:', 'Mus musculus, Saccharomyces cerevisiae, Arabidopsis thaliana, Xenopus tropicalis, Danio rerio, Caenorhabditis elegans', key = 'ex2')
    st.write('Model organisms: Mus musculus, Saccharomyces cerevisiae, Arabidopsis thaliana, Xenopus tropicalis, Danio rerio, Caenorhabditis elegans')
    st.subheader(':green[Step 3. Choose target proteins from the proteome of the target species selected above]')
    st.selectbox('Target proteins:', 'Q9NXK8 Q9Y2M0 Q9NVI1', key = 'ex3')
    st.write('Target proteins (FA pathway): Q9NXK8, Q9Y2M0, Q9NVI1')
    st.subheader(':green[Step 4. Check Model Organism ranking based on protein structure similarity to the target proteins]')
    #st.markdown(f":orange[**The top model organism for your study is: Mus musculus**]")
    st.markdown(
            f"""
            <div style='
                text-align: center; 
                color: orange; 
                font-size: 2rem;           /* Much bigger size */
                font-weight: 800;          /* Bold weight */
                margin-bottom: 20px;       /* Space below */
                text-shadow: 2px 2px 0 #000000; /* Black shadow for readability */
            '>
                The top model organism for your study is: Mus musculus mouse
            </div>
            """,
            unsafe_allow_html=True
        )
    st.image("data/subset/figure_model/ranking2.png", caption="Organisms Ranking based on bitscore (the higher, the better)", use_container_width=True)

    row2 = st.columns(3)
    emojis = [":mouse2:", ":tropical_fish:", ":seedling:"]
    titles = ["House Mouse", "Zebrafish", "Arabidopsis thaliana"]
    ranks = ["1st", "2nd", "3rd"]
    bitscores = ["14527", "3894", "1463"]
    score_1= ["2627", "365", "156"]
    score_2= ["4330", "3025", "1239"]
    score_3= ["7570", "504", "68"]
    proteins1= ["Q9EPX5", "Q6PBJ7", "Q9FM90"]
    proteins2= ["Q69ZT1", "Q1LWH4", "Q5XVJ4"]
    proteins3= ["Q8K368", "B3DJT4", "Q9LY34"]
    # Define the custom heights for the *visible* title containers
    tile_heights = [170, 170, 170] 
    for col, emoji, title, bitscore, height, rank, score1, score2, score3, protein1, protein2, protein3 in zip(row2, emojis, titles, bitscores, tile_heights, ranks, score_1, score_2, score_3, proteins1, proteins2, proteins3):
        with col:
            #with st.container(key="my_blue_container")
            # 1. Create the Visible Title Container with Custom Height
            tile = st.container(height=height)
            # Use a large markdown header for the title inside the container
            tile.markdown(f"## {rank} {emoji} {title}")
            
            with st.expander(":green[Step 5. Check structural alignments]"):
                # Content that appears when the user clicks/expands
                if title == "House Mouse":
                    st.write(f"**Total Bitscore** human vs {title}: **{bitscore}**")
                    #st.markdown("<br>", unsafe_allow_html=True)
                    st.write("---")
                    st.write(f"**Bitscore** human Q9NXK8 (cyan) vs {title} {protein1} (yellow): **{score1}**")
                    alignment = open("data/subset/figure_model/mouse_Q9NXK8.mp4", "rb")
                    bytes = alignment.read()
                    st.video(bytes, autoplay=True, loop=True, muted=True)
                    st.write("---")
                    st.write(f"**Bitscore** human Q9Y2M0 (cyan) vs {title} {protein2} (yellow): **{score2}**")
                    alignment2 = open("data/subset/figure_model/mouse_Q9Y2M0.mp4", "rb")
                    bytes2 = alignment2.read()
                    st.video(bytes2, autoplay=True, loop=True, muted=True)
                    st.write("---")
                    st.write(f"**Bitscore** human Q9NVI1 (cyan) vs {title} {protein3} (yellow): **{score3}**")
                    alignment3 = open("data/subset/figure_model/mouse_Q9NVI1.mp4", "rb")
                    bytes3 = alignment3.read()
                    st.video(bytes3, autoplay=True, loop=True, muted=True)
                elif title == "Zebrafish":
                    st.write(f"**Total Bitscore** human vs {title}: **{bitscore}**")
                    st.write("---")
                    st.write(f"**Bitscore** human Q9NXK8 (cyan) vs {title} {protein1} (yellow): **{score1}**")
                    #st.markdown("<br>", unsafe_allow_html=True)
                    alignment = open("data/subset/figure_model/fish_Q9NXK8.mp4", "rb")
                    bytes = alignment.read()
                    st.video(bytes, autoplay=True, loop=True, muted=True)
                    #st.markdown("<br>", unsafe_allow_html=True)
                    st.write("---")
                    st.write(f"**Bitscore** human Q9Y2M0 (cyan) vs {title} {protein2} (yellow): **{score2}**")
                    #st.markdown("<br>", unsafe_allow_html=True)
                    alignment2 = open("data/subset/figure_model/fish_Q9Y2M0.mp4", "rb")
                    bytes2 = alignment2.read()
                    st.video(bytes2, autoplay=True, loop=True, muted=True)
                    st.write("---")
                    st.write(f"**Bitscore** human Q9NVI1 (cyan) vs {title} {protein3} (yellow): **{score3}**")
                    #st.markdown("<br>", unsafe_allow_html=True)
                    alignment3 = open("data/subset/figure_model/fish_Q9NVI1.mp4", "rb")
                    bytes3 = alignment3.read()
                    st.video(bytes3, autoplay=True, loop=True, muted=True)
                elif title == "Arabidopsis thaliana":
                    st.write(f"**Total Bitscore** human vs {title}: **{bitscore}**")
                    st.write("---")
                    st.write(f"**Bitscore** human Q9NXK8 (cyan) vs {title} {protein1} (yellow): **{score1}**")
                    alignment = open("data/subset/figure_model/plant_Q9NXK8.mp4", "rb")
                    bytes = alignment.read()
                    st.video(bytes, autoplay=True, loop=True, muted=True)
                    st.write("---")
                    st.write(f"**Bitscore** human Q9Y2M0 (cyan) vs {title} {protein2} (yellow): **{score2}**")
                    alignment2 = open("data/subset/figure_model/plant_Q9Y2M0.mp4", "rb")
                    bytes2 = alignment2.read()
                    st.video(bytes2, autoplay=True, loop=True, muted=True)
                    st.write("---")
                    st.write(f"**Bitscore** human Q9NVI1 (cyan) vs {title} {protein3} (yellow): **{score3}**")
                    alignment3 = open("data/subset/figure_model/plant_Q9NVI1.mp4", "rb")
                    bytes3 = alignment3.read()
                    st.video(bytes3, autoplay=True, loop=True, muted=True)

#st.subheader("Start your analysis by clicking on 'Rank Model Organisms' !")
#st.warning("Please note that this analysis might take longer the more structural matching have to performed for all pathway proteins. It is possible to register your email in the sidebar, so that your results will be sent to you as soon as they are ready!")
st.warning("Please note: Analysis time scales with the number of structural matches required. To receive your results automatically upon completion, you may register your email address in the sidebar.")

###------sidebar sign up--------------------
#user sign up for updates and analysis results
st.sidebar.subheader('Model Organism Ranking', divider="violet")
with st.sidebar.expander("Help: How to use the Model Organism Ranking tool"):
    st.write("Start it by clicking on 'Rank Model Organisms'!")
    st.markdown("""
    * Your target species (i.e. Homo sapiens)
    * The proteins in the pathway you are studying (note that their UniProt ID has to be found in the proteome of the target species)
    * The available model organisms you could choose from
    """)
    st.write("This will result in a ranking of the model organisms  based on the structural similarity (bitscore, the higher the better) of the input proteins compared to the target species")
    st.write("If you wish to visualise the structural alignments of your input proteins, turn on 'Show protein structure alignments'.")

st.sidebar.subheader('Receive Results', divider="rainbow")
email_input = st.sidebar.text_input('Your email', key="email_input") 
if st.sidebar.button('Receive results'):
    user_email = st.session_state.email_input.strip()
    if not user_email:
        st.error("Please enter an email address.")
    elif not is_valid_email(user_email):
        st.sidebar.error("Please enter a valid email address (e.g., user@example.com)")
    else:
        if save_email(user_email):
            st.balloons()
            st.sidebar.success(f"Thank you! **{user_email}** is where you will receive updates regarding SHE!")


################analysis -------------------------------------------
proteomes_path = 'data/tree/selected_proteomes_above_25_formatted_organisms.csv'
new_proteomes_path =  'data/core/merged_evo_rate.csv'
proteomes_df = pd.read_csv(proteomes_path)

st.markdown("""
<style>
/* Target the button inside the Streamlit button container */
div.stButton > button {
    width: 100%;             /* Fill the column width */
    height: 60px;            /* Make it taller/bigger */
    font-size: 20px;         /* Larger text */
    font-weight: bold;
    border-radius: 12px;     /* Rounded corners */
    border: 2px solid #4ade80; /* Thin green border default */
    background-color: #1f2937; /* Dark background (matches your dark theme) */
    color: white;
    transition: all 0.3s ease-in-out; /* Smooth animation */
}

/* HOVER STATE: The "Green Aura" */
div.stButton > button:hover {
    background-color: #1f2937; 
    border-color: #22c55e;   /* Brighter green border */
    color: #4ade80;          /* Green text */
    
    /* The Aura Effect: Box Shadow */
    box-shadow: 0 0 20px rgba(74, 222, 128, 0.7), 
                0 0 10px rgba(74, 222, 128, 0.5);
                
    transform: scale(1.02);  /* Slightly grow the button */
}

/* ACTIVE STATE: When clicked */
div.stButton > button:active {
    box-shadow: none;
    transform: scale(0.98);
}
</style>
""", unsafe_allow_html=True)

# 1. Initialize a "state" variable to remember if the analysis is active
if 'ranking_active' not in st.session_state:
    st.session_state['ranking_active'] = False

# Initialize state for the inner "Start Analysis" button
if 'analysis_submitted' not in st.session_state:
    st.session_state['analysis_submitted'] = False

left, middle, right = st.columns([0.5, 2, 0.5])
with middle:
    # This button will pick up the CSS styles defined above
    if st.button('Rank Model Organisms', width='stretch'):
        st.session_state['ranking_active'] = True
# 3. Check the STATE variable, not the button itself
if st.session_state['ranking_active']:
#if middle.button('Rank Model Organisms'):
    user_email = st.session_state.email_input.strip()
    if not user_email:
        random_num = random.randint(100000, 999999)
        user = f"user_{random_num}"
    else: 
        user=user_email
    st.title('Find the best model organism for your study according to the structural phylogenetic tree!')

    # choose target -------------------------------------------------------------------------------------
    #target: default human
    options = proteomes_df['Organism_formatted'].tolist()
    target_search_term = 'homo_sapiens_human'
    try:
        # Find the index regardless of how it's capitalized in the CSV
        default_index = [x.lower() for x in options].index(target_search_term)
    except ValueError:
        default_index = 0
    # 1. .replace('_', ' ') -> "Homo sapiens Human"
    # 2. .capitalize()      -> "Homo sapiens human" (First letter Big, rest small)
    target = st.selectbox(
        'Which organism would you like to study (i.e: Homo sapiens human)', 
        options, 
        index=default_index, 
        key="target", 
        format_func=lambda x: x.replace('_', ' ').capitalize()
    )

    # --- MODEL ORGANISM SELECTION ---
    filtered_model_options = [opt for opt in options if opt != target]
    # 1. Define Callback Functions
    def add_model_organism():
        st.session_state.model_count += 1 

    def remove_model_organism():
        # Prevent removing the last remaining box
        if st.session_state.model_count > 1:
            st.session_state.model_count -= 1

    # 2. Initialize Session State
    if 'model_count' not in st.session_state:
        st.session_state.model_count = 1 
        
    model_organisms = []

    # 3. Render Selectboxes Loop
    for i in range(1, st.session_state.model_count + 1):
        selection = st.selectbox(
            f'Which model organism #{i} would you like to add to the comparison?',
            #proteomes_df['Organism_formatted'],
            filtered_model_options,
            index=None,  # box starts empty
            placeholder="Select an organism...", # Text shown when empty
            format_func=lambda x: x.replace('_', ' ').capitalize(),
            key=f"model_{i}"
        )
        
        # Only add to the list if the user actually selected something
        if selection:
            model_organisms.append(selection)

    # 4. Add/Remove Buttons (Side by Side)
    #col_add, col_remove = st.columns([1.5, 1.5])
    col_l, col_add, col_rem, col_r = st.columns([1, 1, 1, 1])
    with col_add:
        st.button(
            "➕ Add Another Model", 
            on_click=add_model_organism,
            width='stretch'
        )
    with col_rem:
        # Only show remove button if there is more than 1 box
        if st.session_state.model_count > 1:
            st.button(
                "➖ Remove Last Model", 
                on_click=remove_model_organism,
                width='stretch'
            )

    #choose proteins -----------------------------------------------------------------------
    protein_list = st.text_input(
    'Which proteins are you studying in your target organism? (each UniProt ID separated by a space, i.e: Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616). Please note that the protein IDs have to be found in the reference proteome of the target organism selected above.', key="protein_list")
    st.write("---")
    # Ratios: [Spacer, Button, Spacer]
    c_left, c_run, c_right = st.columns([2, 1.5, 2])
    # We use a callback or direct check to set the submission state
    if c_run.button("Start Ranking Analysis", type="primary", width='stretch'):
        st.session_state['analysis_submitted'] = True
    if st.session_state['analysis_submitted']:
        #start ranking-------------------------------------------------------------------------
        if (protein_list == '') or not model_organisms:
            st.warning('Please input the UniProt ID of the protein(s) in your pathway of interest, the target organism to study, and the model organisms you can choose from.')
        else:
            protein_list = protein_list.replace(',', '').upper()
            st.write('Selected Protein(s) in your study: ', protein_list)
            
            # Handle Target Name
            if target == '':
                target_species = 'Homo_sapiens_Human'
            else:
                target_species = target

            status_placeholder = st.empty()
            status_placeholder.info(f"a. Ranking your model organisms...")

            # --- CALL CACHED FUNCTION ---
            # This runs immediately if cached: #st.chache_data()
            df_final, comparison_df, error_msg = prepare_ranking_data2(target_species, model_organisms, protein_list)

            if error_msg:
                print(error_msg)
                #st.session_state['ranking_active'] = False
                st.stop() # Stops execution safely

            # --- SAVE FILES ---
            # Side effects (saving files) should happen outside the cached function
            csv_out = "results/model_trees/" + user + "/"
            os.makedirs(csv_out, exist_ok=True)
            
            txt_path = csv_out + "model_ranking.txt"
            df_final.to_csv(txt_path, index=False, sep='\t')
            
            ranking_df = comparison_df.copy()
            if not ranking_df.empty:
                # Only rename if we actually have columns
                ranking_df.columns = [f'{target}_Protein', 'Structural_ortholog', 'Bit_score', 'Model_Species']
                txt_path2 = csv_out + "protein_hits_ranking.txt"
                ranking_df.to_csv(txt_path2, index=False, sep='\t')
            else:
                # If empty, skipping saving the hits file and warn the user
                st.warning("No structural orthologs were found.")
                # Depending on your preference, you might want to stop here 
                # because the plots below will likely be empty too.
                if df_final['Score'].max() == 0:
                    st.stop()

            # --- DISPLAY LOGIC ---
            st.write("--- Organism Ranking ---")
            st.write('Selected Target Organism: ', target_species)
            st.write('Selected Model Organism(s): ', ', '.join(model_organisms))
            
            top_model_organism = df_final.iloc[0]['Organism']
            # Extract just Genus species (e.g., "Mus musculus")
            top_model_organism = ' '.join(top_model_organism.split('_', maxsplit=2)[:2])
            
            st.markdown(
                f"""
                <div style='
                    text-align: center; 
                    color: orange; 
                    font-size: 2rem;           /* Much bigger size */
                    font-weight: 800;          /* Bold weight */
                    margin-bottom: 20px;       /* Space below */
                    text-shadow: 2px 2px 0 #000000; /* Black shadow for readability */
                '>
                    The top model organism for your study is: {top_model_organism}
                </div>
                """,
                unsafe_allow_html=True
            )
            status_placeholder.empty()
            
            # 1. Prepare Data for Plotting
            df_final.sort_values('Rank', ascending=True, inplace=True)
            df_final['Display_Organism'] = df_final['Organism'].str.split('_', n=2).str[:2].str.join('_')
            df_final['Display_Organism'] = df_final['Display_Organism'].str.replace('_', ' ')
            
            comparison_df['Model_Protein'] = comparison_df['proteinB'].str.split('-').str[1]
            comparison_df['Model_Organism'] = comparison_df['Compared_Species'].str.split('_', n=2).str[:2].str.join('_')
            comparison_df['Model_Organism'] = comparison_df['Model_Organism'].str.replace('_', ' ')
            
            # Create the color grouping based on proteinA
            comparison_df['Color_Group'] = comparison_df['proteinA'].str.split('-').str[1]
            
            all_species = comparison_df['Model_Organism'].unique().tolist()
            target_name = ' '.join(target_species.split('_', maxsplit=2)[:2])

            # 2. Radial Bar Chart (Organism Ranking)
            ordered_species = df_final['Display_Organism'].tolist()
            fixed_cols = ['Rank', 'Organism', 'Score', 'Display_Organism']
            uniprot_id_cols = [col for col in df_final.columns if col not in fixed_cols]
            
            # Initialize hover_data dictionary
            hover_data_config = {'Score': ':.1f', 'Rank': True, 'Display_Organism': False}
            for col in uniprot_id_cols:
                hover_data_config[col] = True 

            fig = px.line_polar(
                df_final,
                r='Score',                 
                theta='Display_Organism',  
                color_discrete_sequence=['#00CC96'],
                line_close=True,
                title='Model Organism Ranking by Bit score',
                hover_data=hover_data_config,
                start_angle=0,
                direction="clockwise", 
                template='plotly_dark'
            )

            fig.update_traces(
                mode='lines+markers', 
                fill='toself',        
                opacity=0.7,      
                marker=dict(size=8)
            )
            
            fig.update_layout(
                polar=dict(
                    radialaxis=dict(visible=True, showticklabels=True, showline=False, angle=90, tickangle=90),
                    angularaxis=dict(categoryorder='array', categoryarray=ordered_species)
                ),
                coloraxis_showscale=False,
            )
            st.plotly_chart(fig, use_container_width=True)

            # 3. Radar Chart (Protein Hits)
            radar_fig = px.line_polar(
                comparison_df,
                r='bit_score',         
                theta='Model_Organism', 
                color='Color_Group',    
                line_close=True,
                title=f'Protein Hit Bitscores by Species compared to {target_name}',
                labels={'bit_score': 'Bit Score', 'Model_Organism': 'Model Organism', 'Color_Group': 'Target Protein (UniProt ID)'},
                hover_data={'Model_Protein': True, 'bit_score': ':.0f', 'Compared_Species': False},
                start_angle=0,
                direction="clockwise", 
                template='plotly_dark'
            )

            radar_fig.update_traces(
                mode='lines+markers', 
                fill='toself',       
                opacity=0.5,      
                marker=dict(size=5)  
            )

            radar_fig.update_polars(
                radialaxis_showgrid=True,
                angularaxis_showgrid=True,
                bgcolor='#000000',
                angularaxis_categoryorder='array', 
                angularaxis_categoryarray=ordered_species,
                radialaxis=dict(
                    visible=True, showticklabels=True, showline=False, angle=90, tickangle=90,
                    range=[0, comparison_df['bit_score'].max() * 1.1] 
                )
            )

            # Extract color map for later use
            color_map = {}
            for trace in radar_fig.data:
                color = trace.line.color 
                if trace.name not in color_map and color:
                    color_map[trace.name] = color
            # st.plotly_chart(radar_fig, width='stretch')

            # --- STRUCTURE ALIGNMENT ---
            if st.toggle(':green[Show protein structure alignments]'):
                st.subheader(f'Structural comparison between {target_name} and model organisms')
                
                unique_target_proteins = comparison_df['proteinA'].unique()
                selected_protein_A = st.radio(
                    "Select the **Target protein** to visualise:",
                    unique_target_proteins,
                    horizontal=True
                )
                
                selected_protein_df = comparison_df[comparison_df['proteinA'] == selected_protein_A]
                
                if selected_protein_df.empty:
                    st.warning(f"No structural comparison hit found for {selected_protein_A}.")
                else:
                    selected_group = selected_protein_df['Color_Group'].iloc[0]
                    selected_protein_df=selected_protein_df.sort_values(by='bit_score', ascending=False)
                    
                    # --- Single Protein Radar Chart ---
                    radar_fig_alone = px.line_polar(
                        selected_protein_df,  
                        r='bit_score',
                        theta='Model_Organism',
                        color='Color_Group',
                        color_discrete_map={selected_group: color_map.get(selected_group, '#FFFFFF')},
                        line_close=True,
                        title=f'Structural hits for {selected_protein_A} from {target_name}',
                        labels={'bit_score': 'Bit Score', 'Model_Organism': 'Model Organism', 'Color_Group': 'Target Protein'},
                        hover_data={'Model_Protein': True, 'bit_score': ':.0f', 'Compared_Species': False},
                        start_angle=0,
                        direction="clockwise",
                        template='plotly_dark'
                    )

                    radar_fig_alone.update_traces(mode='lines+markers', fill='toself', opacity=0.8, marker=dict(size=5))
                    radar_fig_alone.update_polars(
                        radialaxis_showgrid=True, angularaxis_showgrid=True, bgcolor='#000000',
                        angularaxis_categoryorder='array', angularaxis_categoryarray=ordered_species,
                        radialaxis=dict(visible=True, showticklabels=True, showline=False, angle=90, tickangle=90,
                                        range=[0, selected_protein_df['bit_score'].max() * 1.1])
                    )
                    st.plotly_chart(radar_fig_alone, use_container_width=True)
                    
                    # --- STRUCTURE LOOP ---
                    # Pre-load Mapping/PDB info for the target ONCE
                    target_uniprot_id = selected_protein_A.split('-')[1]
                    mapping_file = 'data/core/model/extracted_and_renamed.csv'
                    proteome_map = load_mapping_proteome(mapping_file)
                    proteome_target_id = proteome_map.get(target_species)
                    
                    pdb_target = get_pdb2(target_uniprot_id, proteome_target_id)
                    
                    if not pdb_target:
                        st.error(f"Could not retrieve PDB file for target {target_uniprot_id}")
                    else:
                        for index, row in selected_protein_df.iterrows():
                            hit_proteinB = row['proteinB']
                            hit_species = row['Compared_Species']
                            target_species_clean = target_species.replace('_', ' ')
                            hit_species_clean = hit_species.replace('_', ' ')
                            st.divider()
                            st.write(f"### **{hit_species_clean}** (Bit Score: {row['bit_score']:.0f})")

                            # Get Model PDB
                            model_uniprot_id = hit_proteinB.split('-')[1]
                            proteome_model_id = proteome_map.get(hit_species)
                            pdb_model = get_pdb2(model_uniprot_id, proteome_model_id)

                            if not pdb_model:
                                st.warning(f"PDB file not found for {hit_species}.")
                                continue

                            # --- ALIGNMENT & SVD ---
                            st.write(f"Structural Overlap {target_species_clean} {target_uniprot_id} vs. {hit_species_clean} {model_uniprot_id}")
                            aln_data = run_foldseek_alignment_only(query_path=pdb_model, target_path=pdb_target)

                            #if aln_data:
                            if aln_data is not None and not aln_data.empty:
                                out_file_name = f"foldseek_{target_uniprot_id}_vs_{model_uniprot_id}.pdb"
                                out_file = os.path.join(csv_out, out_file_name)
                                
                                target_fmt = get_format_from_path(pdb_target)
                                model_fmt = get_format_from_path(pdb_model)
                                
                                success, rmsd = calculate_svd_and_merge(
                                    target_file=pdb_target, target_fmt=target_fmt, 
                                    query_file=pdb_model, query_fmt=model_fmt, 
                                    aln_data=aln_data, 
                                    out_name=out_file
                                )

                                if success:
                                    # --- USE THE FRAGMENT ---
                                    render_interactive_structure(
                                        key_suffix=index, # Unique key per row ensures widgets don't conflict
                                        out_file=out_file, 
                                        target_name=selected_protein_A,
                                        target_species=target_species_clean,
                                        query_name=hit_proteinB, 
                                        model_species=hit_species_clean,
                                        tm_score=aln_data.get('tm', 0.0), 
                                        rmsd=rmsd
                                    )
                                else:
                                    st.warning("SVD superposition failed.")
                            else:
                                st.warning("Foldseek could not align the structures. Maybe too many disordered regions?")
                                st.warning("Try it directly on Foldseek search page:")
                                #st.page_link("https://search.foldseek.com/search", label="Foldseek search", icon="🌎")
                                fold_path = "data/logos/foldseek_logo.png" # Replace with your png path
                                fold_b64 = get_base64_image(fold_path)
                                st.markdown(f"""
                                    <a href="https://search.foldseek.com/search" target="_blank" style="
                                        text-decoration: none; 
                                        color: inherit; 
                                        display: inline-flex; 
                                        align-items: center; 
                                        gap: 10px; 
                                        background-color: #f0f2f6;
                                        padding: 10px 15px; 
                                        border-radius: 8px; 
                                        transition: 0.3s;
                                        border: 1px solid rgba(0,0,0,0.1);
                                    "
                                    onmouseover="this.style.backgroundColor='#e0e2e6';"
                                    onmouseout="this.style.backgroundColor='#f0f2f6';"
                                    >
                                        <img src="data:image/png;base64,{fold_b64}" style="
                                            width: 20px; 
                                            height: 20px; 
                                            object-fit: contain;
                                        ">
                                        
                                        <span style="font-weight: 600;">Foldseek search</span>
                                    </a>
                                """, unsafe_allow_html=True)

            # --- DOWNLOAD BUTTON LOGIC ---
            is_random_user = False
            if not user_email:
                is_random_user = True
            user_dir = os.path.join("results", "model_trees", user)
            target_folder_path = user_dir
            user_name = user.split('@')[0]
            zip_file_name = f"results_{user_name}" 
            output_base_dir = "temp_archives" 
            os.makedirs(output_base_dir, exist_ok=True) 
            shutil.make_archive(
                base_name=os.path.join(output_base_dir, zip_file_name),
                format='zip',
                root_dir=os.path.dirname(target_folder_path), 
                base_dir=os.path.basename(target_folder_path)
            )
            zip_file_path = os.path.join(output_base_dir, f"{zip_file_name}.zip")
            with open(zip_file_path, "rb") as file:
                st.sidebar.download_button(
                    label="Download results now",
                    data=file, 
                    file_name=f"{zip_file_name}.zip",
                    mime="application/zip",
                    icon=":material/download:",
                )
        
    st.write("---")
    if st.button("Close Analysis"):
        st.session_state['ranking_active'] = False
        st.session_state['analysis_submitted'] = False # Reset internal state too
        if is_random_user and 'user_dir' in locals():
            try:
                if os.path.exists(zip_file_path): os.remove(zip_file_path)
                if os.path.exists(user_dir): shutil.rmtree(user_dir)
            except Exception:
                pass
        st.rerun()
    if st.session_state['ranking_active'] == False:
        # Cleanup if module is closed externally
        if 'is_random_user' in locals() and is_random_user:
            try:
                if os.path.exists(zip_file_path): os.remove(zip_file_path)
                if os.path.exists(user_dir): shutil.rmtree(user_dir)
            except Exception:
                pass

# --- CITATION
st.subheader('Cite', divider="rainbow")
st.write('Preprint:  Qiuzhen Li, Diandra Daumiller, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')
