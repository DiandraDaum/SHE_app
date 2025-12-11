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

#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
# save user email
EMAIL_FILE = 'data/email_results.csv'
def save_email(email):
    """Saves the email and current timestamp to the CSV file. 
    It ensures the directory exists before attempting to write."""
    try:
        # Ensure the directory exists (necessary for 'data/emails_results.csv')
        os.makedirs(os.path.dirname(EMAIL_FILE), exist_ok=True)
        
        # Check if the file exists to determine if we need to write headers
        file_exists = os.path.exists(EMAIL_FILE)
        
        # Use 'a' (append mode) and 'newline=""' to prevent blank rows
        with open(EMAIL_FILE, 'a', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Write header only if the file did not exist before opening
            if not file_exists:
                writer.writerow(['email', 'timestamp'])
            
            # Write the new email entry
            timestamp = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
            writer.writerow([email, timestamp])
        return True
    except Exception as e:
        st.error(f"An error occurred while saving: {e}")
        return False

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
        st.error(f"Error reading mapping file: {e}")
        sys.exit(1)

def load_mapping_proteome(mapping_file):
    try:
        mapping_df = pd.read_csv(mapping_file)
        return dict(zip(mapping_df['species_name'], mapping_df['proteome_id']))
    except Exception as e:
        st.error(f"Error reading mapping file: {e}")
        sys.exit(1)

def find_comparison_csv(target_dir, compare_dir_id):
    """Finds the Human_vs_Species or Species_vs_Human CSV file."""
    search_pattern = os.path.join(target_dir, '**', f'*vs{compare_dir_id}_top_hits.csv')
    csv_files = glob.glob(search_pattern, recursive=True)
    
    if not csv_files:
        search_pattern_alt = os.path.join(target_dir, '**', f'*vs{compare_dir_id}*.csv')
        csv_files = glob.glob(search_pattern_alt, recursive=True)

    return csv_files[0] if csv_files else None

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
        st.warning(f"Error processing {csv_file}: {e}")
        
    # Convert map to list (Vector) ensures order is identical for all species
    score_vector = [score_map[uid] for uid in uniprot_ids]
    best_hits_df = pd.DataFrame(best_hit_data)
    
    return total_score, presence_map, score_vector, best_hits_df

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
        st.warning(f"CSV file {csv_file} is empty. Skipping.")
    except Exception as e:
        st.error(f"Error processing CSV {csv_file}: {e}")
        
    return total_score

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
                st.warning(f"Could not fully read {csv_file} for self-score: {e}")
                continue
        
        if max_score_for_this_protein == 0.0:
            st.warning(f"  > No score found for {uniprot_id} in any comparison file.")
        
        #st.info(f"  > Max score (baseline) for {uniprot_id} is: {max_score_for_this_protein}")
        total_self_score += max_score_for_this_protein

    return total_self_score

def generate_distance_matrix(species_scores_dict, uniprot_ids, out_dir):
    """
    Generates dmat.npy and ids_in_dmat.csv for the bootstrap script.
    """
    #st.info("Generating distance matrix ...")
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(species_scores_dict, orient='index', columns=uniprot_ids)
    
    # Normalize: log2(x + 1) to handle the wide range of bit scores
    df_log = np.log2(df + 1)
    
    # Calculate Distance (Euclidean distance between score profiles)
    dists = pdist(df_log.values, metric='euclidean')
    dmat = squareform(dists)
    
    # Save dmat.npy
    dmat_path = os.path.join(out_dir, 'dmat.npy')
    np.save(dmat_path, dmat)
    # FIX 1: Changed 'outdir' to 'out_dir' and ensured correct path construction
    dmat_out = dmat_path
    
    # Save ids_in_dmat.csv
    ids_df = pd.DataFrame({'id': df.index})
    ids_path = os.path.join(out_dir, 'ids_in_dmat.csv')
    ids_df.to_csv(ids_path, index=False)
    # FIX 2: Changed 'outdir' to 'out_dir' and ensured correct path construction
    ids_out = ids_path

    return dmat_out, ids_out


# make tree from distance matrix -------------------------------------------
def bootstrap_from_distance_matrix(distance_matrix):
    """Bootstrap by resampling cells in the distance matrix"""
    n = len(distance_matrix)
    bootstrap_matrix = copy.deepcopy(distance_matrix)

    col_pos = np.random.choice(range(n), n, replace=True)
    sampled_vals = np.random.choice(distance_matrix.flatten(), n, replace=True)
    bootstrap_matrix[range(n), col_pos] = sampled_vals

    #Make symmetric
    bootstrap_matrix = np.triu(bootstrap_matrix)
    bootstrap_matrix += bootstrap_matrix.T

    # Ensure diagonal remains zero
    np.fill_diagonal(bootstrap_matrix, 0)

    return bootstrap_matrix

def load_ids(ids_file):
    """Load sequence IDs from file with flexible format handling"""
    try:
        # First try reading as numpy file
        if ids_file.endswith('.npy'):
            return np.load(ids_file).tolist()

        # Try reading as text file
        with open(ids_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        # Try different formats
        if len(lines) == 0:
            raise ValueError("Empty file")

        # Try tab-delimited format
        if '\t' in lines[0]:
            ids = [line.split('\t')[-1] for line in lines]
        else:
            # Assume one ID per line
            ids = lines

        return ids

    except Exception as e:
        print(f"Error reading IDs file: {e}")
        print("File should contain one ID per line or be tab-delimited")
        sys.exit(1)

def parse_biotite_tree(biotite_tree, ids):
    """Convert biotite tree to Bio.Phylo tree"""
    newick = biotite_tree.to_newick(ids)
    return Phylo.read(StringIO(newick), 'newick')

def analyze_bootstrap_support(original_tree, bootstrap_trees):
    """Calculate bootstrap support values"""
    support_tree = Consensus.get_support(original_tree, bootstrap_trees)


    return support_tree

def get_support_statistics(tree):
    """Calculate statistics about bootstrap support values"""
    support_values = [clade.confidence for clade in tree.find_clades()]
    clade_names = [clade.name for clade in tree.find_clades()]

    support_df = pd.DataFrame()
    support_df['clade'] = clade_names
    support_df['confidence'] = support_values

    if not support_values:
        return None

    support_values =  support_df.confidence.dropna()
    stats = {
        'mean_support': np.mean(support_values),
        'median_support': np.median(support_values),
        'min_support': np.min(support_values),
        'max_support': np.max(support_values),
        'strong_support_nodes': len([v for v in support_values if v >= 70]),
        'moderate_support_nodes': len([v for v in support_values if 50 <= v < 70]),
        'weak_support_nodes': len([v for v in support_values if v < 50]),
        'total_nodes': len(support_values)
    }
    return stats, support_df


def load_tree(tree_path):
    # --- 1. Load the Tree ---
    try:
        tree = Tree(tree_path, format=1) 
    except Exception as e:
        st.error(f"Error loading tree: {e}")

    # --- 2. Define a Tree Style (Optional) ---
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "r"  # 'r' for rectangular, 'c' for circular
    ts.arc_start = -180 # 0 degrees
    ts.arc_span = 180 # 180 degrees

    # --- 3. Render the Tree to a File-like Object (SVG is best) ---
    # We use an in-memory string buffer to avoid writing a temp file
    svg_buffer = StringIO()
    tree.render(svg_buffer, w=600, h=400, units="px", tree_style=ts)
    
    # --- 4. Display the SVG in Streamlit ---
    # Streamlit can display SVG content directly as HTML
    st.markdown(svg_buffer.getvalue(), unsafe_allow_html=True)
    
    return tree


def render_tree(tree_path, user):
    # Load the Tree: i.e. tree = Tree(main_tree_file)
    model_tree = Tree(tree_path, format=1) #format=0
    ts = TreeStyle()
    ts.mode = "r"   
    ts.show_leaf_name = True #true or false
    ts.show_branch_length = False
    ts.force_topology = True #not bad
    # Render tree
    model_tree_out="results/model_trees/" + user + "/model_tree.png"
    model_tree.render(
    model_tree_out,
    tree_style=ts,
    #w=183,
    units="px"  
    )
    return model_tree_out


# plotly_chart---------------------------------------------------------------------
def get_tree_coordinates(tree):
    # This function traverses the tree and assigns x, y coordinates
    
    # 1. Assign Y-coordinates (depth/level)
    # The leaves will be at y=0, and the root will be at the max depth
    leaves = tree.get_leaves()
    max_depth = len(leaves)
    
    for i, leaf in enumerate(leaves):
        leaf.y = i
    
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.x = node.dist # X is the accumulated branch length
        else:
            # For internal nodes, Y is the mean of its children's Y coordinates
            node.y = sum(child.y for child in node.children) / len(node.children)
            # X is the max branch length of its children + its own branch length
            node.x = max(child.x for child in node.children) + node.dist
            
    # Normalize X to be distance from root
    root_x = tree.x
    for node in tree.traverse("preorder"):
        node.x = root_x - node.x

    return tree

def create_plotly_data(tree):
    # This collects all the branch line segments and node scatter points
    branch_lines_x = []
    branch_lines_y = []
    node_coords = {'x': [], 'y': [], 'text': [], 'name': []}

    for node in tree.traverse():
        # Add node coordinates for the scatter plot
        node_coords['x'].append(node.x)
        node_coords['y'].append(node.y)
        node_coords['name'].append(node.name if node.name else f"Internal Node {node.x}")
        
        # Create custom hover text
        hover_text = f"**Node:** {node.name if node.name else 'Internal Node'}<br>"
        if not node.is_root():
            hover_text += f"**Dist:** {node.dist:.3f}"
        node_coords['text'].append(hover_text)
        
        # Create line segments for branches
        if not node.is_root():
            parent = node.up
            # 1. Horizontal segment (from child X to parent X at child Y)
            # Extends the X list: [Child_X, Parent_X, None]
            branch_lines_x.extend([node.x, parent.x, None])
            # Extends the Y list: [Child_Y, Child_Y, None]
            branch_lines_y.extend([node.y, node.y, None])
            
            # 2. Vertical segment (from child Y to parent Y at Parent X)
            # This segment connects the horizontal line to the parent node's vertical position
            # Extends the X list: [Parent_X, Parent_X, None]
            branch_lines_x.extend([parent.x, parent.x, None])
            # Extends the Y list: [Child_Y, Parent_Y, None]
            branch_lines_y.extend([node.y, parent.y, None])
            
    return branch_lines_x, branch_lines_y, node_coords

def display_interactive_tree(tree_path):
    
    try:
        # Load the tree from the Newick file
        tree = Tree(tree_path, format=1)
        
        # Calculate coordinates
        tree = get_tree_coordinates(tree)
        
        # Generate Plotly data lists
        branch_lines_x, branch_lines_y, node_coords = create_plotly_data(tree)
        
    except Exception as e:
        st.error(f"Error processing tree file: {e}")
        return

    # --- Create the Plotly Figure ---
    # 1. Branch Traces (Lines)
    trace_branches = go.Scatter(
        x=branch_lines_x, y=branch_lines_y,
        mode='lines',
        line=dict(color='black', width=1),
        hoverinfo='none', # Branches don't need hover
        name='Branches'
    )

    # 2. Node Traces (Scatter Points)
    trace_nodes = go.Scatter(
        x=node_coords['x'], y=node_coords['y'],
        mode='markers+text',
        marker=dict(size=5, color='red'),
        text=node_coords['name'], # Leaf/Node name displayed next to the point
        textposition="middle right",
        textfont=dict(size=10),
        hovertext=node_coords['text'], # Custom hover text
        hoverinfo='text',
        name='Nodes'
    )

    # 3. Create the figure and layout
    fig = go.Figure(data=[trace_branches, trace_nodes])

    fig.update_layout(
        title='Interactive Tree for model organism ranking',
        xaxis=dict(showgrid=False, title='Branch length'),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False), # Hide y-axis labels
        showlegend=False,
        hovermode="closest",
        plot_bgcolor='white',
        margin=dict(l=20, r=20, t=40, b=20)
    )
    
    # --- Display in Streamlit ---
    st.plotly_chart(fig, use_container_width=True)




# reteive protein pdbs
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
        st.error("CRITICAL ERROR: 'foldseek_local' not found. Please run the 'wget' commands in the terminal to download it.")
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
            status_placeholder.warning("Foldseek ran but found **NO STRUCTURAL HITS**.")
            return None
            
        df = pd.read_csv(temp_m8, sep="\t", names=["qstart", "tstart", "qaln", "taln", "tm"])
        
        if df.empty:
            status_placeholder.warning("Foldseek returned an empty alignment table.")
            return None
            
        #st.success("Foldseek alignment completed successfully!")
        status_placeholder.empty()
        return df.iloc[0]
        
    except subprocess.CalledProcessError as e:
        status_placeholder.error(f"FOLDSEEK FAILED.\n\nStderr: {e.stderr}")
        return None
    except Exception as e:
        status_placeholder.error(f"An unexpected error occurred: {e}")
        return None

def get_format_from_path(filepath):
    """Simple helper to determine PDB/CIF format from file extension."""
    if filepath.lower().endswith(('.cif', '.mmcif')):
        return "cif"
    elif filepath.lower().endswith('.pdb'):
        return "pdb"
    return None

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
    
def render_3dmol_alignment(color, structure_path, target_name, query_name, tm_score, rmsd, spin_on=True):
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
        * **Chain legend:** **{target_name}** (**cyan**) vs. **{query_name}** (**Yellow**)
        * **plDDT legend: >90 (dark blue), 81-90 (cyan), 71-80 (green), 61-70 (yellow), 51-60 (orange), <51 (red)**
    """)

    # 5. Render the view
    try:
        showmol(view, height=510, width=800)
    except Exception as e:
         st.error(f"Error rendering 3Dmol viewer with stmol: {e}")


# --- Function to generate the Mol* HTML viewer ---

def render_3dmol(afdb_download_url: str):
    """
    Fetches the CIF data and renders it using py3Dmol via stmol.
    """
    st.info("Fetching protein structure from AlphaFold database...")
    cif_data = None
    
    try:
        # 1. Fetch the CIF data
        cif_response = requests.get(afdb_download_url, timeout=30)
        cif_response.raise_for_status() 
        cif_data = cif_response.text
        
        if not cif_data or not isinstance(cif_data, str):
             st.error("Fetched data is empty or not a string. Cannot render structure.")
             return
             
    except requests.exceptions.RequestException as e:
        st.error(f"Network error during CIF download: {e}")
        return

    # 2. Initialize the py3Dmol viewer object
    # DO NOT call .render()
    view = py3Dmol.view(width=800, height=500, js='https://3dmol.org/build/3Dmol.js')

    # 3. Add the model
    view.addModel(cif_data, 'cif')

    # 4. Set the coloring scheme (lDDT for AlphaFold)
    # Official Hex Codes (approximate):
    # Very High (>90): Dark Blue (#0053D6)
    # High (70-90): Light Blue/Cyan (#64CCEB)
    # Low (50-70): Yellow (#FFDB00)
    # Very Low (<50): Red/Orange (#FF7D45)
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
    #view.render()
    #view.show()
    
    # 5. Render the view using the stmol component
    # The showmol function handles the Streamlit component embedding correctly.
    try:
        showmol(view, height=510, width=800)
    except Exception as e:
         st.error(f"Error rendering 3Dmol viewer with stmol: {e}")

#---------------------------------------

##MAIN##
st.set_page_config(
    page_title="SHE: Analyse",
    page_icon= "data/tree/logo_tree.png", 
    layout="centered"
)


st.title('SHE: Structural History of Eukarya')
#st.image('data/tree/model.png',ccaption='Protein-protein structure comparison to rank the best model organism for your target organism')
st.write("Discover the best model organism to use for your study based on the proteins involved in your pathway")
st.title('Model Organisms search')
st.subheader('Target organism: Homo sapiens')
st.subheader('Target proteins: Q9NXK8, Q9NVI1')
# Create two columns
col1, col2 = st.columns(2)

# Place the first image in the first column
with col1:
    st.image("data/subset/figure_model/ranking.png", caption="Organisms Ranking", use_container_width=True)

# Place the second image in the second column
with col2:
    st.image("data/subset/figure_model/proteins.png", caption="Protein bitscores", use_container_width=True)



COLOR_MOUSE = "#A6EFFF"  # Light Blue
COLOR_FISH = "#EEF4AE"   # Light Green/Turquoise
COLOR_PLANT = "#F8D7B2"  # Light Orange/Peach

css = """
.st-key-my_blue_container{
    background: rgba(100, 100, 200, 0.3);
}

.st-key-my_blue_container .st-key-make_it_pink img{
    background: linear-gradient(90deg,blue, blue, blue);
}
"""
st.html(f"<style>{css}</style>")

row2 = st.columns(3)
emojis = [":mouse2:", ":tropical_fish:", ":seedling:"]
titles = ["House Mouse", "Zebrafish", "Arabidopsis thaliana"]
ranks = ["#1", "#2", "#3"]
bitscores = ["10197", "869", "224"]
score_1= ["2627", "365", "156"]
proteins1= ["Q9EPX5", "Q6PBJ7", "Q9FM90"]
proteins2= ["Q8K368", "B3DJT4", "Q9LY34"]
score_2= ["7570", "504", "68"]
# Define the custom heights for the *visible* title containers
tile_heights = [120, 100, 80] 
for col, emoji, title, bitscore, height, rank, score1, score2, protein1, protein2 in zip(row2, emojis, titles, bitscores, tile_heights, ranks, score_1, score_2, proteins1, proteins2):
    with col:
        #with st.container(key="my_blue_container")
        # 1. Create the Visible Title Container with Custom Height
        tile = st.container(height=height)
        # Use a large markdown header for the title inside the container
        tile.markdown(f"## {rank} {emoji} {title}")
        
        with st.expander("Click for structural alignment details"):
            # Content that appears when the user clicks/expands
            st.write(f"**Total Bitscore** human vs {title}: **{bitscore}**")
            st.write("---")
            st.write(f"**Bitscore** human Q9NXK8 vs {title} {protein1}: **{score1}**")
            if title == "House Mouse":
                alignment = open("data/subset/figure_model/mouse_Q9NXK8.mp4", "rb")
                bytes = alignment.read()
                st.video(bytes, autoplay=True, loop=True, muted=True)
                st.write("---")
                st.write(f"**Bitscore** human Q9NVI1 vs {title} {protein2}: **{score2}**")
                alignment2 = open("data/subset/figure_model/mouse_Q9NVI1.mp4", "rb")
                bytes2 = alignment2.read()
                st.video(bytes2, autoplay=True, loop=True, muted=True)
            elif title == "Zebrafish":
                alignment = open("data/subset/figure_model/fish_Q9NXK8.mp4", "rb")
                bytes = alignment.read()
                st.video(bytes, autoplay=True, loop=True, muted=True)
                st.write("---")
                st.write(f"**Bitscore** human Q9NVI1 vs {title} {protein2}: **{score2}**")
                alignment2 = open("data/subset/figure_model/fish_Q9NVI1.mp4", "rb")
                bytes2 = alignment2.read()
                st.video(bytes2, autoplay=True, loop=True, muted=True)
            elif title == "Arabidopsis thaliana":
                alignment = open("data/subset/figure_model/plant_Q9NXK8.mp4", "rb")
                bytes = alignment.read()
                st.video(bytes, autoplay=True, loop=True, muted=True)
                st.write("---")
                st.write(f"**Bitscore** human Q9NVI1 vs {title} {protein2}: **{score2}**")
                alignment2 = open("data/subset/figure_model/plant_Q9NVI1.mp4", "rb")
                bytes2 = alignment2.read()
                st.video(bytes2, autoplay=True, loop=True, muted=True)



st.write("Start your analysis by checking the 'Find best Model Organism' box in the sidebar. You will be able to choose your final target species, proteins and available model organisms. This will result in a ranking of the model organisms you have available based on the structural similarity of the proteins in the pathway you are studying compared to the target species.")
st.write("If you wish to visualise the structural alignments of your input proteins, turn on 'Show protein structure alignments'.")

st.warning("Please note that this analysis might require some time. It is possible to register your email in the sidebar, so that your results will be sent to you as soon as they are ready!")



###------sidebar sign up--------------------
#user sign up for updates and analysis results
st.sidebar.subheader('Results', divider="rainbow")
st.sidebar.write('Sign up to receive your results')
st.sidebar.text_input('Your email here:', key="email") #strip against whitespace


# Submission Button in the main content area
if st.sidebar.button('Receive Results'):
    # Retrieve the value from session state using the key "email" (updated from "user")
    user_email = st.session_state.email.strip()
    if not user_email:
        # Only check remaining: must not be empty
        st.error("Please enter an email address.")
    else:
        # If we reach here, the email is non-empty and we save it.
        if save_email(user_email):
            # Updated success message
            st.sidebar.success(f"Thank you! Your email: **{user_email}** .")

################analysis -------------------------------------------
proteomes_path = 'data/tree/selected_proteomes_above_25_formatted_organisms.csv'
new_proteomes_path =  'data/core/merged_evo_rate.csv'
proteomes_df = pd.read_csv(proteomes_path)

st.sidebar.subheader('6. Model Organism Search', divider="violet")
if st.sidebar.checkbox('Find best Model organism'):
    user_email = st.session_state.email.strip()
    if not user_email:
        random_num = random.randint(100000, 999999)
        user = f"user_{random_num}"
    else: 
        user=user_email
    st.title('Find the best model organism for your study according to the structural phylogenetic tree!')
    target = st.text_input(
    'Which organism would you like to study (empty = Default: Homo_sapiens_Human)', key="target")
    
    def add_model_organism():
        st.session_state.model_count += 1 # Correctly increments the count
        
    if 'model_count' not in st.session_state:
        st.session_state.model_count = 1 # Start with 1 mandatory model
        
    model_organisms = []
    # Loop renders the correct number of selectboxes based on model_count
    for i in range(1, st.session_state.model_count + 1):
        # Assuming proteomes_df is defined elsewhere
        selection = st.selectbox(
            f'Which model organism ({i}) do you have available (i.e: Saccharomyces cerevisiae, Arabidopsis thaliana, Drosophila melanogaster, Xenopus tropicalis, Danio rerio, Caenorhabditis elegans, Mus musculus)?',
            proteomes_df['Organism_formatted'],
            key=f"model_{i}"
        )
        model_organisms.append(selection)
    # Use the defined function in on_click
    st.button(
        " Add Another Model Organism", 
        key="add_model_btn", 
        on_click=add_model_organism
    )

    protein_list = st.text_input(
    'Which proteins are you studying in your target organism? (each Uniprot ID separated by a space, i.e: Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616)', key="protein_list")

    if (protein_list == '') or not model_organisms:
            st.warning('Please input the Uniprot ID of the protein(s) to study, the organism to study (Defaul to Human), and the model organisms you can to choose from, to find the best model organism for you study.')
    else:
        st.write('Selected Protein(s) to study: ', protein_list)
        #st.write('Selected Model Organism(s): ', model) # 'Danio_rerio_Zebrafish_Brachydanio_rerio Caenorhabditis_elegans Mus_musculus_Mouse'
        st.write('Selected Model Organism(s): ', ', '.join(model_organisms))
        if target == '':
            target_species='Homo_sapiens_Human'
            st.write('Selected Target Organism: ', target_species) #"Homo_sapiens_Human"
        else:
            target_species = target
            st.write('Selected Target Organism: ', target_species) # not "Homo_sapiens_Human"

        status_placeholder = st.empty()
        # Display the loading message in the placeholder
        status_placeholder.info(f"a. Ranking your model organisms...")
        #st.info('a. Ranking your model organisms...')
        uniprot_ids=protein_list.split() # "Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616"
        #compare_species=model.split()
        compare_species = model_organisms
        
        root_dir='data/decompressed/' #350 gb
        
        # --- 1. Load Mapping ---
        mapping='data/core/model/extracted_and_renamed.csv'    #Path to CSV file with columns id,proteome_id,species_name
        species_map = load_mapping(mapping)
        
        # Get the directory ID for the target
        #target_dir_id = species_map[target_species]
        target_dir_id = species_map.get(target_species)
        if not target_dir_id:
            st.error("Target species not found in mapping.")
            sys.exit(1)
        target_dir = os.path.join(root_dir, target_dir_id)
        
        # Storage
        ranking_scores = {}
        itol_data = {}
        tree_data = {} 
        all_best_hits_dfs = [] #new list for homologous names
        
        # 1. Calculate Human (Target) Data first
        # This ensures Human is in the tree!
        human_score, human_vector = get_human_self_vector(target_dir, uniprot_ids)
        ranking_scores[target_species] = human_score
        tree_data[target_species] = human_vector
        itol_data[target_species] = {uid: 1 for uid in uniprot_ids} # Human always has the proteins

        # 2. Calculate Comparison Species
        for species in compare_species:
            if species not in species_map: continue
            
            compare_dir_id = species_map[species]
            #st.info(f"Processing {species}...")
            csv_file = find_comparison_csv(target_dir, compare_dir_id)
            search_col = 'proteinA'
            
            if not csv_file:
                reverse_dir = os.path.join(root_dir, compare_dir_id)
                if os.path.isdir(reverse_dir):
                    csv_file = find_comparison_csv(reverse_dir, target_dir_id)
                    search_col = 'proteinB'

            if csv_file:
                score, presence, vector, best_hits_df = process_orthologs(csv_file, uniprot_ids, search_col)

                if search_col == 'proteinB':
                    # Create a temporary column to hold the original proteinA data
                    temp_col = best_hits_df['proteinA'].copy()
                    
                    # Move the original proteinB data (the target) to proteinA
                    best_hits_df['proteinA'] = best_hits_df['proteinB']
                    
                    # Move the temporary data (the hit) to proteinB
                    best_hits_df['proteinB'] = temp_col
                
                # Add a column to identify the species for later filtering
                best_hits_df['Compared_Species'] = species 
                all_best_hits_dfs.append(best_hits_df)

                ranking_scores[species] = score
                itol_data[species] = presence
                tree_data[species] = vector
            else:
                ranking_scores[species] = 0.0
                itol_data[species] = {uid: -1 for uid in uniprot_ids}
                tree_data[species] = [0.0] * len(uniprot_ids)

        if all_best_hits_dfs:
            # Use pandas.concat to merge the list of DataFrames into one large DataFrame
            comparison_df = pd.concat(all_best_hits_dfs, ignore_index=True)
        else:
            # Handle the case where no comparisons were successfully processed
            st.warning("No ortholog comparison data was generated.")
            comparison_df = pd.DataFrame() # Create an empty DataFrame
            

        df_rank = pd.DataFrame.from_dict(ranking_scores, orient='index', columns=['Score'])
        df_target = df_rank.loc[[target_species]]
        df_other = df_rank.drop(target_species)
        # Sort by Score in descending order (highest score is Rank 2)
        df_other_sorted = df_other.sort_values(by='Score', ascending=False)
        # The target species is guaranteed to be the first row (index 0)
        df_rank_sorted = pd.concat([df_target, df_other_sorted])
        df_rank = df_rank_sorted.reset_index()
        df_rank.columns = ['Organism', 'Score']
        presence_list = []
        for sp, pres in itol_data.items():
            # Start the dictionary for this row with the organism name
            row_data = {'Organism': sp}
            
            # Iterate through all UniProt IDs
            for uid in uniprot_ids:
                # Get the presence value (default to 0 if not found)
                # Apply the -1/1 mapping: 1 if present, -1 if absent/missing
                presence_val = 1 if pres.get(uid, 0) == 1 else 0 #else 0 or -1
                
                # Add a new key-value pair to the row_data dict: Key=UniProt ID, Value=Presence (-1 or 1)
                row_data[uid] = presence_val
            
            # Append the complete row data to the list
            presence_list.append(row_data)

        df_presence = pd.DataFrame(presence_list)
        # Filter target species out
        target_species_row = df_rank[df_rank['Organism'] == target_species]
        df_rank_filtered = df_rank[df_rank['Organism'] != target_species].copy()
        df_rank_filtered.reset_index(drop=True, inplace=True)
        df_rank_filtered.insert(0, 'Rank', df_rank_filtered.index + 1)
        df_final = pd.merge(df_rank_filtered, df_presence, on='Organism', how='left')
        csv_out = "results/model_trees/" + user + "/"
        os.makedirs(csv_out, exist_ok=True)
        csv_path = csv_out + "model_ranking.txt"
        df_final.to_csv(csv_path, index=False, sep='\t')
        csv_path2 = csv_out + "protein_hits_ranking.txt"
        comparison_df.to_csv(csv_path2, index=False, sep='\t')
        # 6. Display the Final Table in Streamlit
        #st.write("\n" + "="*70)
        st.write("--- Organism Ranking ---")
        #st.write(f"Target: {target_species}")
        #st.write(df_final.to_html(index=False), unsafe_allow_html=True)
        top_model_organism = df_final.iloc[0]['Organism']
        top_model_organism = top_model_organism.split('_', maxsplit=2)[:2]
        top_model_organism = ' '.join(top_model_organism)
        #st.write(f"**The top model organism for your study is: {top_model_organism}**")
        st.markdown(f":orange[**The top model organism for your study is: {top_model_organism}**]")
        status_placeholder.empty()
        # show bar chart
        df_final['Display_Organism'] = df_final['Organism'].str.split('_', n=2).str[:2].str.join('_')
        df_final['Display_Organism'] = df_final['Display_Organism'].str.replace('_', ' ')
        comparison_df['Model_Protein'] = comparison_df['proteinB'].str.split('-').str[1]
        comparison_df['Model_Organism'] = comparison_df['Compared_Species'].str.split('_', n=2).str[:2].str.join('_')
        comparison_df['Model_Organism'] = comparison_df['Model_Organism'].str.replace('_', ' ')
        # 2. Create the color grouping based on proteinA (color parameter)
        comparison_df['Color_Group'] = comparison_df['proteinA'].str.split('-').str[1]
        all_species = comparison_df['Model_Organism'].unique().tolist()
        target_name = target_species.split('_', maxsplit=2)[:2]
        target_name = ' '.join(target_name)

        # add presence
        fixed_cols = ['Rank', 'Organism', 'Score', 'Display_Organism']
        uniprot_id_cols = [col for col in df_final.columns if col not in fixed_cols]
        #Initialize the hover_data dictionary for px.bar
        hover_data_config = {'Score': ':.1f', 'Rank': True, 'Display_Organism': False}
        
        # Add all UniProt ID columns to the dictionary
        for col in uniprot_id_cols:
            hover_data_config[col] = True # Setting to True makes the column appear with its value

        #Create the horizontal bar chart figure, passing the complete hover_data_config
        # Create the Radial Bar Chart (Polar Bar)
        fig = px.line_polar(
            df_final,
            r='Score',                  # Length of the bar (formerly x)
            theta='Display_Organism',   # Category around the circle (formerly y)
            color_discrete_sequence=['#00CC96'],
            line_close=True,
            title='Model Organism Ranking by Bit score',
            hover_data=hover_data_config,
            # distinct start angle makes it easier to read
            start_angle=0,
            direction="clockwise", 
            template='plotly_dark' # Optional: change direction if preferred
        )

        fig.update_traces(
            mode='lines+markers', # Ensures markers and lines are drawn
            fill='toself',        # Ensures the area inside the closed line is filled (the polyshape)
            opacity=0.7,      # OPTIONAL: Add opacity to see overlapping shapes
            marker=dict(size=8))   # OPTIONAL: Make the markers visible
        
        # Update layout for Polar specifics
        fig.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    showticklabels=True, 
                    showline=False,
                    angle=90,   # Position the radial axis labels so they don't overlap data
                    tickangle=90 # Rotate text for readability
                ),
                angularaxis=dict(
                    # Apply the sorting here
                    categoryorder='array',
                    categoryarray=all_species 
                )
            ),
            coloraxis_showscale=False,
            # Optional: Add a hole to make it a "Donut" style radial plot for better label readability
            # polar_hole=0.1 
        )
        st.plotly_chart(fig, use_container_width=True)
        #radar plot---------------------------------------------------------
        # 1. Create the cleaner display name for the proteinB (for hover data)

        # 3. Create a unique list of all Compared_Species names to define the angular axes (theta)
        #all_species = comparison_df['Compared_Species'].unique().tolist()

        # --- Radar Chart Generation ---
        radar_fig = px.line_polar(
            comparison_df,
            r='bit_score',         # Radial axis variable (score magnitude)
            theta='Model_Organism', # Angular axis uses the species name
            color='Color_Group',     # Color based on the Target Protein (proteinA)
            #color_discrete_sequence=px.colors.qualitative.Set1, #Plasma, Dark2, Vidid
            line_close=True,
            title=f'Protein Hit Bitscores by Species compared to {target_name}',
            labels={
                'bit_score': 'Bit Score',
                'Model_Organism': 'Model Organism',
                'Color_Group': 'Target Protein (UniProt ID)'
            },
            # <<< NEW: Add Display_Protein to hover data
            hover_data={
                'Model_Protein': True,  # Show the UniProt ID of the hit protein (proteinB)
                'bit_score': ':.0f',      # Show the bitscore clearly
                'Compared_Species': False # Already on the axis, no need to repeat
            },
            start_angle=0,
            direction="clockwise", 
            template='plotly_dark'
        )

        # Update the layout for better readability
        #radar_fig.update_traces(fill='toself')
        radar_fig.update_traces(
            mode='lines+markers', # Ensures markers and lines are drawn
            fill='toself',        # Ensures the area inside the closed line is filled (the polyshape)
            opacity=0.5,      # OPTIONAL: Add opacity to see overlapping shapes
            marker=dict(size=5)   # OPTIONAL: Make the markers visible
        )

        # Customize the polar axis
        radar_fig.update_polars(
            radialaxis_showgrid=True,
            angularaxis_showgrid=True,
            bgcolor='#000000',
            # Ensure the angular axis uses the full list of unique species names
            angularaxis_categoryorder='array', 
            angularaxis_categoryarray=all_species,
            radialaxis=dict(
                visible=True,
                showticklabels=True,
                showline=False,
                angle=90,   # Position the radial axis labels so they don't overlap data
                tickangle=90, # Rotate text for readability
                # Set range based on the overall max score
                range=[0, comparison_df['bit_score'].max() * 1.1] 
            )
        )

        # Render the radar chart in Streamlit
        #st.plotly_chart(radar_fig, use_container_width=True)


        color_map = {}
        for trace in radar_fig.data:
            # trace.name is the Color_Group (the protein ID part)
            # trace.line.color is the actual color (e.g., 'rgba(31, 119, 180, 1)')
            color = trace.line.color 
            
            if trace.name not in color_map and color:
                color_map[trace.name] = color

        # show protein structures --------------------------------------------------------------------------------------
        if st.sidebar.toggle('Show protein structure alignments'):
            st.subheader(f'Structural comparison of input proteins between {target_name} and model organisms')
            #comparison_df: proteinA, proteinB, bit_Score, Compared_Species
            #target_species: usually human (protein A)
            unique_target_proteins = comparison_df['proteinA'].unique()
            selected_protein_A = st.radio(
                "Select the **Target protein** to visualise:",
                unique_target_proteins,
                horizontal=True
            )
            selected_protein_df = comparison_df[comparison_df['proteinA'] == selected_protein_A]
            if selected_protein_df.empty:
                st.warning(f"No comparison hits found for {selected_protein_A}.")
            else:
                selected_group = selected_protein_df['Color_Group'].iloc[0]
                # Extract the human UniProt ID (AF-Q13616-F1-model_v4 -> Q13616)
                #plotly radial chart
                radar_fig_alone = px.line_polar(
                    selected_protein_df,  # <<< Use the FILTERED DataFrame
                    r='bit_score',
                    theta='Model_Organism',
                    color='Color_Group', # REMOVED: Only one protein is being plotted
                    color_discrete_map={selected_group: color_map[selected_group]},
                    line_close=True,
                    title=f'Structural hits for {selected_protein_A} from {target_name}',
                    labels={
                        'bit_score': 'Bit Score',
                        'Model_Organism': 'Model Organism',
                        'Color_Group': 'Target Protein (UniProt ID)'
                    },
                    hover_data={
                        'Model_Protein': True,
                        'bit_score': ':.0f',
                        'Compared_Species': False
                    },
                    start_angle=0,
                    direction="clockwise",
                    template='plotly_dark'
                )

                # Update the layout for better readability
                radar_fig_alone.update_traces(
                    mode='lines+markers',
                    fill='toself',
                    # Since there is only one trace, we can use a higher opacity for a cleaner look
                    opacity=0.8,
                    marker=dict(size=5)
                )

                # Customize the polar axis (keeping the dark mode customizations)
                radar_fig_alone.update_polars(
                    radialaxis_showgrid=True,
                    angularaxis_showgrid=True,
                    bgcolor='#000000', # Black background for the plot interior
                    angularaxis_categoryorder='array',
                    angularaxis_categoryarray=all_species,
                    radialaxis=dict(
                        visible=True,
                        showticklabels=True,
                        showline=False,
                        angle=90,   # Position the radial axis labels so they don't overlap data
                        tickangle=90,
                        range=[0, selected_protein_df['bit_score'].max() * 1.1] # Use max of the *filtered* data
                    )
                )
                # Render the radar chart in Streamlit
                st.plotly_chart(radar_fig_alone, use_container_width=True)

                try:
                    target_uniprot_id = selected_protein_A.split('-')[1]
                    mapping='data/core/model/extracted_and_renamed.csv'    #Path to CSV file with columns id,proteome_id,species_name
                    proteome_map = load_mapping_proteome(mapping)
                    proteome_target_id = proteome_map.get(target_species) #fix here to get proteome id
                    pdb_target = get_pdb2(target_uniprot_id, proteome_target_id)
                    #st.info(f"info: target_uniprot_id {target_uniprot_id}, proteome_target_id {proteome_target_id}, target species {target_species}, structure {pdb_target}.")
                except Exception as e:
                    st.error(f"Parsing Error: {e}, for target uniprot id {target_uniprot_id}")
                    pdb_target = None
                
                # Check if the target (human) structure can be loaded
                if not pdb_target:
                    st.error(f"Could not retrieve PDB file for the target protein {target_uniprot_id} from proteome {proteome_target_id}, structure {pdb_target}.")
                
                # --- 3. Iterate over the hits in all model organisms ---
                #st.write(f"Structural protein Hits for: **{target_uniprot_id}** from {target_name} ({proteome_target_id})")
                
                # We now iterate over the rows of the filtered DataFrame
                for index, row in selected_protein_df.iterrows():
                    
                    # Variables for the current comparison
                    hit_proteinB = row['proteinB']
                    hit_species = row['Compared_Species']
                    
                    st.markdown("---")
                    st.write(f"### **{hit_species}** (Bit Score: {row['bit_score']:.0f})")
                    
                    # 3.1 Get PDB paths for the model organism hit
                    try:
                        model_uniprot_id = hit_proteinB.split('-')[1]
                        proteome_model_id = proteome_map.get(hit_species)
                        pdb_model = get_pdb2(model_uniprot_id, proteome_model_id)
                    except Exception:
                        st.error(f"Could not parse UniProt ID from {hit_proteinB}.")
                        pdb_model = None

                    if not pdb_model:
                        st.warning(f"Skipping structure comparison: PDB file missing for {hit_species}.")
                        continue

                    #Structural Alignment and Superposition ---
                    st.write(f"Structural Overlap {target_species} {target_uniprot_id} vs. {hit_species} {model_uniprot_id}")
                    alignment_key = f'alignment_toggle_{index}'
                    # Ensure PDB/CIF formats are correctly identified for the SVD/Foldseek functions
                    target_fmt = get_format_from_path(pdb_target)
                    model_fmt = get_format_from_path(pdb_model)

                    # 1. Get Alignment MAP from Foldseek (Model is Query, Target is Target)
                    aln_data = run_foldseek_alignment_only(query_path=pdb_model, target_path=pdb_target)
                    #st.info(f"Alignment data: {aln_data}")

                    if aln_data is not None:
                        
                        # Define output path for the merged PDB file
                        out_file_name = f"foldseek_{target_uniprot_id}_vs_{model_uniprot_id}.pdb"
                        out_file = os.path.join(csv_out, out_file_name)
                        
                        # 2. Calculate Superposition locally using SVD
                        # Target (Fixed) is selected_protein_A, Query (Moving) is hit_proteinB
                        success, rmsd = calculate_svd_and_merge(
                            target_file=pdb_target, target_fmt=target_fmt, 
                            query_file=pdb_model, query_fmt=model_fmt, 
                            aln_data=aln_data, 
                            out_name=out_file
                        )
                        
                        if success:
                            # 3. Render the Merged Structure
                            alignment_color = st.radio(
                                "Color by:",
                                ('chain', 'plDDT'),
                                index=0, # 'plDDT' is the default
                                horizontal=True,
                                key=alignment_key)
                            enable_spin = st.checkbox(
                                "Spin Structure", 
                                value=True, 
                                key=f"spin_{alignment_key}"
                            )
                            render_3dmol_alignment(alignment_color,
                                structure_path=out_file, 
                                target_name=selected_protein_A, 
                                query_name=hit_proteinB, 
                                tm_score=aln_data.get('tm', 0.0), 
                                rmsd=rmsd,
                                spin_on=enable_spin
                            )
                        else:
                            st.warning(f"Could not perform SVD superposition for {hit_species}. Check console for details.")
                            
                    else:
                        st.warning("Foldseek could not find a reliable alignment between the two structures.")     

        # download button
        is_random_user = False
        if not user_email:
            is_random_user = True
        user_dir = os.path.join("results", "model_trees", user)
        os.makedirs(user_dir, exist_ok=True) 
        target_folder_path=user_dir
        user_name=user.split('@')[0]
        zip_file_name = f"results_{user_name}" 
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
            st.sidebar.download_button(
                label="Download results now",
                data=file, # The actual byte content of the zip file
                file_name=f"{zip_file_name}.zip", # The name the user will see
                mime="application/zip", # Critical: Set the correct MIME type for zip files
                icon=":material/download:",
            )
        # Clean up the generated zip file after use
        os.remove(zip_file_path)
        if is_random_user:
            shutil.rmtree(user_dir)
            #temp_fs="tmp_fs"
            #shutil.rmtree(temp_fs)