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
from Bio import Phylo
from Bio.Phylo import Consensus
from biotite.sequence.phylo import neighbor_joining
from io import StringIO
import copy
import shutil
from pathlib import Path
from itolapi import Itol, ItolExport


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

# find model organisms --------------------------------------------------------------------------------------------------------------------------
@st.cache_data
def load_mapping(mapping_file):
    try:
        mapping_df = pd.read_csv(mapping_file)
        return dict(zip(mapping_df['species_name'], mapping_df['id']))
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
                
                # For iTOL
                if best_score > 0:
                    presence_map[uniprot_id] = 1
                
                # For Tree (Store actual score)
                score_map[uniprot_id] = best_score

    except Exception as e:
        st.warning(f"Error processing {csv_file}: {e}")
        
    # Convert map to list (Vector) ensures order is identical for all species
    score_vector = [score_map[uid] for uid in uniprot_ids]
    
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

def load_tree(tree_path, user):
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

def load_itol(current_dir):
    # Create the Itol class
    test = Itol()
    # Set the tree file
    tree = current_dir + 'model_tree.newick'
    test.add_file(tree)
    #test.add_file(current_dir + 'colors_tol.txt')
    #test.add_file(current_dir + 'labels.txt')
    #test.add_file(current_dir + 'ranges.txt')
    # Add parameters
    test.params['APIkey'] = 'P8puHywKybkQsfX1lKbxzg'
    test.params['treeName'] = 'Model Organisms Tree'
    # Check parameters
    test.print_variables()
    # Submit the tree
    print('')
    print((
        'Uploading the tree.  This may take some time depending on how large the '
        'tree is and how much load there is on the itol server'
    ))
    good_upload = test.upload()
    if not good_upload:
        print('There was an error:' + test.comm.upload_output)
        sys.exit(1)
    # Read the tree ID
    print('Tree ID: ' + test.comm.tree_id)

    # Read the iTOL API return statement
    print('iTOL output: ' + test.comm.upload_output)

    # Website to be redirected to iTOL tree
    print('Tree Web Page URL: ' + test.get_webpage())

    # Warnings associated with the upload
    print('Warnings: ', test.comm.warnings)

    # Export a pre-made tree to svg
    itol_exporter = ItolExport()
    itol_exporter.set_export_param_value('tree', 'model_tree')
    itol_exporter.set_export_param_value('format', 'svg')
    itol_exporter.set_export_param_value('datasetList', 'dataset1')
    # itol_exporter.export(Path('example_pdf.pdf'))
    # print('exported tree to ',export_location)

    # Export the tree above to svg
    print('Exporting to svg')
    itol_exporter = test.get_itol_export()
    export_location = current_dir + 'model_tree.svg'
    itol_exporter.set_export_param_value('format', 'svg')
    itol_exporter.set_export_param_value('datasetList', 'dataset1')
    itol_exporter.export(export_location)
    print('exported tree to ', export_location)




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
st.image(
    'data/tree/model.png',
    caption='Protein-protein structure comparison to rank the best model organism for your target organism'
)

st.write("Here, you can look for the best model organism to use for your study based on the proteins involved in your pathway. Click on 'Find best Model Organism' in the sidebar to start you analysis!")
st.write("Start by checking the 'Find best Model Organism' box in the sidebar. This will result in a ranking of your model organisms.")
st.write("It is also possible to build a structural phylogenetic tree from the input proteins. To do so, switch on 'Build structural tree'.")
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

st.sidebar.subheader('6. Model Organism Search', divider="violet")
if st.sidebar.checkbox('Find best Model organism'):
    user_email = st.session_state.email.strip()
    if not user_email:
        user="user1"
    else: 
        user=user_email
    st.title('Find the best model organism for your study according to the structural phylogenetic tree!')
    target = st.text_input(
    'Which organism would you like to study (empty = Default: Homo_sapiens_Human)', key="target")
    protein_list = st.text_input(
    'Which proteins are you studying in your target organism? (each Uniprot ID separated by a space, i.e: Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616)', key="protein_list")
    model = st.text_input(
    'Which model organisms do you have available (each species separated by a space i.e: Saccharomyces_cerevisiae_strain_ATCC_204508_S288c_Baker_s_yeast Arabidopsis_thaliana_Mouse_ear_cress Drosophila_melanogaster_Fruit_fly Xenopus_tropicalis_Western_clawed_frog_Silurana_tropicalis Danio_rerio_Zebrafish_Brachydanio_rerio Caenorhabditis_elegans Mus_musculus_Mouse)?', key="model")
    if (protein_list == '') | (model == ''):
        st.warning('Please input the Uniprot ID of the protein(s) to study, the organism to study (Defaul to Human), and the model organisms you can to choose from, to find the best model organism for you study.')
    else:
        st.write('Selected Protein(s) to study: ', protein_list)
        st.write('Selected Model Organism(s): ', model) # 'Danio_rerio_Zebrafish_Brachydanio_rerio Caenorhabditis_elegans Mus_musculus_Mouse'
        if target == '':
            target_species='Homo_sapiens_Human'
            st.write('Selected Target Organism: ', target_species) #"Homo_sapiens_Human"
        else:
            target_species = target
            st.write('Selected Target Organism: ', target_species) # not "Homo_sapiens_Human"

        st.info('a. Ranking your model organisms...')
        uniprot_ids=protein_list.split() # "Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616"
        compare_species=model.split()
        
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
                score, presence, vector = process_orthologs(csv_file, uniprot_ids, search_col)
                ranking_scores[species] = score
                itol_data[species] = presence
                tree_data[species] = vector
            else:
                ranking_scores[species] = 0.0
                itol_data[species] = {uid: -1 for uid in uniprot_ids}
                tree_data[species] = [0.0] * len(uniprot_ids)

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
        df_rank.insert(0, 'Rank', df_rank.index + 1)
        df_final = pd.merge(df_rank, df_presence, on='Organism', how='left')
        # Ssave sf as csv
        csv_out="results/model_trees/" + user + "/organisms_ranking.csv"
        df_final.to_csv(csv_out, index=False)
        # 6. Display the Final Table in Streamlit
        st.write("\n" + "="*70)
        st.write("--- Organism Ranking ---")
        #st.write(df_final)
        st.write(df_final.to_html(index=False), unsafe_allow_html=True)
        top_model_organism = df_final.iloc[1]['Organism']
        st.write(f"🎉 **The top model organism for your study is: {top_model_organism}**")
        st.write("="*70)

        #make tree-----------------------------------------------------------------------------------------------------------------------------
        if st.sidebar.toggle('Build structural evolutionary tree'):
            st.subheader('Structural tree for your proteins in the input model organisms')
            st.info('b. Making structural phylogenetic tree from your input proteins...please wait some minutes')
            # distance matrix
            matrix_out="results/model_trees/" + user + "/"
            if not os.path.exists(matrix_out):
                os.makedirs(matrix_out) #mkdir

            dmat_out, ids_out = generate_distance_matrix(tree_data, uniprot_ids, matrix_out)
            distance_matrix=dmat_out
            ids=ids_out
            n_bootstraps=100
            #st.info(f"id matrix file: {distance_matrix}, and ids: {ids}")

            #make tree----------------------------------
            distance_matrix = np.load(distance_matrix)
            ids = pd.read_csv(ids)['id'].values #Read from csv - keep the format consistens
            # Verify number of IDs matches matrix dimensions
            if len(ids) != len(distance_matrix):
                st.error(f"Error: Number of IDs ({len(ids)}) doesn't match matrix dimensions ({len(distance_matrix)})")
                sys.exit(1)
            #print(f"Loaded {len(ids)} IDs and {distance_matrix.shape} distance matrix")
            trees = []
            bootstrap_matrices = []
            #print("Creating original tree...")
            original_tree = neighbor_joining(distance_matrix)
            trees.append(parse_biotite_tree(original_tree, ids))
            #print(f"Creating {n_bootstraps-1} bootstrap trees...")
            for _ in tqdm(range(n_bootstraps-1)):
                dmat = bootstrap_from_distance_matrix(distance_matrix)
                bootstrap_matrices.append(dmat)
                tree = neighbor_joining(dmat)
                trees.append(parse_biotite_tree(tree, ids))

            #print("Calculating bootstrap support values...")
            original_tree = trees[0]
            bootstrap_replicates = trees[1:]
            supported_tree = analyze_bootstrap_support(original_tree, bootstrap_replicates)

            #print("Calculating support statistics...")
            stats, support_df = get_support_statistics(supported_tree)
            outdir="results/model_trees/" + user + "/"
            support_df.to_csv(outdir+'support_per_clade.csv', index=None)
            #print('Saved support metrics to', args.outdir+'support_per_clade.csv')
            #print("Saving results...")
            # Save bootstrap matrices
            np.save(os.path.join(outdir, 'bootstrap_distance_matrices.npy'),
                    np.array(bootstrap_matrices))
            # Save supported tree
            Phylo.write(supported_tree,
                        os.path.join(outdir, 'bootstrap_tree.newick'),
                        'newick')
            # Save statistics
            with open(os.path.join(outdir, 'bootstrap_stats.txt'), 'w') as f:
                f.write("Bootstrap Support Statistics:\n")
                for stat, value in stats.items():
                    if 'support' in stat:
                        f.write(f"{stat}: {value:.1f}%\n")
                    else:
                        f.write(f"{stat}: {value}\n")

                f.write("\nInterpretation Guide:\n")
                f.write("- Support values ≥70%: Strong support\n")
                f.write("- Support values 50-69%: Moderate support\n")
                f.write("- Support values <50%: Weak support\n")

            # load tree ----------------------------------------------
            model_tree=outdir + "/"
            #svg_model = load_tree(model_tree, user)
            #st.image(svg_model, caption='Best model organism tree', width='stretch')
            tree = load_itol(model_tree)

        # show protein structures --------------------------------------------------------------------------------------
        if st.sidebar.toggle('Show protein structure alignments'):
            st.subheader('Structural comparison of your proteins between the input target and model organisms')
            st.warning('This function is going to be released soon')
            st.info('c. Loading the protein-protein comparisons ...')
            """
            #make df with target_protein_id and top structural homologous model_organism_protein_id1, model_organism_protein_id2, model_organism_protein_id3, etc
            uniprot_list=pd.DataFrame(protein_list, orient='index', columns=['target_protein'])
            target_protein = st.sidebar.selectbox(
            'Which target protein do you want to see the structure comparison with your model organisms of?',
            uniprot_list['target_protein'],
            key = 'protein_selector')
            target_uniprot=target_protein
            for model_protein_id in columns
                model_protein_path=""   #path to compressed pdb comparison
                # render pdb with mol*
                render_3dmol(model_protein_path) # render each model protein homologous
            """

        user_dir = os.path.join("results", "model_trees", user)
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
      