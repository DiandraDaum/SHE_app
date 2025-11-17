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
import os
from collections import defaultdict
from tqdm import tqdm
from datetime import datetime
import glob

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
    """Loads the mapping file and creates a {species_name: directory_id} dict."""
    try:
        mapping_df = pd.read_csv(mapping_file)
        
        # [FIX 1] We only need the 'id' column, which has the full directory name
        # (e.g., "UP000000589_10090"). This is what we will use to find filenames.
        species_map = dict(zip(mapping_df['species_name'], mapping_df['id']))
        
        return species_map
    except FileNotFoundError:
        st.error(f"Error: Mapping file not found at {mapping_file}")
        sys.exit(1)
    except KeyError as e:
        st.error(f"Error: Mapping file must have columns 'species_name' and 'id'. Missing: {e}")
        sys.exit(1)

def find_comparison_csv(target_dir, compare_dir_id):
    """Finds the specific ...vs_COMPARE_ID_top_hits.csv file within the target dir."""
    # [FIX 1] We now search using the full 'dir_id' (e.g., UP000000589_10090)
    search_pattern = os.path.join(target_dir, '**', f'*vs{compare_dir_id}_top_hits.csv')
    csv_files = glob.glob(search_pattern, recursive=True)
    
    if not csv_files:
        # Fallback search if the naming is slightly different
        search_pattern_alt = os.path.join(target_dir, '**', f'*vs{compare_dir_id}*.csv')
        csv_files = glob.glob(search_pattern_alt, recursive=True)

    if not csv_files:
        st.warning(f"No CSV file found in {target_dir} for comparison with {compare_dir_id}")
        return None
    
    return csv_files[0] # Return the first match

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
            else:
                st.info(f"  > No hit found for {uniprot_id} in {csv_file}")

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
    st.info(f"Calculating self-scores for target {target_dir_id}...")
    
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
        
        st.info(f"  > Max score (baseline) for {uniprot_id} is: {max_score_for_this_protein}")
        total_self_score += max_score_for_this_protein

    return total_self_score

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
st.write("Please note that this analysis might require some time. It is possible to register your email in the sidebar, so that your results will be sent to you as soon as they are ready!")

st.write("Start by checking the 'Find best Model Organism' box in the sidebar. This will result in a ranking of your model organisms.")
st.write("It is also possible to build a structural phylogenetic tree from the input proteins. To do so, switch on 'Build structural tree'.")
st.write("If you wish to visualise the structural alignments of your input proteins, turn on 'Show protein structure alignments'.")

###------sidebar sign up--------------------
#user sign up for updates and analysis results
st.sidebar.subheader('Results', divider="rainbow")
st.sidebar.write('Sign up to receive your results by email')
st.sidebar.text_input('Your email', key="email") #strip against whitespace


# Submission Button in the main content area
if st.sidebar. button('Receive Results'):
    # Retrieve the value from session state using the key "email" (updated from "user")
    user_email = st.session_state.email.strip()
    if not user_email:
        # Only check remaining: must not be empty
        st.error("Please enter an email address.")
    else:
        # If we reach here, the email is non-empty and we save it.
        if save_email(user_email):
            # Updated success message
            st.sidebar.success(f"Thank you! **{user_email}** is where you will receive your results.")
            # Clear the input field after successful submission (updated key)
            #st.session_state.email = ""


################analysis -------------------------------------------
#@st.cache_data
st.sidebar.subheader('6. Model Organism Search', divider="violet")
if st.sidebar.checkbox('Find best Model organism'):
    st.title('Find the best model organism for your study according to the structural phylogenetic tree!')
    protein_list = st.sidebar.text_input(
    'Which proteins are you studying? (each Uniprot ID separated by a space, i.e: Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616)', key="protein")
    target = st.sidebar.text_input(
    'Which organism would you like to study (Default: Homo_sapiens_Human)', key="target")
    model = st.sidebar.text_input(
    'Which model organisms do you have available (i.e: Danio_rerio_Zebrafish_Brachydanio_rerio Caenorhabditis_elegans Mus_musculus_Mouse)?', key="model")
    st.write('Selected Protein(s) to study: ', protein_list)
    st.write('Selected Target Organism: ', target) #"Homo_sapiens_Human"
    st.write('Selected Model Organism(s): ', model) # 'Danio_rerio_Zebrafish_Brachydanio_rerio Caenorhabditis_elegans Mus_musculus_Mouse'
    if (protein_list == '') | (model == ''):
        st.warning('Please input the Uniprot ID of the protein(s) to study, the organism to study (Defaul to Human), and the model organisms you can to choose from, to find the best model organism for you study.')
    else:
        st.info('a. Ranking your model organisms...')
        if target == '':
            target_species='Homo_sapiens_Human'
        else:
            target_species = target

        uniprot_ids=protein_list # "Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616"
        compare_species=model
        root_dir='data/decompressed/' #350 gb
        
        # --- 1. Load Mapping ---
        mapping='data/core/model/extracted_and_renamed.csv'    #Path to CSV file with columns id,proteome_id,species_name
        species_map = load_mapping(mapping)
        
        # Get the directory ID for the target (Human)
        if target_species not in species_map:
            st.error(f"Target species {target_species} not in mapping file.")
            sys.exit(1)

        # [FIX 1] Get the directory ID for the target
        target_dir_id = species_map[target_species]
        target_dir = os.path.join(root_dir, target_dir_id)
        
        if not os.path.isdir(target_dir):
            st.error(f"Target directory not found: {target_dir}")
            sys.exit(1)
            
        st.info(f"Target directory set to: {target_dir}")
        
        # --- 2. Score & Rank Organisms ---
        pathway_scores = defaultdict(float)
        
        for species in compare_species:
            if species not in species_map:
                st.warning(f"'{species}' not in mapping file. Skipping.")
                continue
                
            compare_dir_id = species_map[species]
            st.info(f"--- Scoring '{species}' (vs {compare_dir_id}) ---")
            
            # [FIX] Implement bi-directional search
            csv_file = find_comparison_csv(target_dir, compare_dir_id)
            search_col = 'proteinA' # Default: search proteinA for the UniProt ID
            
            if csv_file is None:
                # If Human-vs-Mouse not found, try Mouse-vs-Human
                st.warning(f"Could not find {target_dir_id} vs {compare_dir_id}. Trying reverse search...")
                reverse_target_dir = os.path.join(root_dir, compare_dir_id)
                
                if not os.path.isdir(reverse_target_dir):
                    st.warning(f"Reverse directory {reverse_target_dir} not found. Skipping {species}.")
                    continue
                    
                csv_file = find_comparison_csv(reverse_target_dir, target_dir_id)
                search_col = 'proteinB' # In reverse file, Human proteins are in column B
            
            if csv_file:
                # Pass the correct column to search
                total_score = sum_scores_for_orthologs(csv_file, uniprot_ids, search_col)
                pathway_scores[species] = total_score
            else:
                st.error(f"Could not find any comparison file for {species}. Score set to 0.")
                pathway_scores[species] = 0.0
            
        # --- 3. Get Self-Score for Baseline ---
        target_self_score = get_self_scores(target_dir, target_dir_id, uniprot_ids)
        pathway_scores[target_species] = target_self_score

        # --- 4. Present Final Ranking ---
        if not pathway_scores:
            st.error("No scores were calculated. Check organism names and paths.")
            sys.exit(1)
            
        scores_df = pd.DataFrame.from_dict(pathway_scores, orient='index', columns=['Total_Pathway_Score'])
        scores_df = scores_df.sort_values(by='Total_Pathway_Score', ascending=False)
        
        # Handle the case where max score is 0 to avoid NaN
        max_score = scores_df['Total_Pathway_Score'].max()
        if max_score > 0:
            scores_df['Normalized_Score'] = scores_df['Total_Pathway_Score'] / max_score
        else:
            scores_df['Normalized_Score'] = 0.0
        
        print("\n" + "="*30)
        print("--- Model Organism Ranking ---")
        print(f"Based on structural similarity for UniProt IDs: {', '.join(uniprot_ids)}")
        print(scores_df.to_string(formatters={'Normalized_Score': '{:,.2f}'.format}))
        print("="*30)
        
    if st.sidebar.toggle('Build structural tree'):
        st.warning('This function is going to be released soon')')
        st.info('b. Making structural phylogenetic tree from your input proteins...')

    if st.sidebar.toggle('Show protein structure alignments'):
        st.warning('This function is going to be released soon')
        st.info('c. Loading the protein-protein comparisons ...')
        
