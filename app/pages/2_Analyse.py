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

from collections import defaultdict
from tqdm import tqdm

#Code to analyse:
# model organism search
#gene set search --> tree build

EMAIL_FILE = 'data/email_results.csv'
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
@st.cache_data
# save user email
def save_email(email):
    """Saves the email and current timestamp to the CSV file. 
    It ensures the directory exists before attempting to write."""
    try:
        # Ensure the directory exists (necessary for 'data/email_results.csv')
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
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            writer.writerow([email, timestamp])
        return True
    except Exception as e:
        st.error(f"An error occurred while saving: {e}")
        return False

# find model organisms -------------------
def load_data(families_file, mapping_file):
    """
    Loads families and mapping files.
    REMOVED: Loading hits_file. We will do that later.
    """
    print("Loading families and mapping files...")
    try:
        families_df = pd.read_csv(families_file)
        mapping_df = pd.read_csv(mapping_file)
    except FileNotFoundError as e:
        print(f"Error: File not found. {e}")
        sys.exit(1)
    except pd.errors.EmptyDataError as e:
        print(f"Error: File is empty. {e}")
        sys.exit(1)

    # Create species_name -> proteome_id map
    try:
        species_to_proteome = dict(zip(mapping_df['species_name'], mapping_df['proteome_id']))
    except KeyError:
        print(f"Error: Mapping file {mapping_file} must have 'species_name' and 'proteome_id' columns.")
        sys.exit(1)
    
    print("Families and mapping data loaded.")
    return families_df, species_to_proteome

def load_relevant_hits(hits_file, target_proteins):
    """
    NEW FUNCTION: Streams the giant hits_file and only keeps
    hits relevant to our target proteins.
    """
    print(f"\nStreaming {hits_file} to find relevant hits...")
    
    # Get the set of full protein IDs we care about
    target_ids = {data['full_id'] for data in target_proteins.values()}
    
    if not target_ids:
        print("  > No target protein IDs found. Returning empty hits DataFrame.")
        return pd.DataFrame(columns=['protein1', 'protein2', 'bit_score'])

    chunksize = 10_000_000
    relevant_chunks = []

    try:
        reader = pd.read_csv(hits_file, chunksize=chunksize, 
                           usecols=["protein1", "protein2", "bit_score"])
        
        for chunk in tqdm(reader, desc="Filtering hits file"):
            # Keep rows where EITHER protein1 OR protein2 is in our target set
            relevant_chunk = chunk[
                (chunk['protein1'].isin(target_ids)) | 
                (chunk['protein2'].isin(target_ids))
            ]
            
            if not relevant_chunk.empty:
                relevant_chunks.append(relevant_chunk)

    except FileNotFoundError:
        print(f"  > ERROR: Hits file not found at {hits_file}")
        sys.exit(1)
    except Exception as e:
        print(f"  > ERROR reading hits file: {e}")
        sys.exit(1)

    if not relevant_chunks:
        print("  > No hits found for any target proteins.")
        return pd.DataFrame(columns=['protein1', 'protein2', 'bit_score'])

    # Combine all the small, relevant chunks into one DataFrame
    relevant_hits_df = pd.concat(relevant_chunks)
    print(f"  > Finished. Found {len(relevant_hits_df)} relevant hits.")
    return relevant_hits_df

def find_target_proteins_by_uniprot(uniprot_ids, target_species, species_to_proteome, families_df):
    """
    Finds the full protein IDs and family IDs for the target UniProt IDs
    within a *specific* target species.
    (This function is UNCHANGED)
    """
    print(f"\nFinding target proteins for '{target_species}'...")
    
    if target_species not in species_to_proteome:
        print(f"  > ERROR: Target species '{target_species}' not found in mapping file. Check spelling.")
        sys.exit(1)
        
    target_proteome_id = species_to_proteome[target_species]
    print(f"  > Target Proteome ID: {target_proteome_id}")
    
    target_proteins = {} # uniprot_id -> {full_id, family_id}
    
    for uniprot_id in uniprot_ids:
        # Precise search: Must contain *both* the proteome ID and the UniProt ID
        matches = families_df[
            families_df['protein'].str.contains(str(target_proteome_id), na=False) &
            families_df['protein'].str.contains(f"-{uniprot_id}-", na=False)
        ]
        
        # Fallback for different naming conventions
        if matches.empty:
            matches = families_df[
                families_df['protein'].str.contains(str(target_proteome_id), na=False) &
                families_df['protein'].str.contains(uniprot_id, na=False)
            ]

        if matches.empty:
            print(f"  > WARNING: Could not find '{uniprot_id}' in proteome '{target_proteome_id}'. Skipping this ID.")
            continue
        
        if len(matches) > 1:
            print(f"  > Warning: Found {len(matches)} matches for '{uniprot_id}'. Using first one.")
            
        full_id = matches['protein'].values[0]
        family_id = matches['family_id'].values[0]
        target_proteins[uniprot_id] = {'full_id': full_id, 'family_id': family_id}
        print(f"  > Found '{uniprot_id}': {full_id} (Family: {family_id})")
        
    return target_proteins

def find_best_hit_in_organism(target_protein_id, target_family_id, compare_proteome_id, families_df, hits_df):
    """
    Finds the single best bit_score for a target protein in a specific organism.
    (This function is UNCHANGED, but now receives the *small* relevant hits_df)
    """
    
    # 1. Get all members of the target's family
    family_member_ids = set(families_df[families_df['family_id'] == target_family_id]['protein'])
    
    # 2. Filter family members to *only* those in the organism we are comparing
    compare_organism_proteins = {
        p for p in family_member_ids if str(p).startswith(str(compare_proteome_id))
    }
    
    if not compare_organism_proteins:
        return 0.0 # No proteins from this family exist in the compare organism

    # 3. Find all hits between the target and this organism's proteins
    # This query is now fast because hits_df is small
    hits_v1 = hits_df[
        (hits_df['protein1'] == target_protein_id) & 
        (hits_df['protein2'].isin(compare_organism_proteins))
    ]
    hits_v2 = hits_df[
        (hits_df['protein2'] == target_protein_id) & 
        (hits_df['protein1'].isin(compare_organism_proteins))
    ]
    
    all_hits = pd.concat([hits_v1, hits_v2])
    
    if all_hits.empty:
        # This can happen if the only hits were below the threshold used
        # to create all_hits.csv
        return 0.0 
        
    # 4. Return the single best bit_score
    return all_hits['bit_score'].max()

#---------------------------------------

##MAIN##
st.set_page_config(
    page_title="SHE: Analyse",
    page_icon= "data/tree/logo_tree.png", 
    layout="centered"
)



st.title('SHE: Structural History of Eukarya')
st.subheader('Strucutral Phylogenetic Tree')
st.image(
    'data/tree/All/Intro.png',
    caption='Strucutral Phylogenetic Tree of Eukarya'
)

st.write("Here, you can look for the best model organism to use for your study based on the proteins involved in your pathway. Click on 'Find best Model Organism' in the sidebar to start you analysis!")
st.write("What about studying the evolution of your proteins of interest? Click on 'Study Protein Evolution' in the sidebar and start building your tree!")
st.write("Please note that these analysis might require some time. It is possible to register your email in the sidebar, so that your results will be sent to you as soon as they are ready!")


###------sidebar sign up--------------------
#user sign up for updates and analysis results
st.sidebar.subheader('Results')
st.sidebar.write('Sign up to receive your results by email')
st.sidebar.text_input('Your email', key="email") #strip against whitespace


# Submission Button in the main content area
if st.sidebar. button('Receive Results'):
    # Retrieve the value from session state using the key "email" (updated from "user")
    user_email = st.session_state.email.strip()
    
    if not user_email:
        # Only check remaining: must not be empty
        st.error("Please enter an email address.")
        
        # If we reach here, the email is non-empty and we save it.
        if save_email(user_email):
            # Updated success message
            st.success(f"Thank you! **{user_email}** is where you will receive your results.")
            # Clear the input field after successful submission (updated key)
            st.session_state.email = ""


################analysis -------------------------------------------
st.sidebar.subheader('Best model organism')
if st.sidebar.checkbox('Find best Model organism'):
    st.title('Find the best model organism for your study according to the structural phylogenetic tree!')
    protein_list = st.sidebar.text_input(
    'Which proteins are you studying? (each Uniprot ID separated by a space, i.e: Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616)', key="protein")
    target = st.sidebar.text_input(
    'Which organism would you like to study (Default: Homo_sapiens_Human)', key="target")
    model = st.sidebar.text_input(
    'Which model organisms do you have available?', key="model")
    st.write('Selected Protein(s) to study: ', protein_list)
    st.write('Selected Target Organism: ', target) #"Homo_sapiens_Human"
    st.write('Selected Model Organism(s): ', model) # 'Danio_rerio_Zebrafish_Brachydanio_rerio Caenorhabditis_elegans Mus_musculus_Mouse'
    if (protein_list == '') | (model == ''):
        st.warning('Please input the Uniprot ID of the protein(s) to study, the organism to study (Defaul to Human), and the model organisms you can to choose from, to find the best model organism for you study.')
    else:
        if target == '':
            target_species='Homo_sapiens_Human'
        else:
            target_species = target
        """
        uniprot_ids=protein_list # "Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616"
        families=   #Path to the protein_families.csv file: protein_families_173_without0_DSU_OUTPUT.csv
        mapping='data/core/model/extracted_and_renamed.csv'    #Path to CSV file with 'proteome_id' and 'species_name' columns
        
        families_df, species_to_proteome = load_data(families, mapping)
        target_proteins = find_target_proteins_by_uniprot(uniprot_ids, target_species, species_to_proteome, families_df)
        if not target_proteins:
            st.warning("No target proteins were found in the selected model organisms")
        hits=   #Path to the all_hits.csv file (raw homologous pairs): all_homologous_173_without0_pairs_from_decompressed_files.csv
        hits_df = load_relevant_hits(hits, target_proteins)
        
        ####print score
        print("\nRanking model organisms...")
        pathway_scores = defaultdict(float)
        compare = model #List of model organisms to compare (e.g., 'Mus_musculus_Mouse')
        species_to_score = list(compare)
        if target_species not in species_to_score:
            species_to_score.append(target_species)
        for species in species_to_score:
            if species not in species_to_proteome:
                st.warning(f" Skipping'{species}' because canno tbe mapped.")
                continue

            compare_proteome_id = species_to_proteome[species]
            print(f"  > Scoring '{species}'...")
            total_pathway_score = 0.0
            for uniprot_id, data in target_proteins.items():
                # This function now receives the small, pre-filtered hits_df
                best_gene_score = find_best_hit_in_organism(
                    data['full_id'],
                    data['family_id'],
                    compare_proteome_id,
                    families_df,
                    hits_df # <-- This is the small DataFrame!
                )
                
                if species == target_species and best_gene_score == 0.0:
                    # This query is also fast now
                    max_score1 = hits_df[hits_df['protein1'] == data['full_id']]['bit_score'].max()
                    max_score2 = hits_df[hits_df['protein2'] == data['full_id']]['bit_score'].max()
                    best_gene_score = max(np.nan_to_num(max_score1), np.nan_to_num(max_score2))
                    print(f"    - (Self-score for {uniprot_id} set to: {best_gene_score})")

                total_pathway_score += best_gene_score
                
            pathway_scores[species] = total_pathway_score
        
        if not pathway_scores:
            st.warning("Error: No scores were calculated. Check organism names.")
        
        scores_df = pd.DataFrame.from_dict(pathway_scores, orient='index', columns=['Total_Pathway_Score'])
        scores_df = scores_df.sort_values(by='Total_Pathway_Score', ascending=False)
        
        scores_df['Normalized_Score'] = scores_df['Total_Pathway_Score'] / scores_df['Total_Pathway_Score'].max()
        
        st.success("\n--- Model Organism Ranking ---")
        st.write(f"Based on structural similarity for UniProt IDs: {', '.join(uniprot_ids)}")
        st.write(scores_df.to_string(formatters={'Normalized_Score': '{:,.2f}'.format}))"""


st.sidebar.subheader('Structural evolution of a set of protein')
if st.sidebar.checkbox('Study Protein evolution'):
    st.title('Build the structural phylogenetic tree from your input proteins!')
    prots = st.sidebar.text_input(
    'Which proteins are you studying?', key="prots")
    st.write('Selected Proteins: ', prots)
    if (prots == ''):
        st.warning('Please input the Uniprot ID of the protein(s) you want to build the structural evolutionary tree from.')
    else:
        pass