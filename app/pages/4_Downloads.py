import streamlit as st
import pandas as pd
import base64
import os
import glob
import shutil
import foldcomp as fc
from pathlib import Path
import requests

# --- CONFIGURATION ---
NAV_ITEMS = [
    {"label": "SHE", "icon": "🏠", "path": "app.py"},
    {"label": "Model Organisms Ranking", "icon": "🐁", "path": "pages/1_Model_Organism_Ranking.py"},
    {"label": "Explore SHE tree", "icon": "🌳", "path": "pages/2_Explore_Tree.py"},
    {"label": "Interactive SHE Tree", "icon": "🔎", "path": "pages/3_Interactive_Tree.py"},
    {"label": "Downloads", "icon": ":material/download:", "path": "pages/4_Downloads.py"},
]

# --- CACHED FUNCTIONS ---
@st.cache_data
def get_base64_image(image_path):
    """Reads an image from disk and converts to Base64 (Cached)."""
    if not os.path.exists(image_path):
        return ""
    with open(image_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode()

@st.cache_data
def load_index_file(file_path):
    """Loads the large index CSV for structure mapping. Vital for performance."""
    if os.path.exists(file_path):
        return pd.read_csv(file_path)
    return None

@st.cache_data
def get_file_content(file_path):
    """Reads a text file (like .newick) into memory."""
    if os.path.exists(file_path):
        with open(file_path, "rb") as f:
            return f.read()
    return None

# hadles input protein mismatches
#@st.cache_data(show_spinner=False)
def get_uniprot_suggestion(invalid_id, target_species_proteome):
    """
    Queries UniProt to find out why an ID is invalid and suggests a fix.
    """
    # Clean inputs
    target_proteome_id = target_species_proteome.strip().upper()
    base_url = "https://rest.uniprot.org/uniprotkb"

    try:
        # 1. Look up the ID
        response = requests.get(f"{base_url}/{invalid_id}")
        
        if response.status_code == 404:
            return None, "ID does not exist in UniProt: please check the correctness of your UniProt ID."
        
        if response.status_code != 200: 
            return None, "UniProt API error."
            
        data = response.json()
        
        # 2. Check if ID is Outdated
        primary_accession = data.get('primaryAccession')
        if primary_accession and primary_accession != invalid_id:
            return primary_accession, f"The ID **{invalid_id}** is outdated. The current primary ID is **{primary_accession}**."

        # Extract useful metadata
        input_gene = None
        if 'genes' in data and len(data['genes']) > 0:
            input_gene = data['genes'][0].get('geneName', {}).get('value', '')
            
        org_sci_name = data.get('organism', {}).get('scientificName', '')
        input_tax_id = data.get('organism', {}).get('taxonId')

        # 3. Check Proteome/Species Mismatch
        found_proteomes = []
        if 'uniProtKBCrossReferences' in data:
            for ref in data['uniProtKBCrossReferences']:
                if ref.get('database') == 'Proteomes':
                    found_proteomes.append(ref.get('id'))
        
        # If the target proteome is NOT in the protein's list, we investigate why.
        if target_proteome_id not in found_proteomes:
            
            # Check if it's actually a Species Mismatch or just a missing Reference Tag
            is_genuine_mismatch = True 
            
            if input_gene:
                # Search for the gene in the TARGET Proteome
                # We request 'organism_id' to verify if the species matches our input
                query = f"gene_exact:{input_gene} AND proteome:{target_proteome_id} AND reviewed:true"
                search_url = f"{base_url}/search?query={query}&fields=accession,organism_id,organism_name&size=1"
                
                search_resp = requests.get(search_url)
                if search_resp.status_code == 200:
                    results = search_resp.json().get('results', [])
                    if results:
                        best_match = results[0]['primaryAccession']
                        match_tax_id = results[0].get('organism', {}).get('taxonId')
                        
                        # --- CRITICAL FIX ---
                        # If the TaxID of the result (Target Proteome) matches the TaxID of the input (User ID),
                        # THEN IT IS THE SAME SPECIES.
                        # It is NOT a mismatch. It is likely just Unreviewed.
                        # We allow the code to continue to Step 4.
                        if match_tax_id == input_tax_id:
                            is_genuine_mismatch = False
                        else:
                            # TaxIDs are different (e.g., Mouse vs Human). This is a real mismatch.
                            return best_match, f"The ID **{invalid_id}** belongs to *{org_sci_name}*. Did you mean the ID **{best_match}** (Gene: {input_gene}) found in your target proteome **{target_proteome_id}**?"

            # If we confirmed it's a different species, return error.
            # If we confirmed it's the SAME species (is_genuine_mismatch=False), we SKIP this block.
            if is_genuine_mismatch:
                 # If we didn't find a gene match to verify, we default to assuming mismatch 
                 # UNLESS we can check Step 4. 
                 # However, if it's truly the same species, Step 4 will likely catch it.
                 # To be safe, if we haven't proven it's the same species, we warn.
                 # But for Q99870, we proved it above, so we skip this return.
                 
                 # Only return if we actually found a different species or really can't map it.
                 if input_gene: # If we had a gene but failed to match TaxID or failed search
                     # If the search returned nothing, we can't prove match, so we might return error or fallback.
                     # Let's fallback to Step 4 if we are unsure.
                     pass 
                 else:
                     actual_prot = found_proteomes[0] if found_proteomes else "Unknown"
                     return None, f"The ID **{invalid_id}** belongs to *{org_sci_name}* (Proteome: {actual_prot}), not your target **{target_proteome_id}**."

        # 4. Check if it is Unreviewed (TrEMBL)
        # Q99870 will now reach here because we skipped the error above!
        entry_type = data.get('entryType', '') 
        
        if "unreviewed" in entry_type.lower() and input_gene and input_tax_id:
            query = f"gene_exact:{input_gene} AND organism_id:{input_tax_id} AND reviewed:true"
            search_url = f"{base_url}/search?query={query}&fields=accession&size=1"
            
            s_resp = requests.get(search_url)
            if s_resp.status_code == 200:
                results = s_resp.json().get('results', [])
                if results:
                    reviewed_id = results[0]['primaryAccession']
                    if reviewed_id != invalid_id:
                        return reviewed_id, f"The ID **{invalid_id}** is an unreviewed (TrEMBL) entry. The reviewed (Swiss-Prot) entry for gene **{input_gene}** is **{reviewed_id}**."

        # 5. Default Fallback
        return None, "ID exists but was probably not reviewed or included in the reference proteome used for AFDB v4."

    except Exception as e:
        return None, f"Could not connect to UniProt: {e}"
    
# --- HELPER FUNCTIONS ---

def zip_directory(source_path, output_name):
    """Zips a directory and returns the path to the zip file."""
    output_base_dir = "temp_archives"
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Create the zip archive
    # root_dir is the parent of the folder we want to zip
    # base_dir is the folder itself relative to root
    archive_path = shutil.make_archive(
        base_name=os.path.join(output_base_dir, output_name),
        format='zip',
        root_dir=os.path.dirname(source_path),
        base_dir=os.path.basename(source_path)
    )
    return archive_path

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="SHE: Downloads",
    page_icon="data/tree/logo_tree.png",
    layout="centered"
)

# --- GLOBAL STYLES ---
st.markdown("""
    <style>
    /* 1. Header Cleanup */
    [data-testid="stHeaderActionElements"] { display: none !important; }
    h1, h2, h3, h4, h5, h6 { padding-top: 0rem !important; }

    /* 2. Banner Styles */
    .banner-container {
        background-size: cover; background-position: center;
        border-radius: 15px; padding: 100px 20px;
        text-align: center; position: relative; margin-bottom: 20px;
    }
    .banner-overlay {
        background-color: rgba(0, 0, 0, 0.4);
        position: absolute; top: 0; left: 0; width: 100%; height: 100%;
        border-radius: 15px; z-index: 1;
    }
    .banner-text {
        position: relative; z-index: 2; color: white;
        transition: all 0.3s ease-in-out; text-shadow: 2px 2px 5px rgba(0,0,0,0.8);
    }
    .banner-text:hover {
        text-shadow: 0 0 10px #fff, 0 0 20px #fff, 0 0 30px #00e6b8ff;
        transform: scale(1.02);
    }

    /* 3. Expander Styling */
    [data-testid="stExpander"] {
        border: 1px solid rgba(250, 250, 250, 0.2);
        border-radius: 12px;
        transition: all 0.3s ease-in-out;
    }
    [data-testid="stExpander"]:hover {
        border-color: #22c55e;
        box-shadow: 0 0 20px rgba(74, 222, 128, 0.4);
    }
    [data-testid="stExpander"]:hover summary { color: #4ade80 !important; }
    [data-testid="stExpander"] summary svg { fill: #9ca3af !important; transition: fill 0.3s; }
    [data-testid="stExpander"]:hover summary svg { fill: #4ade80 !important; }
    
    /* 4. Align Images in Columns */
    [data-testid="stImage"] {
        display: flex; align-items: flex-start !important; justify-content: center; height: 100%;
    }
    img { margin-top: 0 !important; }
    </style>
""", unsafe_allow_html=True)

# --- SIDEBAR ---
with st.sidebar:
    st.title("Browse")
    for item in NAV_ITEMS:
        st.page_link(item["path"], label=item["label"], icon=item["icon"])
    st.divider()
    
    st.subheader('Feedback', divider="yellow")
    st.write("Do you like it?")
    sentiment_mapping = [":material/thumb_down:", ":material/thumb_up:"]
    selected = st.feedback("thumbs")
    if selected is not None:
        st.markdown(f"You selected: {sentiment_mapping[selected]}")

# --- MAIN HEADER ---
st.markdown("""
    <h2 style='
        text-align: center; font-size: 3rem; font-weight: 800; 
        background: linear-gradient(to right, #00FF7F, #40E0D0); 
        -webkit-background-clip: text; -webkit-text-fill-color: transparent; 
        color: transparent;'>
        Downloads
    </h2>
""", unsafe_allow_html=True)

img_b64 = get_base64_image('data/logos/SHE.png')
st.markdown(
    f"""
    <div class="banner-container" style='background-image: url("data:image/png;base64,{img_b64}");'>
        <div class="banner-overlay"></div>
        <div class="banner-text">
            <h2>Download data and results from our analysis</h2>
        </div>
    </div>
    """,
    unsafe_allow_html=True
)

st.write('Please select the type of file you want to download below.')

# --- SECTION 1: BIT SCORES ---
with st.expander('Alignment scores by Organism'):
    st.subheader('Protein-protein alignment bitscores by Organism')
    st.caption('Download Foldseek proteome-proteome comparison bitscores (.csv)')
    proteome_bit = st.text_input('Which proteome: UniProt ID? (e.g. UP000005640)', key="proteome_bit")
    if proteome_bit:
        proteome_bit = proteome_bit.strip().upper()
        st.page_link(f"https://www.uniprot.org/proteomes/{proteome_bit}", label="UniProt Proteomes", icon="🌎")
        # ACTION BUTTON: Prevents heavy operation from running on every keystroke
        if st.button("Prepare Download", key="btn_bit"):
            root_dir_bitscores = 'data/decompressed'
            file_pattern = os.path.join(root_dir_bitscores, f"{proteome_bit}_*")
            matching_paths = glob.glob(file_pattern)

            if matching_paths:
                with st.spinner("Zipping files... This may take a moment."):
                    try:
                        zip_path = zip_directory(matching_paths[0], f"{proteome_bit}_bitscores")
                        st.snow()
                        with open(zip_path, "rb") as f:
                            st.download_button(
                                label=f"Download bitscores for {proteome_bit} (.zip)",
                                data=f,
                                file_name=f"{proteome_bit}_bitscores.zip",
                                mime="application/zip",
                                icon=":material/download:",
                            )
                        # Cleanup done automatically on next run or via OS, 
                        # strictly deleting here might break the download stream if not careful.
                    except Exception as e:
                        print(f"Error creating zip file: {e}")
                        st.stop()
            else:
                st.error(f"Proteome {proteome_bit} not found.")

# --- SECTION 2: FOLDCOMP RESULTS ---
with st.expander('Complete Foldcomp results by Proteome'):
    st.subheader('Foldcomp results per Proteome')
    st.caption('See https://github.com/steineggerlab/foldcomp for usage instructions.')
    proteome_fc = st.text_input('Which proteome: UniProt ID (ie. UP000005640)?', key="proteome_fc")
    if proteome_fc:
        proteome_fc = proteome_fc.strip().upper()
        st.page_link(f"https://www.uniprot.org/proteomes/{proteome_fc}", label="UniProt Proteomes", icon="🌎")
        if st.button("Prepare Download", key="btn_fc"):
            subset_dir = 'data/subset/subset'
            proteome_dir_pattern = os.path.join(subset_dir, f"{proteome_fc}_*")
            matching_paths = glob.glob(proteome_dir_pattern)

            if matching_paths:
                with st.spinner("Compressing folder..."):
                    try:
                        zip_path = zip_directory(matching_paths[0], f"{proteome_fc}_foldseek")
                        st.snow()
                        with open(zip_path, "rb") as f:
                            st.download_button(
                                label=f"Download Foldseek results for {proteome_fc} (.zip)",
                                data=f,
                                file_name=f"{proteome_fc}_foldseek.zip",
                                mime="application/zip",
                                icon=":material/download:",
                            )
                    except Exception as e:
                        print(f"Error: {e}")
            else:
                st.error(f"Proteome {proteome_fc} not found.")

# --- SECTION 3: STRUCTURE PREDICTIONS ---
with st.expander('Protein Structure Predictions'):
    st.subheader('Protein structure predictions by protein')
    st.caption('Download predicted protein structures (.pdb) from AFDB v4')
    
    col_search = st.columns(2)
    with col_search[0]:
        proteome_struct = st.text_input('Proteome ID (e.g. UP000000226)', key="proteome_struct")
        proteome_struct = proteome_struct.strip().upper()
        st.page_link(f"https://www.uniprot.org/proteomes/{proteome_struct}", label="View in UniProt Proteomes", icon="🌎")
    with col_search[1]:
        protein_struct = st.text_input('Protein ID (e.g. V7ATE2)', key="protein_struct")
        protein_struct = protein_struct.strip().upper()
        st.page_link(f"https://www.uniprot.org/uniprotkb/{protein_struct}/entry", label="View in UniProtKB", icon="🌎")
    if proteome_struct and protein_struct:
        if st.button("Prepare Download", key="btn_struct"):
            index_dir = 'data/subset/index'
            subset_dir = 'data/subset/subset'
            
            # 1. Find Index
            proteome_index_pattern = os.path.join(index_dir, f"{proteome_struct}_*_fc_id_mapping.csv")
            matching_proteome = glob.glob(proteome_index_pattern)
            
            if not matching_proteome:
                st.error(f"Proteome {proteome_struct} not found. Please check the input Proteome ID.")
            else:
                # 2. Load Index (Cached)
                with st.spinner("Searching index..."):
                    df_index = load_index_file(matching_proteome[0])
                    af_protein_id = f"AF-{protein_struct}-F1-model_v4"
                    
                    if df_index is not None:
                        match = df_index[df_index['afdb_id'] == af_protein_id]
                        
                        if match.empty:
                            suggested_id, reason_msg = get_uniprot_suggestion(protein_struct, proteome_struct)
                            if suggested_id:
                                st.error(f"Protein {protein_struct} not found: {reason_msg}")
                            else:
                                st.error(f"Protein {protein_struct} not found in proteome {proteome_struct}: {reason_msg}")
                                st.page_link(f"https://www.uniprot.org/uniprotkb/{af_protein_id}/entry", label=f"Check {af_protein_id} existence in {proteome_struct} proteome and potential updated reviewed ID on UniProt", icon="🌎")
                        else:
                            # 3. Extract from Foldcomp
                            retrieved_index = match['fc_db_index'].iloc[0]
                            
                            # Find DB
                            wild_proteome_dir = os.path.join(subset_dir, f"{proteome_struct}_*")
                            proteome_dirs = glob.glob(wild_proteome_dir)
                            
                            if proteome_dirs:
                                proteome_dir = proteome_dirs[0]
                                full_proteome_id = Path(proteome_dir).name
                                db_path = os.path.join(proteome_dir, full_proteome_id)
                                
                                try:
                                    with fc.open(db_path) as db:
                                        target_name, pdb_data = db[retrieved_index]
                                    # Provide download directly from memory
                                    st.snow()
                                    st.download_button(
                                        label=f"Download {target_name}.pdb",
                                        data=pdb_data,
                                        file_name=f"{target_name}.pdb",
                                        mime="application/octet-stream",
                                        icon=":material/download:",
                                    )
                                except Exception as e:
                                    print(f"Error reading database: {e}")
                            else:
                                st.error("Database files not found for {proteome_struct}.")


# --- SECTION 4: TREE DOWNLOADS ---
with st.expander('Tree by Clade'):
    st.subheader('Tree (newick file) by Clade')
    st.caption('Upload these files to iTOL (https://itol.embl.de/) for visualization.')
    
    kingdoms = ["Eukarya", "Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"]
    
    # Create columns for a grid layout of buttons
    t_cols = st.columns(3)
    
    for i, k in enumerate(kingdoms):
        file_name = f"kingdoms_{k}_tree.newick"
        file_path = os.path.join("data", "tree", k, file_name)
        
        # Read file into memory once
        file_content = get_file_content(file_path)
        
        if file_content:
            with t_cols[i % 3]:
                label = "Download Eukarya" if k == "Eukarya" else f"Download {k}"
                st.download_button(
                    label=label,
                    data=file_content,
                    file_name=f"{k}_tree.newick",
                    mime="text/plain", # Newick is text, not png
                    icon=":material/download:",
                    width='stretch'
                )

# --- FOOTER ---
st.markdown("<br><hr>", unsafe_allow_html=True)
st.subheader('Cite', divider="rainbow")
#st.write('Preprint:  Qiuzhen Li, Diandra Daumiller, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')
st.write('Qiuzhen Li, Diandra Daumiller, Patrick Bryant: The structural History of Eukarya. bioRxiv 2026.02.04.703859; doi: https://doi.org/10.64898/2026.02.04.703859')

# Logo Footer
col1, col2, col3, col4 = st.columns(4)
with col2:
    st.markdown("<br>", unsafe_allow_html=True)
    st.image("data/logos/scilifelab.png", width='content')
with col3:
    st.image("data/logos/su.png", width='content')