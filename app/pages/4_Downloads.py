import streamlit as st
import pandas as pd
import base64
import os
import glob
import shutil
import foldcomp as fc
from pathlib import Path

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
                        st.error(f"Error creating zip file: {e}")
            else:
                st.error(f"Proteome {proteome_bit} not found.")

# --- SECTION 2: FOLDCOMP RESULTS ---
with st.expander('Complete Foldcomp results by Proteome'):
    st.subheader('Foldcomp results per Proteome')
    st.caption('See https://github.com/steineggerlab/foldcomp for usage instructions.')
    
    proteome_fc = st.text_input('Which proteome: UniProt ID?', key="proteome_fc")
    
    if proteome_fc:
        st.page_link(f"https://www.uniprot.org/proteomes/{proteome_fc}", label="UniProt Proteomes", icon="🌎")
        
        if st.button("Prepare Download", key="btn_fc"):
            subset_dir = 'data/subset/subset'
            proteome_dir_pattern = os.path.join(subset_dir, f"{proteome_fc}_*")
            matching_paths = glob.glob(proteome_dir_pattern)

            if matching_paths:
                with st.spinner("Compressing folder..."):
                    try:
                        zip_path = zip_directory(matching_paths[0], f"{proteome_fc}_foldseek")
                        with open(zip_path, "rb") as f:
                            st.download_button(
                                label=f"Download Foldseek results for {proteome_fc} (.zip)",
                                data=f,
                                file_name=f"{proteome_fc}_foldseek.zip",
                                mime="application/zip",
                                icon=":material/download:",
                            )
                    except Exception as e:
                        st.error(f"Error: {e}")
            else:
                st.error(f"Proteome {proteome_fc} not found.")

# --- SECTION 3: STRUCTURE PREDICTIONS ---
with st.expander('Protein Structure Predictions'):
    st.subheader('Protein structure predictions by protein')
    st.caption('Download predicted protein structures (.pdb) from AFDB v4')
    
    col_search = st.columns(2)
    with col_search[0]:
        proteome_struct = st.text_input('Proteome ID (e.g. UP000000226)', key="proteome_struct")
    with col_search[1]:
        protein_struct = st.text_input('Protein ID (e.g. V7ATE2)', key="protein_struct")

    if proteome_struct and protein_struct:
        st.page_link(f"https://www.uniprot.org/uniprotkb/{protein_struct}/entry", label="View in UniProtKB", icon="🌎")
        
        if st.button("Find Structure", key="btn_struct"):
            index_dir = 'data/subset/index'
            subset_dir = 'data/subset/subset'
            
            # 1. Find Index
            proteome_index_pattern = os.path.join(index_dir, f"{proteome_struct}_*_fc_id_mapping.csv")
            matching_proteome = glob.glob(proteome_index_pattern)
            
            if not matching_proteome:
                st.error(f"Proteome {proteome_struct} not found.")
            else:
                # 2. Load Index (Cached)
                with st.spinner("Searching index..."):
                    df_index = load_index_file(matching_proteome[0])
                    af_protein_id = f"AF-{protein_struct}-F1-model_v4"
                    
                    if df_index is not None:
                        match = df_index[df_index['afdb_id'] == af_protein_id]
                        
                        if match.empty:
                            st.error(f"Structure {af_protein_id} not found in this proteome.")
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
                                    st.download_button(
                                        label=f"Download {target_name}.pdb",
                                        data=pdb_data,
                                        file_name=f"{target_name}.pdb",
                                        mime="application/octet-stream",
                                        icon=":material/download:",
                                    )
                                    st.snow()
                                except Exception as e:
                                    st.error(f"Error reading database: {e}")
                            else:
                                st.error("Database files missing.")


# --- SECTION 4: TREE DOWNLOADS ---
with st.expander('Tree by Clade'):
    st.subheader('Tree (newick file) by Clade')
    st.caption('Upload these files to iTOL (https://itol.embl.de/) for visualization.')
    
    kingdoms = ["All", "Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"]
    
    # Create columns for a grid layout of buttons
    t_cols = st.columns(3)
    
    for i, k in enumerate(kingdoms):
        file_name = f"kingdoms_{k}_tree.newick"
        file_path = os.path.join("data", "tree", k, file_name)
        
        # Read file into memory once
        file_content = get_file_content(file_path)
        
        if file_content:
            with t_cols[i % 3]:
                label = "Download Eukarya" if k == "All" else f"Download {k}"
                st.download_button(
                    label=label,
                    data=file_content,
                    file_name=f"{k}_tree.newick",
                    mime="text/plain", # Newick is text, not png
                    icon=":material/download:",
                    use_container_width=True
                )

# --- FOOTER ---
st.markdown("<br><hr>", unsafe_allow_html=True)
st.subheader('Cite', divider="rainbow")
st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')

# Logo Footer
col1, col2, col3, col4 = st.columns(4)
with col2:
    st.markdown("<br>", unsafe_allow_html=True)
    st.image("data/logos/scilifelab.png", use_container_width=False)
with col3:
    st.image("data/logos/su.png", use_container_width=False)