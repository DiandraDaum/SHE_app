import os
import streamlit as st
import pandas as pd
import base64

# --- CONFIGURATION ---
NAV_ITEMS = [
    {"label": "SHE", "icon": "🏠", "path": "app.py"},
    {"label": "Model Organisms Ranking", "icon": "🐁", "path": "pages/1_Model_Organism_Ranking.py"},
    {"label": "Explore SHE tree", "icon": "🌳", "path": "pages/2_Explore_Tree.py"},
    {"label": "Interactive SHE Tree", "icon": "🔎", "path": "pages/3_Interactive_Tree.py"},
    {"label": "Downloads", "icon": ":material/download:", "path": "pages/4_Downloads.py"},
]

# --- CACHED FUNCTIONS (PERFORMANCE) ---

@st.cache_data
def kingdom_color(kingdom):
    """Returns the hex color associated with a specific kingdom."""
    colors = {
        'Opisthokonta': "#bfdbfe",
        'Viridiplantae': "#bbf7d0",
        'Amoebozoa': "#fecdd3",
        'Sar': "#f3f4f6",
        'Discoba': "#e9d5ff",
        'Rhodophyta': "#c3f8f1",
        'Cryptophyceae': "#f6fad5",
        'Apusozoa': "#ddd6fe"
    }
    return colors.get(kingdom, "#FAF8F8")

@st.cache_data
def get_base64_image(image_path):
    """Reads an image from disk and converts to Base64 (Cached)."""
    if not os.path.exists(image_path):
        return ""
    with open(image_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode()

@st.cache_data
def load_csv_data(file_path):
    """Loads CSV data into a DataFrame (Cached)."""
    if os.path.exists(file_path):
        return pd.read_csv(file_path)
    return None

# --- PAGE SETUP ---
st.set_page_config(
    page_title="SHE: Explore",
    page_icon="data/tree/logo_tree.png",
    layout="centered"
)

# --- GLOBAL STYLES ---
st.markdown("""
    <style>
    /* 1. Header Cleanup */
    [data-testid="stHeaderActionElements"] { display: none !important; }
    h1, h2, h3, h4, h5, h6 { padding-top: 0rem !important; }

    /* 2. Banner & Text Styles */
    .banner-container {
        background-size: cover;
        background-position: center;
        border-radius: 15px;
        padding: 100px 20px;
        text-align: center;
        position: relative;
        margin-bottom: 20px;
    }
    .banner-overlay {
        background-color: rgba(0, 0, 0, 0.4);
        position: absolute; top: 0; left: 0; width: 100%; height: 100%;
        border-radius: 15px; z-index: 1;
    }
    .banner-text {
        position: relative; z-index: 2; color: white;
        transition: all 0.3s ease-in-out;
        text-shadow: 2px 2px 5px rgba(0,0,0,0.8);
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
    </style>
""", unsafe_allow_html=True)

# --- SIDEBAR ---
with st.sidebar:
    st.title("Browse")
    for item in NAV_ITEMS:
        st.page_link(item["path"], label=item["label"], icon=item["icon"])
    st.divider()

    # Data for selection
    df_kingdoms = ["Eukarya", "Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"]
    
    st.subheader('Visualization Options')
    layout = st.selectbox('Tree Layout', ['Circular', 'Linear'])
    labels = st.toggle('Labels', value=False)
    
    st.subheader('1. Visualise Clades', divider=True)
    kingdom = st.selectbox(
        'Which major Taxon? (default: Eukarya)',
        df_kingdoms,
        key='kingdom_selector'
    )

    # Dynamic Style for the Selectbox
    color = kingdom_color(kingdom)
    st.markdown(f"""
        <style>
        div[data-baseweb="select"] > div:first-child {{
            background-color: {color};
            border-radius: 0.5rem;
        }}
        div[data-baseweb="select"] > div:first-child div,
        div[data-baseweb="select"] > div:first-child span,
        div[data-baseweb="select"] svg {{
            color: black !important;
            fill: black !important;
            -webkit-text-fill-color: black !important;
            font-weight: 500;
        }}
        </style>
    """, unsafe_allow_html=True)

# --- MAIN PAGE HEADER ---
st.markdown("""
    <h2 style='
        text-align: center; 
        font-size: 3rem; 
        font-weight: 800; 
        background: linear-gradient(to right, #00FF7F, #40E0D0); 
        -webkit-background-clip: text; 
        -webkit-text-fill-color: transparent; 
        color: transparent;'>
        Explore SHE tree
    </h2>
""", unsafe_allow_html=True)

# Banner Image
img_b64 = get_base64_image('data/logos/SHE.png')
st.markdown(
    f"""
    <div class="banner-container" style='background-image: url("data:image/png;base64,{img_b64}");'>
        <div class="banner-overlay"></div>
        <div class="banner-text">
            <h2>Discover the eukaryotic structural phylogeny and evolutionary rates</h2>
        </div>
    </div>
    """,
    unsafe_allow_html=True
)

# Info Tiles
row1 = st.columns(2)
with row1[0].expander("**1. Explore Clades**"):
    st.write(':blue[Start exploring the new structural phylogenetic tree by selecting a major taxon in the sidebar!]')
    st.write(':blue[You can also choose the tree visualisation style in the sidebar.]')

with row1[1].expander("**2. Proteome evolutionary rate**"):
    st.write(':green[Compare the innovative structural evolutionary rate of different proteomes across eukarya! Use the "Visualise Evolutionary Rate" section in the sidebar.]')


# --- SECTION 1: KINGDOM TREES ---
if 'show_tree' not in st.session_state:
    st.session_state.show_tree = True

button_label = '⬇️ Hide main Tree' if st.session_state.show_tree else '⬆️ Show main Tree'

if st.button(button_label):
    st.session_state.show_tree = not st.session_state.show_tree

if st.session_state.show_tree:
    st.title(f'Clade tree: {kingdom}')
    
    # Wiki Link Logic
    wiki_suffix = f"{kingdom}_supergroup" if kingdom == "Sar" else ("Eukaryote" if kingdom == "Eukarya" else kingdom)
    st.page_link(f"https://en.wikipedia.org/wiki/{wiki_suffix}", label="Wiki", icon="🌎")

    # Construct Path Efficiently
    suffix = "_labels.svg" if labels else ".svg"
    # Filename structure: kingdoms_Amoebozoa_Circular_tree1000.svg
    filename = f"kingdoms_{kingdom}_{layout}_tree1000{suffix}"
    kingdom_tree_path = os.path.join("data", "tree", kingdom, filename)

    if os.path.exists(kingdom_tree_path):
        st.image(kingdom_tree_path, caption=f"{kingdom} Phylogenetic Tree", width='stretch')
    else:
        st.warning(f"Tree visualization not found: {filename}")

# Branch Rates Data (Optimized with Cache)
if st.sidebar.toggle('Show species branch rate in kingdom'):
    csv_filename = f"kingdoms_{kingdom}_branch_rates.csv"
    csv_path = os.path.join("data", "tree", kingdom, csv_filename)
    
    # Use the cached function
    df_rates = load_csv_data(csv_path)
    
    if df_rates is not None:
        st.dataframe(df_rates)
    else:
        st.info("Branch rate data not available for this kingdom.")

# --- SECTION 2: EVOLUTIONARY RATE ---
st.sidebar.divider()
st.sidebar.subheader('2. Evolutionary rate', divider="green")

if st.sidebar.toggle('See tree with evolutionary rate'):
    st.title('Evolutionary rate')
    st.caption(f"Showing {kingdom} ({layout}) - Yellow: High rate, Purple: Low rate")

    # Construct Path
    suffix_evo = "_labels.svg" if labels else ".svg"
    # Filename structure: kingdoms_Amoebozoa_Circular_tree_evorate.svg
    filename_evo = f"kingdoms_{kingdom}_{layout}_tree_evorate{suffix_evo}"
    evo_tree_path = os.path.join("data", "tree", kingdom, filename_evo)

    if os.path.exists(evo_tree_path):
        st.image(evo_tree_path, caption=f'Evolutionary Rate ({layout})', width='stretch')
    else:
        st.warning(f"Evolutionary rate tree not found for {kingdom} in {layout} layout.")

st.sidebar.divider()
st.sidebar.page_link("pages/4_Downloads.py", label="Download trees", icon="📥")

# --- FOOTER ---
st.markdown("<br><hr>", unsafe_allow_html=True)
st.subheader('Cite', divider="rainbow")
st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')