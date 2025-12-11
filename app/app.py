import streamlit as st
import pandas as pd
import numpy as np
import time
import csv
from ete3 import Tree, TreeStyle, NodeStyle
import sys
import pathlib
import base64
import textwrap
import os
from datetime import datetime

# --- CONFIGURATION & CONSTANTS ---
EMAIL_FILE = 'data/email_updates.csv'
VIDEO_PATH_BG = "data/tree/bg.mp4" 

# Define navigation items for the custom header
NAV_ITEMS = [
    {"label": "Home", "icon": "🏠", "path": "app.py"},
    {"label": "Tree Viewer", "icon": "🌳", "path": "pages/1_Explore_Tree.py"},
    {"label": "Tree Viewer", "icon": "🍃", "path": "pages/2_Interactive_Tree.py"},
    {"label": "Model Organisms", "icon": "🐁", "path": "pages/3_Model_organism_ranking.py"},
    {"label": "Downloads", "icon": "⬇️", "path": "pages/4_Downloads.py"},
]

# Define features for the main grid section
FEATURES = [
    {"icon": "🌳", "title": "Structural Phylogeny", "description": "Explore the structural phylogenetic tree for the main subkingdoms including Opisthokonta, Viridiplantae, SAR, Amoebozoa, Discoba, and Rhodophyta.", "link_path": "pages/1_Explore_Tree.py", "demo": "data/subset/figure_model/plant_Q9NXK8.mp4"}, #"demo": "data/demos/demo1.mp4"
    {"icon": "📈", "title": "Evolutionary Rate", "description": "Visualize colored branches by evolutionary rate to identify fast and slow-evolving proteomes.", "link_path": "pages/1_Explore_Tree.py", "demo": "data/subset/figure_model/mouse_Q9NXK8.mp4"},
    {"icon": "🍄", "title": "Taxon Search", "description": "Zoom on specific taxa or species to discover new evolutionary relationships.", "link_path": "pages/2_Interactive_Tree.py", "demo": "data/subset/figure_model/fish_Q9NXK8.mp4"},
    {"icon": "🔎", "title": "Protein Structural Homologous", "description": "Input UniProt IDs to analyze presence/absence of structural homologs across all eukaryotic organisms.", "link_path": "pages/2_Interactive_Tree.py", "demo": "data/subset/figure_model/plant_Q9NVI1.mp4"}, # Reusing Model Organisms page for search
    {"icon": "🐁", "title": "Model Organism Ranking", "description": "Select proteins of interest to find the best model organisms where they are most conserved.", "link_path": "pages/3_Model_organism_ranking.py", "demo": "data/subset/figure_model/mouse_Q9NVI1.mp4"},
    {"icon": "⬇️", "title": "Data Downloads", "description": "Download FoldSeek comparison bitscores, Newick files, and predicted protein structures from each proteome.", "link_path": "pages/4_Downloads.py", "demo": "data/subset/figure_model/fish_Q9NVI1.mp4"},
]

# Define methodology steps
STEPS = [
    {"number": "01", "title": "Proteome Collection", "description": "1,542 complete eukaryotic proteomes from UniProt with predicted structures."},
    {"number": "02", "title": "Structural Alignment", "description": "All-vs-all protein alignment using FoldSeek for rapid structure comparison."},
    {"number": "03", "title": "Distance Matrix", "description": "Calculate proteome similarity from sum of bit scores between organisms."},
    {"number": "04", "title": "Tree Construction", "description": "Phylogenetic inference using neighbour-joining with bootstrapping."},
]


# FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
def save_email(email):
    """Saves the email and current timestamp to the CSV file."""
    try:
        os.makedirs(os.path.dirname(EMAIL_FILE), exist_ok=True)
        file_exists = os.path.exists(EMAIL_FILE)
        
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

def get_base64_of_bin_file(bin_file):
    """Encodes a binary file (like a video) to base64 for embedding in HTML."""
    try:
        with open(bin_file, 'rb') as f:
            data = f.read()
        return base64.b64encode(data).decode()
    except FileNotFoundError:
        return None

@st.cache_data
def embed_video_background(video_path):
    """Embeds the video and custom CSS for styling the landing page."""
    video_base64 = get_base64_of_bin_file(video_path)
    if video_base64 is None:
        video_url = "" 
    else:
        video_url = f"data:video/mp4;base64,{video_base64}"

    st.markdown(
        f"""
        <style>
        /* BASE SETUP */
        .stApp {{ background: none !important; }}
        .stApp > header, .stApp > div:not(.video-container) {{
            z-index: 10; position: relative; background: none !important; 
        }}
        .video-background {{
            position: fixed; top: 0; left: 0; min-width: 100%; min-height: 100%;
            width: auto; height: auto; z-index: -100; overflow: hidden; object-fit: cover;
        }}
        
        /* THEME COLORS */
        :root {{
            --primary-color: #10b981; 
            --secondary-color: #3b82f6; 
            --foreground-color: #ffffff;
            --muted-foreground-color: #a0a0a0;
            --background-color: #0d0d0d;
        }}
        
        .text-gradient {{
            background: linear-gradient(90deg, var(--primary-color), var(--secondary-color)); 
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            color: transparent; 
        }}

        .glass-panel {{
            background-color: rgba(255, 255, 255, 0.1); 
            backdrop-filter: blur(8px); 
            -webkit-backdrop-filter: blur(8px); 
            border: 1px solid rgba(255, 255, 255, 0.2);
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }}
        
        /* METHODOLOGY STYLES */
        .method-step {{ text-align: center; position: relative; padding: 0 1rem; }}
        .method-number-bubble {{
            width: 4rem; height: 4rem; border-radius: 50%;
            display: flex; align-items: center; justify-content: center; margin: 0 auto 1rem; 
            position: relative; z-index: 10;
        }}
        .method-number {{ font-size: 1.25rem; font-weight: 700; }}
        .method-line {{
            position: absolute; top: 32px; left: 50%; right: 50%; height: 2px;
            background: linear-gradient(to right, var(--primary-color), var(--secondary-color));
        }}
        @media (min-width: 768px) {{
            .method-line {{ left: 8px; right: 8px; }}
            .method-step:not(:first-child) .method-line {{ left: -50%; }}
            .method-step:not(:last-child) .method-line {{ right: -50%; }}
        }}
        </style>
        
        <div class="video-container">
            <video class="video-background" autoplay muted loop>
                <source src="{video_url}" type="video/mp4">
            </video>
        </div>
        """,
        unsafe_allow_html=True
    )

# --- APP EXECUTION ---

# 1. Set Page Config
st.set_page_config(
    page_title="SHE: Structural History of Eukarya",
    page_icon= "data/tree/logo_tree.png", 
    layout="wide"
)

# 2. Add Styles for Dashboard (Applied Globally)
st.markdown("""
<style>
    /* --- NEW DASHBOARD CARD STYLES --- */
    .dashboard-card {
        background-color: #111827; /* Dark Navy/Black */
        border: 1px solid rgba(255, 255, 255, 0.1);
        border-radius: 16px;
        padding: 24px;
        height: 100%;
        display: flex;
        flex-direction: column;
        justify-content: space-between;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.5);
        transition: all 0.3s ease;
    }
    .dashboard-card:hover {
        border-color: #10b981; /* Green hover effect */
        box-shadow: 0 10px 15px -3px rgba(16, 185, 129, 0.2);
    }
    
    .card-header h3 {
        font-size: 1.4rem;
        font-weight: 700;
        color: #ffffff;
        margin-bottom: 0.5rem;
    }
    
    .card-description {
        color: #9ca3af; /* Gray-400 */
        font-size: 0.95rem;
        margin-bottom: 1.5rem;
        line-height: 1.5;
    }
    
    /* Feature Icon Styling */
    .card-icon {
        font-size: 2rem;
        margin-bottom: 1rem;
        display: inline-block;
        background: rgba(255,255,255,0.05);
        padding: 10px;
        border-radius: 12px;
    }

    /* VISUAL: PROGRESS BARS (Model Organisms) */
    .prog-row { display: flex; align-items: center; margin-bottom: 12px; font-size: 0.85rem; }
    .prog-label { width: 120px; color: #e5e7eb; }
    .prog-track { flex-grow: 1; height: 8px; background: #374151; border-radius: 4px; overflow: hidden; margin: 0 12px; }
    .prog-fill { height: 100%; border-radius: 4px; }
    .prog-val { width: 30px; text-align: right; color: #9ca3af; }

    /* VISUAL: HEATMAP (Evolutionary Rate) */
    .heatmap-grid { display: grid; grid-template-columns: repeat(5, 1fr); gap: 6px; margin-top: 1rem; margin-bottom: 1rem; }
    .heatmap-cell { height: 28px; border-radius: 4px; opacity: 0.9; }
    
    /* VISUAL: HISTOGRAM (Homologs) */
    .histo-container { display: flex; align-items: flex-end; gap: 8px; height: 60px; margin-top: 10px; }
    .histo-bar-box { flex: 1; display: flex; flex-direction: column; align-items: center; gap: 4px; }
    .histo-bar { width: 100%; border-radius: 4px 4px 0 0; opacity: 0.9; }
    .histo-label { font-size: 0.7rem; color: #6b7280; }
    .histo-val { font-size: 0.8rem; font-weight: bold; color: #fff; }
</style>
""", unsafe_allow_html=True)


# 3. Main Title & Tree Video
st.markdown("""
    <h1 style='text-align: center; color: green; font-size: 4rem; font-weight: bold;'>
        SHE: Structural History of Eukarya
    </h1>
""", unsafe_allow_html=True)
st.markdown("<h2 style='text-align: center; color: white;'>Beyond traditional phylogenetics — leverage proteome-wide structural comparisons for groundbreaking insights", unsafe_allow_html=True)

st.sidebar.success("Find more analysis options in the panels above.")

try:
    video_tree = open("data/tree/final_tree.mp4", "rb")
    tree_bytes = video_tree.read()
    st.video(tree_bytes, autoplay=True, loop=True, muted=True)
except FileNotFoundError:
    st.warning("Tree video not found (data/tree/final_tree.mp4)")

# 4. Main Content Container
with st.container():

    # --- TOP DASHBOARD SECTION ---
    st.markdown("<br>", unsafe_allow_html=True)

    # Row 1: Model Organisms (Left) & Tree/Evolution (Right)
    dashboard_row1 = st.columns([1, 1], gap="medium")

    # --- CARD 1: MODEL ORGANISM RANKING ---
    with dashboard_row1[0]:
        st.markdown("""
        <div class="dashboard-card">
            <div class="card-header">
                <h3>Model Organism Ranking</h3>
                <p class="card-description">
                    Select your proteins of interest and discover the optimal model organism 
                    where they are most structurally conserved.
                </p>
                <div style="background: #1f2937; padding: 15px; border-radius: 12px; margin-bottom: 15px;">
                    <div style="color: #6b7280; font-size: 0.75rem; margin-bottom: 10px;">Conservation Score</div>
                    <div class="prog-row">
                        <div class="prog-label">D. melanogaster</div>
                        <div class="prog-track"><div class="prog-fill" style="width: 94%; background: #38bdf8;"></div></div>
                        <div class="prog-val">94%</div>
                    </div>
                    <div class="prog-row">
                        <div class="prog-label">C. elegans</div>
                        <div class="prog-track"><div class="prog-fill" style="width: 87%; background: #38bdf8;"></div></div>
                        <div class="prog-val">87%</div>
                    </div>
                    <div class="prog-row">
                        <div class="prog-label">S. cerevisiae</div>
                        <div class="prog-track"><div class="prog-fill" style="width: 82%; background: #38bdf8;"></div></div>
                        <div class="prog-val">82%</div>
                    </div>
                    <div class="prog-row">
                        <div class="prog-label">A. thaliana</div>
                        <div class="prog-track"><div class="prog-fill" style="width: 71%; background: #4ade80;"></div></div>
                        <div class="prog-val">71%</div>
                    </div>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        st.page_link("pages/2_Model_organisms.py", label="Find Best Model Organism", icon="🐁", use_container_width=True)

    # --- CARD 2: TREE & EVOLUTIONARY RATE ---
    with dashboard_row1[1]:
        st.markdown("""
        <div class="dashboard-card">
            <div class="card-header">
                <h3>Tree & Evolutionary Rate</h3>
                <p class="card-description">
                    Identify fast and slow-evolving lineages by visualizing evolutionary rates across the tree.
                </p>
                <div style="background: #1f2937; padding: 15px; border-radius: 12px; margin-bottom: 15px;">
                   <div style="display:flex; justify-content:space-between; margin-bottom:5px;">
                        <span style="color: #6b7280; font-size: 0.75rem;">Evolutionary Rate</span>
                        <span style="color: #6b7280; font-size: 0.75rem;">Slow <span style="display:inline-block; width:40px; height:4px; background:linear-gradient(90deg, #4ade80, #ef4444); vertical-align:middle; margin:0 5px;"></span> Fast</span>
                   </div>
                   <div class="heatmap-grid">
                        <div class="heatmap-cell" style="background:#4ade80;"></div> <div class="heatmap-cell" style="background:#4ade80;"></div>
                        <div class="heatmap-cell" style="background:#a3e635;"></div> <div class="heatmap-cell" style="background:#facc15;"></div> <div class="heatmap-cell" style="background:#ea580c;"></div> <div class="heatmap-cell" style="background:#4ade80;"></div>
                        <div class="heatmap-cell" style="background:#a3e635;"></div>
                        <div class="heatmap-cell" style="background:#facc15;"></div>
                        <div class="heatmap-cell" style="background:#facc15;"></div>
                        <div class="heatmap-cell" style="background:#4ade80;"></div>
                        <div class="heatmap-cell" style="background:#4ade80;"></div>
                        <div class="heatmap-cell" style="background:#ea580c;"></div>
                        <div class="heatmap-cell" style="background:#ef4444;"></div> <div class="heatmap-cell" style="background:#ef4444;"></div>
                        <div class="heatmap-cell" style="background:#ea580c;"></div>
                   </div>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        st.page_link("pages/1_Explore_Tree.py", label="Explore Evolutionary Rates", icon="🌳", use_container_width=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # --- CARD 3: HOMOLOGS (Full Width) ---
    st.markdown("""
    <div class="dashboard-card" style="flex-direction: row; display: block;">
        <div style="display: flex; flex-wrap: wrap; align-items: center;">
            <div style="flex: 1.5; padding-right: 20px; min-width: 300px;">
                <h3 style="color:white; margin:0; font-size: 1.5rem; font-weight: 700;">Protein Homolog Discovery</h3>
                <p style="color:#9ca3af; margin-top:10px; margin-bottom: 20px;">
                    Enter any UniProt ID to instantly discover structural homologs across all 24.3 million proteins.
                    Map the evolutionary distribution of your protein family.
                </p>
            </div>
            <div style="flex: 1; min-width: 250px;">
                <div style="background: #1f2937; padding: 15px; border-radius: 12px;">
                    <div style="color: #6b7280; font-size: 0.75rem;">Structural Homolog Distribution</div>
                    <div class="histo-container">
                        <div class="histo-bar-box">
                            <div class="histo-val">342</div>
                            <div class="histo-bar" style="height: 100%; background: #38bdf8;"></div>
                            <div class="histo-label">Opis</div>
                        </div>
                        <div class="histo-bar-box">
                            <div class="histo-val">156</div>
                            <div class="histo-bar" style="height: 45%; background: #4ade80;"></div>
                            <div class="histo-label">Viri</div>
                        </div>
                        <div class="histo-bar-box">
                            <div class="histo-val">89</div>
                            <div class="histo-bar" style="height: 25%; background: #facc15;"></div>
                            <div class="histo-label">SAR</div>
                        </div>
                        <div class="histo-bar-box">
                            <div class="histo-val">23</div>
                            <div class="histo-bar" style="height: 10%; background: #c084fc;"></div>
                            <div class="histo-label">Amoe</div>
                        </div>
                         <div class="histo-bar-box">
                            <div class="histo-val">12</div>
                            <div class="histo-bar" style="height: 5%; background: #fb923c;"></div>
                            <div class="histo-label">Disc</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    st.page_link("pages/3_Downloads.py", label="Search Homologs & Download Data", icon="🔎", use_container_width=True)

    st.markdown("<br><br>", unsafe_allow_html=True) 

    # --- A. STATISTICS TILES SECTION ---
    css = """
        .green-tile-container {
            background: rgba(16, 185, 129, 0.1); 
            border: 1px solid rgba(16, 185, 129, 0.2);
            min-height: 100px; 
            padding: 15px;
            border-radius: 12px; 
            display: flex;
            flex-direction: column;
            justify-content: center;
            text-align: center;
        }
        .green-tile-container h2 { margin: 0; color: #10b981; font-weight: 800; }
        .green-tile-container p { margin: 0; color: #a0a0a0; font-size: 0.9rem; }
        """
    st.html(f"<style>{css}</style>")

    row1 = st.columns(4)
    titles = ["Organisms", "Proteins", "Kingdoms", "Powered By"]
    values = ["1542", "24.3M", "6", "FoldSeek"]

    for i, (col, value, title) in enumerate(zip(row1, values, titles)):
        with col:
            st.markdown(
                f"""
                <div class="green-tile-container">
                    <h2>{value}</h2>
                    <p>{title}</p>
                </div>
                """,
                unsafe_allow_html=True
            )

    # --- B. FEATURES GRID SECTION (Updated to Match Dashboard Style) ---
    st.write("---") 
    st.markdown("""
        <div style="text-align: center; max-width: 42rem; margin: 0 auto 3rem;">
            <h2 style="font-size: 2rem; font-weight: 700; color: green; margin-bottom: 1rem;">
                Discover What Makes SHE Different
            </h2>
            <p style="color: var(--muted-foreground-color);">
                Explore the complete structural phylogeny of eukaryotes with our interactive analysis platform.
            </p>
        </div>
    """, unsafe_allow_html=True)

    # Features Grid using the new 'dashboard-card' style
    cols = st.columns(3)
    
    for i, feature in enumerate(FEATURES):
        col_index = i % 3
        with cols[col_index]:
            # Use 'dashboard-card' here instead of 'glass-panel'
            st.markdown(f"""
                <div class="dashboard-card" style="min-height: 280px; margin-bottom: 20px;">
                    <div class="card-header">
                        <div class="card-icon">{feature["icon"]}</div>
                        <h3>{feature["title"]}</h3>
                        <p class="card-description">{feature["description"]}</p>
                    </div>
                </div>
            """, unsafe_allow_html=True)
            
            # Place the link directly under the visual card
            st.page_link(feature["link_path"], label=f"Go to: {feature['title']}", icon="🔗")
        
            with st.popover(f"▶️ Watch Demo", use_container_width=True):
                st.markdown(f"**Demo for {feature['title']}**")
                if "demo" in feature:
                    try:
                        with open(feature["demo"], "rb") as video_file:
                            video_bytes = video_file.read()
                        st.video(video_bytes, autoplay=True, loop=True, muted=True)
                    except FileNotFoundError:
                        st.warning(f"Video file not found: {feature.get('demo')}")
                else:
                    st.warning("No video demo available for this feature yet.")
                st.markdown("<br>", unsafe_allow_html=True)

    st.markdown("<hr style='margin: 4rem 0;'>", unsafe_allow_html=True) 


    # --- C. METHODOLOGY SECTION ---
    st.markdown("""
        <div style="text-align: center; max-width: 42rem; margin: 0 auto 3rem;">
            <h2 style="font-size: 2rem; font-weight: 700; color: green; margin-bottom: 1rem;">
                How It Works
            </h2>
            <p style="color: var(--muted-foreground-color);">
                Our novel methodology compares entire proteomes based on 3D protein structure rather than sequence alone.
            </p>
        </div>
    """, unsafe_allow_html=True)
    
    step_cols = st.columns(4)
    for i, step in enumerate(STEPS):
        with step_cols[i]:
            st.markdown(f"""
                <div class="method-step">
                    <div class="glass-panel method-number-bubble">
                        <span class="method-number text-gradient">{step["number"]}</span>
                    </div>
                    <div class="step-content">
                        <h3 class="step-title">{step["title"]}</h3>
                        <p style="font-size: 0.875rem; color: var(--muted-foreground-color);">{step["description"]}</p>
                    </div>
                </div>
            """, unsafe_allow_html=True)


    st.markdown("<hr style='margin: 4rem 0;'>", unsafe_allow_html=True) 

    # --- D. CITATION AND UPDATES ---
    st.subheader('Cite', divider="rainbow")
    st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')

    st.sidebar.subheader('Updates', divider="rainbow")
    st.sidebar.write('Sign up for updates!')
    email = st.sidebar.text_input('Your email', key="email_input") 
    if st.sidebar.button('Subscribe'):
        user_email = st.session_state.email_input.strip()
        if not user_email:
            st.error("Please enter an email address.")
        else:
            if save_email(user_email):
                st.balloons()
                st.sidebar.success(f"Thank you! **{user_email}** is where you will receive updates regarding SHE!")