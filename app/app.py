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
# NOTE: Placeholder path. Replace with your actual video file path if it exists.
VIDEO_PATH_BG = "data/tree/bg.mp4" 

# Define navigation items for the custom header
NAV_ITEMS = [
    {"label": "Home", "icon": "🏠", "path": "app.py"},
    {"label": "Tree Viewer", "icon": "🌳", "path": "pages/1_Explore_Tree.py"},
    {"label": "Model Organisms", "icon": "🐁", "path": "pages/2_Model_organisms.py"},
    {"label": "Downloads", "icon": "⬇️", "path": "pages/3_Downloads.py"},
]

# Define features for the main grid section
FEATURES = [
    {"icon": "🌳", "title": "Interactive Phylogeny", "description": "Explore the structural phylogenetic tree for the main subkingdoms including Opisthokonta, Viridiplantae, SAR, Amoebozoa, Discoba, and Rhodophyta.", "link_path": "pages/1_Explore_Tree.py"},
    {"icon": "🔬", "title": "Evolutionary Rate", "description": "Visualize colored branches by evolutionary rate to identify fast and slow-evolving proteomes.", "link_path": "pages/1_Explore_Tree.py"},
    {"icon": "🧬", "title": "Taxon Filtering", "description": "Zoom on specific taxa or species to discover new evolutionary relationships.", "link_path": "pages/1_Explore_Tree.py"},
    {"icon": "🔎", "title": "Protein Search", "description": "Input UniProt IDs to analyze presence/absence of structural homologs across all eukaryotic organisms.", "link_path": "pages/1_Explore_Tree.py"}, # Reusing Model Organisms page for search
    {"icon": "🐁", "title": "Model Organism Ranking", "description": "Select proteins of interest to find the best model organisms where they are most conserved.", "link_path": "pages/2_Model_organisms.py"},
    {"icon": "⬇️", "title": "Data Downloads", "description": "Download FoldSeek comparison bitscores, Newick files, and predicted protein structures from each proteome.", "link_path": "pages/3_Downloads.py"},
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
    """Saves the email and current timestamp to the CSV file. 
    It ensures the directory exists before attempting to write."""
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
        # Don't fail the app if the video isn't there, just show a warning
        st.warning(f"Warning: Background video file not found at path: {video_path}. Using plain background.")
        video_url = "" # Set to empty if file is missing
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
        
        /* 1. CUSTOM THEME COLORS & EFFECTS */
        :root {{
            --primary-color: #10b981; /* Tailwind green-500 */
            --secondary-color: #3b82f6; /* Tailwind blue-500 */
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
        
        /* 2. FIXED HEADER STYLES */
        .fixed-header {{
            position: fixed; top: 0; left: 0; right: 0; z-index: 999;
            height: 64px; display: flex; align-items: center; justify-content: space-between;
            padding: 0 1rem;
        }}
        .logo-container {{ display: flex; align-items: center; gap: 0.75rem; }}
        .logo-title {{ font-size: 1.125rem; font-weight: 700; line-height: 1; }}
        .logo-subtitle {{ font-size: 0.625rem; color: var(--muted-foreground-color); margin-top: -0.25rem; line-height: 1; }}
        .nav-desktop {{ display: none; }}
        @media (min-width: 768px) {{ .nav-desktop {{ display: flex; align-items: center; gap: 0.25rem; }} }}
        .nav-button, .nav-button-active {{
            display: flex; align-items: center; gap: 0.5rem; padding: 0.25rem 0.75rem; border-radius: 0.375rem; 
            font-size: 0.875rem; font-weight: 500; cursor: pointer; transition: background-color 0.15s;
            text-decoration: none; color: #e0e0e0;
        }}
        .nav-button:hover {{ background-color: rgba(255, 255, 255, 0.1); color: var(--foreground-color); }}
        .nav-button-active {{ background-color: rgba(255, 255, 255, 0.2); color: var(--primary-color); font-weight: 600; }}
        
        /* 3. MAIN CONTENT OFFSET & HERO STYLES */
        .main-content-offset {{ margin-top: 64px; padding-top: 5rem; }}
        .hero-title {{ 
            font-size: 4rem; line-height: 1; font-weight: 800; text-align: center; margin-bottom: 2rem;
            color: var(--foreground-color);
        }}
        @media (min-width: 768px) {{ .hero-title {{ font-size: 5rem; }} }}
        .hero-subtitle {{ color: var(--muted-foreground-color); text-align: center; max-width: 48rem; margin: 0 auto 3rem; font-size: 1.125rem; }}

        .hero-stats-box {{
            padding: 1rem; border-radius: 0.75rem; text-align: center;
        }}
        .stat-value {{ font-size: 1.5rem; font-weight: 700; }}
        
        /* 4. FEATURE GRID STYLES */
        .feature-grid {{
            display: grid; grid-template-columns: repeat(1, minmax(0, 1fr)); gap: 1.5rem;
        }}
        @media (min-width: 768px) {{ .feature-grid {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }} }}
        @media (min-width: 1024px) {{ .feature-grid {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }} }}
        
        .feature-tile {{
            border-radius: 1rem; padding: 1.5rem; min-height: 220px;
            display: flex; flex-direction: column; 
            transition: all 0.3s ease;
        }}
        .feature-tile:hover {{
            border-color: var(--primary-color);
            box-shadow: 0 0 10px rgba(16, 185, 129, 0.3); /* subtle glow on hover */
        }}
        .feature-icon-box {{
            width: 3rem; height: 3rem; border-radius: 0.75rem; 
            display: flex; align-items: center; justify-content: center; margin-bottom: 1rem;
            background-color: rgba(16, 185, 129, 0.2); /* primary/10 */
        }}
        .feature-title {{ 
            font-size: 1.125rem; font-weight: 600; color: var(--foreground-color); margin-bottom: 0.5rem;
        }}
        .feature-description {{
            font-size: 0.875rem; color: var(--muted-foreground-color); line-height: 1.5;
        }}
        
        /* 5. METHODOLOGY STYLES */
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
        .step-title {{ font-size: 1.125rem; font-weight: 600; color: var(--foreground-color); margin-bottom: 0.5rem; }}

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

# --- APP EXECUTION ---

# 1. Set Page Config
st.set_page_config(
    page_title="SHE: Structural History of Eukarya",
    page_icon= "data/tree/logo_tree.png", 
    layout="wide" # Use wide layout for the rich landing page design
)

# 2. Embed Background Video & Custom CSS
#embed_video_background(VIDEO_PATH_BG)

st.title('SHE: Structural History of Eukarya')

st.sidebar.success("Find more analysis options in the panels above.")

# 3. Custom Fixed Header (replicated from previous response)
def create_custom_header(nav_items):
    nav_html_items = ""
    for item in nav_items:
        is_active = item["label"] == "Home"
        button_class = "nav-button-active" if is_active else "nav-button"
        # The path links won't navigate correctly in the fixed header, but provide visual structure
        nav_html_items += f"""
        <a href="#" class="{button_class}">
            {item["icon"]} {item["label"]}
        </a>
        """

    st.markdown(
        f"""
        <header class="glass-panel fixed-header">
            <div class="logo-container">
                <span style="color:var(--primary-color); font-size: 32px;">🌳</span>
                <div class="logo-text">
                    <span class="logo-title text-gradient">SHE</span>
                    <span class="logo-subtitle">Structural History of Eukarya</span>
                </div>
            </div>
            <nav class="nav-desktop">
                {nav_html_items}
            </nav>
        </header>
        """,
        unsafe_allow_html=True
    )

#create_custom_header(NAV_ITEMS)
video_tree = open("data/tree/final_tree.mp4", "rb")
tree_bytes = video_tree.read()
st.video(tree_bytes, autoplay=True, loop=True, muted=True)
# 4. Main Content Container (starts below the fixed header)
with st.container():

    # CTA Buttons (using Streamlit components for functionality)
    cta_col1, cta_col2, cta_col3 = st.columns([1, 1, 1])
    with cta_col1:
        st.page_link("pages/1_Explore_Tree.py", label="🌳 Explore the Tree", icon="👉", help="explore the **structural phylogenetic tree** for all Eukarya, a kingdom, a taxon, a species or a set of proteins")
    with cta_col2:
        st.page_link("pages/2_Model_organisms.py", label="🐁 Model Organisms search", icon="👉", help="find the **best model organism** for your study based on the structural comparison of the protein(s) in your pathway of interest.")
    with cta_col3:
        st.page_link("pages/3_Downloads.py", label="⬇️ Download Data", icon="👉", help="Access all raw alignment scores and files.")

    st.markdown("<br><br>", unsafe_allow_html=True) # Spacer
    # --- A. STATISTICS TILES SECTION ---

    # Define your custom CSS with the GREEN color
    css = """
        .blue-tile-container {
            /* ⭐ CHANGE MADE HERE: Using a light green RGBA value */
            background: rgba(100, 200, 100, 0.3); 
            min-height: 120px; 
            padding: 10px;
            border-radius: 8px; 
            
            display: flex;
            flex-direction: column;
            justify-content: center;
        }
        
        /* Ensure markdown header (##) spacing is minimal inside the tile */
        .blue-tile-container h2 {
            margin: 0; 
            padding: 0;
            line-height: 1.2;
        }
        
        /* Style the paragraph text (the title) */
        .blue-tile-container p {
            margin: 0;
            padding: 0;
            font-size: 16px; 
            opacity: 0.8;
        }
        """
    st.html(f"<style>{css}</style>")

    # --- Application Logic ---

    row1 = st.columns(4)
    titles = ["Organisms", "Proteins", "Kingdoms", "Powered By"]
    values = ["1542", "24.3M", "6", "FoldSeek"]

    for i, (col, value, title) in enumerate(zip(row1, values, titles)):
        with col:
            unique_key = f"my_green_container_{i}" # Changed key name for clarity
            
            with st.container(key=unique_key):
                
                # Single st.markdown call for the entire tile.
                st.markdown(
                    f"""
                    <div class="blue-tile-container">
                        <h2>{value}</h2>
                        <p>{title}</p>
                    </div>
                    """,
                    unsafe_allow_html=True
                )


    # --- B. FEATURES SECTION ---

    st.markdown("""
        <div style="text-align: center; max-width: 42rem; margin: 0 auto 3rem;">
            <h2 style="font-size: 2rem; font-weight: 700; color: var(--foreground-color); margin-bottom: 1rem;">
                Powerful Visualization Tools
            </h2>
            <p style="color: var(--muted-foreground-color);">
                Explore the complete structural phylogeny of eukaryotes with our interactive analysis platform.
            </p>
        </div>
    """, unsafe_allow_html=True)

    # Features Grid (Using custom HTML/CSS for tiles, but st.page_link for navigation below the tile)
    
    # We use a wrapper to use the custom grid CSS class
    st.markdown('<div class="feature-grid">', unsafe_allow_html=True)
    
    # Use st.columns inside the markdown wrapper for the content flow
    # Since we have 6 features, we can use 3 columns and iterate through them twice
    # This is a Streamlit-specific compromise for grid layout
    
    # Create the columns outside the loop
    cols = st.columns(3)
    
    for i, feature in enumerate(FEATURES):
        col_index = i % 3  # Distribute tiles across the 3 columns
        with cols[col_index]:
            st.markdown(f"""
                <div class="glass-panel feature-tile">
                    <div class="feature-icon-box" style="color: var(--primary-color);">
                        {feature["icon"]}
                    </div>
                    <h3 class="feature-title">{feature["title"]}</h3>
                    <p class="feature-description">{feature["description"]}</p>
                </div>
            """, unsafe_allow_html=True)
            
            # CORRECTED LINE: Using single quotes for 'title' dictionary access
            st.page_link(feature["link_path"], label=f"Go to {feature['title'].split()[0]} Page", icon="🔗")
            st.markdown("<br>", unsafe_allow_html=True) # Small space between tiles
            
    st.markdown('</div>', unsafe_allow_html=True)
    
    st.markdown("<hr style='margin: 4rem 0;'>", unsafe_allow_html=True) # Separator


    # --- C. METHODOLOGY SECTION ---

    st.markdown("""
        <div style="text-align: center; max-width: 42rem; margin: 0 auto 3rem;">
            <h2 style="font-size: 2rem; font-weight: 700; color: var(--foreground-color); margin-bottom: 1rem;">
                How It Works
            </h2>
            <p style="color: var(--muted-foreground-color);">
                Our novel methodology compares entire proteomes based on 3D protein structure rather than sequence alone.
            </p>
        </div>
    """, unsafe_allow_html=True)
    
    # Methodology Steps (Using 4 columns)
    step_cols = st.columns(4)
    
    for i, step in enumerate(STEPS):
        with step_cols[i]:
            # The line requires complex CSS positioning. We will use the number bubble and a centered layout.
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


    st.markdown("<hr style='margin: 4rem 0;'>", unsafe_allow_html=True) # Separator


    # --- D. CITATION AND UPDATES (Original content, slightly restyled) ---

    st.subheader('Cite', divider="rainbow")
    st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')

    st.sidebar.subheader('Updates', divider="rainbow")
    st.sidebar.write('Sign up for updates!')
    # Use a unique key for the text input
    email = st.sidebar.text_input('Your email', key="email_input") 
    if st.sidebar.button('Subscribe'):
        # Correctly reference the value using st.session_state
        user_email = st.session_state.email_input.strip()
        if not user_email:
            # Only check remaining: must not be empty
            st.error("Please enter an email address.")
        else:
            # If we reach here, the email is non-empty and we save it.
            if save_email(user_email):
                st.balloons()
                # Updated success message
                st.sidebar.success(f"Thank you! **{user_email}** is where you will receive updates regarding SHE!")