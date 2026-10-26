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
import re

# --- CONFIGURATION & CONSTANTS ---
EMAIL_FILE = 'data/email_updates.csv'
VIDEO_PATH_BG = "data/tree/bg.mp4" 

# Define navigation items for the custom header
NAV_ITEMS = [
    {"label": "SHE", "icon": "🏠", "path": "app.py"},
    {"label": "Model Organism Ranking", "icon": "🐁", "path": "pages/1_Model_Organism_Ranking.py"},
    {"label": "Explore SHE tree", "icon": "🌳", "path": "pages/2_Explore_Tree.py"},
    {"label": "Interactive SHE Tree", "icon": "🔎", "path": "pages/3_Interactive_Tree.py"},
    {"label": "Downloads", "icon": ":material/download:", "path": "pages/4_Downloads.py"},
]

# Define features for the main grid section
FEATURES = [
    {"icon": "🌳", "title": "Structural Phylogeny", "description": "Explore the structural phylogenetic tree for the main eukaryotic clades.", "link_path": "pages/2_Explore_Tree.py", "demo": "data/demos/she.mp4"}, #"demo": "data/demos/demo1.mp4"
    {"icon": "📈", "title": "Evolutionary Rate", "description": "Visualize colored branches by evolutionary rate to identify fast and slow-evolving proteomes.", "link_path": "pages/2_Explore_Tree.py", "demo": "data/demos/evorate.mp4"},
    {"icon": "🍄", "title": "Taxon Interactive Search", "description": "Zoom on specific taxa or species to discover new evolutionary relationships.", "link_path": "pages/3_Interactive_Tree.py", "demo": "data/demos/taxon.mp4"},
    {"icon": "🔎", "title": "Structural Homologs Search", "description": "Analyze the presence/absence of protein structural orthologs across all eukaryotic organisms.", "link_path": "pages/3_Interactive_Tree.py", "demo": "data/demos/homolog.mp4"}, # Reusing Model Organisms page for search
    {"icon": "🐁", "title": "Model Organism Ranking", "description": "Find the best model organism match to study your biological pathway with.", "link_path": "pages/1_Model_Organism_Ranking.py", "demo": "data/demos/model.mp4"},
    {"icon": "📥", "title": "Data Downloads", "description": "Download the trees, the structural alignment scores and the predicted protein structures.", "link_path": "pages/4_Downloads.py", "demo": "data/demos/demo6.mp4"},
]

# Define methodology steps
STEPS = [
    {"number": "01", "title": "Proteome Collection", "description": "1,542 complete eukaryotic proteomes from UniProt with predicted protein structures."},
    {"number": "02", "title": "Structural Alignment", "description": "Foldseek all-vs-all protein structure comparison at the proteome level - 296,000B alignments."},
    {"number": "03", "title": "Distance Matrix", "description": " Calculate proteome similarity from sum of bit scores between organisms."},
    {"number": "04", "title": "Tree Construction", "description": "Phylogenetic inference using neighbour-joining with bootstrapping."},
]


# FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
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

#@st.cache_data
@st.cache_data(show_spinner=False)
def get_base64_of_bin_file(bin_file):
    """Encodes a binary file (like a video) to base64 for embedding in HTML."""
    try:
        with open(bin_file, 'rb') as f:
            data = f.read()
        return base64.b64encode(data).decode()
    except FileNotFoundError:
        return None

@st.cache_data(show_spinner=False)
def get_video_bytes(video_path):
    """
    Reads a video file into memory. 
    Cached: Prevents re-reading 6 video files on every page refresh.
    """
    if not os.path.exists(video_path):
        return None
    with open(video_path, "rb") as f:
        return f.read()

#@st.cache_data
@st.cache_data(show_spinner=False)
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

#@st.cache_data
@st.cache_data(show_spinner=False)
def get_base64_image(image_path):
    """Converts a local file to a base64 string for HTML embedding."""
    if not os.path.exists(image_path):
        return "" # Return empty if file not found to prevent errors
    with open(image_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode()

#@st.cache_data
@st.fragment
def set_background(image_path):
    """
    Adds a background image to the main Streamlit app container.
    """
    if not os.path.exists(image_path):
        st.error(f"Background image not found at: {image_path}")
        return
    with open(image_path, "rb") as f:
        data = f.read()
    # Convert image to Base64
    #encoded_image = base64.b64encode(data).decode()
    encoded_image = get_base64_image(image_path)
    # Inject CSS with the Base64 image
    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url("data:image/png;base64,{encoded_image}");
            background-size: cover;        /* Scale image to cover entire screen */
            background-position: center;   /* Center the image */
            background-repeat: no-repeat;  /* Do not repeat */
            background-attachment: fixed;  /* Keep background fixed when scrolling */
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

# --- APP EXECUTION ---
#sidebar navigation
with st.sidebar:
    st.title("Browse")
    for item in NAV_ITEMS:
        st.page_link(item["path"], label=item["label"], icon=item["icon"])
st.sidebar.divider()
st.markdown("""
<style>
    /* Force all standard Streamlit text to use your readable style */
    .stMarkdown p, .stMarkdown li {
        color: #ffffff !important;
        font-weight: 600 !important;
        text-shadow: 2px 2px 0 #000000, 0 0 10px rgba(0,0,0,0.8);
        line-height: 1.6;
        font-size: 1.1rem;
    }
</style>
""", unsafe_allow_html=True)
# 1. Set Page Config
st.set_page_config(
    page_title="SHE: Structural History of Eukarya",
    page_icon= "data/tree/logo_tree.png", 
    layout="wide"
)
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

# 2. Add Styles for Dashboard (Applied Globally)
st.markdown("""
<style>
    /* --- DASHBOARD CARD STYLES --- */
    .dashboard-card {
        background-color: #111827;
        border: 1px solid rgba(255, 255, 255, 0.1);
        border-radius: 16px;
        padding: 24px;
        
        /* Forces cards to match height */
        min-height: 500px; 
        
        display: flex;
        flex-direction: column;
        justify-content: space-between;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.5);
        transition: all 0.3s ease;
    }
    
    /* Hover effect now applies when hovering the LINK too */
    a:hover .dashboard-card {
        border-color: #10b981;
        box-shadow: 0 10px 15px -3px rgba(16, 185, 129, 0.2);
    }
    
    .card-header h3 {
        font-size: 1.4rem;
        font-weight: 700;
        color: #ffffff;
        margin-bottom: 0.5rem;
    }
    
    .card-description {
        color: #9ca3af;
        font-size: 0.95rem;
        margin-bottom: 1.5rem;
        line-height: 1.5;
        min-height: 80px; 
        flex-grow: 1; 
    }

    /* --- IMAGE STYLING --- */
    .card-visual-img {
        width: 100%;
        height: 220px;       
        object-fit: cover;   
        object-position: center;
        border-radius: 8px;
        border: 1px solid rgba(255,255,255,0.1);
        opacity: 0.9;
        display: block;
    }
    
    /* --- LINK STYLING  (replace heatmap)--- */
    /* Removes the default blue underline from the <a> tag */
    a.card-link {
        text-decoration: none;
        color: inherit;
        display: block;
        height: 100%; /* Ensures the click area fills the column */
    }
</style>
""", unsafe_allow_html=True)


# 3. Main Title & Tree Video
st.sidebar.success("Find more analysis options in the panels above.")
#video_tree = open("data/tree/final_tree.mp4", "rb")
#tree_bytes = video_tree.read()
#st.video(tree_bytes, autoplay=True, loop=True, muted=True)

# 1. Define image path
image_path = 'data/logos/interactive.png'
set_background(image_path)
if os.path.exists(image_path):
    img_b64 = get_base64_image(image_path)
# 3. Render the Banner with Overlay Text
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
            <h1 style='
                font-size: 5rem; 
                font-weight: 800; 
                margin-bottom: 1rem;
                line-height: 1.1;
                position: relative;
                z-index: 10;
                text-shadow: none !important;
                background: linear-gradient(to right, #00FF7F, #e9f5ef); 
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                filter: brightness(1.2) drop-shadow(4px 4px 0px rgba(0,0,0,1));
                '>
                Structural History of Eukarya
            </h1>
            <h3>(SHE)</h3>
            <h2>Beyond traditional phylogenetics</h2>
            <h3>Leverage proteome-wide structural comparisons for groundbreaking insights</h3>
        </div>
    </div>
    """,
    unsafe_allow_html=True
)
##40E0D0
img_ranking = get_base64_image("data/subset/figure_model/ranking2.png") 
img_evolution = get_base64_image("data/logos/evorate.png")
img_homologs = get_base64_image("data/logos/homolog.png")

with st.container():
    st.markdown("<br>", unsafe_allow_html=True)

    # Row 1: Model Organisms (Left) & Tree/Evolution (Right)
    dashboard_row1 = st.columns([1, 1], gap="medium")

    # --- CARD 1 ---
    with dashboard_row1[0]:
        # Wrapped in <a> tag pointing to the URL slug
        st.markdown(f"""
        <a href="Model_Organism_Ranking" target="_self" class="card-link">
            <div class="dashboard-card">
                <div class="card-header">
                    <h3>🐁 Model Organism Ranking</h3>
                    <p class="card-description">
                        Find the best model organism for your study by combined structural matching for all proteins in your pathway.
                    </p>
                </div>
                <div style="background: #1f2937; padding: 10px; border-radius: 12px;">
                    <img src="data:image/png;base64,{img_ranking}" class="card-visual-img">
                </div>
            </div>
        </a>
        """, unsafe_allow_html=True)
        # Keep the st.page_link as a secondary navigation option
        #st.page_link("pages/1_Model_Organism_Ranking.py", label="Find Best Model Organism", icon="🐁", width=True)

    # --- CARD 2 ---
    with dashboard_row1[1]:
        st.markdown(f"""
        <a href="Explore_Tree" target="_self" class="card-link">
            <div class="dashboard-card">
                <div class="card-header">
                    <h3>🌳 SHE Tree & Evolutionary Rates</h3>
                    <p class="card-description">
                        Explore the structural phylogeny of Eukarya and identify fast/slow-evolving proteomes by visualizing the evolutionary rates across the tree.
                    </p>
                </div>
                <div style="background: #1f2937; padding: 10px; border-radius: 12px;">
                   <img src="data:image/png;base64,{img_evolution}" class="card-visual-img">
                </div>
            </div>
        </a>
        """, unsafe_allow_html=True)
        #st.page_link("pages/2_Explore_Tree.py", label="Explore Evolutionary Rates", icon="📈", width=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # --- CARD 3 (Horizontal Layout) ---
    st.markdown(f"""
    <a href="Interactive_Tree" target="_self" class="card-link">
        <div class="dashboard-card" style="flex-direction: row; display: flex; align-items: center; min-height: auto;">
            <div style="flex: 1.5; padding-right: 20px; min-width: 300px;">
                <h3 style="color:white; margin:0; font-size: 1.5rem; font-weight: 700;">🔎 Interactive SHE tree</h3>
                <p style="color:#9ca3af; margin-top:10px; margin-bottom: 20px;">
                    Zoom on specific taxa and discover protein structural homologs across all 24.3 million proteins.
                    Map the evolutionary distribution of your protein family.
                </p>
            </div>
            <div style="flex: 1; min-width: 250px;">
                <div style="background: #1f2937; padding: 10px; border-radius: 12px;">
                    <img src="data:image/png;base64,{img_homologs}" class="card-visual-img" style="height: 180px;">
                </div>
            </div>
        </div>
    </a>
    """, unsafe_allow_html=True)
    #st.page_link("pages/3_Interactive_Tree.py", label="Search Homologs", icon="🔎", width=True)
    st.markdown("<br>", unsafe_allow_html=True)
    # --- A. STATISTICS TILES SECTION ---
    css = """
        /* 1. Base Style */
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
            
            /* Smooth animation for the hover effect */
            transition: all 0.3s ease-in-out; 
        }
        
        .green-tile-container h2 { margin: 0; color: #10b981; font-weight: 800; }
        .green-tile-container p { margin: 0; color: #a0a0a0; font-size: 0.9rem; }

        /* 2. HOVER STATE: The "Green Aura" */
        .green-tile-container:hover {
            border-color: #22c55e;   /* Brighter green border */
            
            /* The Aura Effect: Green Glow */
            box-shadow: 0 0 20px rgba(74, 222, 128, 0.4), 
                        0 0 10px rgba(74, 222, 128, 0.2);
            
            /* Optional: Slight pop effect */
            transform: scale(1.03);
        }
        """
    st.html(f"<style>{css}</style>")



    row1 = st.columns(4)
    titles = ["Organisms", "Proteins", "Structure Alignments", "Powered By"]
    values = ["1,541", "24.3M", "296T", "Foldseek"]

    for i, (col, value, title) in enumerate(zip(row1, values, titles)):
        with col:
            if title == "Powered By":
                logo_path = "data/logos/foldseek_logo.png" 
                logo_b64 = get_base64_image(logo_path)
                # --- SPECIAL FOLDSEEK CARD ---
                st.markdown(
                    f"""
                    <a href="https://github.com/steineggerlab/foldseek" target="_blank" style="text-decoration: none; color: inherit;">
                        <div class="green-tile-container" style="
                            background-image: url('data:image/png;base64,{logo_b64}'); 
                            background-size: cover; 
                            background-position: center;
                            background-blend-mode: overlay;
                            background-color: rgba(0, 50, 0, 0.7); /* Dark green tint to make text readable */
                            cursor: pointer;
                            ">
                            <h2 style="color: white; text-shadow: 2px 2px 4px #000;">{value}</h2>
                            <p style="color: white; font-weight: bold; text-shadow: 1px 1px 2px #000;">{title}</p>
                        </div>
                    </a>
                    """,
                    unsafe_allow_html=True
                )
            else:
                # --- STANDARD CARDS (Semi-Transparent Green) ---
                st.markdown(
                    f"""
                    <div class="green-tile-container" style="
                        background-color: rgba(0, 50, 0, 0.80); 
                        color: white;
                        backdrop-filter: blur(4px); /* Optional: blurs the background behind the card */
                        ">
                        <h2 style="color: white; margin: 0; text-shadow: 1px 1px 2px black;">{value}</h2>
                        <p style="color: rgba(255, 255, 255, 0.9); margin: 0; text-shadow: 1px 1px 2px black;">{title}</p>
                    </div>
                    """,
                    unsafe_allow_html=True
                )

    # --- B. FEATURES GRID SECTION (Updated to Match Dashboard Style) ---
    st.write("---") 
    st.markdown("""
    <div style="text-align: center; max-width: 42rem; margin: 0 auto 3rem;">
        <h2 style='
            font-size: 2rem; 
            font-weight: 800; 
            margin-bottom: 1rem;
            /* Gradient text needs a filter drop-shadow, not text-shadow */
            background: linear-gradient(to right, #00FF7F, #40E0D0);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            filter: drop-shadow(3px 3px 0px rgba(0,0,0,0.8));
            '>
            Discover What Makes SHE Different
        </h2>
        <p class="readable-text" style="font-size: 1.2rem;">
            Explore the complete structural phylogeny of eukaryotes with our interactive analysis platform.
        </p>
    </div>
""", unsafe_allow_html=True)

    # Features Grid using the new 'dashboard-card' style
    cols = st.columns(3)

    for i, feature in enumerate(FEATURES):
        col_index = i % 3
        with cols[col_index]:
            
            # 1. URL Slug Logic
            file_name = os.path.basename(feature["link_path"]) 
            slug_raw = os.path.splitext(file_name)[0] 
            if slug_raw[0].isdigit() and "_" in slug_raw:
                slug = slug_raw.split("_", 1)[1]
            else:
                slug = slug_raw

            # 2. Render Card
            st.markdown(f"""
                <a href="{slug}" target="_self" class="card-link">
                    <div class="dashboard-card" style="min-height: 200px; margin-bottom: 20px;">
                        <div class="card-header">
                            <div class="card-icon">{feature["icon"]}</div>
                            <h3>{feature["title"]}</h3>
                            <p class="card-description">{feature["description"]}</p>
                        </div>
                    </div>
                </a>
            """, unsafe_allow_html=True)
            
            # 3. The Demo Button (Optimized)
            with st.popover(f"▶️ Watch Demo", width='stretch'):
                st.markdown(f"**Demo for {feature['title']}**")
                if "demo" in feature and feature["demo"]:
                    video_bytes = get_video_bytes(feature["demo"])
                    if video_bytes:
                        st.video(video_bytes, autoplay=True, loop=True, muted=True)
                else:
                    st.warning("No video demo available for this option yet.")
                st.markdown("<br>", unsafe_allow_html=True)

    st.markdown("<hr style='margin: 4rem 0;'>", unsafe_allow_html=True)


#updates
st.sidebar.subheader('Updates', divider="rainbow")
st.sidebar.write('Sign up for updates!')
email = st.sidebar.text_input('Your email', key="email_input") 
if st.sidebar.button('Subscribe'):
    user_email = st.session_state.email_input.strip()
    if not user_email:
        st.error("Please enter an email address.")
    elif not is_valid_email(user_email):
        st.sidebar.error("Please enter a valid email address (e.g., user@example.com).")
    else:
        if save_email(user_email):
            st.balloons()
            st.sidebar.success(f"Thank you! **{user_email}** is where you will receive updates regarding SHE!")


# --- C. METHODOLOGY SECTION ---
# Add text-shadow: 0 2px 4px rgba(0,0,0,0.6);
# 1. Define a reusable style for maximum readability
st.markdown("""
<style>
    .readable-text {
        color: #ffffff !important;
        font-weight: 600 !important; /* Semi-bold */
        /* The "Double Shadow" technique: 
           Layer 1: A hard dark drop shadow 
           Layer 2: A soft dark glow around the letters */
        text-shadow: 2px 2px 0 #000000, 0 0 10px rgba(0,0,0,0.8);
        line-height: 1.5;
    }
</style>
""", unsafe_allow_html=True)

# 2. The Header Section
st.markdown("""
    <div style="text-align: center; max-width: 42rem; margin: 0 auto 3rem;">
        <h2 style='
            font-size: 2rem; 
            font-weight: 800; 
            margin-bottom: 1rem;
            /* Gradient text needs a filter drop-shadow, not text-shadow */
            background: linear-gradient(to right, #00FF7F, #40E0D0);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            filter: drop-shadow(3px 3px 0px rgba(0,0,0,0.8));
            '>
            How SHE Works
        </h2>
        <p class="readable-text" style="font-size: 1.2rem;">
            SHE compares entire proteomes across all of Eukarya based on 3D protein structures - capturing relationships beyond what sequence and single gene analysis can provide.
        </p>
    </div>
""", unsafe_allow_html=True)

# 3. The Steps Loop (FIXED)
step_cols = st.columns(4)
for i, step in enumerate(STEPS):
    with step_cols[i]:
        # NOTICE: The HTML strings below are pushed to the left to avoid being seen as code blocks
        st.markdown(f"""
<div class="method-step">
    <div class="glass-panel method-number-bubble">
        <span class="method-number text-gradient">{step["number"]}</span>
    </div>
    <div class="step-content">
        <h3 class="step-title" style="color: white; text-shadow: 3px 3px 0 #000; font-weight: 900;">
            {step["title"]}
        </h3>
        <p style="color: white !important; font-weight: 700 !important; font-size: 1rem; text-shadow: 2px 2px 0 #000000; background-color: transparent; line-height: 1.5;">
            {step["description"]}
        </p>
    </div>
</div>
""", unsafe_allow_html=True)
        
st.write("---")

# --- D. CITATION AND UPDATES ---
st.markdown("""
    <h3 style='
        color: white;
        text-shadow: 3px 3px 0 #000000;  /* The hard black shadow */
        font-weight: 800;
        font-size: 1.75rem;
        margin-bottom: 0px;             /* Remove gap between text and line */
    '>
        Cite
    </h3>
""", unsafe_allow_html=True)
st.subheader('', divider="rainbow")
#st.write('Preprint:  Qiuzhen Li, Diandra Daumiller, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')
st.write('Qiuzhen Li, Diandra Daumiller, Patrick Bryant: The structural History of Eukarya. bioRxiv 2026.02.04.703859; doi: https://doi.org/10.64898/2026.02.04.703859')
#align figs
st.markdown("""
    <style>
    [data-testid="stImage"] {
        display: flex;
        align-items: flex-start !important;
        justify-content: center;
        height: 100%;
    }
    img {
        margin-top: 0 !important;
    }
    </style>
""", unsafe_allow_html=True)
colm1, colm2, colm3, colm4 = st.columns(4)
with colm2:
    st.markdown("<br>", unsafe_allow_html=True)
    st.image("data/logos/scilifelab.png", width='content')
with colm3:
    st.image("data/logos/su.png", width='content')