import os
import streamlit as st
import pandas as pd
import numpy as np
import csv
import sys
import pathlib
import requests
from io import StringIO
import copy
import glob
import math

#github_path = '/' #change to your path
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
@st.cache_data
# for coloring the selectbox
def kingdom_color(kingdom):
    if kingdom == 'Opisthokonta':
        return "#bfdbfe"
    elif kingdom == 'Viridiplantae':
        return "#bbf7d0"
    elif kingdom == 'Amoebozoa':
        return "#fecdd3"
    elif kingdom == 'Sar':
        return "#f3f4f6"
    elif kingdom == 'Discoba':
        return "#e9d5ff"
    elif kingdom == 'Rhodophyta':
        return "#c3f8f1"
    elif kingdom == 'Cryptophyceae':
        return "#f6fad5"
    elif kingdom == 'Apusozoa':
        return "#ddd6fe"
    else:
        # This will be used for internal nodes ('Unknown') or leaves with no match
        return "#FAF8F8"
    

#### intro -----------------------------------------------------------------------------------------------------------------------------------
st.set_page_config(
    page_title="SHE: Explore",
    page_icon= "data/tree/logo_tree.png", 
    layout="wide"
)
#st.title(':green[SHE: Structural History of Eukarya]')

st.subheader('Structural Phylogenetic Tree')
st.markdown("""
<style>
/* 1. Target the Expander Container (The 'details' HTML tag) */
[data-testid="stExpander"] {
    border: 1px solid rgba(250, 250, 250, 0.2); /* Subtle default border */
    border-radius: 12px;                        /* Rounded corners */
    transition: all 0.3s ease-in-out;           /* Smooth animation */
}

/* 2. HOVER STATE: The "Green Aura" */
[data-testid="stExpander"]:hover {
    border-color: #22c55e;   /* Brighter green border */
    
    /* The Aura Effect: Green Glow */
    box-shadow: 0 0 20px rgba(74, 222, 128, 0.4), 
                0 0 10px rgba(74, 222, 128, 0.2);
}

/* 3. OPTIONAL: Make the Title Green on Hover too */
[data-testid="stExpander"]:hover summary {
    color: #4ade80 !important;
}

/* 4. OPTIONAL: Customize the internal icon color */
[data-testid="stExpander"] summary svg {
    fill: #9ca3af !important; /* Default gray */
    transition: fill 0.3s ease-in-out;
}
[data-testid="stExpander"]:hover summary svg {
    fill: #4ade80 !important; /* Green on hover */
}
</style>
""", unsafe_allow_html=True)
row1 = st.columns(2)
#emojis = [":deciduous_tree:", ":mouse2:", ":material/download:"]
titles = ["1. Explore Clades", "2. Proteome evolutionary rate"]

# Create the tiles and their content
for col, content in zip(row1,  titles):
    #with col:   tile = st.container(height=120)   tile.title(f"{emoji} {content}")
    with col.expander(f"**{content}**"):
        # --- kingdoms ---
        if content == "1. Explore Clades":
            st.write(':blue[Start exploring the new structural phylogenetic tree by selecting a major taxon in the sidebar!]')
            st.write(':blue[You can also choose the tree visualisation style in the sidebar.]')
            
        # --- evorate ---
        elif content == "2. Proteome evolutionary rate":
            st.write(':green[Compare the innovative structural evolutionary rate of different proteomes across eukarya! Use the "Visualise Evolutionary Rate" section in the sidebar.]')
                


# Add sidebar ---------------------------------------------------------------------------------------------------------------------------
#data
df = pd.DataFrame({
    'kingdom': ["Eukarya", "Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"],
    'number_species': [1000000, 10, 20, 30, 40, 50, 70] # change to actual number of species
    })

st.sidebar.subheader('Visualization Options')
layout = st.sidebar.selectbox(
    'Tree Layout',
    ['Circular', 'Linear']
)

#add labels
if st.sidebar.toggle('Labels'):
    labels = True
else:
    labels = False
    
# Choose Taxon
st.sidebar.subheader('1. Visualise Clades', divider=True)
kingdom = st.sidebar.selectbox(
    'Which kingdom? (default: Eukarya)',
    df['kingdom'],
    key = 'kingdom_selector')

# color for the currently selected kingdom
color = kingdom_color(kingdom)
# CSS to style the selectbox based on the current selection.
st.markdown(
    f"""
    <style>
    /* 1. Target the background of the selectbox */
    div[data-baseweb="select"] > div:first-child {{
        background-color: {color};
        border-radius: 0.5rem;
    }}
    
    /* 2. Force the text (selected value) to be BLACK */
    div[data-baseweb="select"] > div:first-child div,
    div[data-baseweb="select"] > div:first-child span {{
        color: black !important;
        -webkit-text-fill-color: black !important; /* Necessary for some browsers */
        font-weight: 500; /* Optional: make it slightly bolder for readability */
    }}
    
    /* 3. Force the dropdown arrow icon to be BLACK */
    div[data-baseweb="select"] svg {{
        fill: black !important;
        color: black !important;
    }}
    </style>
    """,
    unsafe_allow_html=True
)


############################Visualization options --------------------------------------------------------

#show kingdom trees
# 1. Initialize session state for visibility
if 'show_tree' not in st.session_state:
    st.session_state.show_tree = True
if st.session_state.show_tree:
    button_label = '⬇️ Hide Kingdom Tree'
else:
    button_label = '⬆️ Show Kingdom Tree'
if st.button(button_label):
    # This line toggles the visibility state (True -> False, False -> True)
    st.session_state.show_tree = not st.session_state.show_tree


# Display the kingdom tree image --------------------------------------------------------------
if st.session_state.show_tree:
    st.title('Clade tree: ' + kingdom)
    'Selected Clade: ', kingdom
    'Selected layout: ', layout
    if labels == True:
        'Labels: ON'
    else:
        'Labels: OFF'
    #kingdom_tree, kingdom_tree_out = load_tree_1000(kingdom, layout)
    if labels == True:
        kingdom_tree_out= 'data/tree/' + kingdom + '/kingdoms_' + kingdom + "_" + layout + "_tree1000_labels.svg"
    else:
        kingdom_tree_out= 'data/tree/' + kingdom + '/kingdoms_' + kingdom + "_" + layout + "_tree1000.svg"
    if kingdom == "Sar":
        what = kingdom + "_supergroup"
        st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
    elif kingdom == "Eukarya":
        st.page_link("https://en.wikipedia.org/wiki/Eukaryote", label="Wiki", icon="🌎")
    else: 
        what = kingdom
        st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
    st.image(
        kingdom_tree_out, # main_tree_out
        caption= kingdom + ' Phylogenetic Tree',
        width='stretch'
    )
if st.sidebar.toggle('Show species branch rate in kingdom'):
    kingdom_csv = 'data/tree/' + kingdom + '/kingdoms_' + kingdom + "_branch_rates.csv" #i.e kingdoms_amoebozoa_branch_rates.csv
    kingdom_csv = pd.read_csv(kingdom_csv)
    kingdom_csv

# show evorate tree color ------------------------------------------------------------------------
st.sidebar.subheader('2. Evolutionary rate', divider="green")
if st.sidebar.toggle('See tree with evolutionary rate'):
    st.title('Evolutionary rate:')
    'Selected Kingdom: ', kingdom
    'Selected layout: ', layout
    tree_svg = pathlib.Path('data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Linear_tree_evorate.svg")
    if tree_svg.exists():
        if layout == 'Circular':
            if labels==True:
                tree_c='data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Circular_tree_evorate_labels.svg"
                st.image(
                    tree_c,
                    caption= 'Circular Tree with evolutionary rate; color legend: yellow high, purple low',
                    width='stretch')
            elif labels==False:
                tree_c='data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Circular_tree_evorate.svg"
                st.image(
                    tree_c,
                    caption= 'Circular Tree with evolutionary rate; color legend: yellow high, purple low',
                    width='stretch')
        elif layout=='Linear':
            if labels==True:
                tree_r = 'data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Linear_tree_evorate_labels.svg"
                st.image(
                    tree_r,
                    caption= 'Linear Tree with evolutionary rate; color legend: yellow high, purple low',
                    width='stretch')
            elif labels==False:
                tree_r = 'data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Linear_tree_evorate.svg"
                st.image(
                    tree_r,
                    caption= 'Linear Tree with evolutionary rate; color legend: yellow high, purple low',
                    width='stretch')

st.sidebar.page_link("pages/4_Downloads.py", label="Download trees", icon="⬇️")
# --- CITATION
st.subheader('Cite', divider="rainbow")
st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')
