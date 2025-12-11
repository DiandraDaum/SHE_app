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
        return "#080808"
    

#### intro -----------------------------------------------------------------------------------------------------------------------------------
st.set_page_config(
    page_title="SHE: Explore",
    page_icon= "data/tree/logo_tree.png", 
    layout="centered"
)
st.title('SHE: Structural History of Eukarya')

st.subheader('Structural Phylogenetic Tree')

row1 = st.columns(2)
#emojis = [":deciduous_tree:", ":mouse2:", ":material/download:"]
titles = ["1. Explore Kingdoms", "2. Proteome evolutionary rate"]

# Create the tiles and their content
for col, content in zip(row1,  titles):
    #with col:   tile = st.container(height=120)   tile.title(f"{emoji} {content}")
    with col.expander(f"**{content}**"):
        # --- kingdoms ---
        if content == "1. Explore Kingdoms":
            st.write(':blue[Start exploring the new structural phylogenetic tree by selecting a kingdom in the sidebar!]')
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
st.sidebar.subheader('1. Visualise Kingdoms', divider=True)
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
    /* Target the selectbox input area */
    div[data-baseweb="select"] > div:first-child {{
        background-color: {color};
        border-radius: 0.5rem; /* Optional: for a nicer box look */
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
    st.title('Kingdom tree: ' + kingdom)
    'Selected Kingdom: ', kingdom
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