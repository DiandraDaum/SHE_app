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

# Make download button for all comparisons
# Traffic analysis code



#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
#@st.cache_data

st.set_page_config(
    page_title="SHE: Downloads",
    page_icon= "data/tree/logo_tree.png", 
    layout="centered"
)

st.title('SHE: Structural History of Eukarya')
#-----------
st.sidebar.write("Do you like it?")
sentiment_mapping = [":material/thumb_down:", ":material/thumb_up:"]
selected = st.sidebar.feedback("thumbs")
if selected is not None:
    st.sidebar.markdown(f"You selected: {sentiment_mapping[selected]}")

st.subheader('Downloads')
st.write('Here it is possible to download the data (structures and trees) from our analysis.')
st.write('Please, select the type of file you want to download from the sidebar.')

st.sidebar.subheader('Downloads', divider="yellow")
if st.sidebar.checkbox('Protein-protein alignments by Organism'):
    st.subheader('Protein-protein alignments by Organism')
    st.write('Here it is possible to download the Foldseek proteome-proteome comparisons by organism')


if st.sidebar.checkbox('Tree by Kingdom'):
    st.subheader('Tree (newick file) by Kingdom')
    st.write('Here it is possible to download the structural newick tree for each kigdom. To open it, try uploading it on iTOL (https://itol.embl.de/upload.cgi)!')
    df = pd.DataFrame({
        'kingdom': ["All", "Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"],
        'number_species': [10, 20, 30, 40, 50, 70, 80] # change to actual number of species
        })

    for kingdom in df['kingdom']:
        file = "data/tree/" + kingdom + "/" + "kingdoms_" + kingdom + "_tree.newick"
        file_name= kingdom + "_tree.newick"
        label1="Download Tree " + kingdom
        label2="Download Tree Eukarya"
        if kingdom == "All":
            with open(file, "rb") as file:
                st.download_button(
                label=label2,
                data=file,
                file_name=file_name,
                mime="image/png",
                icon=":material/download:",
            )
        else:
            #print(f"Downloading file: {file}")
            with open(file, "rb") as file:
                st.download_button(
                    label=label1,
                    data=file,
                    file_name=file_name,
                    mime="image/png",
                    icon=":material/download:",
                )




