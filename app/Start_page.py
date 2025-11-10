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

#github_path = '/' #change to your path
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
EMAIL_FILE = 'data/email_updates.csv'
def save_email(email):
    """Saves the email and current timestamp to the CSV file. 
    It ensures the directory exists before attempting to write."""
    try:
        # Ensure the directory exists (necessary for 'data/emails_results.csv')
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
#### intro -----------------------------------------------------------------------------------------------------------------------------------
st.set_page_config(
    page_title="SHE: Start page",
    page_icon= "data/tree/logo_tree.png", 
    layout="centered"
)

st.title('SHE: Structural History of Eukarya')

st.sidebar.success("Find more analysis options in the panels above.")


st.subheader('Strucutral Phylogenetic Tree')
#, main_tree_out = load_main_tree()
# Display the main_tree image
st.image(
    'data/tree/All/Intro.png',
    caption='Strucutral Phylogenetic Tree of Eukarya'
)

st.write('Comparative genomics is hindered by the varying evolutionary rates of genes and genomes.')
st.write('Therefore, ribosomal RNA genes have been the primary focus for constructing phylogenetic relationships, as they are more stable.')
st.write('This approach ignores most of the evolutionary changes in an organism arising across its genome.')
st.write('In addition, sequence studies capture genetic changes that may take place, but these do not necessarily impact the resulting protein structures.')
st.write('New developments in protein structure prediction and alignment enable comparisons on the structural level, capturing more distant and detailed evolutionary relationships.')
st.write('However, these comparisons are still constrained to the single protein level.')
st.write('To capture a more complete picture of evolution, we here construct evolutionary relationships based on predictions of all protein structures in an organism, obtaining proteome-wide comparisons.')
st.write('We define a similarity metric between proteomes based on the structural bit score of aligned protein regions and compute the evolutionary history for all Eukarya with available proteomes in UniProt (2.96x1014 structure alignments between 24.3 million proteins).')
st.write('This provides a comprehensive and detailed picture of proteome structural similarity and evolutionary rates across eukaryotes, highlighting the importance of structural novelty at the proteome level for organism complexity, with humans as a clear outlier.')


st.write('Start exploring the new structural tree by selecting the panels in at top of the sidebar!')
st.write('1. Explore: check the structrual phylogenetic tree for all Eukarya, a kingdom, a taxon or a species.')
st.write('2. Analyse: find the best model organism for your study or track the structural evolution of a set of proteins.')
st.write('3. Download: retrieve all the alignments and tree data.')
st.write('You will be able to change the tree visualisation style in the sidebar.')


st.subheader('Cite')
st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: The structural History of Eukarya')

st.subheader('Updates')
st.write('Sign up for updates!')
email = st.text_input('Your email', key="email") #strip against whitespace
# Submission Button in the main content area
# Submission Button in the main content area
if st.button('Subscribe'):
    # Retrieve the value from session state using the key "email" (updated from "user")
    user_email = st.session_state.email.strip()
    
    if not user_email:
        # Only check remaining: must not be empty
        st.error("Please enter an email address.")
    else:
        # No email format validation is performed here.
        
        # If we reach here, the email is non-empty and we save it.
        if save_email(user_email):
            # Updated success message
            st.success(f"Thank you! **{user_email}** is where you will receive updates about SHE.")
            # Clear the input field after successful submission (updated key)
            st.session_state.email = ""


st.subheader('Demo')
st.write('Check up our demo to start explorying SHE')
#video_path = "data/tree/app_demo.mp4"
video_path="https://rd0b3c972.serve.scilifelab.se/files/project-vol/Start_page/app_demo.mp4"

try:
    # 1. Check if the file exists using os.path.exists
    if not os.path.exists(video_path):
        # Raise an exception if the file doesn't exist
        raise FileNotFoundError(f"The file '{video_path}' was not found.")

    # 2. If the file exists, open and read the bytes
    with open(video_path, "rb") as video_intro:
        video_bytes = video_intro.read()
    
    # 3. Display the video in Streamlit
    st.video(video_bytes)
    
except FileNotFoundError as e:
    # Handle the case where the file is not found
    st.warning(f"Demo video could not be loaded: {e}")
    st.info("The application will run without the video.")