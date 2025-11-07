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
def save_email(email):
    """Saves the email and current timestamp to the CSV file."""
    try:
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
    page_title="SHE",
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

st.subheader('Demo')
st.write('Check up our demo to start explorying SHE')
video_intro = open("data/tree/app_demo.mp4", "rb")
video_bytes = video_file.read()
st.video(video_bytes)

st.subheader('Updates')
st.write('Sign up for updates!')
email = st.sidebar.text_input('Your email', key="email") #strip against whitespace

EMAIL_FILE = 'data/emails_updates.csv'
# Submission Button in the main content area
if st.button('Subscribe'):
    # Retrieve the value from session state using the key "user"
    user_email = st.session_state.email.strip()
    
    if not user_email:
        st.error("Please enter an email address before subscribing.")
    elif not is_valid_email(user_email):
        st.error("Please enter a valid email format (e.g., user@example.com).")
    else:
        # Check if email is already in the CSV (optional, but good practice)
        try:
            with open(EMAIL_FILE, 'r', newline='', encoding='utf-8') as f:
                reader = csv.reader(f)
                # Skip header row if present
                next(reader, None) 
                
                # Simple check for duplicates (case-insensitive)
                existing_emails = [row[0].lower() for row in reader]
                if user_email.lower() in existing_emails:
                    st.warning(f"The email **{user_email}** is already subscribed!")
                    # Stop processing here
                    st.stop()
        except FileNotFoundError:
            # File doesn't exist yet, which is fine
            pass 
        except Exception as e:
            # Handle potential file reading errors
            st.error(f"Could not check for duplicates: {e}")
            st.stop()

        # If we reach here, the email is new and valid
        if save_email(user_email):
            st.success(f"Thank you! **{user_email}** has been added to the updates list.")
            # Clear the input field after successful submission
            st.session_state.user = ""






st.subheader('Cite')
st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: The structural History of Eukarya')