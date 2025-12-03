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
from datetime import datetime
import base64

# --- App Content (will appear on top of the video) ---
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
            timestamp = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
            writer.writerow([email, timestamp])
        return True
    except Exception as e:
        st.error(f"An error occurred while saving: {e}")
        return False


# background video
def get_base64_of_bin_file(bin_file):
    """
    Reads a binary file and returns its Base64 encoded string.
    """
    try:
        with open(bin_file, 'rb') as f:
            data = f.read()
        return base64.b64encode(data).decode()
    except FileNotFoundError:
        st.error(f"Error: Video file not found at path: {bin_file}. Please check the path.")
        return None

def embed_video_background(video_path):
    """
    Embeds a local MP4 video as a fixed background using Base64 encoding.
    """
    video_base64 = get_base64_of_bin_file(video_path)
    
    if video_base64 is None:
        return # Stop execution if file wasn't found

    # 1. Construct the Base64 URL for the video source
    video_url = f"data:video/mp4;base64,{video_base64}"

    # 2. Inject the HTML/CSS for the background video
    st.markdown(
        f"""
        <style>
        /* This hides all Streamlit content (which we will put back later) */
        .stApp > header, .stApp > div:not(.video-container) {{
            z-index: 10; /* Bring content above the video */
            position: relative; 
            background: none; /* Make content background transparent */
        }}

        /* Video element to cover the entire screen */
        .video-background {{
            position: fixed;
            top: 0;
            left: 0;
            min-width: 100%;
            min-height: 100%;
            width: auto;
            height: auto;
            z-index: -100; /* Puts the video behind all other content */
            overflow: hidden;
        }}
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

video_path_bg = "data/tree/bg.mp4" 
# Call the function to set the background
embed_video_background(video_path_bg)
#st.sidebar.video(tree_bytes, autoplay=True, loop=True, muted=True)

# htlm video
def embed_html_background(html_file_path):
    """
    Reads a local HTML file and embeds its content as a fixed, full-screen background 
    behind the Streamlit app content.
    """
    try:
        with open(html_file_path, 'r', encoding='utf-8') as f:
            html_content = f.read()

    except FileNotFoundError:
        st.error(f"Error: HTML file not found at path: {html_file_path}. Please check the path.")
        return 

    # Inject the HTML/CSS
    st.markdown(
        f"""
        <style>
        /* CSS styles here... */
        .stApp {{ background: none; }}
        .stApp > header, .stApp > div:not(.html-container) {{
            z-index: 10;
            position: relative; 
            background: none;
        }}
        .html-background {{
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            z-index: -100;
            overflow: hidden;
        }}
        </style>
        
        <div class="html-container">
            <div class="html-background">
                {html_content}
            </div>
        </div>
        """,
        unsafe_allow_html=True
    )


def embed_html_background_iframe(html_file_path):
    """
    Embeds a full HTML file (with its own <!DOCTYPE>, <html>, etc.) 
    as a fixed, full-screen background using a robust iframe.
    
    NOTE: The HTML file must be accessible via Streamlit's static path.
    """
    # Use the HTML file path directly as the iframe source (src)
    st.markdown(
        f"""
        <style>
        /* 1. Ensure Streamlit content is on top and transparent */
        .stApp > header, .stApp > div:not(.iframe-container) {{
            z-index: 10;
            position: relative; 
            background: none;
        }}

        /* 2. CSS for the Fixed Background Layer using an iframe */
        .iframe-background {{
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            border: none; 
            z-index: -100; /* Puts it behind all Streamlit content */
            overflow: hidden;
            pointer-events: none; /* Allows mouse clicks to reach Streamlit widgets on top */
        }}
        </style>
        
        <div class="iframe-container">
            <iframe class="iframe-background" src="{html_file_path}">
            </iframe>
        </div>
        """,
        unsafe_allow_html=True
    )

#html_path_bg = "background.html" 
#embed_html_background_iframe(html_path_bg)


#### intro -----------------------------------------------------------------------------------------------------------------------------------
st.set_page_config(
    page_title="SHE: Start page",
    page_icon= "data/tree/logo_tree.png", 
    layout="centered"
)



st.title('SHE: Structural History of Eukarya')

st.sidebar.success("Find more analysis options in the panels above.")

st.subheader('Structural Phylogenetic Tree from proteome-wide structure alignments')
#st.image('data/tree/All/Intro.png', caption='Strucutral Phylogenetic Tree of Eukarya')


st.write('Classical phylogenetic trees constructed based on sequence alignments of a small set of genes (i.e. rRNAs) are hindered by limitations.')
st.write('To capture more distant and detailed picture of evolution, here we constructed evolutionary relationships based on the comparison of all predicted protein structures in an organism (2.96x1014 structure alignments between 24.3 million proteins).')
st.write('The Structural History of Eukarya (SHE) provides a comprehensive picture of proteome structural similarity and evolutionary rates across eukaryotes (1542 species).')

st.write('Start exploring the new proteome-wide structural phylogenetic tree by selecting the panels in the sidebar!')
row1 = st.columns(3)
emojis = [":tree:", ":mouse:", ":download:"]
titles = ["#1", "#2", "#3"]
for col, emoji, content in zip(row1, emojis, titles):
    tile = col.container(height=120)
    tile.title(f"{emoji} {content}")
    if content == "#1":
        st.write('Explore Tree: check the structuRal phylogenetic tree for all Eukarya, a kingdom, a taxon, a species or a set of proteins. You will be able to change the tree visualisation style in the sidebar.')
        st.page_link("pages/1_Explore_Tree.py", label="Page1", icon=":tree:")
        #demo here
        video_path_demo1 = "data/tree/app_demo3.mp4"
        if os.path.exists(video_path_demo1):
            with open(video_path_demo1, "rb") as video_intro:
                video_bytes = video_intro.read()
    elif content== "#2":
        st.write('Model Organisms: find the best model organism for your study based on the structural comparison of the protein(s) in your pathway of interest.')
        st.page_link("pages/2_Model_organisms.py", label="Page2", icon=":mouse:")
        #demo here
        video_path_demo2 = "data/tree/app_demo2.mp4"
        if os.path.exists(video_path_demo2):
            with open(video_path_demo2, "rb") as video_intro:
                video_bytes = video_intro.read()
    elif content == "#3":
        st.write('Download: retrieve all the structural alignments scores, protein predictions and tree files.')
        st.page_link("pages/3_Downloads.py", label="Page3", icon=":material/download:")
        #demo here
        video_path_demo3 = "data/tree/app_demo3.mp4"
        if os.path.exists(video_path_demo3):
            with open(video_path_demo3, "rb") as video_intro:
                video_bytes = video_intro.read()

st.sidebar.subheader('Updates', divider="rainbow")
st.sidebar.write('Sign up for updates!')
email = st.sidebar.text_input('Your email', key="email") #strip against whitespace
if st.sidebar. button('Subscribe'):
    # Retrieve the value from session state using the key "email" (updated from "user")
    user_email = st.session_state.email.strip()
    if not user_email:
        # Only check remaining: must not be empty
        st.error("Please enter an email address.")
    else:
        # If we reach here, the email is non-empty and we save it.
        if save_email(user_email):
            st.balloons()
            # Updated success message
            st.sidebar.success(f"Thank you! **{user_email}** is where you will receive updates regarding SHE!")
            # Clear the input field after successful submission (updated key)
            #st.session_state.email = ""


video_tree = open("data/tree/final_tree.mp4", "rb")
tree_bytes = video_tree.read()
st.sidebar.video(tree_bytes, autoplay=True, loop=True, muted=True)


st.subheader('Demo')
st.write('Check up our demo to start exploring SHE')
video_path_demo = "data/tree/app_demo.mp4"

try:
    # Check if the file exists using os.path.exists
    if os.path.exists(video_path_demo): # local path
        # If the file exists, open and read the bytes
        with open(video_path_demo, "rb") as video_intro:
            video_bytes = video_intro.read()
    
    st.video(video_bytes)
except:
    # in case the file is not found
    st.warning(f"Demo video could not be loaded")



st.subheader('Cite')
st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')