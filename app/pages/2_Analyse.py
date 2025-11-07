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

#Code to analyse:
# model organism search
#gene set search --> tree build


#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
@st.cache_data
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


###------
#user sign up for updates and analysis results
st.subheader('Results')
st.write('Sign up to receive your results by email')
email = st.sidebar.text_input('Your email', key="email") #strip against whitespace

EMAIL_FILE = 'data/emails_results.csv'
# Submission Button in the main content area
if st.button('Subscribe'):
    # Retrieve the value from session state using the key "user"
    user_email = st.session_state.email.strip()
    
    if not user_email:
        st.error("Please enter an email address ")
    elif not is_valid_email(user_email):
        st.error("Please enter a valid email format (e.g., name.surname@example.com).")
    else:
        # Check if email is already in the CSV (optional, but good practice)
        try:
            with open(EMAIL_FILE, 'r', newline='', encoding='utf-8') as f:
                reader = csv.reader(f)
                # Skip header row if present
                next(reader, None) 

        # If we reach here, the email is new and valid
        if save_email(user_email):
            st.success(f"Thank you! **{user_email}** is where you will receive your results.")
            # Clear the input field after successful submission
            st.session_state.user = ""


