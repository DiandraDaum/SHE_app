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

EMAIL_FILE = 'data/email_results.csv'
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
#@st.cache_data
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

###------
#user sign up for updates and analysis results
st.sidebar.subheader('Results')
st.sidebar.write('Sign up to receive your results by email')
st.sidebar.text_input('Your email', key="email") #strip against whitespace


# Submission Button in the main content area
if st.sidebar. button('Receive Results'):
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
            st.success(f"Thank you! **{user_email}** is where you will receive your results.")
            # Clear the input field after successful submission (updated key)
            st.session_state.email = ""


