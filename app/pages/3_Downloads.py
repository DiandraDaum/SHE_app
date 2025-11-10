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

st.subheader('Downloads')

st.write('Here it is possible to download the Foldseek proteome-proteome comparisons by organism')



#-----------
st.sidebar.write("Do you like it?")
sentiment_mapping = [":material/thumb_down:", ":material/thumb_up:"]
selected = st.sidebar.feedback("thumbs")
if selected is not None:
    st.sidebar.markdown(f"You selected: {sentiment_mapping[selected]}")
