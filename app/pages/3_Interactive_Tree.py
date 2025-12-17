import os
import streamlit as st
import pandas as pd
import numpy as np
import time
import csv
from ete3 import Tree, TreeStyle, NodeStyle
import sys
import pathlib
import base64
import time
import requests
from urllib.parse import urlparse
import py3Dmol
from stmol import showmol
from stmol import *
from IPython.display import HTML
import plotly.graph_objects as go
import plotly.express as px
from scipy.spatial.distance import pdist, squareform
from Bio import Phylo
from Bio.Phylo import Consensus
from biotite.sequence.phylo import neighbor_joining
from io import StringIO
import copy
import shutil
import glob
from pathlib import Path
import foldcomp as fc
import random
from collections import defaultdict
from tqdm import tqdm
import math

# Define navigation items for the custom header
NAV_ITEMS = [
    {"label": "SHE", "icon": "🏠", "path": "app.py"},
    {"label": "Model Organisms Ranking", "icon": "🐁", "path": "pages/1_Model_Organism_Ranking.py"},
    {"label": "Explore SHE tree", "icon": "🌳", "path": "pages/2_Explore_Tree.py"},
    {"label": "Interactive SHE Tree", "icon": "🔎", "path": "pages/3_Interactive_Tree.py"},
    {"label": "Downloads", "icon": ":material/download:", "path": "pages/4_Downloads.py"},
]
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
#@st.cache_data
# for coloring the selectbox
@st.cache_data
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
    
#for find species
@st.cache_data
def extract_closest_species(newick_file, input_species, num_neighbors): #species
    """
    Extracts a subtree containing a target species and its N closest relatives
    based on phylogenetic distance.

    Args:
        newick_file (str): Path to the input Newick tree file.
        target_species (str): The name of the target species to find neighbors for.
                              (e.g., 'Homo_sapiens'). This is case-sensitive.
        num_neighbors (int): The number of closest species to extract. default = 10?
        output_file (str): Path to save the resulting subtree in Newick format.
    """
    # Step 1: Load the main tree.
    newick_file = 'data/tree/Eukarya/kingdoms_Eukarya_tree.newick'
    print(f"Loading tree from '{newick_file}'...")
    try:
        main_tree = Tree(newick_file, format=1)
    except FileNotFoundError:
        print(f"Error: Newick file not found at '{newick_file}'")
        return
    except Exception as e:
        print(f"An error occurred while loading or parsing the Newick tree: {e}")
        return

    # Step 2: Find the target species leaf node in the tree.
    # The leaf name might contain extra info (e.g., 'Homo_sapiens:proteinID'),
    # so we search for a leaf name that starts with the target species name.
    target_node = None
    all_leaf_names = main_tree.get_leaf_names()
    target_species = input_species # Assume it is the species name by default
    # Check if the input looks like an ID (UniProt starts with 'UP', Organism ID is digit-only)

    # --- Step 3: Find the target species leaf node in the tree ---

    #target species could be leaf name, UniProt ID or proteome ID
    # Find the full leaf name that corresponds to the target species
    #target_species = input_species
    matching_full_names = [name for name in all_leaf_names if name.startswith(target_species)] 
    
    if not matching_full_names:
        print(f"Error: Target species '{target_species}' could not be found in the tree.")
        return
        
    # In case of multiple partial matches, we take the first one.
    # This assumes names like 'Homo_sapiens' and 'Homo_sapiens_neanderthalensis' are distinct enough.
    full_target_name = matching_full_names[0]
    target_node = main_tree.get_leaves_by_name(full_target_name)[0]
    print(f"Found target node: '{full_target_name}'")

    # Step 3: Calculate the distance from the target to all other leaves.
    distances = []
    print("Calculating distances to all other species...")
    for leaf in main_tree.iter_leaves():
        # No need to calculate distance to self, but doing so is harmless.
        distance = target_node.get_distance(leaf)
        distances.append((distance, leaf.name))

    # Step 4: Sort species by distance and select the closest N neighbors.
    distances.sort(key=lambda x: x[0])

    # The first item will be the target itself (distance 0). We need it plus N neighbors.
    total_to_keep = num_neighbors + 1
    
    if len(distances) < total_to_keep:
        print(f"Warning: The tree has only {len(distances)-1} other species, which is less than the requested {num_neighbors}.")
        print("The output tree will contain all species.")
        leaves_to_keep = [name for dist, name in distances]
    else:
        # Get the names of the target and its N closest neighbors.
        closest_species = distances[:total_to_keep]
        leaves_to_keep = [name for dist, name in closest_species]

    print(f"\nTop {num_neighbors} closest species to '{full_target_name}':")
    for dist, name in distances[1:total_to_keep]:
        print(f"- {name} (distance: {dist:.4f})")

    # Step 5: Prune the tree and save the result.
    main_tree.prune(leaves_to_keep, preserve_branch_length=True)
        # Step 3: Save the new, smaller subtree.
    if len(main_tree) > 0:
        output_path = 'results/species_trees/'+ target_species +'_subtree.newick'
        main_tree.write(outfile=output_path, format=1)
        print(f"\nSuccess! Extracted subtree with {len(leaves_to_keep)} species saved to '{output_path}'")
    else:
        print("The resulting subtree has no leaves. No file was written.")
    
    return output_path, target_species

@st.cache_data
def load_species_tree(species, layout, species_out):
# File paths
    #species_out =  'SHE/species_trees/' + species +'_subtree.newick'
    # Load the Tree: i.e. tree = Tree(main_tree_file)
    species_tree = Tree(species_out, format=1) #format=0
    svg_path =  'results/species_trees/'+ species +'_subtree.svg' 

    # Configure the TreeStyle
    if layout == 'Circular':
        ts = TreeStyle()
        ts.mode = "c"          # Draw tree in circular mode otherwise "p"
        ts.rotation=210
        ts.show_leaf_name = True #true or false
        ts.show_branch_length = True
        ts.scale = None #big circles with dotted lines
        ts.optimal_scale_level = "full" #mid or full
        ts.force_topology = True #not bad
        # Render tree
        species_tree.render(
        svg_path,
        tree_style=ts,
        w=500,
        units="px" 
        )
    elif layout == 'Linear':
        ts = TreeStyle()
        ts.mode = "r"   
        ts.show_leaf_name = True #true or false
        ts.show_branch_length = False
        ts.force_topology = True #not bad
        # Render tree
        species_tree.render(
        svg_path,
        tree_style=ts,
        #w=183,
        units="px"  
        )
    
    return svg_path

#for find taxon
@st.cache_data
def extract_subtree(csv_file, newick_file, taxonomic_group): #taxon
    """
    Extracts a subtree from a larger phylogenetic tree based on a taxonomic group.

    Args:
        csv_file (str): Path to the CSV file containing organism names and their
                        taxonomic lineage.
        newick_file (str): Path to the Newick file of the main phylogenetic tree.
        taxonomic_group (str): The name of the taxonomic group to extract (e.g., 'Primates').
                               This is case-sensitive.
    """
    # Step 1: Find all organisms belonging to the specified taxonomic group from the CSV.
    output_file = 'results/extracted_trees/' + taxonomic_group +'_subtree.newick' #Path to save the resulting subtree in Newick format.
    organisms_to_keep = set()
    print(f"Searching for organisms in the '{taxonomic_group}' group...")
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            # Check for required columns in the CSV header
            if 'Organism_formatted' not in reader.fieldnames or 'Taxonomic lineage' not in reader.fieldnames:
                print("Error: CSV file must contain 'Organism_formatted' and 'Taxonomic lineage' columns.")
                return

            for row in reader:
                # This checks if the taxonomic group string is part of the full lineage string.
                if taxonomic_group in row['Taxonomic lineage']:
                    organisms_to_keep.add(row['Organism_formatted'])
    except FileNotFoundError:
        print(f"Error: CSV file not found at '{csv_file}'")
        return
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")
        return

    if not organisms_to_keep:
        print(f"No organisms found for the taxonomic group '{taxonomic_group}'. No tree will be generated.")
        return

    print(f"Found {len(organisms_to_keep)} organisms in the '{taxonomic_group}' group.")

    # Step 2: Load the main tree and prune it to keep only the desired organisms.
    print(f"Loading tree from '{newick_file}'...")
    try:
        main_tree = Tree(newick_file, format=1)
    except FileNotFoundError:
        print(f"Error: Newick file not found at '{newick_file}'")
        return
    except Exception as e: # Catches potential ete3 parsing errors
        print(f"An error occurred while loading or parsing the Newick tree: {e}")
        return
        
    # The 'prune' method from ete3 efficiently keeps only the specified leaves.
    # We must match the organism names from the CSV to the full leaf names in the tree,
    # which might contain extra info (e.g., 'Organism_Name:ProteinID').
    
    # Create a map from the simple organism name to the full leaf name in the tree.
    all_leaf_names = main_tree.get_leaf_names()
    org_name_map = {name.split(':')[0] if ':' in name else name: name for name in all_leaf_names}
    
    leaves_in_tree_to_keep = []
    for org in organisms_to_keep:
        if org in org_name_map:
            leaves_in_tree_to_keep.append(org_name_map[org])

    if not leaves_in_tree_to_keep:
        print(f"Warning: None of the {len(organisms_to_keep)} organisms from the CSV were found as leaves in the tree file.")
        print("Please check for inconsistencies in organism names between your files.")
        return

    print(f"Pruning tree to keep {len(leaves_in_tree_to_keep)} matching leaves.")
    main_tree.prune(leaves_in_tree_to_keep, preserve_branch_length=True)

    # Step 3: Save the new, smaller subtree.
    if len(main_tree) > 0:
        main_tree.write(outfile=output_file, format=1)
        print(f"\nSuccess! Extracted subtree saved to '{output_file}'")
    else:
        print("The resulting subtree has no leaves. No file was written.")
    return output_file

@st.cache_data
def load_extracted_tree(taxon, layout, taxon_out):
# File paths
    taxon = taxon
    #taxon_tree_path = tree
    taxon_tree_path = 'results/extracted_trees/' + taxon +'_subtree.newick'
    taxon_tree_out = taxon_out

    # Load the Tree: i.e. tree = Tree(main_tree_file)
    taxon_tree = Tree(taxon_tree_path, format=1) #format=0

    # Configure the TreeStyle
    if layout == 'Circular':
        ts = TreeStyle()
        ts.mode = "c"          # Draw tree in circular mode otherwise "p"
        ts.rotation=210
        ts.show_leaf_name = True #true or false
        ts.show_branch_length = True
        ts.scale = None #big circles with dotted lines
        ts.optimal_scale_level = "full" #mid or full
        ts.force_topology = True #not bad
        # Render tree
        taxon_tree.render(
        taxon_tree_out,
        tree_style=ts,
        w=500,
        units="px" 
        )
    elif layout == 'Linear':
        ts = TreeStyle()
        ts.mode = "r"   
        ts.show_leaf_name = True #true or false
        ts.show_branch_length = False
        ts.force_topology = True #not bad
        # Render tree
        taxon_tree.render(
        taxon_tree_out,
        tree_style=ts,
        #w=183,
        units="px"  
        )
    
    return taxon_tree_out

# for find protein: color tree if a protein is present in species
def get_species_id(protein_name):
    """Robustly parses a protein name to find a UniProt proteome ID."""
    if isinstance(protein_name, str) and protein_name.startswith("UP"):
        return protein_name.split('-')[0].split('_')[0]
    return None

def protein_in_family(protein, families_csv_path):
    """
    Find the homologs family of a protein UniProt ID.
    """
    families_df = pd.read_csv(families_csv_path)
    #new uniport id column
    #families_df['uniprot_id'] = families_df['protein'].split('-')[2].split('-')[0] #if UP000054538_930991-AF-A0A0D0DNZ3-F1-model_v4, keep A0A0D0DNZ3
    families_df['uniprot_id'] = families_df['protein'].astype(str).str.split('-').str[2] #if UP000054538_930991-AF-A0A0D0DNZ3-F1-model_v4, keep A0A0D0DNZ3
    #print(families_df['uniprot_id'].head())
    # Look for a direct match in the identified column
    match_protein = families_df[families_df['uniprot_id'] == protein]
    if not match_protein.empty:
        family = match_protein['family_id'].iloc[0].strip()
        print(f"Mapped Protein UniProt ID '{protein}' to family: {family}")
        return family, families_df
    else:
        st.warning(f"This protein does not seem to have a structural homolog...")
        st.stop()
        sys.exit(1)



def get_tree_coordinates(tree):
    # This function traverses the tree and assigns x, y coordinates
    
    # 1. Assign Y-coordinates (depth/level)
    # The leaves will be at y=0, and the root will be at the max depth
    leaves = tree.get_leaves()
    max_depth = len(leaves)
    
    for i, leaf in enumerate(leaves):
        leaf.y = i
    
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.x = node.dist # X is the accumulated branch length
        else:
            # For internal nodes, Y is the mean of its children's Y coordinates
            node.y = sum(child.y for child in node.children) / len(node.children)
            # X is the max branch length of its children + its own branch length
            node.x = max(child.x for child in node.children) + node.dist
            
    # Normalize X to be distance from root
    root_x = tree.x
    for node in tree.traverse("preorder"):
        node.x = root_x - node.x

    return tree

@st.cache_data
def create_plotly_data(tree):
    # This collects all the branch line segments and node scatter points
    branch_lines_x = []
    branch_lines_y = []
    node_coords = {'x': [], 'y': [], 'text': [], 'name': []}

    for node in tree.traverse():
        # Add node coordinates for the scatter plot
        node_coords['x'].append(node.x)
        node_coords['y'].append(node.y)
        node_coords['name'].append(node.name if node.name else f"Internal Node {node.x}")
        
        # Create custom hover text
        hover_text = f"**Node:** {node.name if node.name else 'Internal Node'}<br>"
        if not node.is_root():
            hover_text += f"**Dist:** {node.dist:.3f}"
        node_coords['text'].append(hover_text)
        
        # Create line segments for branches
        if not node.is_root():
            parent = node.up
            # 1. Horizontal segment (from child X to parent X at child Y)
            # Extends the X list: [Child_X, Parent_X, None]
            branch_lines_x.extend([node.x, parent.x, None])
            # Extends the Y list: [Child_Y, Child_Y, None]
            branch_lines_y.extend([node.y, node.y, None])
            
            # 2. Vertical segment (from child Y to parent Y at Parent X)
            # This segment connects the horizontal line to the parent node's vertical position
            # Extends the X list: [Parent_X, Parent_X, None]
            branch_lines_x.extend([parent.x, parent.x, None])
            # Extends the Y list: [Child_Y, Parent_Y, None]
            branch_lines_y.extend([node.y, parent.y, None])
            
    return branch_lines_x, branch_lines_y, node_coords

def display_interactive_tree(tree_path):
    
    try:
        # Load the tree from the Newick file
        tree = Tree(tree_path, format=1)
        
        # Calculate coordinates
        tree = get_tree_coordinates(tree)
        
        # Generate Plotly data lists
        branch_lines_x, branch_lines_y, node_coords = create_plotly_data(tree)
        
    except Exception as e:
        print(f"Error processing tree file: {e}")
        return

    # --- Create the Plotly Figure ---
    #proteomes_path = 'data/tree/selected_proteomes_above_25_formatted_organisms.csv'
    proteomes_path = 'data/core/merged_evo_rate.csv'
    proteomes_df = pd.read_csv(proteomes_path)
    proteomes_df.columns = proteomes_df.columns.str.strip().str.replace(' ', '_')
    proteomes_df['Kingdom'] = proteomes_df['Taxonomic_lineage'].str.split(', ').str[1]
    organism_to_kingdom = proteomes_df.set_index('Organism_formatted')['Kingdom'].to_dict()
    # Generate the list of colors for the nodes based on the names in node_coords
    node_colors = []
    for node_name in node_coords['name']:
        # Look up the kingdom. Use 'Unknown' if the name isn't found 
        # (usually for internal tree nodes)
        kingdom = organism_to_kingdom.get(node_name, 'Unknown') 
        
        # Determine the final color using the function
        color = kingdom_color(kingdom)
        node_colors.append(color)

    # 1. Branch Traces (Lines)
    trace_branches = go.Scatter(
        x=branch_lines_x, y=branch_lines_y,
        mode='lines',
        line=dict(color='black', width=1),
        hoverinfo='none', # Branches don't need hover
        name='Branches'
    )

    # 2. Node Traces (Scatter Points)
    trace_nodes = go.Scatter(
        x=node_coords['x'], y=node_coords['y'],
        mode='markers+text',
        #marker=dict(size=5, color='red'),
        marker=dict(
            size=10, 
            color=node_colors, 
            line=dict(width=1, color='Black') # Add a border for visibility
        ),
        text=node_coords['name'], # Leaf/Node name displayed next to the point
        textposition="middle right",
        textfont=dict(size=10),
        hovertext=node_coords['text'], # Custom hover text
        hoverinfo='text',
        name='Nodes'
    )

    legend_traces = []
    unique_kingdoms = sorted(proteomes_df['Kingdom'].unique()) # Get all unique kingdoms

    for kingdom in unique_kingdoms:
        color = kingdom_color(kingdom)
        # Dummy trace for the legend entry
        legend_traces.append(
            go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(size=10, color=color, line=dict(width=1, color='Black')),
                name=kingdom,
                legendgroup='kingdoms'
            )
        )

    # 3. Create the figure and layout
    fig = go.Figure(data=[trace_branches, trace_nodes])

    fig.update_layout(
        title='Interactive Tree for evolution visualisation',
        xaxis=dict(showgrid=False, title='Branch length'),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False), # Hide y-axis labels
        showlegend=True, #False
        hovermode="closest",
        plot_bgcolor='white',
        margin=dict(l=20, r=20, t=40, b=20),
        height=800
    )
    st.plotly_chart(fig, use_container_width=True)



# new plotly circular
def polar_coordinates_transform_c(tree):
    """
    Transforms the Linear (X, Y) coordinates of the tree nodes 
    into polar (Radius, Angle) coordinates for a Circular layout.
    """
    # Max Y coordinate (used for normalization to 360 degrees)
    max_y = max(node.y for node in tree.get_leaves())
    
    for node in tree.traverse():
        # 1. Radius (r) is the X-coordinate (accumulated branch length from root)
        r = node.x
        
        # 2. Angle (theta) is based on the Y-coordinate (normalized leaf order)
        # Normalize Y to 0-360 degrees (0 to 2*pi radians)
        theta_rad = (node.y / max_y) * 2 * math.pi
        
        # Store polar coordinates on the node
        node.r = r
        node.theta_rad = theta_rad
        
        # Store Cartesian conversion for plotting (optional, but often helpful)
        node.x_polar = r * math.cos(theta_rad)
        node.y_polar = r * math.sin(theta_rad)
    return tree


def create_plotly_data_c(tree, layout):
    branch_lines_x = []
    branch_lines_y = []
    node_coords = {'x': [], 'y': [], 'text': [], 'name': []}

    # Determine which coordinates to use
    if layout == 'Circular':
        x_attr = 'x_polar'
        y_attr = 'y_polar'
    else: # 'Linear'
        x_attr = 'x'
        y_attr = 'y'

    for node in tree.traverse():
        
        # Get coordinates dynamically
        node_x = getattr(node, x_attr)
        node_y = getattr(node, y_attr)
        
        # Add node coordinates for the scatter plot
        node_coords['x'].append(node_x)
        node_coords['y'].append(node_y)
        # ... (rest of node_coords dictionary population remains the same) ...
        node_coords['name'].append(node.name if node.name else f"Internal Node {node_x}")
        hover_text = f"**Node:** {node.name if node.name else 'Internal Node'}<br>"
        if not node.is_root():
            hover_text += f"**Dist:** {node.dist:.3f}"
        node_coords['text'].append(hover_text)
        
        # Create line segments for branches
        if not node.is_root():
            parent = node.up
            parent_x = getattr(parent, x_attr)
            parent_y = getattr(parent, y_attr)
            
            if layout == 'Linear':
                # 1. Horizontal segment (Child_X to Parent_X at Child_Y)
                branch_lines_x.extend([node_x, parent_x, None])
                branch_lines_y.extend([node_y, node_y, None])
                
                # 2. Vertical segment (Child_Y to Parent_Y at Parent_X)
                branch_lines_x.extend([parent_x, parent_x, None])
                branch_lines_y.extend([node_y, parent_y, None])
            
            elif layout == 'Circular':
                # For Circular, we connect the child directly to the parent
                # The curvature will depend on the projection, but the core 
                # segments are the straight lines between nodes in the new space.
                branch_lines_x.extend([node_x, parent_x, None])
                branch_lines_y.extend([node_y, parent_y, None])
            
    return branch_lines_x, branch_lines_y, node_coords


# dark mode for species and taxa trees
def display_interactive_tree_c_zoom(tree_path, name):
    
    try:
        tree = Tree(tree_path, format=1)
        tree = get_tree_coordinates(tree) 
        num_species = len(tree.get_leaves())
    except Exception as e:
        print(f"Error processing tree file: {e}")
        return

    taxon=name
    # --- 1. Layout Selection ---
    available_layouts = ['Linear']
    MIN_SPECIES_FOR_CIRCULAR = 400
    if num_species >= MIN_SPECIES_FOR_CIRCULAR:
        available_layouts.append('Circular')
        
    col_layout, col_label = st.columns(2)
    with col_layout:
        tree_layout = st.radio("Select Tree Layout:", available_layouts, horizontal=True)
    with col_label:
        show_labels_option = st.radio("Show Labels:", ('Yes', 'No'), index=0, horizontal=True)
    
    layout_key = tree_layout
    show_labels = (show_labels_option == 'Yes')
    
    # --- 2. Coordinate Calculation ---
    try:
        if layout_key == 'Circular':
            tree = polar_coordinates_transform_c(tree) 
        branch_lines_x, branch_lines_y, node_coords = create_plotly_data_c(tree, layout=layout_key)
    except Exception as e:
        print(f"Error generating coordinates: {e}")
        return
        
    # --- 3. Proteomes & Coloring ---
    try:
        #proteomes_path = 'data/tree/selected_proteomes_above_25_formatted_organisms.csv'
        proteomes_path = 'data/core/merged_evo_rate.csv'
        proteomes_df = pd.read_csv(proteomes_path)
        proteomes_df.columns = proteomes_df.columns.str.strip().str.replace(' ', '_')
        proteomes_df['Kingdom'] = proteomes_df['Taxonomic_lineage'].str.split(', ').str[1]
        organism_to_kingdom = proteomes_df.set_index('Organism_formatted')['Kingdom'].to_dict()
        unique_kingdoms = sorted(proteomes_df['Kingdom'].unique())
        
        node_colors = []
        for node_name in node_coords['name']:
            kingdom = organism_to_kingdom.get(node_name, 'Unknown') 
            color = kingdom_color(kingdom)
            node_colors.append(color)
    except:
        node_colors = ['#1f77b4'] * len(node_coords['x'])
        unique_kingdoms = ['Default']
    
    # --- 4. NEW: Find Node Feature ---
    # Get all names
    raw_node_names = list(node_coords['name'])
    
    # FILTER: Keep only names that contain at least one letter (a-z, A-Z)
    # This removes "9606", "58.56", "100", etc.
    filtered_node_names = [
        name for name in raw_node_names 
        if name and any(c.isalpha() for c in str(name))
    ]
    
    # Use the filtered list for the dropdown
    search_node = st.selectbox(
        "🔍 Find a Species:",
        options=["None"] + sorted(filtered_node_names), 
        index=0,
        key="node_finder"
    )

    # --- 5. Define Traces ---
    trace_branches = go.Scatter(
        x=branch_lines_x, y=branch_lines_y,
        mode='lines',
        line=dict(color='lightgrey', width=1),
        hoverinfo='none', 
        name='Branches'
    )
    
    node_text = node_coords['name'] if show_labels else [''] * len(node_coords['x'])

    trace_nodes = go.Scatter(
        x=node_coords['x'], y=node_coords['y'],
        mode='markers+text',
        marker=dict(size=8, color=node_colors, line=dict(width=1, color='white')),
        text=node_text, 
        textposition="middle right",
        textfont=dict(size=14, color='white'),
        hovertext=node_coords['text'],
        hoverinfo='text',
        name='Nodes'
    )

    data_traces = [trace_branches, trace_nodes]

    # --- 6. Add Highlight Trace if Search is Active ---
    if search_node != "None":
        try:
            # Find the index of the selected node
            idx = node_coords['name'].index(search_node)
            target_x = node_coords['x'][idx]
            target_y = node_coords['y'][idx]
            
            # Create a "Highlight" trace (Big Red Star)
            trace_highlight = go.Scatter(
                x=[target_x], y=[target_y],
                mode='markers+text',
                marker=dict(
                    size=20, 
                    color='#FF0000', # Red
                    symbol='star',   # Star shape
                    line=dict(width=2, color='white')
                ),
                text=[search_node],
                textposition="top center",
                textfont=dict(size=16, color='#FF0000', family="Arial Black"),
                name='Selected'
            )
            data_traces.append(trace_highlight)
            
            st.success(f"📍 Found: {search_node}")
            
        except ValueError:
            st.warning("Node not found in coordinates.")

    # --- 7. Legend Traces ---
    for kingdom in unique_kingdoms:
        color = kingdom_color(kingdom)
        data_traces.append(
            go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(size=10, color=color, line=dict(width=1, color='white')),
                name=kingdom,
                legendgroup='kingdoms'
            )
        )

    # --- 8. Layout ---
    fig = go.Figure(data=data_traces)
    
    layout_params = dict(
        title=f'{taxon} Tree ({num_species} Species)', # tree title
        height=800, 
        showlegend=True,
        hovermode="closest",
        template="plotly_dark",     
        plot_bgcolor='black',      
        paper_bgcolor='black',      
        font=dict(color='white'),
        margin=dict(l=20, r=20, t=40, b=20),
        autosize=True
    )
    
    # Adjust layout based on type
    if layout_key == 'Linear':
        layout_params.update(xaxis=dict(showgrid=False, title='Branch length'), yaxis=dict(showgrid=False, showticklabels=False))
    else:
        layout_params.update(xaxis=dict(showgrid=False, showticklabels=False), yaxis=dict(showgrid=False, showticklabels=False))

    fig.update_layout(**layout_params)
    st.plotly_chart(fig, use_container_width=True)


##for homologs
def display_protein_presence_tree_plotly(protein, families_csv_path, new_proteomes_path, main_tree_path, layout='Linear'):
    """
    Generates an interactive Plotly tree (Dark Mode) colored by protein presence.
    """
    
    # --- 1. Load Proteome Mapping (Name -> ID) ---
    try:
        new_proteomes_csv = pd.read_csv(new_proteomes_path)
        new_proteomes_csv.columns = new_proteomes_csv.columns.str.strip().str.replace(' ', '_')
        mapping_df = new_proteomes_csv[['Organism_formatted', 'Proteome_Id']].drop_duplicates()
        name_to_id_map = mapping_df.set_index('Organism_formatted')['Proteome_Id'].to_dict()
    except Exception as e:
        print(f"Error loading proteomes: {e}")
        return

    # --- 2. Determine Presence/Absence ---
    # We use your existing helper function logic here
    try:
        family, family_df = protein_in_family(protein, families_csv_path)
        
        if family is None:
            st.warning(f"Protein {protein} not found in families file.")
            return

        # Define helper if missing
        if 'get_species_id' not in globals():
            get_species_id_func = lambda x: x.split('_')[0] 
        else:
            get_species_id_func = get_species_id

        family_df['proteome_id'] = family_df['protein'].apply(get_species_id_func)
        
        # Filter for the relevant family rows
        target_rows = family_df[family_df['family_id'] == family]
        
        # --- NEW LOGIC START ---
        # Create a dictionary: { 'UP000005640': 'P12345, Q99999' }
        # This groups multiple hits for the same species into one string
        proteome_to_uniprot = target_rows.groupby('proteome_id')['protein'].apply(lambda x: ', '.join(x)).to_dict()
        
        # The set of IDs is simply the keys of this new dict
        present_proteomes = set(proteome_to_uniprot.keys())
        
    except Exception as e:
        print(f"Error determining protein family presence: {e}")

    # --- 3. Generate Tree Coordinates (Plotly) ---
    try:
        tree = Tree(main_tree_path, format=1)
        tree = get_tree_coordinates(tree) 
        if layout == 'Circular':
            tree = polar_coordinates_transform_c(tree) 
        branch_x, branch_y, node_coords = create_plotly_data_c(tree, layout=layout) 
    except Exception as e:
        print(f"Error generating coordinates: {e}")
        st.stop()

    # --- 4. Assign Colors to Nodes ---
    node_colors = []
    hover_texts = []
    present_species_info = [] 
    #present_species_names = []
    
    # Colors matching your original style (but adapted for Dark Mode)
    COLOR_PRESENT = "#a0fc97" # Pale Green
    COLOR_ABSENT  = "#a4a4aa" # Grey
    COLOR_EMPTY   = "rgba(0,0,0,0)"
    
    for node_name in node_coords['name']:
            # 1. Check if the node is a valid species in our mapping
            proteome_id = name_to_id_map.get(node_name)
            
            if proteome_id:
                # It IS a species -> Check Presence vs Absence
                if proteome_id in proteome_to_uniprot:
                    # Case A: Present
                    node_colors.append(COLOR_PRESENT)
                    hit_ids = proteome_to_uniprot[proteome_id]
                    hover_texts.append(f"<b>{node_name}</b><br>Status: Present<br>Ortholog(s): {hit_ids}")
                    present_species_info.append(f"{node_name}\t{hit_ids}")
                else:
                    # Case B: Absent (But it is a known species)
                    node_colors.append(COLOR_ABSENT)
                    hover_texts.append(f"<b>{node_name}</b><br>Status: Absent")
            else:
                # Case C: Not a species name (Internal node, bootstrap, or empty)
                node_colors.append(COLOR_EMPTY)
                if node_name:
                    hover_texts.append(f"Node: {node_name}")
                else:
                    hover_texts.append("")

    # --- 5. Build Plotly Figure (Dark Mode) ---
    
    # Trace 1: Branches
    trace_branches = go.Scatter(
        x=branch_x, y=branch_y,
        mode='lines',
        line=dict(color='lightgrey', width=1), # Light grey for dark bg
        hoverinfo='none',
        name='Branches'
    )

    # Trace 2: Nodes
    # Trace Nodes (Partial Update)
    trace_nodes = go.Scatter(
        x=node_coords['x'], y=node_coords['y'],
        mode='markers', 
        marker=dict(size=8, color=node_colors, line=dict(width=1, color='white')),
        #marker=dict(size=8, color=COLOR_EMPTY, line=dict(width=1, color='white')),
        text=node_coords['name'],
        hovertext=hover_texts,
        hoverinfo='text',
        name='Node'
    )
    
    # Trace 3: Legend Entries (Manual)
    legend_traces = [
        go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=10, color=COLOR_PRESENT, line=dict(width=1, color='white')),
            name=f'Present ({protein})'
        ),
        go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=10, color=COLOR_ABSENT, line=dict(width=1, color='white')),
            name='Absent'
        )
    ]

    # Layout Configuration
    fig = go.Figure(data=[trace_branches, trace_nodes] + legend_traces)
    
    layout_params = dict(
        title=f"Structural orthologs of {protein}",
        height=900,
        showlegend=True,
        hovermode="closest",
        template="plotly_dark",     # Dark Theme
        plot_bgcolor='black',       # Pitch black background
        paper_bgcolor='black',
        font=dict(color='white'),
        margin=dict(l=20, r=20, t=40, b=20),
        autosize=True
    )

    if layout == 'Linear':
        layout_params.update(
            xaxis=dict(showgrid=False, title='Branch length', showticklabels=True),
            yaxis=dict(showgrid=False, showticklabels=False)
        )
    else: # Circular
        layout_params.update(
            xaxis=dict(showgrid=False, showticklabels=False, scaleanchor="y", scaleratio=1),
            yaxis=dict(showgrid=False, showticklabels=False)
        )

    fig.update_layout(**layout_params)

    present_species_names = present_species_info
    
    return fig, present_species_names



# --- Function to generate the Mol* HTML viewer ---
@st.cache_data
def render_3dmol(afdb_download_url: str):
    """
    Fetches the CIF data and renders it using py3Dmol via stmol.
    """
    #st.info("Fetching protein structure from AlphaFold database...")
    cif_data = None
    
    try:
        # 1. Fetch the CIF data
        cif_response = requests.get(afdb_download_url, timeout=30)
        #cif_response.raise_for_status() 
        cif_data = cif_response.text
        
        if not cif_data or not isinstance(cif_data, str):
             st.error("Cannot render structure, entry not found in AFDB.")
             
    except requests.exceptions.RequestException as e:
        st.warning("Cannot render structure, entry not found in AFDB.")
        print(f"Network error during CIF download: {e}")

    # 2. Initialize the py3Dmol viewer object
    # DO NOT call .render()
    view = py3Dmol.view(width=800, height=500, js='https://3dmol.org/build/3Dmol.js')
    view.setBackgroundColor('black')
    # 3. Add the model
    view.addModel(cif_data, 'cif')
    # 4. Set the coloring scheme (lDDT for AlphaFold)
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
    view.zoomTo()
    try:
        showmol(view, height=600, width=600)
    except Exception as e:
         st.warning(f"Cannot render with 3Dmol viewer")
         print(f"Error rendering 3Dmol viewer with stmol: {e}")

##### protein evolution set
def load_mapping(mapping_file):
    try:
        mapping_df = pd.read_csv(mapping_file)
        return dict(zip(mapping_df['species_name'], mapping_df['id']))
    except Exception as e:
        print(f"Error reading mapping file: {e}")
        sys.exit(1)

def find_comparison_csv(target_dir, compare_dir_id):
    """Finds the Human_vs_Species or Species_vs_Human CSV file."""
    search_pattern = os.path.join(target_dir, '**', f'*vs{compare_dir_id}_top_hits.csv')
    csv_files = glob.glob(search_pattern, recursive=True)
    
    if not csv_files:
        search_pattern_alt = os.path.join(target_dir, '**', f'*vs{compare_dir_id}*.csv')
        csv_files = glob.glob(search_pattern_alt, recursive=True)

    return csv_files[0] if csv_files else None

def process_orthologs(csv_file, uniprot_ids, search_col='proteinA'):
    """
    Calculates scores for a comparison species.
    Returns: (Total Score, Presence Map, Score Vector)
    """
    total_score = 0.0
    presence_map = {uid: -1 for uid in uniprot_ids}
    score_map = {uid: 0.0 for uid in uniprot_ids} # Default 0.0 if missing

    try:
        df = pd.read_csv(csv_file, header=None, usecols=[0, 1, 4], on_bad_lines='warn')
        df.columns = ['proteinA', 'proteinB', 'bit_score']
        
        for uniprot_id in uniprot_ids:
            # Search logic (your original feature)
            hits = df[df[search_col].str.contains(f"-{uniprot_id}-", na=False)]
            if hits.empty:
                hits = df[df[search_col].str.contains(uniprot_id, na=False)]
            
            if not hits.empty:
                best_score = hits['bit_score'].max()
                total_score += best_score
                
                # For iTOL
                if best_score > 0:
                    presence_map[uniprot_id] = 1
                
                # For Tree (Store actual score)
                score_map[uniprot_id] = best_score

    except Exception as e:
        print(f"Error processing {csv_file}: {e}")
        
    # Convert map to list (Vector) ensures order is identical for all species
    score_vector = [score_map[uid] for uid in uniprot_ids]
    
    return total_score, presence_map, score_vector

def get_human_self_vector(target_dir, uniprot_ids):
    """
    [NEW FEATURE]
    Calculates the 'Perfect' vector for Human by finding the max self-score
    across all files. This allows Human to be plotted on the tree.
    """
    #st.info("Calculating Human self-scores for phylogenetic tree...")
    
    # Initialize with 0.0
    self_scores_map = {uid: 0.0 for uid in uniprot_ids}
    
    search_pattern = os.path.join(target_dir, '**', f'*vs*_top_hits.csv')
    all_csv_files = glob.glob(search_pattern, recursive=True)
    
    if not all_csv_files:
        return 0.0, [0.0]*len(uniprot_ids)

    # Optimization: We just need to fill the map with the max score found anywhere
    for csv_file in all_csv_files:
        try:
            # Read in chunks to save memory
            df_chunk = pd.read_csv(csv_file, header=None, usecols=[0, 4], on_bad_lines='warn',
                                 chunksize=100000, dtype={0: 'str', 4: 'float'})
            for chunk in df_chunk:
                chunk.columns = ['proteinA', 'bit_score']
                
                for uniprot_id in uniprot_ids:
                    # Check if we already have a high score (heuristic optimization)
                    # but to be safe and match original logic, we scan.
                    hits = chunk[chunk['proteinA'].str.contains(f"-{uniprot_id}-", na=False)]
                    if hits.empty:
                        hits = chunk[chunk['proteinA'].str.contains(uniprot_id, na=False)]
                    
                    if not hits.empty:
                        local_max = hits['bit_score'].max()
                        if local_max > self_scores_map[uniprot_id]:
                            self_scores_map[uniprot_id] = local_max
        except Exception as e:
            print(f"Exception in get_human_self_vector:{e}")
            
    # Create the vector and total sum
    human_vector = [self_scores_map[uid] for uid in uniprot_ids]
    total_self_score = sum(human_vector)
    
    return total_self_score, human_vector

def sum_scores_for_orthologs(csv_file, uniprot_ids, search_col='proteinA'):
    """
    Loads a single CSV file and sums the *best hit* bit scores
    for each target UniProt ID.
    
    [FIX] Added 'search_col' to specify which column to search
    for the UniProt ID (for bi-directional search).
    """
    total_score = 0.0
    try:
        df = pd.read_csv(csv_file, header=None, usecols=[0, 1, 4], on_bad_lines='warn')
        df.columns = ['proteinA', 'proteinB', 'bit_score']
        
        for uniprot_id in uniprot_ids:
            # [FIX] Use the dynamic 'search_col' argument
            # This is 'proteinA' for a forward search
            # and 'proteinB' for a reverse search.
            hits = df[df[search_col].str.contains(f"-{uniprot_id}-", na=False)]
            
            if hits.empty:
                hits = df[df[search_col].str.contains(uniprot_id, na=False)]
            
            if not hits.empty:
                best_score = hits['bit_score'].max()
                total_score += best_score
            #else:
                #st.info(f"  > No hit found for {uniprot_id} in {csv_file}")

    except pd.errors.EmptyDataError:
        print(f"CSV file {csv_file} is empty. Skipping.")
    except Exception as e:
        print(f"Error processing CSV {csv_file}: {e}")
        
    return total_score

def get_self_scores(target_dir, target_dir_id, uniprot_ids):
    """
    Calculates the "self-score" for the target species, which will be
    the 1.0 baseline.
    
    [FIX] This function no longer looks for a "self-CSV".
    Instead, it scans all available CSVs to find the
    absolute max score for each target protein, which serves
    as the robust self-score baseline.
    """
    #st.info(f"Calculating self-scores for target {target_dir_id}...")
    
    # Find ALL comparison CSVs in the target directory
    search_pattern = os.path.join(target_dir, '**', f'*vs*_top_hits.csv')
    all_csv_files = glob.glob(search_pattern, recursive=True)
    
    if not all_csv_files:
        st.error(f"No comparison CSVs found at all in {target_dir}. Cannot calculate self-scores.")
        return 0.0

    total_self_score = 0.0
    
    for uniprot_id in uniprot_ids:
        max_score_for_this_protein = 0.0
        
        # Search all CSVs for the highest score for this one protein
        for csv_file in all_csv_files:
            try:
                # Load only the proteinA and bit_score columns
                df_chunk = pd.read_csv(csv_file, header=None, usecols=[0, 4], on_bad_lines='warn',
                                       chunksize=100000, dtype={0: 'str', 4: 'float'})
                
                for chunk in df_chunk:
                    chunk.columns = ['proteinA', 'bit_score']
                    
                    # Find hits for the UniProt ID
                    hits = chunk[chunk['proteinA'].str.contains(f"-{uniprot_id}-", na=False)]
                    if hits.empty:
                        hits = chunk[chunk['proteinA'].str.contains(uniprot_id, na=False)]
                    
                    if not hits.empty:
                        local_max = hits['bit_score'].max()
                        if local_max > max_score_for_this_protein:
                            max_score_for_this_protein = local_max
                            
            except pd.errors.EmptyDataError:
                # This file is empty, skip it
                continue
            except Exception as e:
                # Log other errors (e.g., file corrupted) but continue
                print(f"Could not fully read {csv_file} for self-score: {e}")
                continue
        
        if max_score_for_this_protein == 0.0:
            st.warning(f"  > No score found for {uniprot_id} in any comparison file.")
        
        #st.info(f"  > Max score (baseline) for {uniprot_id} is: {max_score_for_this_protein}")
        total_self_score += max_score_for_this_protein

    return total_self_score

def generate_distance_matrix(species_scores_dict, uniprot_ids, out_dir):
    """
    Generates dmat.npy and ids_in_dmat.csv for the bootstrap script.
    """
    #st.info("Generating distance matrix ...")
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(species_scores_dict, orient='index', columns=uniprot_ids)
    
    # Normalize: log2(x + 1) to handle the wide range of bit scores
    df_log = np.log2(df + 1)
    
    # Calculate Distance (Euclidean distance between score profiles)
    dists = pdist(df_log.values, metric='euclidean')
    dmat = squareform(dists)
    
    # Save dmat.npy
    dmat_path = os.path.join(out_dir, 'dmat.npy')
    np.save(dmat_path, dmat)
    # FIX 1: Changed 'outdir' to 'out_dir' and ensured correct path construction
    dmat_out = dmat_path
    
    # Save ids_in_dmat.csv
    ids_df = pd.DataFrame({'id': df.index})
    ids_path = os.path.join(out_dir, 'ids_in_dmat.csv')
    ids_df.to_csv(ids_path, index=False)
    # FIX 2: Changed 'outdir' to 'out_dir' and ensured correct path construction
    ids_out = ids_path

    return dmat_out, ids_out


# make tree from distance matrix -------------------------------------------

def bootstrap_from_distance_matrix(distance_matrix):
    """Bootstrap by resampling cells in the distance matrix"""
    n = len(distance_matrix)
    bootstrap_matrix = copy.deepcopy(distance_matrix)

    col_pos = np.random.choice(range(n), n, replace=True)
    sampled_vals = np.random.choice(distance_matrix.flatten(), n, replace=True)
    bootstrap_matrix[range(n), col_pos] = sampled_vals

    #Make symmetric
    bootstrap_matrix = np.triu(bootstrap_matrix)
    bootstrap_matrix += bootstrap_matrix.T

    # Ensure diagonal remains zero
    np.fill_diagonal(bootstrap_matrix, 0)

    return bootstrap_matrix

def load_ids(ids_file):
    """Load sequence IDs from file with flexible format handling"""
    try:
        # First try reading as numpy file
        if ids_file.endswith('.npy'):
            return np.load(ids_file).tolist()

        # Try reading as text file
        with open(ids_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        # Try different formats
        if len(lines) == 0:
            raise ValueError("Empty file")

        # Try tab-delimited format
        if '\t' in lines[0]:
            ids = [line.split('\t')[-1] for line in lines]
        else:
            # Assume one ID per line
            ids = lines

        return ids

    except Exception as e:
        print(f"Error reading IDs file: {e}")
        print("File should contain one ID per line or be tab-delimited")
        sys.exit(1)

def parse_biotite_tree(biotite_tree, ids):
    """Convert biotite tree to Bio.Phylo tree"""
    newick = biotite_tree.to_newick(ids)
    return Phylo.read(StringIO(newick), 'newick')

def analyze_bootstrap_support(original_tree, bootstrap_trees):
    """Calculate bootstrap support values"""
    support_tree = Consensus.get_support(original_tree, bootstrap_trees)


    return support_tree

def get_support_statistics(tree):
    """Calculate statistics about bootstrap support values"""
    support_values = [clade.confidence for clade in tree.find_clades()]
    clade_names = [clade.name for clade in tree.find_clades()]

    support_df = pd.DataFrame()
    support_df['clade'] = clade_names
    support_df['confidence'] = support_values

    if not support_values:
        return None

    support_values =  support_df.confidence.dropna()
    stats = {
        'mean_support': np.mean(support_values),
        'median_support': np.median(support_values),
        'min_support': np.min(support_values),
        'max_support': np.max(support_values),
        'strong_support_nodes': len([v for v in support_values if v >= 70]),
        'moderate_support_nodes': len([v for v in support_values if 50 <= v < 70]),
        'weak_support_nodes': len([v for v in support_values if v < 50]),
        'total_nodes': len(support_values)
    }
    return stats, support_df

@st.cache_data
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
    encoded_image = base64.b64encode(data).decode()

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

@st.cache_data
def get_base64_image(image_path):
    """Converts a local file to a base64 string for HTML embedding."""
    if not os.path.exists(image_path):
        return "" # Return empty if file not found to prevent errors
    with open(image_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode()

#### intro -----------------------------------------------------------------------------------------------------------------------------------
st.set_page_config(
    page_title="SHE: Analyse",
    page_icon= "data/tree/logo_tree.png", 
    layout="wide"
)

#sidebar navigation
with st.sidebar:
    st.title("Browse")
    for item in NAV_ITEMS:
        st.page_link(item["path"], label=item["label"], icon=item["icon"])
st.sidebar.divider()

#set_background('data/logos/interactive.png')
#st.title(':green[SHE: Structural History of Eukarya]')
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
#st.title('Interactive Structural Phylogenetic Tree')
st.markdown("""
    <h2 style='
        text-align: center; 
        font-size: 3rem; 
        font-weight: bold; 
        /* Gradient Background */
        background: linear-gradient(to right, #00FF7F, #40E0D0); 
        /* Clip background to text */
        -webkit-background-clip: text; 
        -webkit-text-fill-color: transparent; 
        background-clip: text; 
        color: transparent;'>
        Interactive Structural Phylogenetic Tree
    </h2>
""", unsafe_allow_html=True)
image_path = 'data/logos/interactive.png'
# 2. Check if file exists and encode
if os.path.exists(image_path):
    img_b64 = get_base64_image(image_path)

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
            <h2>Discover interactive analysis options from the SHE tree</h2>
        </div>
    </div>
    """,
    unsafe_allow_html=True
)

#st.write("Discover our interactive analysis option of the eukaryotic structural tree")
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
row2 = st.columns(3)
titles = ["1. Zoom on Species or Taxon", "2. Protein structural homologs", "3. Structural evolution of a set of proteins"]
# Create the tiles and their content
for col, content in zip(row2, titles):
    #with col:   tile = st.container(height=120)   tile.title(f"{emoji} {content}")
    with col.expander(f"**{content}**"):
        # --- zoom ---
        if content == "1. Zoom on Species or Taxon":
            st.write(':orange[Interested in discovering the new structural phylogenetics of a specific species? Use the "Find Species" option.]')
            st.write(':orange[Interested in discovering the new structural phylogenetics of a taxon? Use the "Find Taxon" option.]')
            st.checkbox('Find Species', key="find_species_box")
            st.checkbox('Find Taxon', key="find_taxon_box")

        # --- Explore Tree Content ---
        elif content == "2. Protein structural homologs":
            st.write(':red[Interested in discovering in which organisms a certain protein can be found? Use the "Find protein structural homologs" option.]')
            st.checkbox('Find Protein Structural Homologs', key="homologs")

        # --- Model Organism Search Content ---
        elif content == "3. Structural evolution of a set of proteins":
            st.write(":red[What about studying the evolution of your proteins of interest? Click on 'Study the structural evolution set of proteins' and start building your tree!]")
            st.checkbox('Study the Structural Evolution of a Set of Proteins', key='set')



# Add sidebar ---------------------------------------------------------------------------------------------------------------------------
#data
df = pd.DataFrame({
    'kingdom': ["Eukarya", "Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"],
    'number_species': [1000000, 10, 20, 30, 40, 50, 70] # change to actual number of species
    })

##############Extract taxonomic group ------------------------------------------------
main_tree_path = 'data/tree/Eukarya/kingdoms_Eukarya_tree.newick'
#proteomes_path = 'data/tree/selected_proteomes_above_25_formatted_organisms.csv'
new_proteomes_path =  'data/core/merged_evo_rate.csv'
proteomes_df = pd.read_csv(new_proteomes_path) # old: proteomes_path

# find species ------------------------------------------------------------------------------------------------------
st.sidebar.subheader('1. Zoom on Taxonomic Group', divider="orange")
with st.sidebar.expander("Help: How to find the structural tree for a species or a taxon"):
    st.write("Select if you want to find:")
    st.markdown("""
    * a species (i.e. Homo sapiens)
    * a taxa (i.e. Primates)
    """)
    st.write("The resulting structural phylogenetic tree is going to be interactive and colored based on 6 subkingdoms.")
    st.write("For the species search, it is possible to choose the number of related neighbor species in the tree by using the appropriate slide.")
    

#if st.checkbox('Find Species'):
if st.session_state.get("find_species_box"):
    st.title('Species subtree:')
    species = ''
    renamed_species = ''
    #
    options = proteomes_df['Organism_formatted'].tolist()
    target_search_term = 'homo_sapiens_human'
    try:
        # Find the index regardless of how it's capitalized in the CSV
        default_index = [x.lower() for x in options].index(target_search_term)
    except ValueError:
        default_index = 0
    # 1. .replace('_', ' ') -> "Homo sapiens Human"
    # 2. .capitalize()      -> "Homo sapiens human" (First letter Big, rest small)
    species = st.selectbox(
        'Which organism would you like to study (i.e: Homo sapiens human)', 
        options, 
        index=default_index, 
        key="target", 
        format_func=lambda x: x.replace('_', ' ').capitalize()
    )
    #species = st.selectbox('Which species (name - Homo_sapiens)?', proteomes_df['Organism_formatted'], key = 'species_selector')
    if species == '':
        st.warning('Please enter species (i.e. Proteome ID - i.e. UP000005640 or Taxon ID - i.e. 9606 or name - Homo_sapiens) to extract a subtree.' \
            'Find all species names in "Show species branch rate in Kindom Eukarya"')
    else:
        what = species
        what = what.split('_')
        what= what[:2]
        what= ' '.join(what)
        st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
        num_neighbors = st.slider("How many close species?", 2, 1541, 10)
        #num_neighbors = 10
        species_out, renamed_species = extract_closest_species(main_tree_path, species, num_neighbors)
        #handle different input types
        if renamed_species != species:
            st.write('Selected Species: ',renamed_species)
            what=renamed_species
            what = what.split('_')
            what= what[:2]
            what= ' '.join(what)
            st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
        #species = renamed_species
        st.write('Close species in the tree share similar protein structures at the proteome level.')
        display_interactive_tree_c_zoom(species_out, renamed_species)


# find taxon --------------------------------------------------------------------------------------------------
#if st.checkbox('Find Taxon'): #it does not handle two new taxa in the same session
if st.session_state.get("find_taxon_box"):
    st.title('Taxa subtree:')
    taxon = st.text_input(
    'Which taxon?', key="taxon")
    st.write('Selected Taxon: ', taxon)
    st.write('Close species in the tree share similar protein structures at the proteome level.')
    taxon_path = pathlib.Path('results/extracted_trees/' + taxon +'_subtree.svg')
    if taxon == '':
        st.warning('Please enter a taxon (i.e. Primates) to extract a subtree')
    else:
        what = taxon
        st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
        taxon_tree = extract_subtree(new_proteomes_path, main_tree_path, taxon)
        display_interactive_tree_c_zoom(taxon_tree, taxon)
            


############### Protein tree extraction -----------------------------------------------------------------------------------------
st.sidebar.subheader('2. Structural Homologs', divider="red")
with st.sidebar.expander("Help: How to Find Protein Structural Homologs"):
    st.write("Digit the protein (UniProt ID) you want to find the orthologs of.")
    st.write("The resulting structural phylogenetic tree will show you all the species sharing a protein structural homolog of your input protein.")

#if st.checkbox('Find Protein Structural Homologs'):
if st.session_state.get("homologs"):
    protein = ''
    st.title('Tree of all organisms sharing a protein structure:')
    protein = st.text_input(
    'Which protein (UniProt ID)?', key="protein")
    st.write('Selected Protein: ', protein)
    if protein == '':
            st.warning('Please enter protein uniprot ID i.e. P05023 or A0A0P8ZLX7')
    else:
        # uniprot link
        st.page_link("https://www.uniprot.org/uniprotkb/" + protein + "/entry", label="UniProt", icon="🌎")
        # --- Configuration AFDB  ---
        try: 
            UNIPROT_ID = protein 
            AFDB_API_URL = (
                "https://alphafold.ebi.ac.uk/api/prediction/" + UNIPROT_ID + "?")
                #f"{protein}?") #https://alphafold.ebi.ac.uk/api/prediction/Q5VSL9
            st.info(f"Loading predicted structure **{UNIPROT_ID}** from AlphaFold DB...")
        except:
            st.warning(f"**{UNIPROT_ID}** not found in AlphaFold DB...")
        # Use requests.get() to query the API and retrieve the JSON metadata
        response = requests.get(AFDB_API_URL, timeout=59)
        # Check for general request errors (e.g., 404, 500)
        #response.raise_for_status() 
        # Parse the JSON response
        data = response.json()
        # The API returns a list of prediction objects (usually one)
        if data and isinstance(data, list) and data[0].get('cifUrl'):
            # Extract the direct download link from the JSON response
            AFDB_DOWNLOAD_URL = data[0]['cifUrl']
            # Render the Mol* viewer
            render_3dmol(AFDB_DOWNLOAD_URL)
            #st.success(f"Structure loaded successfully (AFDB ID: **{AFDB_ID}**).")
            st.caption(f"Source: [AlphaFold DB - {UNIPROT_ID}]({AFDB_DOWNLOAD_URL})")
        else:
            # Handle cases where the API returns a valid response but no structure data
            st.warning(f"Structure data for {UNIPROT_ID} was **not found** in the AFDB.")
        
        # get tree
        layout = st.radio(
            "Select Tree Layout:",
            ('Linear', 'Circular'),
            index=0,
            horizontal=True
        )
        #protein_svg = pathlib.Path('results/protein_trees/' + protein +'_presence_tree_' + layout +'.svg')
        for protein in [protein]:
            families_csv_path = 'data/core/protein_families_238.csv' #i.e. ./data/core/protein_families_238.csv
            fig, present_species = display_protein_presence_tree_plotly(
                protein=protein, 
                families_csv_path=families_csv_path, 
                new_proteomes_path=new_proteomes_path, 
                main_tree_path=main_tree_path, 
                layout=layout
            )
            # Render the returned figure
            if fig:
                st.plotly_chart(fig, use_container_width=True)
                
                # 3. Handle the Species List and Download Button
                if present_species:
                    # Display a success message with the count
                    st.markdown(f"**Found {len(present_species)} species with structural orthologs.**")
                    
                    # Convert list to a string (one species per line)
                    list_data = "\n".join(present_species)
                    
                    # Create the download button
                    st.download_button(
                        label=f"📥 Download list of species with structural orthologs for {protein}",
                        data=list_data,
                        file_name=f"{protein}_structural_orthologs_species.txt",
                        mime="text/plain",
                        help="Download a text file containing the names of all species colored green in the tree."
                    )
                else:
                    st.info("No structural orthologs found in the tree.")


st.sidebar.subheader('3. Protein Set Structural Evolution', divider="red")
with st.sidebar.expander("Help: how to build the Structural Tree from a Set of Proteins"):
    st.write("Input the UniProt ID of the proteins to build the tree from.")
    st.write("A structural phylogenetic tree is going to be build based on your input protein.")

#if st.checkbox('Study the Structural Evolution of a Set of Proteins'):
if st.session_state.get("set"):
    proteins = ''
    st.title('Tree built from your protein structures:')
    st.warning('Please not that this function is still under development and might take some time. Contact us for feedback :)')
    proteins = st.text_input(
    'Which proteins (each UniProt ID separated by a space, i.e: Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616)?', key="proteins")
    st.write('Selected Proteins: ', proteins)
    if proteins == '':
            st.warning('Please enter a list of protein UniProt ID split by a space i.e. "P05023 A0A0P8ZLX7"')
    else:
        random_num = random.randint(100000, 999999)
        user = f"user_{random_num}"
        uniprot_ids=proteins.split()
        target_species='Homo_sapiens_Human'
        organisms_list=[]
        for species in proteomes_df["Organism_formatted"]:
            organisms_list.append(species)
        compare_species=organisms_list #all 1542 species
        root_dir='data/decompressed/' #350 gb
        # --- 1. Load Mapping ---
        mapping='data/core/model/extracted_and_renamed.csv'    #Path to CSV file with columns id,proteome_id,species_name
        species_map = load_mapping(mapping)
        
        # Get the directory ID for the target
        #target_dir_id = species_map[target_species]
        target_dir_id = species_map.get(target_species)
        if not target_dir_id:
            st.error("Target species not found in mapping.")
            sys.exit(1)
        target_dir = os.path.join(root_dir, target_dir_id)
        
        # Storage
        ranking_scores = {}
        itol_data = {}
        tree_data = {} 
        
        # 1. Calculate Human (Target) Data first
        # This ensures Human is in the tree!
        human_score, human_vector = get_human_self_vector(target_dir, uniprot_ids)
        ranking_scores[target_species] = human_score
        tree_data[target_species] = human_vector
        itol_data[target_species] = {uid: 1 for uid in uniprot_ids} # Human always has the proteins

        # 2. Calculate Comparison Species
        for species in compare_species:
            if species not in species_map: continue
            
            compare_dir_id = species_map[species]
            #st.info(f"Processing {species}...")
            
            csv_file = find_comparison_csv(target_dir, compare_dir_id)
            search_col = 'proteinA'
            
            if not csv_file:
                reverse_dir = os.path.join(root_dir, compare_dir_id)
                if os.path.isdir(reverse_dir):
                    csv_file = find_comparison_csv(reverse_dir, target_dir_id)
                    search_col = 'proteinB'

            if csv_file:
                score, presence, vector = process_orthologs(csv_file, uniprot_ids, search_col)
                ranking_scores[species] = score
                itol_data[species] = presence
                tree_data[species] = vector
            else:
                ranking_scores[species] = 0.0
                itol_data[species] = {uid: -1 for uid in uniprot_ids}
                tree_data[species] = [0.0] * len(uniprot_ids)

        df_rank = pd.DataFrame.from_dict(ranking_scores, orient='index', columns=['Score'])
        # Sort by Score in descending order (highest score is Rank 2)
        df_rank_sorted = df_rank.sort_values(by='Score', ascending=False)
        # The target species is guaranteed to be the first row (index 0)
        df_rank = df_rank_sorted.reset_index()
        df_rank.columns = ['Organism', 'Score']
        presence_list = []
        for sp, pres in itol_data.items():
            # Start the dictionary for this row with the organism name
            row_data = {'Organism': sp}
            
            # Iterate through all UniProt IDs
            for uid in uniprot_ids:
                # Get the presence value (default to 0 if not found)
                # Apply the -1/1 mapping: 1 if present, -1 if absent/missing
                presence_val = 1 if pres.get(uid, 0) == 1 else 0 #else 0 or -1
                
                # Add a new key-value pair to the row_data dict: Key=UniProt ID, Value=Presence (-1 or 1)
                row_data[uid] = presence_val
            
            # Append the complete row data to the list
            presence_list.append(row_data)

        df_presence = pd.DataFrame(presence_list)
        # Filter target species out
        target_species_row = df_rank[df_rank['Organism'] == target_species]
        df_rank_filtered = df_rank[df_rank['Organism'] != target_species].copy()
        df_rank_filtered.reset_index(drop=True, inplace=True)
        df_rank_filtered.insert(0, 'Rank', df_rank_filtered.index + 1)
        df_final = pd.merge(df_rank_filtered, df_presence, on='Organism', how='left')
        csv_out = "results/proteins_evolution/" + user + "/"
        if not os.path.exists(csv_out):
            os.makedirs(csv_out) #mkdir
        csv_path = csv_out + "organisms_ranking.csv"
        df_final.to_csv(csv_path, index=False)
        
        # distance matrix---------------------------------------
        matrix_out=csv_out
        distance_matrix, ids = generate_distance_matrix(tree_data, uniprot_ids, matrix_out)
        n_bootstraps=100
        #make tree----------------------------------
        distance_matrix = np.load(distance_matrix)
        ids = pd.read_csv(ids)['id'].values #Read from csv - keep the format consistens
        # Verify number of IDs matches matrix dimensions
        if len(ids) != len(distance_matrix):
            st.error(f"Error: Number of IDs ({len(ids)}) doesn't match matrix dimensions ({len(distance_matrix)})")
            sys.exit(1)
        #print(f"Loaded {len(ids)} IDs and {distance_matrix.shape} distance matrix")
        trees = []
        bootstrap_matrices = []
        #print("Creating original tree...")
        original_tree = neighbor_joining(distance_matrix)
        trees.append(parse_biotite_tree(original_tree, ids))
        #print(f"Creating {n_bootstraps-1} bootstrap trees...")
        for _ in tqdm(range(n_bootstraps-1)):
            dmat = bootstrap_from_distance_matrix(distance_matrix)
            bootstrap_matrices.append(dmat)
            tree = neighbor_joining(dmat)
            trees.append(parse_biotite_tree(tree, ids))

        #print("Calculating bootstrap support values...")
        original_tree = trees[0]
        bootstrap_replicates = trees[1:]
        supported_tree = analyze_bootstrap_support(original_tree, bootstrap_replicates)

        #print("Calculating support statistics...")
        stats, support_df = get_support_statistics(supported_tree)
        outdir=matrix_out
        support_df.to_csv(outdir+'support_per_clade.csv', index=None)
        #print('Saved support metrics to', args.outdir+'support_per_clade.csv')
        #print("Saving results...")
        # Save bootstrap matrices
        np.save(os.path.join(outdir, 'bootstrap_distance_matrices.npy'),
                np.array(bootstrap_matrices))
        # Save supported tree
        Phylo.write(supported_tree,
                    os.path.join(outdir, 'bootstrap_tree.newick'),
                    'newick')
        # Save statistics
        with open(os.path.join(outdir, 'bootstrap_stats.txt'), 'w') as f:
            f.write("Bootstrap Support Statistics:\n")
            for stat, value in stats.items():
                if 'support' in stat:
                    f.write(f"{stat}: {value:.1f}%\n")
                else:
                    f.write(f"{stat}: {value}\n")

            #f.write("\nInterpretation Guide:\n")
            #f.write("- Support values ≥70%: Strong support\n")
            #f.write("- Support values 50-69%: Moderate support\n")
            #f.write("- Support values <50%: Weak support\n")
        # load tree ----------------------------------------------
        evo_tree= outdir + "bootstrap_tree.newick"
        st.info("Displaying evolutionary tree of orthologs...")
        display_interactive_tree(evo_tree)
        # download button
        is_random_user = True
        user_dir = os.path.join("results", "proteins_evolution", user)
        os.makedirs(user_dir, exist_ok=True) 
        target_folder_path=user_dir
        user_name=user.split('@')[0]
        zip_file_name = f"results_{user_name}" 
        output_base_dir = "temp_archives" 
        os.makedirs(output_base_dir, exist_ok=True) 
        shutil.make_archive(
            base_name=os.path.join(output_base_dir, zip_file_name),
            format='zip',
            root_dir=os.path.dirname(target_folder_path), # Start zipping from the parent directory
            base_dir=os.path.basename(target_folder_path)  # Only include the target folder itself
        )
        # Get the full path to the newly created ZIP file and make download availabel
        zip_file_path = os.path.join(output_base_dir, f"{zip_file_name}.zip")
        with open(zip_file_path, "rb") as file:
            st.sidebar.download_button(
                label="Download results now",
                data=file, # The actual byte content of the zip file
                file_name=f"{zip_file_name}.zip", # The name the user will see
                mime="application/zip", # Critical: Set the correct MIME type for zip files
                icon=":material/download:",
            )
        # Clean up the generated zip file after use
        os.remove(zip_file_path)
        if is_random_user:
            shutil.rmtree(user_dir)

# --- CITATION
st.subheader('Cite', divider="rainbow")
st.write('Preprint:  Qiuzhen Li, Daumiller Diandra, Patrick Bryant: <The structural History of Eukarya. bioRxiv ; doi: >')
