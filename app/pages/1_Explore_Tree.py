import os
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
import streamlit as st
import pandas as pd
import numpy as np
import time
import csv
from ete3 import Tree, TreeStyle, NodeStyle
import sys
import pathlib
import base64
import textwrap
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
        return "#ffffff"
    
#for find species
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
    newick_file = 'data/tree/All/kingdoms_All_tree.newick'
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

    #target species could be leaf name, Uniprot ID or proteome ID
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
    Find the homologous family of a protein uniprot ID.
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
        print(f"Mapped Protein Uniprot ID '{protein}' to family: {family}")
        return family, families_df
    else:
        st.warning(f"This protein does not seem to have a structural homologous...")

def protein_tree_color(protein, families_csv_path, new_proteomes_path, main_tree):
    """
    Given a protein uniprot ID, find its family and color the tree of all organisms sharing that protein.
    Input: 
    protein: uniprot ID (i.e. P12345)
    families_csv_path: path to the protein families csv file (i.e. ./data/core/protein_families_238.csv)
    new_proteomes_path: path to the updated proteomes csv file (i.e. ./ data/core/merged_evo_rate.csv)
    main_tree: path to the main tree newick file (i.e. ./data/tree/All/kingdoms_All_tree.newick)

    Output:
    presence_tree_svg: path to the svg file of the tree with colored leaves
    """

    # 1 Load the names file and create the name map
    st.info("Generating tree: this might require some seconds, please wait")
    print(f"--- Step 1: Loading Name Map from '{new_proteomes_path}' ---")
    try:
        new_proteomes_csv = pd.read_csv(new_proteomes_path)
        new_proteomes_csv.columns = new_proteomes_csv.columns.str.strip().str.replace(' ', '_')
    except FileNotFoundError:
        print(f"FATAL ERROR: The names file '{new_proteomes_path}' was not found.")
        sys.exit(1)
    if 'ID' not in new_proteomes_csv.columns or 'Organism_formatted' not in new_proteomes_csv.columns:
        print(f"FATAL ERROR: The required columns 'ID' and/or 'Organism_formatted' were not found in '{new_proteomes_csv}'.")
        sys.exit(1)
    mapping_df = new_proteomes_csv[['ID', 'Organism_formatted']].drop_duplicates().set_index('ID')
    id_to_name_map = mapping_df['Organism_formatted'].to_dict()
    print(f"Successfully created a name map for {len(id_to_name_map)} unique IDs.")
    
    # 2 from input protein find family
    family, family_df = protein_in_family(protein, families_csv_path)
    if family is None:
        st.warning(f"Protein uniprot ID '{protein}' not found in mapping data.")
        return
    if len(family_df) == 0:
        print("FATAL ERROR: No proteins with 'UP...' IDs were found in the data file.")
        sys.exit(1)
    
    print(f"--- Step 2: Preparing data for '{family}' from '{families_csv_path}' ---")

    family_df['proteome_id'] = family_df['protein'].apply(get_species_id) #UP000054538_930991-AF-A0A0D0DNZ3-F1-model_v4 -> UP000054538
    family_df.dropna(subset=['proteome_id'], inplace=True)
    #print(family_df['proteome_id'].head())
    #0    UP000244409
    #1    UP000002149
    #Name: proteome_id, dtype: object
    

    presence_absence = pd.crosstab(family_df['family_id'], family_df['proteome_id']) > 0 # if protein present in species proteome True else False
    if family not in presence_absence.index:
        print(f"FATAL ERROR: Family '{family}' not found in the filtered dataset.")
        sys.exit(1)
    family_pattern = presence_absence.loc[family] #series with species_id as index and True/False as values
    #print(family_pattern.head())
    #proteome_id
    #UP000000226    False    
    #UP000000304     True
    #Name: fam_2011, dtype: bool

     # 3. Load the tree to get all leaf names
    print(f"--- Step 3: Loading Tree from '{main_tree}' ---")
    try:
        t = Tree(main_tree)
    except FileNotFoundError:
        print(f"FATAL ERROR: The tree file '{main_tree}' was not found.")
        sys.exit(1)
    all_leaf_ids = t.get_leaf_names()
    print(f"Found {len(all_leaf_ids)} leaves in the tree.")
    #print(all_leaf_ids)
    #species = [leaf.name for leaf in t]

    # 4. Render tree with colored leaves
    print(f"--- Step 4: Rendering Tree with Colored Leaves ---")
    # Set dashed blue lines in all leaves
    nst1 = NodeStyle()
    nst1["bgcolor"] = "#a0fc97" #Present
    nst2 = NodeStyle()
    nst2["bgcolor"] = "#a4a4aa" #Not present

    #set colors for each leaf
    for leaf_id in all_leaf_ids:
        match_row = new_proteomes_csv[new_proteomes_csv['Organism_formatted'] == leaf_id]
        species_id =  match_row["Proteome_Id"].values[0] if not match_row.empty else None
        is_present = family_pattern.get(species_id, False) #default False if species_id not found in family_pattern
        for leaf in t:
            if leaf.name == leaf_id:
                #print(f"Leaf_id: {leaf_id}, Species_id: {species_id}, Matching leaf: {leaf.name}, In family: {is_present}")
                if is_present:
                    leaf.set_style(nst1)
                else:
                    leaf.set_style(nst2)
        
    #Define the output paths for both layouts
    circular_path = 'results/protein_trees/' + protein + '_presence_tree_Circular.svg'
    linear_path = 'results/protein_trees/' + protein + '_presence_tree_Linear.svg'

    # --- Render Circular Tree ---
    # Configure the TreeStyle for Circular layout
    ts_circular = TreeStyle()
    ts_circular.mode = "c"          # Draw tree in circular mode
    ts_circular.rotation = 210
    ts_circular.show_leaf_name = True
    ts_circular.show_branch_length = False
    ts_circular.scale = None
    ts_circular.optimal_scale_level = "full"
    ts_circular.force_topology = True
    # Render Circular tree
    t.render(
        circular_path,
        tree_style=ts_circular,
        w=1000,
        units="px" 
    )

    # --- Render Rectangular (Linear) Tree ---
    # Configure the TreeStyle for Rectangular (Linear) layout
    ts_linear = TreeStyle()
    ts_linear.mode = "r"          # Draw tree in rectangular mode
    ts_linear.show_leaf_name = True
    ts_linear.show_branch_length = False
    ts_linear.force_topology = True
    # Render Rectangular (Linear) tree
    t.render(
        linear_path,
        tree_style=ts_linear,
        w=500,
        units="px"  
    )
    
    return circular_path, linear_path

# color tree based on evorate
def load_evotree(kingdom): 
    # File paths
    evo_rate =  'data/core/merged_evo_rate.csv'
    evorate_df = pd.read_csv(evo_rate)
    evorate_df.columns = evorate_df.columns.str.strip().str.replace(' ', '_')
    tree_path = 'data/tree/'  # i.e. ./SHE/tree/Amoebozoa
    tree_file = tree_path + kingdom + '/kingdoms_' + kingdom+ '_tree.newick' # i.e. ./SHE/tree/Amoebozoa/kingdoms_Amoebozoa_tree.newick
    tree = Tree(tree_file) #format=0


    # define color per kingdom
    nst1 = NodeStyle()
    nst1["bgcolor"] = "#0d0887"
    nst2 = NodeStyle()
    nst2["bgcolor"] = "#4b0091" 
    nst3 = NodeStyle()
    nst3["bgcolor"] = "#7901a1" 
    nst4 = NodeStyle()
    nst4["bgcolor"] = "#a52a9e" 
    nst5 = NodeStyle()
    nst5["bgcolor"] = "#cb4786" 
    nst6 = NodeStyle()
    nst6["bgcolor"] = "#e4715f" 
    nst7 = NodeStyle()
    nst7["bgcolor"] = "#f69c3a" 
    nst8 = NodeStyle()
    nst8["bgcolor"] = "#feca1c" 
    nst9 = NodeStyle()
    nst9["bgcolor"] = "#f5f829" 
    nst10 = NodeStyle()
    nst10["bgcolor"] = "#f0f921" 


    #extract kingdom leaves only
    tree_leaf_names = {leaf.name for leaf in tree.iter_leaves()}
    evorate_kingdom_df = evorate_df[
    evorate_df['Organism_formatted'].isin(tree_leaf_names)
    ].copy() # work on a separate DataFrame for kingdom only

    # Ensure 'Rate' is numeric 
    evorate_kingdom_df['Rate'] = pd.to_numeric(evorate_kingdom_df['Rate'], errors='coerce')
    evorate_kingdom_df.dropna(subset=['Rate'], inplace=True)

    #split df into 10 groups based on evorate
    evorate_kingdom_df['Rate_Group'] = pd.qcut(
        evorate_kingdom_df['Rate'],
        q=10,
        labels=False,
        duplicates='drop' # Handles identical rate values
    ) + 1

    # Convert to a Series mapping organism name to its rate group for quick lookup
    evo_group_map = evorate_kingdom_df.set_index('Organism_formatted')['Rate_Group']

    #select color styles
    styles = [nst1, nst2, nst3, nst4, nst5, nst6, nst7, nst8, nst9, nst10]

    #loop through each leaf to color it based on evorate
    for leaf in tree:
        organism_name = leaf.name
        # Use the map to get the group number directly
        group_num = evo_group_map.get(organism_name)
        if group_num is not None:
            # group_num will be an integer from 1 to 10.
            # Use it to index the styles list (group_num - 1)
            if 1 <= group_num <= 10:
                style_to_apply = styles[group_num - 1]
                leaf.set_style(style_to_apply)
            else:
                print(f"Warning: Unexpected group number {group_num} for leaf '{organism_name}'")
        else:
            # organisms not found in the evorate_df
            print(f"Warning: No matching evorate found for leaf '{organism_name}'")

    # Render circular layout
    tree_out_c= tree_path + kingdom+  '/kingdoms_' + kingdom + "_Circular_tree_evorate.svg" #png or pdf or svg
    ts_c = TreeStyle()
    ts_c.mode = "c"          # Draw tree in circular mode otherwise "p"
    ts_c.rotation=210
    ts_c.show_leaf_name = False #true or false
    ts_c.show_branch_length = False
    ts_c.scale = None #big circles with dotted lines
    ts_c.optimal_scale_level = "full" #mid or full
    ts_c.force_topology = True #not bad
    ts_c.root_opening_factor = 1 #new
    #ts_c.layout_fn = labels
    # Render tree
    tree.render(
    tree_out_c,
    tree_style=ts_c,
    w=1000,
    units="px" 
    )

    #linear layout
    tree_out_r= tree_path + kingdom +  '/kingdoms_' + kingdom + "_Linear_tree_evorate.svg" #png or pdf or svg
    ts_r = TreeStyle()
    ts_r.mode = "r"   
    ts_r.show_leaf_name = True #true or false
    ts_r.show_branch_length = False
    ts_r.force_topology = True #not bad
    ts_r.root_opening_factor = 1 #new
    #ts_r.layout_fn = labels
    # Render tree
    tree.render(
    tree_out_r,
    tree_style=ts_r,
    w=500,
    units="px" 
    )
    
    return tree_out_c, tree_out_r

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
        st.error(f"Error processing tree file: {e}")
        return

    # --- Create the Plotly Figure ---
    proteomes_path = 'data/tree/selected_proteomes_above_25_formatted_organisms.csv'
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
        margin=dict(l=20, r=20, t=40, b=20)
    )
    
    # --- Display in Streamlit ---
    st.plotly_chart(fig, use_container_width=True)

# --- Function to generate the Mol* HTML viewer ---

def render_3dmol(afdb_download_url: str):
    """
    Fetches the CIF data and renders it using py3Dmol via stmol.
    """
    st.info("Fetching protein structure from AlphaFold database...")
    cif_data = None
    
    try:
        # 1. Fetch the CIF data
        cif_response = requests.get(afdb_download_url, timeout=30)
        cif_response.raise_for_status() 
        cif_data = cif_response.text
        
        if not cif_data or not isinstance(cif_data, str):
             st.error("Fetched data is empty or not a string. Cannot render structure.")
             return
             
    except requests.exceptions.RequestException as e:
        st.error(f"Network error during CIF download: {e}")
        return

    # 2. Initialize the py3Dmol viewer object
    # DO NOT call .render()
    view = py3Dmol.view(width=800, height=500, js='https://3dmol.org/build/3Dmol.js')

    # 3. Add the model
    view.addModel(cif_data, 'cif')

    # 4. Set the coloring scheme (lDDT for AlphaFold)
    # Official Hex Codes (approximate):
    # Very High (>90): Dark Blue (#0053D6)
    # High (70-90): Light Blue/Cyan (#64CCEB)
    # Low (50-70): Yellow (#FFDB00)
    # Very Low (<50): Red/Orange (#FF7D45)
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
    #view.render()
    #view.show()
    
    # 5. Render the view using the stmol component
    # The showmol function handles the Streamlit component embedding correctly.
    try:
        showmol(view, height=510, width=800)
    except Exception as e:
         st.error(f"Error rendering 3Dmol viewer with stmol: {e}")

##### protein evolution set
def load_mapping(mapping_file):
    try:
        mapping_df = pd.read_csv(mapping_file)
        return dict(zip(mapping_df['species_name'], mapping_df['id']))
    except Exception as e:
        st.error(f"Error reading mapping file: {e}")
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
        st.warning(f"Error processing {csv_file}: {e}")
        
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
        except Exception:
            continue
            
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
        st.warning(f"CSV file {csv_file} is empty. Skipping.")
    except Exception as e:
        st.error(f"Error processing CSV {csv_file}: {e}")
        
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
                st.warning(f"Could not fully read {csv_file} for self-score: {e}")
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



#### intro -----------------------------------------------------------------------------------------------------------------------------------
st.set_page_config(
    page_title="SHE: Explore",
    page_icon= "data/tree/logo_tree.png", 
    layout="centered"
)
st.title('SHE: Structural History of Eukarya')

st.subheader('Structural Phylogenetic Tree')
#, main_tree_out = load_main_tree()
# Display the main_tree image
st.image(
    'data/tree/precomputed.png',
    caption='Strucutral Phylogenetic Tree of Eukarya'
)


#data
df = pd.DataFrame({
    'kingdom': ["All", "Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"],
    'number_species': [1000000, 10, 20, 30, 40, 50, 70] # change to actual number of species
    })



# Add sidebar ---------------------------------------------------------------------------------------------------------------------------
# Choose Taxon
st.sidebar.subheader('1. Kingdoms', divider=True)
kingdom = st.sidebar.selectbox(
    'Which kingdom? (default: All)',
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

st.write('Close species in the tree share similar protein structures at the proteome level.')
st.write('Start exploring the new structural tree by selecting a kingdom in the sidebar!')
st.write('You can also choose the tree visualisation style in the sidebar.')
st.write('Interested in discovering the new structural phylogenetics of specific species or taxa? Use the "Extract Taxonomic Group" section in the sidebar.')
st.write('Interested in discovering in which organisms a certain protein can be found? Use the "Look for a protein" section in the sidebar.')
st.write("What about studying the evolution of your proteins of interest? Click on 'Protein Evolution > Study the structural evolution set of proteins' in the sidebar and start building your tree!")


############################Visualization options --------------------------------------------------------
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
    elif kingdom == "All":
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
            tree_c='data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Circular_tree_evorate.svg"
            st.image(
                tree_c,
                caption= 'Circular Tree with evolutionary rate; color legend: yellow high, purple low',
                width='stretch')
        else:
            tree_r = 'data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Linear_tree_evorate.svg"
            st.image(
                tree_r,
                caption= 'Linear Tree with evolutionary rate; color legend: yellow high, purple low',
                width='stretch')
    else:
        tree_c, tree_r = load_evotree(kingdom)
        if layout == 'Circular':
            st.image(
                tree_c,
                caption= 'Circular Tree with evolutionary rate; color legend: yellow high, purple low',
                width='stretch')
        else:
            st.image(
                tree_r,
                caption= 'Linear Tree with evolutionary rate; color legend: yellow high, purple low',
                width='stretch')


##############Extract taxonomic group ------------------------------------------------
main_tree_path = 'data/tree/All/kingdoms_All_tree.newick'
proteomes_path = 'data/tree/selected_proteomes_above_25_formatted_organisms.csv'
new_proteomes_path =  'data/core/merged_evo_rate.csv'
proteomes_df = pd.read_csv(proteomes_path)

# find species ------------------------------------------------------------------------------------------------------
st.sidebar.subheader('3. Zoom on Taxonomic Group', divider="orange")
if st.sidebar.checkbox('Find species'):
    st.title('Species subtree:')
    species = ''
    renamed_species = ''
    species = st.sidebar.selectbox(
    'Which species (name - Homo_sapiens)?', proteomes_df['Organism_formatted'], key = 'species_selector')
    if species == '':
        st.warning('Please enter species (i.e. Proteome ID - i.e. UP000005640 or Taxon ID - i.e. 9606 or name - Homo_sapiens) to extract a subtree.' \
            'Find all species names in "Show species branch rate in Kindom All"')
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
        display_interactive_tree(species_out)


# find taxon --------------------------------------------------------------------------------------------------
if st.sidebar.checkbox('Find taxon'): #it does not handle two new taxa in the same session
    st.title('Taxa subtree:')
    taxon = st.sidebar.text_input(
    'Which taxon?', key="taxon")
    st.write('Selected Taxon: ', taxon)
    st.write('Close species in the tree share similar protein structures at the proteome level.')
    taxon_path = pathlib.Path('results/extracted_trees/' + taxon +'_subtree.svg')
    if taxon == '':
        st.warning('Please enter a taxon (i.e. Primates) to extract a subtree')
    else:
        what = taxon
        st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
        taxon_tree = extract_subtree(proteomes_path, main_tree_path, taxon)
        display_interactive_tree(taxon_tree)
            


############### Protein tree extraction -----------------------------------------------------------------------------------------
st.sidebar.subheader('Protein Evolution', divider="red")
if st.sidebar.checkbox('4.Find protein structural homologous'):
    protein = ''
    st.title('Tree of all organisms sharing a protein structure:')
    protein = st.sidebar.text_input(
    'Which protein (Uniprot ID)?', key="protein")
    st.write('Selected Protein: ', protein)
    if protein == '':
            st.warning('Please enter protein uniprot ID i.e. P05023 or A0A0P8ZLX7')
    else:
        # uniprot link
        st.page_link("https://www.uniprot.org/uniprotkb/" + protein + "/entry", label="UniProt", icon="🌎")
        #AFDB link
        #st.page_link("https://alphafold.ebi.ac.uk/search/text/" + protein, label="AFDB", icon="🌎")
        # --- Configuration AFDB  ---
        UNIPROT_ID = protein 
        AFDB_API_URL = (
            "https://alphafold.ebi.ac.uk/api/prediction/" + UNIPROT_ID + "?")
            #f"{protein}?") #https://alphafold.ebi.ac.uk/api/prediction/Q5VSL9
        st.markdown(f"Loading predicted structure **{UNIPROT_ID}** from AlphaFold DB...")

        # Use requests.get() to query the API and retrieve the JSON metadata
        response = requests.get(AFDB_API_URL, timeout=59)
        # Check for general request errors (e.g., 404, 500)
        response.raise_for_status() 
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
            st.warning(f"Error: Structure data for {UNIPROT_ID} was **not found** in the API response. Please make sure to input a protein Uniprot ID.")
        
        # get tree
        protein_svg = pathlib.Path('results/protein_trees/' + protein +'_presence_tree_' + layout +'.svg')
        for protein in [protein]:
            #if protein == '':
                #st.warning('Please enter protein uniprot ID without spaces.')
            #elif protein_svg.exists():
            if protein_svg.exists():
                st.image(
                    protein_svg, #svg file
                    caption= protein + ' Structural Tree: green if protein homologous structures are present in species proteome, grey if absent',
                    width='stretch')
                protein = ''
            else:
                families_csv_path = 'data/core/protein_families_238.csv' #i.e. ./data/core/protein_families_238.csv
                circular_svg, linear_svg = protein_tree_color(protein, families_csv_path, new_proteomes_path, main_tree_path)
                if os.path.exists(circular_svg) and os.path.getsize(linear_svg) > 1000:
                    if layout == 'Circular':
                        st.image(
                            circular_svg, #svg file
                            caption= protein + ' Structural Tree: green if protein homologous structures are present in species proteome, grey if absent',
                            width='stretch')
                        protein = ''
                    elif layout == 'Linear':
                        st.image(
                            linear_svg, #svg file
                            caption= protein + ' Structural Tree: green if protein homologous structures are present in species proteome, grey if absent',
                            width='stretch')
                        protein = ''
                else:
                    st.warning("QApplication is creating an empty svg. Restart the app")

if st.sidebar.checkbox('5.Study the structural evolution set of proteins'):
    proteins = ''
    st.title('Tree built from your protein structures:')
    st.warning('This function is going to be released soon')
    proteins = st.sidebar.text_input(
    'Which proteins (each Uniprot ID separated by a space, i.e: Q9BXW9 Q9NXK8 Q9Y2M0 Q9NVI1 Q13616)?', key="proteins")
    st.write('Selected Proteins: ', proteins)
    if proteins == '':
            st.warning('Please enter a list of protein Uniprot ID split by a space i.e. "P05023 A0A0P8ZLX7"')
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
        
        # distance matrix---------------------------------------
        matrix_out="results/proteins_evolution/" + user + "/"
        if not os.path.exists(matrix_out):
            os.makedirs(matrix_out) #mkdir
        distance_matrix, ids = generate_distance_matrix(tree_data, uniprot_ids, matrix_out)
        n_bootstraps=10
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