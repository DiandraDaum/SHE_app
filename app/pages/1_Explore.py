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
import requests

#github_path = '/' #change to your path
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
@st.cache_data
# for coloring the selectbox
def kingdom_color(kingdom):
    if kingdom == 'Opisthokonta':
        return "#bfdbfe" #Opisthokonta
    elif kingdom == 'Viridiplantae':
        return "#bbf7d0" #Viridiplantae
    elif kingdom == 'Amoebozoa':
        return "#fecdd3" #Amoebozoa
    elif kingdom == 'Sar':
        return "#f3f4f6" #Sar
    elif kingdom == 'Discoba':
        return "#e9d5ff" #Discoba
    elif kingdom == 'Rhodophyta':
        return "#c3f8f1" #Rhodophyta
    elif kingdom == 'Cryptophyceae':
        return "#f6fad5" #Cryptophyceae
    elif kingdom == 'Apusozoa':
        return "#ddd6fe" #Apusozoa
    else:
        return "#ffffff" #default white
    
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


# --- Step 2: Flexible Target Input Handling ---

    target_species = input_species # Assume it is the species name by default

    # Check if the input looks like an ID (UniProt starts with 'UP', Organism ID is digit-only)
    if input_species.startswith('UP'):
        is_proteome_id = True
        is_organism_id = False
    elif [(input_species.isdigit()) and (len(input_species) > 3) ]:
        is_organism_id = True
        is_proteome_id = False
    else:
        print("input_species is already a name")
    
    # Only perform lookup if we suspect the input is an ID AND the mapping data is available
    if [(is_proteome_id == True) | (is_organism_id == True)]:
        #new here
        #load df
        #proteomes_df = pd.DataFrame()
        # Load the mapping DataFrame
        proteomes_df = pd.read_csv('data/tree/selected_proteomes_above_25_formatted_organisms.csv')
        proteomes_df.columns = proteomes_df.columns.str.strip().str.replace(' ', '_')
        
        # Determine the column to search based on the suspected ID type
        if is_proteome_id:
            search_col = "Proteome_Id"
            # Ensure Organism_id in the DataFrame is treated as string for accurate match
            proteomes_df["Proteome_Id"] = proteomes_df["Proteome_Id"].astype(str)
            target_value = input_species
        elif is_organism_id:
            search_col = "Organism_Id"
            # Ensure Organism_id in the DataFrame is treated as string for accurate match
            proteomes_df["Organism_Id"] = proteomes_df["Organism_Id"].astype(str)
            target_value = input_species
        else:
            # Should not happen
            search_col = None
            
        if search_col:
            # Look for a direct match in the identified column
            match_row = proteomes_df[proteomes_df[search_col] == target_value]
            
            if not match_row.empty:
                # Strip whitespace from the resulting name just in case
                target_species = match_row["Organism_formatted"].iloc[0].strip()
                print(f"Mapped ID '{input_species}' ({search_col}) to species: {target_species}")
            else:
                print(f"Warning: ID '{input_species}' not found in mapping data. Treating input as species name.")
                return
        else:
            print(f"Input '{input_species}' is treated directly as the target species name.")

    else:
        print(f"Input '{input_species}' is treated directly as the target species name.")
    
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
        st.warning(f"Sorry, your Protein uniprot ID '{protein}' was not found. Try another one please")
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


# --- Function to generate the Mol* HTML viewer ---
def get_molstar_html(afdb_id: str, download_url: str) -> str:
    """Generates the HTML/JavaScript code to embed the Mol* viewer 
    and load the AlphaFold structure directly from the EBI URL."""
    
    # Ensure this script URL is correct and accessible.
    MOLSTAR_SCRIPT = "https://www.ebi.ac.uk/pdbe/pdb-component-library/js/pdbe-molstar-plugin-3.2.0.js"
    
    # We use Mol*'s simple Viewer component, setting the molecule ID and its URL
    return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Mol* Viewer</title>
        <style>
            #molstar {{
                width: 100%;
                height: 500px; /* Adjust height as needed */
            }}
        </style>
        <script src="{MOLSTAR_SCRIPT}"></script>
    </head>
    <body>
        <div id="molstar"></div>
        <script>
            // Initialize the Mol* plugin
            var viewerInstance = new PDBeMolstarPlugin();

            // Configuration for the viewer
            var options = {{
                customData: {{
                    // Check the download URL insertion
                    url: '{download_url}', 
                    format: 'mmcif',      
                    binary: false         
                }},
                // Check the ID insertion
                moleculeId: '{afdb_id}',
                
                // Display options
                hideControls: true, 
                bgColor: {{ r: 255, g: 255, b: 255 }}, 
                
                customPreset: {{
                    name: 'AlphaFold View',
                    residuePred: true
                }}
            }};

            // Render the viewer into the 'molstar' div
            viewerInstance.render('molstar', options);
        </script>
    </body>
    </html>
    """




#### intro -----------------------------------------------------------------------------------------------------------------------------------

st.title('SHE: Structural History of Eukarya')

st.subheader('Strucutral Phylogenetic Tree')
#, main_tree_out = load_main_tree()
# Display the main_tree image
st.image(
    'data/tree/All/Intro.png',
    caption='Strucutral Phylogenetic Tree of Eukarya'
)


#data
df = pd.DataFrame({
    'kingdom': ["All", "Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"],
    'number_species': [1000000, 10, 20, 30, 40, 50, 70] # change to actual number of species
    })



# Add sidebar ---------------------------------------------------------------------------------------------------------------------------
# Choose Taxon
st.sidebar.subheader('Kingdoms')
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
st.write('Interested in discovering the new strucutural phylogenetics of specific species or taxa? Use the "Extract Taxonomic Group" section in the sidebar.')
st.write('Interested in discovering in which organisms a certain protein can be found? Use the "Look for a protein" section in the sidebar.')


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
if st.sidebar.toggle('See evolutionary rate'):
    st.title('Evolutionary rate:')
    'Selected Kingdom: ', kingdom
    'Selected layout: ', layout
    tree_svg = pathlib.Path('data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Linear_tree_evorate.svg")
    if tree_svg.exists():
        if layout == 'Circular':
            tree_c='data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Circular_tree_evorate.svg"
            st.image(
                tree_c,
                caption= 'Circular Tree with evolutionary rate',
                width='stretch')
        else:
            tree_r = 'data/tree/' + kingdom +  '/kingdoms_' + kingdom + "_Linear_tree_evorate.svg"
            st.image(
                tree_r,
                caption= 'Linear Tree with evolutionary rate',
                width='stretch')
    else:
        tree_c, tree_r = load_evotree(kingdom)
        if layout == 'Circular':
            st.image(
                tree_c,
                caption= 'Circular Tree with evolutionary rate',
                width='stretch')
        else:
            st.image(
                tree_r,
                caption= 'Linear Tree with evolutionary rate',
                width='stretch')


##############Extract taxonomic group ------------------------------------------------
main_tree_path = 'data/tree/All/kingdoms_All_tree.newick'
proteomes_path = 'data/tree/selected_proteomes_above_25_formatted_organisms.csv'
new_proteomes_path =  'data/core/merged_evo_rate.csv'
proteomes_df = pd.read_csv(proteomes_path)

# find species ------------------------------------------------------------------------------------------------------
st.sidebar.subheader('Zoom on Taxonomic Group')
if st.sidebar.checkbox('Find species'):
    st.title('Species subtree:')
    species = ''
    renamed_species = ''
    species = st.sidebar.selectbox(
    'Which species (Proteome ID - i.e. UP000005640 or Taxon ID - i.e. 9606 or name - Homo_sapiens)?',
    proteomes_df['Organism_formatted'],
    key = 'species_selector')
    #species = st.sidebar.text_input('Which species (Proteome ID - i.e. UP000005640 or Taxon ID - i.e. 9606 or name - Homo_sapiens)?', key="species").strip() #strip against whitespace
    st.write('Selected Species: ', species)
    st.write('Close species in the tree share similar protein structures at the proteome level.')
    species_svg = pathlib.Path('results/species_trees/' + species +'_subtree.svg')
    if species == '':
        st.warning('Please enter species (i.e. Proteome ID - i.e. UP000005640 or Taxon ID - i.e. 9606 or name - Homo_sapiens) to extract a subtree.' \
            'Find all species names in "Show species branch rate in Kindom All"')
    elif species_svg.exists():
        # Display
        what = species
        what = what.split('_')
        what= what[:2]
        what= ' '.join(what)
        st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
        st.image(
            species_svg, #extracted tree svg path
            caption= species + ' Phylogenetic Tree',
            width='stretch')
        species = ''
    else:
        with st.spinner(f'Extracting and drawing tree for {species}...'):
            # Call extraction function
            num_neighbors = 10
            species_out, renamed_species = extract_closest_species(main_tree_path, species, num_neighbors)
            #handle different input types
            if renamed_species != species:
                st.write('Selected Species: ',renamed_species)
            what=renamed_species
            what = what.split('_')
            what= what[:2]
            what= ' '.join(what)
            st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
            species_svg2 = pathlib.Path('results/species_trees/' + renamed_species +'_subtree.svg')
            if species_svg2.exists():
                st.image(
                    'results/species_trees/'+ renamed_species +'_subtree.svg', #extracted tree svg path
                    caption= renamed_species + ' Phylogenetic Tree',
                    width='stretch')
                species = ''
                renamed_species = ''
                species_out = ''
            else:
                #species = renamed_species
                if species_out != '':
                    layout = 'Linear' #default
                    svg_path = load_species_tree(renamed_species, layout, species_out) #svg_path = './results/species_trees/'+ species +'_subtree.svg' 
                    # Wait up to 30 seconds for the file to be non-empty (adjust as needed)
                    timeout = 30
                    start_time = time.time()
                    # Display tree
                    if svg_path:
                        if os.path.exists(svg_path) and os.path.getsize(svg_path) > 800: #bytes: 1 svg is usually 24 KB
                            st.image(
                                svg_path, #extracted tree svg path
                                caption= renamed_species + ' Phylogenetic Tree',
                                width='stretch')
                            species = ''
                            renamed_species = ''
                            species_out = ''
                        else:
                            st.write(f'tree is ok, but svg is empty: species {species}, renamed_species {renamed_species}')
                            svg_path2 = load_species_tree(renamed_species, layout, species_out)
                            svg_path = svg_path2
                            st.image(
                                svg_path, #extracted tree svg path
                                caption= renamed_species + ' Phylogenetic Tree',
                                width='stretch')
                            species = ''
                            renamed_species = ''
                            species_out = ''

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
    elif taxon_path.exists():
        # Display the main_tree image
        what = taxon
        st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
        st.image(
            'results/extracted_trees/'+ taxon +'_subtree.svg', #extracted tree svg path
            caption= taxon + ' Phylogenetic Tree',
            width='stretch')
    else:
        if taxon:
            what = taxon
            st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
            extract_subtree(proteomes_path, main_tree_path, taxon)
            #data_load_state = st.text('Retrieving Tree...')
            layout = 'Linear' #default
            taxon_out = 'results/extracted_trees/' + taxon +'_subtree.svg'
            taxon_out = load_extracted_tree(taxon, layout, taxon_out)
            svg_path = 'results/extracted_trees/'+ taxon +'_subtree.svg' 
            # Wait up to 30 seconds for the file to be non-empty (adjust as needed)
            timeout = 60
            start_time = time.time()
            # Wait until the file exists AND has content
            while not (os.path.exists(svg_path) and os.path.getsize(svg_path) > 10) and (time.time() - start_time) < timeout:
                time.sleep(1)
            # Display the main_tree image
            if taxon_out: #only display if the function worked
                if os.path.exists(svg_path) and os.path.getsize(svg_path) > 1000:
                    st.image(
                        taxon_out, #extracted tree svg path
                        caption= taxon + ' Phylogenetic Tree',
                        width='stretch')


############### Protein tree extraction -----------------------------------------------------------------------------------------
st.sidebar.subheader('Look for a protein')
if st.sidebar.checkbox('Find protein'):
    protein = ''
    st.title('Tree of all organisms sharing a protein:')
    protein = st.sidebar.text_input(
    'Which protein?', key="protein")
    st.write('Selected Protein: ', protein)
    # uniprot link
    st.page_link("https://www.uniprot.org/uniprotkb/" + protein + "/entry", label="UniProt", icon="🌎")
    #AFDB link
    st.page_link("https://alphafold.ebi.ac.uk/search/text/" + protein, label="AFDB", icon="🌎")
    # --- Configuration AFDB  ---
    # The protein structure can be identified by its UniProt Accession ID (e.g., HRAS)
    UNIPROT_ID = protein 
    # The AlphaFold DB ID format is typically AF-{UNIPROT_ID}-F1 (Fragment 1)
    #AFDB_ID = f"AF-{UNIPROT_ID}-F1-v6" 

    # The official EMBL-EBI download URL for the mmCIF file
    #AFDB_DOWNLOAD_URL = (
        #"https://alphafold.ebi.ac.uk/entry/"
        #f"{AFDB_ID}.cif")

    AFDB_API_URL = (
        "https://alphafold.ebi.ac.uk/api/prediction/"
        f"{protein}?") #https://alphafold.ebi.ac.uk/api/prediction/Q5VSL9
    st.markdown(f"Loading predicted structure **{protein}** from AlphaFold DB...")

    # Use requests.get() to query the API and retrieve the JSON metadata
    response = requests.get(AFDB_API_URL, timeout=10)
    # Check for general request errors (e.g., 404, 500)
    response.raise_for_status() 
    # Parse the JSON response
    data = response.json()
    # The API returns a list of prediction objects (usually one)
    if data and isinstance(data, list) and data[0].get('cifUrl'):
        
        # Extract the direct download link from the JSON response
        AFDB_DOWNLOAD_URL = data[0]['cifUrl']
        
        # Extract the AFDB ID (e.g., AF-P05023-F1-v4) for the Mol* viewer
        # We can extract the ID from the download URL's path
        # URL looks like: https://alphafold.ebi.ac.uk/files/AF-P05023-F1-v4.cif
        import os 
        from urllib.parse import urlparse
        
        parsed_url = urlparse(AFDB_DOWNLOAD_URL)
        # Gets 'AF-P05023-F1-v4.cif' and removes the extension
        AFDB_ID = os.path.splitext(os.path.basename(parsed_url.path))[0]
        
        # Render the Mol* viewer
        molstar_html = get_molstar_html(AFDB_ID, AFDB_DOWNLOAD_URL)
        
        st.components.v1.html(molstar_html, height=550)
        
        st.success(f"Structure loaded successfully (AFDB ID: **{AFDB_ID}**).")
        st.caption(f"Source: [AlphaFold DB - {UNIPROT_ID}]({AFDB_DOWNLOAD_URL})")

    else:
        # Handle cases where the API returns a valid response but no structure data
        st.error(f"Error: Structure data for {UNIPROT_ID} was **not found** in the API response.")
        

    # get tree
    protein_svg = pathlib.Path('results/protein_trees/' + protein +'_presence_tree_' + layout +'.svg')
    for protein in [protein]:
        if protein == '':
            st.warning('Please enter protein uniprot ID i.e. P05023 or A0A0P8ZLX7')
        elif protein_svg.exists():
            st.image(
                protein_svg, #svg file
                caption= protein + ' Phylogenetic Tree; green if present, grey if absent in species proteome',
                width='stretch')
            protein = ''
        else:
            families_csv_path = 'data/core/protein_families_238.csv' #i.e. ./data/core/protein_families_238.csv
            circular_svg, linear_svg = protein_tree_color(protein, families_csv_path, new_proteomes_path, main_tree_path)
            if layout == 'Circular':
                st.image(
                    circular_svg, #svg file
                    caption= protein + ' Phylogenetic Tree; green if present, grey if absent in species proteome',
                    width='stretch')
                protein = ''
            elif layout == 'Linear':
                st.image(
                    linear_svg, #svg file
                    caption= protein + ' Phylogenetic Tree; green if present, grey if absent in species proteome',
                    width='stretch')
                protein = ''


    
