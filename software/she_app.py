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

# source ~/.venv/bin/activate 
# streamlit run ~/Documents/GitHub/SHE_app/software/she_app.py
github_path = '/Users/diandra/Documents/GitHub/' #change to your path
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
#@st.cache_data

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
    output_file = github_path + 'SHE_app/results/extracted_trees/' + taxonomic_group +'_subtree.newick' #Path to save the resulting subtree in Newick format.
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
    taxon_tree_path = github_path + 'SHE_app/results/extracted_trees/' + taxon +'_subtree.newick'
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
    newick_file = github_path + 'SHE_app/data/tree/All/kingdoms_All_tree.newick'
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
        proteomes_df = pd.read_csv(github_path + 'SHE_app/data/tree/selected_proteomes_above_25_formatted_organisms.csv')
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
        output_path = github_path + 'SHE_app/results/species_trees/'+ target_species +'_subtree.newick'
        main_tree.write(outfile=output_path, format=1)
        print(f"\nSuccess! Extracted subtree with {len(leaves_to_keep)} species saved to '{output_path}'")
    else:
        print("The resulting subtree has no leaves. No file was written.")
    
    return output_path, target_species

def load_species_tree(species, layout, species_out):
# File paths
    #species_out = github_path + 'SHE/species_trees/' + species +'_subtree.newick'
    # Load the Tree: i.e. tree = Tree(main_tree_file)
    species_tree = Tree(species_out, format=1) #format=0
    svg_path = github_path + 'SHE_app/results/species_trees/'+ species +'_subtree.svg' 

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






#for protein tree color
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
        print(f"Warning: Protein ID '{protein}' not found in mapping data.")
        return None
    
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

    family_df['species_id'] = family_df['protein'].apply(get_species_id)
    family_df.dropna(subset=['species_id'], inplace=True)

    presence_absence = pd.crosstab(family_df['family_id'], family_df['species_id']) > 0
    if family not in presence_absence.index:
        print(f"FATAL ERROR: Family '{family}' not found in the filtered dataset.")
        sys.exit(1)
    family_pattern = presence_absence.loc[family]

     # 3. Load the tree to get all leaf names
    print(f"--- Step 3: Loading Tree from '{main_tree}' ---")
    try:
        t = Tree(main_tree)
    except FileNotFoundError:
        print(f"FATAL ERROR: The tree file '{main_tree}' was not found.")
        sys.exit(1)
    #all_leaf_ids = t.get_leaf_names()
    #species = [leaf.name for leaf in t]
    #print(f"Found {len(all_leaf_ids)} leaves in the tree.")

    # 4. Render tree with colored leaves
    print(f"--- Step 4: Rendering Tree with Colored Leaves ---")
    # Set dashed blue lines in all leaves
    nst1 = NodeStyle()
    nst1["bgcolor"] = "#fb370c" #Present
    nst2 = NodeStyle()
    nst2["bgcolor"] = "#6f74f3" #Not present

    #set colors for each leaf
    #for leaf_id in all_leaf_ids:
    for leaf in t:
        output_name = id_to_name_map.get(leaf, leaf)
        species_id = leaf.split('_')[0].strip()
        is_present = family_pattern.get(species_id, False) #output True or False?
        if is_present:
            leaf.set_style(nst1)
        else:
            leaf.set_style(nst2)

    layout= 'Linear' #default
    presence_tree_svg = github_path + 'SHE_app/results/protein_trees/'+ protein +'_presence_tree_' + layout + '.svg'
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
        t.render(
        presence_tree_svg,
        tree_style=ts,
        w=1000,
        units="px" 
        )
    elif layout == 'Linear':
        ts = TreeStyle()
        ts.mode = "r"   
        ts.show_leaf_name = True #true or false
        ts.show_branch_length = False
        ts.force_topology = True #not bad
        # Render tree
        t.render(
        presence_tree_svg,
        tree_style=ts,
        #w=183,
        units="px"  
        )
    return presence_tree_svg


#### intro -----------------------------------------------------------------------------------------------------------------------------------
st.title('SHE: Structural History of Eukarya')

st.subheader('Strucutral Phylogenetic Tree')
#, main_tree_out = load_main_tree()
# Display the main_tree image
st.image(
    github_path + 'SHE_app/data/tree/All/Figure1.png',
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

# 1. Get the color for the currently selected kingdom
color = kingdom_color(kingdom)

# 2. Inject CSS to style the selectbox based on the current selection.
#    This targets the Streamlit element that displays the *currently selected value*.
#    The ID/class names might change slightly with Streamlit updates.
#    '.stSelectbox div[data-baseweb="select"]' is a common target for the selectbox input area.
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
#st.write(f"You selected **{kingdom}**. The box is now colored: `{color}`")

st.write('Close species in the tree share similar protein structures at the proteome level.')
st.write('Start exploring the new structural tree by selecting a kingdom in the sidebar!')
st.write('You can also choose the tree visualisation style in the sidebar.')
st.write('Interested in discovering the new strucutural phylogenetics of specific species or taxon? Use the "Extract Taxonomic Group" section in the sidebar.')
st.write('Interested in discovering in which organisms a certain protein can be found? Use the "Look for a protein" section in the sidebar.')



#Visualization options
st.sidebar.subheader('Visualization Options')
layout = st.sidebar.selectbox(
    'Tree Layout',
    ['Circular', 'Linear']
)


if st.sidebar.toggle('Labels'):
    labels = True
else:
    labels = False



#show trees
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


# Display the main_tree image
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
        kingdom_tree_out= github_path + 'SHE_app/data/tree/' + kingdom + '/kingdoms_' + kingdom + "_" + layout + "_tree1000_labels.svg"
    else:
        kingdom_tree_out= github_path + 'SHE_app/data/tree/' + kingdom + '/kingdoms_' + kingdom + "_" + layout + "_tree1000.svg"
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
    kingdom_csv = github_path + 'SHE_app/data/tree/' + kingdom + '/kingdoms_' + kingdom + "_branch_rates.csv" #i.e kingdoms_amoebozoa_branch_rates.csv
    kingdom_csv = pd.read_csv(kingdom_csv)
    kingdom_csv
    #chart_data = df
    #chart_data

#Extract taxonomic group ------------------------------------------------
main_tree_path = github_path + 'SHE_app/data/tree/All/kingdoms_All_tree.newick'
proteomes_path = github_path + 'SHE_app/data/tree/selected_proteomes_above_25_formatted_organisms.csv'
new_proteomes_path = github_path + 'SHE_app/data/core/merged_evo_rate.csv'
proteomes_df = pd.read_csv(proteomes_path)


st.sidebar.subheader('Extract Taxonomic Group')
if st.sidebar.checkbox('Find species'):
    st.title('Species subtree:')
    species = ''
    renamed_species = ''
    species = st.sidebar.text_input(
    'Which species (Proteome ID - i.e. UP000005640 or Taxon ID - i.e. 9606 or name - Homo_sapiens)?', key="species").strip() #strip against whitespace
    st.write('Selected Species: ', species)
    st.write('Close species in the tree share similar protein structures at the proteome level.')
    species_svg = pathlib.Path(github_path + 'SHE_app/results/species_trees/' + species +'_subtree.svg')
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
            species_svg2 = pathlib.Path(github_path + 'SHE_app/results/species_trees/' + renamed_species +'_subtree.svg')
            if species_svg2.exists():
                st.image(
                    github_path + 'SHE_app/results/species_trees/'+ renamed_species +'_subtree.svg', #extracted tree svg path
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


if st.sidebar.checkbox('Find taxon'): #it does not handle two new taxa in the same session
    st.title('Taxa subtree:')
    taxon = st.sidebar.text_input(
    'Which taxon?', key="taxon")
    st.write('Selected Taxon: ', taxon)
    st.write('Close species in the tree share similar protein structures at the proteome level.')
    taxon_path = pathlib.Path(github_path + 'SHE_app/results/extracted_trees/' + taxon +'_subtree.svg')
    if taxon == '':
        st.warning('Please enter a taxon (i.e. Primates) to extract a subtree')
    elif taxon_path.exists():
        # Display the main_tree image
        what = taxon
        st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
        st.image(
            github_path + 'SHE_app/results/extracted_trees/'+ taxon +'_subtree.svg', #extracted tree svg path
            caption= taxon + ' Phylogenetic Tree',
            width='stretch')
    else:
        if taxon:
            what = taxon
            st.page_link("https://en.wikipedia.org/wiki/" + what, label="Wiki", icon="🌎")
            extract_subtree(proteomes_path, main_tree_path, taxon)
            #data_load_state = st.text('Retrieving Tree...')
            layout = 'Linear' #default
            taxon_out = github_path + 'SHE_app/results/extracted_trees/' + taxon +'_subtree.svg'
            taxon_out = load_extracted_tree(taxon, layout, taxon_out)
            svg_path = github_path + 'SHE_app/results/extracted_trees/'+ taxon +'_subtree.svg' 
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


## Proteins extraction
st.sidebar.subheader('Look for a protein')
if st.sidebar.checkbox('Find protein'):
    protein = ''
    st.title('Tree of all organisms sharing a protein:')
    protein = st.sidebar.text_input(
    'Which protein?', key="protein")
    st.write('Selected Protein: ', protein)
    protein_svg = pathlib.Path(github_path + 'SHE_app/results/protein_trees/' + protein +'_presence_tree_*.svg')
    if protein == '':
        st.warning('Please enter protein uniprot ID i.e. W1PLW2')
    elif protein_svg.exists():
        st.image(
            protein_svg, #svg file
            caption= protein + ' Phylogenetic Tree; red if present, blue if absent in species proteome',
            width='stretch')
    else:
        families_csv_path = github_path + 'SHE_app/data/core/protein_families_238.csv' #i.e. ./data/core/protein_families_238.csv
        presence_tree_svg = protein_tree_color(protein, families_csv_path, new_proteomes_path, main_tree_path)
        st.image(
            presence_tree_svg, #svg file
            caption= protein + ' Phylogenetic Tree; red if present, blue if absent in species proteome',
            width='stretch')

    
