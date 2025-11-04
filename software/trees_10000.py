import pandas as pd
import numpy as np
import time
import csv
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, TextFace, LegendFace, SquareFace
#from phytreeviz import TreeViz
import sys
import pathlib
import base64
import textwrap

# source .venv/bin/activate
# python ./SHE/app/trees_1000.py
github_path = '/Users/diandra/Documents/GitHub/'
#FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------
#@st.cache_data
def load_main_tree():
    # File paths
    main_path = './SHE/tree/All/'
    main_tree_file = main_path + 'kingdoms_All_tree.newick'
    main_tree_out= main_path + "main_tree.svg" #png or pdf or svg

    # Load the Tree: i.e. tree = Tree(main_tree_file)
    main_tree = Tree(main_tree_file) #format=0
    #main_tree = TreeViz(main_tree_file) #not working

    # Configure the TreeStyle
    ts = TreeStyle()
    ts.mode = "c"          # Draw tree in circular mode otherwise "r"
    ts.rotation=210
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.scale = None #big circles with dotted lines
    ts.optimal_scale_level = "full" #mid or full
    ts.force_topology = True #not bad
    #ts.complete_branch_lines_when_necessary = True

    # Render tree
    main_tree.render(
        main_tree_out,
        tree_style=ts,
        w=2000,
        #w=800,
        #h=800,
        units="px"  # Ensure units are pixels for direct control, or w=183, units="mm"
    )
    
    return main_tree, main_tree_out



def labels(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=50) #50 for main tree, 10 for subtrees
        faces.add_face_to_node(N, node, 0, position="aligned") #main tree
        #faces.add_face_to_node(N, node, 0, position="branch-right") #subtrees

def load_main_tree_color():
    # File paths
    main_path = './SHE/tree/All/'
    main_tree_file = main_path + 'kingdoms_All_tree.newick'
    #main_tree_out= main_path + 'kingdoms_All_Linear_tree1000.svg' #png or pdf or svg
    main_tree_out= main_path + 'main_tree_color_legend.svg'

    # Set dashed blue lines in all leaves
    nst1 = NodeStyle()
    nst1["bgcolor"] = "#bfdbfe" #Opisthokonta
    nst2 = NodeStyle()
    nst2["bgcolor"] = "#bbf7d0" #Viridiplantae
    nst3 = NodeStyle()
    nst3["bgcolor"] = "#fecdd3" #Amoebozoa
    nst4 = NodeStyle()
    nst4["bgcolor"] = "#f3f4f6" #Sar
    nst5 = NodeStyle()
    nst5["bgcolor"] = "#e9d5ff" #Discoba
    nst6 = NodeStyle()
    nst6["bgcolor"] = "#c3f8f1" #Rhodophyta
    nst7 = NodeStyle()
    nst7["bgcolor"] = "#f6fad5" #Cryptophyceae
    nst8 = NodeStyle()
    nst8["bgcolor"] = "#ddd6fe" #Apusozoa

    proteomes_path = './SHE/tree/selected_proteomes_above_25_formatted_organisms.csv'
    proteomes_df = pd.read_csv(proteomes_path)
    proteomes_df.columns = proteomes_df.columns.str.strip().str.replace(' ', '_')
    # Create a new 'Kingdom' column
    proteomes_df['Kingdom'] = proteomes_df['Taxonomic_lineage'].str.split(', ').str[1]
    print(proteomes_df['Kingdom'].value_counts()) #1542 ok
    #proteomes_df['Organism_formatted2'] = proteomes_df['Organism_formatted'].str.replace('_', ' ')

    # Load the Tree: i.e. tree = Tree(main_tree_file)
    main_tree = Tree(main_tree_file) #format=0
    species = [leaf.name for leaf in main_tree]
    print(f"Number of species in the main tree: {len(species)}") #
    #print(leaf.name for leaf in main_tree)

    #set colors
    for leaf in main_tree:
        proteome_row = proteomes_df[proteomes_df['Organism_formatted'] == leaf.name]
        #print(f"Leaf: {leaf.name},  \n Proteome Row: {proteome_row}") #ok
        if not proteome_row.empty:
            kingdom = proteome_row['Kingdom'].values[0]
            #print(f"Leaf: {leaf.name}, Kingdom: {kingdom}") #Kingdom: 669    Opisthokonta
            #set colors per kingdom
        if kingdom == 'Ophistokonta':
            leaf.set_style(nst1)
        elif kingdom == 'Viridiplantae':
            leaf.set_style(nst2)
        elif kingdom == 'Amoebozoa':
            leaf.set_style(nst3)
        elif kingdom == 'Sar':
            leaf.set_style(nst4)
        elif kingdom == 'Discoba':
            leaf.set_style(nst5)
        elif kingdom == 'Rhodophyta':
            leaf.set_style(nst6)
        elif kingdom == 'Cryptophyceae':
            leaf.set_style(nst7)
        elif kingdom == 'Apusozoa':
            leaf.set_style(nst8)
        else:
            print(f"Warning: No matching proteome found for leaf '{leaf.name}'") 
                

    # --- Create the LegendFace Object ---
    legend_face = LegendFace(max_height=200, max_width=200)
    # Add a title
    title = TextFace("Kingdoms:", fsize=12, bold=True)
    # Add title to the first column, first row (0, 0)
    legend_face.add_face(title, column=0, row=0)

    # Iterate over the data to add colored squares and labels
    for i, (color, label) in enumerate(LEGEND_DATA):
        # Add a colored square (size=10, no label text)
        square = SquareFace(size=10, color=color, label="")
        # Place in the first column (col=0) of the legend grid, starting from row 1
        legend_face.add_face(square, column=0, row=i + 1)
        
    # Add the text label
    label_face = TextFace(f" {label}")
    # Place in the second column (col=1) of the legend grid
    legend_face.add_face(label_face, column=1, row=i + 1)
    # Configure the TreeStyle
    ts = TreeStyle()
    ts.mode = "r"          # Draw tree in circular mode otherwise "p"
    ts.rotation=210
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.scale = None #big circles with dotted lines
    ts.optimal_scale_level = "full" #mid or full
    ts.force_topology = True #not bad
    ts.root_opening_factor = 0.3
    #ts.layout_fn = labels
    # Assign the LegendFace object to the TreeStyle's legend attribute
    ts.legend = legend_face
    # Set the position of the legend (1 = 'top-right', 2 = 'float', etc.)
    ts.legend_position = 0

    # Render tree
    main_tree.render(
        main_tree_out,
        tree_style=ts,
        w=1000, #1000 for circular
        units="px"  # Ensure units are pixels for direct control, or w=183, units="mm"
    )
    
    return main_tree, main_tree_out


def load_tree_1000(kingdom, layout): #labels true or false
# File paths
    tree_path = './SHE/tree/' + kingdom  # i.e. ./SHE/tree/Amoebozoa
    kingdom_tree_file = tree_path + '/kingdoms_' + kingdom + '_tree.newick' # i.e. ./SHE/tree/Amoebozoa/kingdoms_Amoebozoa_tree.newick
    # Load the Tree: i.e. tree = Tree(main_tree_file)
    kingdom_tree = Tree(kingdom_tree_file) #format=0
    #main_tree = TreeViz(main_tree_file) #not working

    # define color per kingdom
    nst1 = NodeStyle()
    nst1["bgcolor"] = "#bfdbfe" #Opisthokonta
    nst2 = NodeStyle()
    nst2["bgcolor"] = "#bbf7d0" #Viridiplantae
    nst3 = NodeStyle()
    nst3["bgcolor"] = "#fecdd3" #Amoebozoa
    nst4 = NodeStyle()
    nst4["bgcolor"] = "#f3f4f6" #Sar
    nst5 = NodeStyle()
    nst5["bgcolor"] = "#e9d5ff" #Discoba
    nst6 = NodeStyle()
    nst6["bgcolor"] = "#c3f8f1" #Rhodophyta
    nst7 = NodeStyle()
    nst7["bgcolor"] = "#f6fad5" #Cryptophyceae
    nst8 = NodeStyle()
    nst8["bgcolor"] = "#ddd6fe" #Apusozoa

    #set colors per kingdom
    for leaf in kingdom_tree:
            if kingdom == 'Ophistokonta':
                leaf.set_style(nst1)
            elif kingdom == 'Viridiplantae':
                leaf.set_style(nst2)
            elif kingdom == 'Amoebozoa':
                leaf.set_style(nst3)
            elif kingdom == 'Sar':
                leaf.set_style(nst4)
            elif kingdom == 'Discoba':
                leaf.set_style(nst5)
            elif kingdom == 'Rhodophyta':
                leaf.set_style(nst6)
            elif kingdom == 'Cryptophyceae':
                leaf.set_style(nst7)
            elif kingdom == 'Apusozoa':
                leaf.set_style(nst8)
    else:
        print(f"Warning: No matching proteome found for leaf '{leaf.name}'") 

    # Configure the TreeStyle
    if layout == 'Circular':
        kingdom_tree_out= tree_path + '/kingdoms_' + kingdom + "_Circular_tree1000_labels.svg" #png or pdf or svg
        ts = TreeStyle()
        ts.mode = "c"          # Draw tree in circular mode otherwise "p"
        ts.rotation=210
        ts.show_leaf_name = False #true or false
        ts.show_branch_length = False
        ts.scale = None #big circles with dotted lines
        ts.optimal_scale_level = "full" #mid or full
        ts.force_topology = True #not bad
        ts.root_opening_factor = 100 #new
        ts.layout_fn = labels
        # Render tree
        kingdom_tree.render(
        kingdom_tree_out,
        tree_style=ts,
        w=1000,
        units="px"  # Ensure units are pixels for direct control, or w=183, units="mm"
        )
    elif layout == 'Linear':
        kingdom_tree_out= tree_path + '/kingdoms_' + kingdom + "_Linear_tree1000.svg" #png or pdf or svg
        ts = TreeStyle()
        ts.mode = "r"          # Draw tree in circular mode otherwise "c,r"
        ts.show_leaf_name = False #true or false
        ts.show_branch_length = False
        ts.force_topology = True #not bad
        ts.root_opening_factor = 1 #new
        #ts.layout_fn = labels
        # Render tree
        kingdom_tree.render(
        kingdom_tree_out,
        tree_style=ts,
        w=500,
        units="px"  # Ensure units are pixels for direct control, or w=183, units="mm"
        )
    
    return kingdom_tree, kingdom_tree_out



def load_evotree(kingdom): 
    # File paths
    evo_rate = github_path + 'SHE_app/data/core/merged_evo_rate.csv'
    evorate_df = pd.read_csv(evo_rate)
    evorate_df.columns = evorate_df.columns.str.strip().str.replace(' ', '_')
    tree_path = github_path + 'SHE_app/tree/'  # i.e. ./SHE/tree/Amoebozoa
    tree_file = tree_path + '/kingdoms_' + kingdom + '_tree.newick' # i.e. ./SHE/tree/Amoebozoa/kingdoms_Amoebozoa_tree.newick
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
    tree_out_c= tree_path + '/kingdoms_' + kingdom + "_Circular_tree_evorate.svg" #png or pdf or svg
    ts_c = TreeStyle()
    ts_c.mode = "c"          # Draw tree in circular mode otherwise "p"
    ts_c.rotation=210
    ts_c.show_leaf_name = False #true or false
    ts_c.show_branch_length = False
    ts_c.scale = None #big circles with dotted lines
    ts_c.optimal_scale_level = "full" #mid or full
    ts_c.force_topology = True #not bad
    ts_c.root_opening_factor = 100 #new
    #ts_c.layout_fn = labels
    # Render tree
    tree.render(
    tree_out_c,
    tree_style=ts_c,
    w=1000,
    units="px" 
    )

    #linear layout
    tree_out_r= tree_path + '/kingdoms_' + kingdom + "_Linear_tree_evorate.svg" #png or pdf or svg
    ts_r = TreeStyle()
    ts_r.mode = "r"   
    ts_r.show_leaf_name = True #true or false
    ts_r.show_branch_length = False
    ts_r.force_topology = True #not bad
    ts_r.root_opening_factor = 1 #new
    ts_r.layout_fn = labels
    # Render tree
    tree.render(
    tree_out_r,
    tree_style=ts_r,
    w=500,
    units="px" 
    )
    
    return tree_out_c, tree_out_r



#### intro -----------------------------------------------------------------------------------------------------------------------------------
#st.title('SHE: Strucutral History of Eukarya')

#st.subheader('Strucutral Phylogenetic Tree')
#main_tree, main_tree_out = load_main_tree()

# Choose Bootstrap
#st.sidebar.subheader('Bootstrap')
#bootstrap = st.sidebar.slider(
    #'Number of bootstraps (default: 1000)',
    # 0, 2000, 1000)
#'Selected Bootstrap: ', bootstrap

#normalization
#if st.sidebar.checkbox('Normalize branch lengths'):
    #pass # add normalization function
    #data_load_state = st.text('Retrieving Tree...')
    #'Starting a long computation...'


# Add sidebar ---------------------------------------------------------------------------------------------------------------------------

df = pd.DataFrame({
    'kingdom': ["Opisthokonta", "Amoebozoa", "Sar", "Viridiplantae", "Discoba", "Rhodophyta"],
    'number_species': [10, 20, 30, 40, 50, 70] # change to actual number of species
    })


for kingdom in df['kingdom']:
    #create trees
    print("Processing", kingdom)
    #layout = 'Circular' #or Circular
    load_evotree(kingdom)
    #kingdom_tree, kingdom_tree_out = load_tree_1000(kingdom, layout)

#load_main_tree_color()










