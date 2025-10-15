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
            if kingdom == 'Opisthokonta':
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
        
    # --- 2. Define Legend Data ---
    LEGEND_DATA = [
        ("#bfdbfe", "Opisthokonta"),
        ("#bbf7d0", "Viridiplantae"),
        ("#fecdd3", "Amoebozoa"),
        ("#f3f4f6", "Sar"),
        ("#e9d5ff", "Discoba"),
        ("#c3f8f1", "Rhodophyta"),
        ("#f6fad5", "Cryptophyceae"),
        ("#ddd6fe", "Apusozoa"),
    ]
    # --- 3. Create the LegendFace Object ---
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

"""
for kingdom in df['kingdom']:
    #create trees
    print("Processing", kingdom)
    layout = 'Circular' #or Circular
    kingdom_tree, kingdom_tree_out = load_tree_1000(kingdom, layout)"""

load_main_tree_color()

















