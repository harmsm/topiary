"""
NOTE: This is a simplified, pure-Python implementation of the DOWNPASS algorithm
used by pastml for ancestral gap reconstruction. It was derived from the pastml
codebase, but does not include any of the original pastml code. That library is
available at https://github.com/evolbioinfo/pastml. When citing topiary, please
cite the original pastml paper: https://doi.org/10.1093/molbev/msz131. 
"""

import ete3

import numpy as np
import pandas as pd

def _get_most_common_states(state_arrays):
    """
    Get the most common states from a list of boolean state arrays.
    Each array has shape (num_sites, 2).
    """
    counts = np.sum(state_arrays, axis=0)
    max_counts = counts.max(axis=1, keepdims=True)
    return counts == max_counts

def get_ancestral_gaps(alignment_file,tree_file,prediction_method="DOWNPASS"):
    """
    Get ancestral gaps from an alignment and raxml output tree file. Gaps are
    reconstructed by parsimony using the DOWNPASS algorithm as implemented in
    pastml.

    Parameters
    ----------
    alignment_file : str
        phy file used to generate ancestors in RAxML
    tree_file : str
        output tree file with labeled internal nodes
    prediction_method : str, default="DOWNPASS"
        method to reconstruct gaps.

    Returns
    -------
    gap_anc_dict : dict
        dictionary keying internal node names to lists of True (gap), False (no
        gap), and None (gapping unclear) for each site in that ancestor.
    """

    # Read the alignment file
    counter = 0
    leaf_names = []
    with open(alignment_file) as f:
        for line in f:

            # First line
            if counter == 0:
                col = line.split()
                num_taxa = int(col[0])
                num_sites = int(col[1])
                counter += 1

                char_matrix = np.zeros((num_taxa,num_sites),dtype=np.bool_)
                continue

            # Next, blank line
            if line.strip() == "":
                counter += 1
                continue

            # Alternating leaf id and sequence lines
            if counter % 2 == 0:
                leaf_names.append(line.strip())
                counter += 1
            else:
                index = (counter - 3)//2
                char_matrix[index,:] = np.array([c == "-" for c in line.strip()])
                counter += 1

    if prediction_method != "DOWNPASS":
        raise ValueError(f"Only prediction_method='DOWNPASS' is supported natively. Got {prediction_method}.")

    # Load the tree, keeping the internal node names
    tree = ete3.Tree(tree_file,format=1)

    # 1. Initialize states for each node
    # state array: shape (num_sites, 2), where [:, 0] is False state, [:, 1] is True state
    leaf_name_to_index = {name: i for i, name in enumerate(leaf_names)}
    
    for node in tree.traverse():
        if node.is_leaf():
            idx = leaf_name_to_index[node.name]
            node.ps_orig = np.zeros((num_sites, 2), dtype=np.bool_)
            node.ps_orig[:, 0] = ~char_matrix[idx, :]
            node.ps_orig[:, 1] = char_matrix[idx, :]
        else:
            node.ps_orig = np.ones((num_sites, 2), dtype=np.bool_)
        
        node.ps_bu = np.copy(node.ps_orig)
        node.ps_preset = np.copy(node.ps_orig)

    # 2. UPPASS
    for node in tree.traverse('postorder'):
        if not node.is_leaf():
            children_states = _get_most_common_states([child.ps_bu for child in node.children])
            state_intersection = node.ps_bu & children_states
            has_intersection = state_intersection.any(axis=1, keepdims=True)
            node.ps_bu = np.where(has_intersection, state_intersection, node.ps_bu)

    # 3. DOWNPASS
    for node in tree.traverse('preorder'):
        if node.is_root():
            node.ps_up = np.ones((num_sites, 2), dtype=np.bool_)
        else:
            state_arrays = [node.up.ps_up] + [sibling.ps_bu for sibling in node.up.children if sibling != node]
            node.ps_up = _get_most_common_states(state_arrays)
            
        if not node.is_leaf():
            state_arrays = [node.ps_up] + [child.ps_bu for child in node.children]
            down_up_states = _get_most_common_states(state_arrays)
        else:
            down_up_states = node.ps_up
            
        preset_states = node.ps_preset
        state_intersection = down_up_states & preset_states
        has_intersection = state_intersection.any(axis=1, keepdims=True)
        node.ps_final = np.where(has_intersection, state_intersection, preset_states)

    # 4. Extract results
    gap_anc_dict = {}
    for node in tree.traverse('preorder'):
        if node.name in leaf_names:
            continue

        # Convert state array back to list of True, False, None
        final_state = node.ps_final
        
        is_false = final_state[:, 0] & ~final_state[:, 1]
        is_true = ~final_state[:, 0] & final_state[:, 1]
        
        # Build the list
        state_list = []
        for f, t in zip(is_false, is_true):
            if f:
                state_list.append(False)
            elif t:
                state_list.append(True)
            else:
                state_list.append(None)
                
        gap_anc_dict[node.name] = state_list

    return gap_anc_dict
