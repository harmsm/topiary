import pytest
import topiary
import ete4 as ete
import os

from topiary.draw.species import species_tree

def test_species_tree():
    # Create a simple species tree
    tree = "((A:1.0,B:1.0):1.0,C:1.0);"
    T = ete.Tree(tree)
    for n in T.leaves():
        n.add_prop("species", f"Species_{n.name}")
    
    # Test species_tree call
    ret = species_tree(T, return_canvas=True)
    assert ret is not None
    # We can check that the name_dict was constructed correctly by inspecting PrettyTree inner state 
    # if it was exposed, but here we'll just check it runs and returns canvas.

    # Test with output_file
    output_file = "test_species_tree.pdf"
    ret = species_tree(T, output_file=output_file)
    assert os.path.exists(output_file)
    os.remove(output_file)
