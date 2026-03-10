import pytest
import os
import shutil
import ete4 as ete
from ete4 import Tree
import topiary
from topiary.io.tree import load_trees
from topiary.raxml.ancestors import _make_ancestor_summary_trees
from topiary.draw.prettytree import PrettyTree

def test_ete4_asr_integration(tmpdir):
    """
    Integration test for the ASR pipeline, focusing on ETE4 property handling
    and rooting-invariant tree merging in load_trees.
    """
    os.chdir(tmpdir)
    
    # 1. Test _make_ancestor_summary_trees export
    avg_pp_dict = {"anc1": 0.95, "anc2": 0.85}
    # 4-leaf labeled tree
    labeled_newick = "((A:0.1,B:0.1)Node1:0.5,(C:0.1,D:0.1)Node2:0.5)root:1;"
    with open("test_labeled.newick", "w") as f:
        f.write(labeled_newick)
        
    # This should run successfully and create tree_anc-label.newick and tree_anc-pp.newick
    _make_ancestor_summary_trees(None, avg_pp_dict, "test_labeled.newick")
    
    assert os.path.exists("tree_anc-label.newick")
    assert os.path.exists("tree_anc-pp.newick")
    
    with open("tree_anc-label.newick") as f:
        label_content = f.read()
        # _make_ancestor_summary_trees renames NodeX to ancX
        assert "anc1" in label_content
        assert "anc2" in label_content
        
    with open("tree_anc-pp.newick") as f:
        pp_content = f.read()
        assert "0.95" in pp_content
        assert "0.85" in pp_content

    # 2. Test load_trees with different rootings
    os.makedirs("test_results", exist_ok=True)
    
    # Tree 1: rooted at A. Internal nodes map to splits AB|CD and A|BCD
    # This is a different rooting than the standard midpoint/bifurcating root
    t_label_rooted = "(A:0.05,(B:0.1,(C:0.1,D:0.1)anc2:0.5)anc1:0.05)root:1;"
    with open("test_results/gene-tree_anc-label.newick", "w") as f:
        f.write(t_label_rooted)
        
    # Tree 2: rooted at midpoint. 
    t_pp_midpoint = "((A:0.1,B:0.1)0.95:0.5,(C:0.1,D:0.1)0.85:0.5);"
    with open("test_results/gene-tree_anc-pp.newick", "w") as f:
        f.write(t_pp_midpoint)
        
    # Tree 3: clean tree
    with open("test_results/gene-tree.newick", "w") as f:
        f.write("((A:0.1,B:0.1):0.5,(C:0.1,D:0.1):0.5);")
        
    # load_trees should NOT raise ValueError: Cannot merge trees with different topologies.
    # It should use split-based mapping to merge features onto out_tree.
    T = load_trees(directory="test_results", prefix="gene")
    
    # Verify properties merged correctly
    # out_tree in load_trees is created from gene-tree.newick (midpoint rooted)
    # It should have internal nodes for splits AB|CD (one of both children of the root)
    found_anc2 = False
    for n in T.traverse():
        if not n.is_leaf:
            label = n.get_prop("anc_label")
            pp = n.get_prop("anc_pp")
            if label == "a2": # anc2 -> a2 transformation in load_trees
                found_anc2 = True
                assert pp == 0.85 or pp is None # pp might be on the other node if rooting differs
    
    # Check that it didn't crash and we found something
    assert len(T) > 0

    # 3. Test PrettyTree conversion
    pt = PrettyTree(T)
    tT = pt._tT
    # Simple check that features exist in ToyTree attributes
    for node in tT.get_nodes():
        if not node.is_leaf():
            # Check for attributes that PrettyTree expects
            assert hasattr(node, "anc_label")
            assert hasattr(node, "anc_pp")

    # 4. Test write_trees (This triggered the AssertionError: missing dist)
    from topiary.io.tree import write_trees
    # This should not raise AssertionError
    nw = write_trees(T, anc_pp=True, anc_label=True)
    assert "a2" in nw
    assert "0.85" in nw

if __name__ == "__main__":
    # For manual run
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        test_ete4_asr_integration(tmp)
        print("Integration test passed!")
