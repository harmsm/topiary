import pytest
import topiary
import ete4 as ete
import pandas as pd
import numpy as np
import os
from unittest.mock import MagicMock, patch

from topiary.draw.tree import tree

def test_tree(tmpdir):
    
    T = ete.Tree("((A:1.0,B:1.0):1.0,C:1.0);")
    for n in T.traverse():
        if not n.is_leaf:
            n.add_prop("event","D")
            n.add_prop("anc_pp",0.9)
            n.add_prop("bs_support",100.0)
            n.add_prop("anc_label","a1")
    
    df = pd.DataFrame({
        "species":["Species A", "Species B", "Species C"],
        "name":["A", "B", "C"],
        "ott":["ottA", "ottB", "ottC"],
        "uid":["uidA", "uidB", "uidC"],
        "keep":[True, True, True]
    })

    # Mock Supervisor
    mock_supervisor = MagicMock()
    mock_supervisor.output_dir = "mock_output"
    mock_supervisor.tree_prefix = "mock_prefix"
    mock_supervisor.df = df

    # Mock load_trees and PrettyTree
    with patch("topiary.draw.tree.load_trees", return_value=T):
        with patch("topiary.draw.tree.PrettyTree") as MockPrettyTree:
            instance = MockPrettyTree.return_value
            instance.default_size = 10
            instance.plotted_properties = ["event", "bs_support", "anc_pp"]
            
            # Call tree with supervisor
            ret = tree(mock_supervisor, return_canvas=True)
            
            # Verify PrettyTree was initialized
            MockPrettyTree.assert_called()
            
            # Verify drawing methods were called. 
            # It should be called for event, bs_support, and anc_pp
            assert instance.draw_nodes.call_count >= 3
            
            # Test call with string path
            with patch("topiary.draw.tree.Supervisor", return_value=mock_supervisor):
                ret = tree("some_path", return_canvas=True)
                
            # Test with output_file
            output_file = os.path.join(tmpdir, "test_tree.pdf")
            ret = tree(mock_supervisor, output_file=output_file)
            instance.render.assert_called_with(str(output_file))

    # Test coverage of logic for node_color and node_size
    with patch("topiary.draw.tree.load_trees", return_value=T):
        with patch("topiary.draw.tree.PrettyTree") as MockPrettyTree:
            with patch("topiary.draw.tree.construct_sizemap", return_value=(lambda x: 10, None)):
                instance = MockPrettyTree.return_value
                instance.default_size = 10
                instance.plotted_properties = ["event", "bs_support", "anc_pp"]

                # node_color is not None (should bypass specific colors)
                tree(mock_supervisor, node_color="red", node_size=20)
                instance.draw_nodes.assert_any_call(color="red", size=20)

                # node_size as list/array (forcing ValueError/TypeError logic in tree.py)
                node_size_list = [5, 10]
                tree(mock_supervisor, node_size=node_size_list)

                # node_size as dict
                node_size_dict = {"D": 20, "S": 10}
                tree(mock_supervisor, node_size=node_size_dict)

                # Test anc_link_path
                tree(mock_supervisor, anc_link_path="http://example.com")
                instance.draw_node_labels.assert_called()

                # Test event_color=None, bs_color=None, pp_color=None
                instance.draw_nodes.reset_mock()
                tree(mock_supervisor, event_color=None, bs_color=None, pp_color=None,
                     anc_label=True, event_label=True, bs_label=True, pp_label=True)
                instance.draw_nodes.assert_not_called()
                instance.draw_node_labels.assert_called()
