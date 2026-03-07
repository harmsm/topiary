import topiary
import sys
import os

def test_debug_path():
    print(f"DEBUG TOPIARY FILE: {topiary.__file__}")
    print(f"DEBUG SYS PATH: {sys.path}")
    # Read opentree/tree.py and print some lines
    tree_file = os.path.join(os.path.dirname(topiary.__file__),"opentree","tree.py")
    print(f"DEBUG TREE FILE: {tree_file}")
    with open(tree_file,'r') as f:
        lines = f.readlines()
        if len(lines) >= 114:
            print(f"DEBUG LINE 114: {lines[113].strip()}")
