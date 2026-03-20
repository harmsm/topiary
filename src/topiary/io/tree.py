"""
Load a tree into an ete4 tree data structure.
"""

from topiary._private.check import check_bool

import ete4 as ete
from ete4 import Tree
import dendropy as dp

import glob
import os
import re

def read_tree(tree,fmt=None):
    """
    Load a tree into an ete4 tree data structure.

    Parameters
    ----------
    tree : ete4.Tree or dendropy.Tree or str
        some sort of tree. can be an ete4.Tree (returns self), a dendropy Tree
        (converts to newick and drops root), a newick file or a newick string.
    fmt : int or None
        format for reading tree from newick. 0-9 or 100. (See Notes for what
        these mean). If fmt is None, try to parse without a format descriptor,
        then these formats in numerical order.

    Returns
    -------
    tree : ete4.Tree
        an ete4 tree object.

    Notes
    -----
    `fmt` number is read directly by ete4. See their documentation for how these
    are read (https://etetoolkit.github.io/ete/tutorial/tutorial_trees.html).
    As of ETE4.4.0, these numbers mean:

    + 0: flexible with support values
    + 1: flexible with internal node names
    + 2: all branches + leaf names + internal supports
    + 3: all branches + all names
    + 4: leaf branches + leaf names
    + 5: internal and leaf branches + leaf names
    + 6: internal branches + leaf names
    + 7: leaf branches + all names
    + 8: all names
    + 9: leaf names
    + 100: topology only

    """

    # Already an ete4 tree.
    if isinstance(tree,ete.Tree):
        return tree

    # Convert dendropy tree into newick (drop root)
    if issubclass(type(tree),dp.Tree):
        tree = tree.as_string(schema="newick",suppress_rooting=True)

    # If we get here, we need to convert. If fmt is not specified, try to parse
    # without a format string.
    if fmt is None:

        try:
            t = Tree(tree)
        except Exception:

            # Try all possible formats now, in succession
            w = "\n\nCould not parse tree without format string. Going to try different\n"
            w += "formats. Please check output carefully.\n\n"
            print(w)

            formats = list(range(10))
            formats.append(100)

            t = None
            for f in formats:
                try:
                    t = Tree(tree,parser=f)
                    w = f"\n\nSuccessfully parsed tree with format style {f}.\n"
                    w += "Please see ete4 documentation for details.\n\n"
                    print(w)
                    break

                except Exception:
                    continue

            if t is None:
                err = "\n\nCould not parse tree!\n\n"
                raise ValueError(err)

    else:
        # Try a conversion with the specified format
        t = Tree(tree,parser=fmt)

    return t



def _ete4_node_dict(T,rooted=False):
    """
    Create dictionary keying ete4 tree nodes to tuple of descendants.
    """
    all_leaves = set(list(T.leaf_names()))
    ref_leaf = min(all_leaves)

    node_dict = {}
    for node in T.traverse():
        d1 = set(list(node.leaf_names()))
        if d1 == all_leaves:
            split = ("ROOT",)
        else:
            if rooted:
                split = tuple(sorted(list(d1)))
            else:
                if ref_leaf in d1:
                    split = tuple(sorted(list(all_leaves - d1)))
                else:
                    split = tuple(sorted(list(d1)))

        if split not in node_dict:
            node_dict[split] = node
        else:
            # If we have multiple nodes with the same split (e.g. a unary 
            # root edge), pick the one that bifurcates (has more children). 
            if len(node.children) > len(node_dict[split].children):
                node_dict[split] = node

    # Convert back to list of one node to maintain compatibility with 
    # _map_tree_to_tree loop, but ensure it's unique per split. 
    for split in node_dict:
        node_dict[split] = [node_dict[split]]

    return node_dict

def _toytree_node_dict(T,rooted=False):
    """
    Create dictionary keying toytree tree nodes to tuple of descendants.
    """
    all_leaves = set(T.get_tip_labels())
    ref_leaf = min(all_leaves)

    node_dict = {}
    for node in T.get_nodes():
        d1 = set(node.get_leaf_names())
        if d1 == all_leaves:
            split = ("ROOT",)
        else:
            if rooted:
                split = tuple(sorted(list(d1)))
            else:
                if ref_leaf in d1:
                    split = tuple(sorted(list(all_leaves - d1)))
                else:
                    split = tuple(sorted(list(d1)))

        if split not in node_dict:
            node_dict[split] = node
        else:
            # If we have multiple nodes with the same split, pick node
            # with more children. 
            if len(node.children) > len(node_dict[split].children):
                node_dict[split] = node

    # Convert back to list of one node to maintain compatibility with 
    # _map_tree_to_tree loop, but ensure it's unique per split. 
    for split in node_dict:
        node_dict[split] = [node_dict[split]]

    return node_dict


def _map_tree_to_tree(T1,T2,rooted=False):
    """
    Map nodes between one tree and another based on shared descendants.

    Parameters
    ----------
    T1 : ete4.Tree or toytree.tree
        one tree to compare
    T2 : ete4.Tree or toytree.tree
        second tree to compare

    Returns
    -------
    shared_nodes : list
        list of tuples containing shared nodes between the two trees
    T1_nodes : list
        list of nodes from T1 that are not in T2
    T2_nodes : list
        list of nodes from T2 that are not in T1
    """

    # Construct dictionary keying node to tuple of descendants for T1
    if isinstance(T1,ete.Tree):
        T1_node_dict = _ete4_node_dict(T1,rooted=rooted)
    else:
        T1_node_dict = _toytree_node_dict(T1,rooted=rooted)

    # Construct dictionary keying node to tuple of descendants for T2
    if isinstance(T2,ete.Tree):
        T2_node_dict = _ete4_node_dict(T2,rooted=rooted)
    else:
        T2_node_dict = _toytree_node_dict(T2,rooted=rooted)

    T1_nodes = []
    shared_nodes = []
    T2_keys = list(T2_node_dict.keys())
    for key in T1_node_dict:
        if key in T2_node_dict:
            for n1 in T1_node_dict[key]:
                for n2 in T2_node_dict[key]:
                    shared_nodes.append((n1,n2))
            
            try:
                T2_keys.remove(key)
            except ValueError:
                pass
        else:
            for n1 in T1_node_dict[key]:
                T1_nodes.append(n1)

    # If key is left in T2_keys, the node is only in T2
    T2_nodes = []
    for key in T2_keys:
        for n2 in T2_node_dict[key]:
            T2_nodes.append(n2)

    return shared_nodes, T1_nodes, T2_nodes

def _get_trees_from_directory(directory=None,
                               prefix=None,
                               T_clean=None,
                               T_support=None,
                               T_anc_label=None,
                               T_anc_pp=None,
                               T_event=None):
    """
    Load trees from directory and validate consistency.
    """

    # Load trees from the directory
    if directory is not None:

        tree_files = glob.glob(os.path.join(directory,"*.newick"))
        to_path = dict([(os.path.split(t)[-1],t) for t in tree_files])

        if prefix is None:

            num_rec = len([t for t in tree_files if os.path.basename(t).startswith("reconciled")])
            if num_rec > 0:
                prefix = "reconciled"
            else:
                prefix = "gene"

        if T_clean is None:
            try:
                T_clean = Tree(to_path[f"{prefix}-tree.newick"],parser=0)
            except KeyError:
                pass

        if T_support is None:
            try:
                T_support = Tree(to_path[f"{prefix}-tree_supports.newick"],parser=0)
            except KeyError:
                pass

        if T_event is None:
            try:
                T_event = Tree(to_path[f"{prefix}-tree_events.newick"],parser=1)
            except KeyError:
                pass

        if T_anc_label is None:
            try:
                T_anc_label = Tree(to_path[f"{prefix}-tree_anc-label.newick"],parser=1)
            except KeyError:
                pass

        if T_anc_pp is None:
            try:
                T_anc_pp = Tree(to_path[f"{prefix}-tree_anc-pp.newick"],parser=0)
            except KeyError:
                pass

    # This is the order of priority for getting branch lengths from the tree.
    # The output tree will be copied from the first non-None tree in this
    # list.
    T_list = [T_event,T_anc_pp,T_anc_label,T_support,T_clean]

    # Make sure trees were actually passed in. If none were, return None.
    T_list = [T for T in T_list if T is not None]
    if len(T_list) == 0:
        return None, prefix, {}

    # Make sure all trees have the same descendants
    ref_leaves = set(list(T_list[0].leaf_names()))
    for T in T_list[1:]:
        test_leaves = set(list(T.leaf_names()))
        if len(set(ref_leaves) - set(test_leaves)) != 0:
            err = "All trees must have the same leaves.\n"
            raise ValueError(err)
    
    # Bundle trees for return
    trees = {"T_clean":T_clean,
             "T_support":T_support,
             "T_anc_label":T_anc_label,
             "T_anc_pp":T_anc_pp,
             "T_event":T_event}

    return T_list, prefix, trees

def _get_root_node(T, root_on):
    """
    Get the node in T that matches the root_on node from another tree.
    """
    # This matches based on identical descendant sets
    node_dict = _ete4_node_dict(T, rooted=True)
    
    # Try to find a node that has one of the two root splits as its
    # descendants.
    for s in root_on:
        if s is not None and s in node_dict:
            return node_dict[s][0]
            
    return None

def _synchronize_tree_rooting(T_list, prefix):
    """
    Synchronize all trees in T_list to the same root, stashing root props.
    """

    # Identify the tree to use for rooting information. 
    root_template = T_list[0].copy()
    
    # Force the template to be binary rooted if it's not already.
    # We do this by unrooting and re-rooting to midpoint if multifurcating at 
    # the root.
    if len(root_template.children) != 2:
        root_template.unroot()
        try:
            root_template.root.del_prop("support")
            root_template.root.del_prop("dist")
        except (AttributeError, KeyError):
            pass
        root_template.set_midpoint_outgroup()

    # Identify the two splits at the root.
    root_on = []
    for c in root_template.children:
        leaves = list(c.leaf_names())
        leaves.sort()
        root_on.append(tuple(leaves))

    # Synchronize ALL trees in the list to this binary rooting.
    for T in T_list:
        
        # Stash current root properties before we unroot/re-root
        root_props = T.root.props.copy()

        # Expressly stash name and support. These are not always in props. 
        root_name = T.root.name
        root_support = T.root.support
        
        # Strip root node properties as they cause ETE4 assertion 
        # errors during rooting and can lead to duplication if they 
        # end up on a child node after re-rooting. 
        for p in list(T.root.props.keys()):
            try:
                T.root.del_prop(p)
            except (AttributeError, KeyError):
                pass
        
        # Reset name and support to avoid them being moved to children. 
        # Note: ETE4's unroot() asserts that the root does NOT have 
        # branch properties like support or dist.
        T.root.name = ""
        try:
            T.root.del_prop("support")
        except (AttributeError, KeyError):
            pass
        try:
            T.root.del_prop("dist")
        except (AttributeError, KeyError):
            pass
        
        # If the tree is bifurcating at the root, ensure the children have 
        # identical support values. ETE4.unroot() will raise an AssertionError 
        # if the root node has two children with different support values. 
        if len(T.children) == 2:
            s0 = getattr(T.children[0],"support",0.0)
            s1 = getattr(T.children[1],"support",0.0)
            if s0 != s1:
                T.children[0].support = s0
                T.children[1].support = s0

        T.unroot()
        new_outgroup = _get_root_node(T, root_on)
        if new_outgroup is not None:
            try:
                T.set_outgroup(new_outgroup)
            except Exception:
                T.set_midpoint_outgroup()
        else:
            T.set_midpoint_outgroup()
        
        # Ensure binary root
        if len(T.children) != 2:
            T.set_midpoint_outgroup()

        # Re-apply stashed properties to the new root
        T.root.props.update(root_props)
        if root_name:
            T.root.name = root_name
        if root_support is not None:
            T.root.support = root_support

    return T_list, root_on

def _prepare_output_tree(T_list, root_on, prefix):
    """
    Create the reference tree and initialize features.
    """
    # Use first tree in list (highest priority for branch lengths) as the 
    # reference tree to return. 
    out_tree = T_list[0].copy()

    # Clean out_tree internal nodes of all properties that might interfere with 
    # rooting/drawing. 
    for n in out_tree.traverse():
        if not n.is_leaf:
            try:
                n.del_prop("support")
                n.del_prop("dist")
            except (AttributeError, KeyError):
                pass
            n.name = ""

    # Root out_tree precisely to the same split as others if it appears unrooted.
    # An unrooted tree in ETE4 typically has 3+ children at the root node. 
    # A rooted tree is unary (1 child) or binary (2 children).
    if len(out_tree.children) > 2:
        out_tree.unroot()
        new_outgroup = _get_root_node(out_tree, root_on)
        if new_outgroup is not None:
            try:
                out_tree.set_outgroup(new_outgroup)
            except Exception:
                out_tree.set_midpoint_outgroup()
        else:
            out_tree.set_midpoint_outgroup()

    # Ensure binary root
    if len(out_tree.children) != 2:
        out_tree.set_midpoint_outgroup()

    # Initialize empty features
    for n in out_tree.traverse():
        if not n.is_leaf:
            n.add_prop("event",None)
            n.add_prop("anc_pp",None)
            n.add_prop("anc_label",None)
            n.add_prop("bs_support",None)

    return out_tree

def _merge_tree_features(out_tree, T_list, root_on, prefix, 
                        T_clean, T_support, T_anc_label, T_anc_pp, T_event):
    """
    Map features from individual trees onto the output tree.
    """

    # features_to_load maps trees with information to copy (keys) to what
    # feature we should extract from that tree. Values are (feature_in_tree,
    # name_of_feature_in_out_tree, whether_to_allow_root_value).
    features_to_load = {}
    if T_clean is not None:
        features_to_load[T_clean] = [("dist","dist",True),("support","support",False)]

    # Map other features. rooted_ok_val only applied if this is a reconciled
    # tree (where root nodes often have events like Speciation). For ancestors,
    # we generally don't want the root node to have a label. 
    root_ok_val = False
    for T, out_f, in_f, root_ok in [(T_event,"event","name",True),
                                    (T_anc_pp,"anc_pp","support",root_ok_val),
                                    (T_anc_label,"anc_label","name",root_ok_val),
                                    (T_support,"bs_support","support",False)]:
        if T is not None:
            if T not in features_to_load:
                features_to_load[T] = []
            features_to_load[T].append((in_f,out_f,root_ok))

    # Ensure we have dist and support (required for ETE4 write)
    has_dist = False
    has_support = False
    for T in features_to_load:
        for f in features_to_load[T]:
            if f[1] == "dist": has_dist = True
            if f[1] == "support": has_support = True
    
    if not has_dist or not has_support:
        for T in T_list:
            if T is not None:
                if T not in features_to_load:
                    features_to_load[T] = []
                if not has_dist:
                    features_to_load[T].append(("dist","dist",True))
                if not has_support:
                    features_to_load[T].append(("support","support",False))
                break

    stash_values = {}
    for T in features_to_load:

        # Since all have the same root and descendants, tree nodes should be
        # identical between trees and uniquely identified by their descendants
        shared, T_alone, out_alone = _map_tree_to_tree(T,out_tree,rooted=True)
        
        # Pull out_feature for error message if topology mismatch.
        out_feature_for_err = features_to_load[T][0][1]
        
        if len(T_alone) > 0 or len(out_alone) > 0:
            print(f"Topology mismatch during {out_feature_for_err}!")
            print(f"T_alone: {[list(n.leaf_names()) for n in T_alone]}")
            print(f"out_alone: {[list(n.leaf_names()) for n in out_alone]}")
            err = "Cannot merge trees with different topologies.\n"
            raise ValueError(err)

        # Map data from features on input tree to the features on the output
        # tree
        for s in shared:
            in_node = s[0]
            out_node = s[1]

            for in_feature, out_feature, root_allowed in features_to_load[T]:

                value = in_node.get_prop(in_feature)
                if value is None:
                    value = in_node.get_prop(f"_{in_feature}")
                if value is None:
                    value = getattr(in_node,in_feature,None)

                # small hack --> anc to a
                if out_feature == "anc_label" and value is not None:
                    value = re.sub("anc","a",str(value))

                if value is not None:
                    # If root node, pull out value if not allowed on root
                    if out_node.is_root:
                        if not root_allowed:
                            if out_feature == "anc_label" or out_feature == "anc_pp":
                                if value is not None and value != "":
                                    stash_values[out_feature] = value

                            if out_feature in ["dist","support","name"]:
                                try:
                                    out_node.del_prop(out_feature)
                                except (KeyError,AttributeError):
                                    pass
                            else:
                                out_node.add_prop(out_feature,None)
                        else:
                            out_node.add_prop(out_feature,value)
                    else:
                        out_node.add_prop(out_feature,value)

    # If we have stashed values, we need to restore them. They were at the 
    # root of the tree they came from. 
    if stash_values.get("anc_label") or stash_values.get("anc_pp"):
        
        # Heuristic: try to find an internal node that lacks an ancestral label.
        # This is likely the node that corresponds to the original root split.
        # (This keeps the labels on internal nodes rather than the root node 
        # itself, which matches Topiary's historical behavior and ensures 
        # they appear in Newick format 1). 
        target_node = None
        for n in out_tree.traverse():
            if n.is_leaf or n.is_root:
                continue
            
            # Use a very broad check for empty
            val = n.get_prop("anc_label")
            if val is None or str(val).strip() == "" or str(val) == "None":
                target_node = n
                break
        
        # If no internal node found, fallback to root
        if target_node is None:
            target_node = out_tree.root
        
        # Restore stashed values independently
        for out_f in ["anc_label", "anc_pp"]:
            if out_f in stash_values:
                # If anc_label, make sure we do the aX conversion
                val = stash_values[out_f]
                if out_f == "anc_label":
                    val = re.sub("anc","a",str(val))
                
                # Store the value. We allow overwriting None/empty/0.0. 
                current = target_node.get_prop(out_f)
                if current is None or str(current).strip() == "" or current == 0.0 or target_node.is_root:
                    target_node.add_prop(out_f, val)

    return out_tree

def load_trees(directory=None,
               prefix=None,
               T_clean=None,
               T_support=None,
               T_anc_label=None,
               T_anc_pp=None,
               T_event=None):
    """
    Generate an ete4 tree with features 'event', 'anc_pp', 'anc_label',
    and 'bs_support' on internal nodes. This information is read from the input
    ete4 trees or the specified topiary output directory. The tree is rooted
    using T_event. If this tree is not specified, the midpoint root is used.
    Trees are read from the directory first, followed by any ete4 trees
    specified as arguments. (This allows the user to override trees from the
    directory if desired). If no trees are passed in, returns None.

    Warning: this will modify input ete4 trees as it works on the trees
    rather than copies.

    Parameters
    ----------
    directory : str
        output directory from a topiary calculation that has .newick files
        in it. Function will load all trees in that directory.
    prefix : str, optional
        what type of trees to plot from the directory. should be "reconciled"
        or "gene". If None, looks for reconciled trees. If it finds any, these
        prefix = "reconciled"
    T_clean : ete4.Tree, optional
        clean tree (leaf labels and branch lengths, nothing else). Stored as
        {}-tree.newick in output directories.
    T_support : ete4.Tree, optional
        support tree (leaf labels, branch lengths, supports). Stored as
        {}-tree_supports.newick in output directories.
    T_anc_label : ete4.Tree, optional
        ancestor label tree (leaf labels, branch lengths, internal names)
        Stored as {}-tree_anc-label.newick.
    T_anc_pp : ete4.Tree, optional
        ancestor posterior probability tree (leaf labels, branch lengths,
        posterior probabilities as supports) Stored as {}-tree_anc-pp.newick.
    T_event : ete4.Tree, optional
        tree with reconciliation events as internal labels (leaf labels,
        branch lengths, event labels). Stored as reconciled-tree_events.newick

    Returns
    -------
    merged_tree : ete4.Tree or None
        rooted tree with features on internal nodes. Return None if no trees
        are passed in.
    """

    T_list, prefix, trees = _get_trees_from_directory(directory,
                                                      prefix,
                                                      T_clean,
                                                      T_support,
                                                      T_anc_label,
                                                      T_anc_pp,
                                                      T_event)
    
    if T_list is None:
        return None

    # Synchronize rooting across all trees
    T_list, root_on = _synchronize_tree_rooting(T_list,prefix)

    # Prepare output tree
    out_tree = _prepare_output_tree(T_list,root_on,prefix)

    # Merge features from all trees onto out_tree
    out_tree = _merge_tree_features(out_tree,T_list,root_on,prefix,**trees)

    return out_tree

def write_trees(T,
                name_dict=None,
                out_file=None,
                overwrite=False,
                anc_pp=True,
                anc_label=True,
                bs_support=True,
                event=True):
    """
    Write out an ete.Tree as a newick format. This function looks for features
    set by :code:`load_trees` and then writes an individual tree out with each
    feature. The features are :code:`anc_pp`, :code:`anc_label`, :code:`bs_support`,
    and :code:`event`. This will write out trees for any of these features 
    present; not all features need to be in place for this function to work. 

    Parameters
    ----------
    T : ete.TreeNode
        ete tree with information loaded into appropriate features. This is the
        tree returned by :code:`load_trees`. 
    name_dict : dict
        name_dict : dict, optional
        dictionary mapping strings in node.name to more useful names. (Can be
        generated using :code:`topiary.draw.core.create_name_dict`). If not 
        specified, trees are written out with uid as tip names
    out_file : str, optional
        output file. If defined, write the newick string the file.
    overwrite : bool, default=False
        whether or not to overwrite an existing file
    anc_pp : bool, default=True
        whether or not to write a tree with anc_pp as support values
    anc_label : bool, default=True
        whether or not to write a tree with anc_label as internal node names
    bs_support : bool, default=True
        whether or not to write a tree with bs_support as support values
    event : bool, default=True
        whether or not to write a tree with events as internal node names

    Returns
    -------
    tree : str
        Newick string representation of the output tree(s)
    """
    
    # --------------------------------------------------------------------------
    # Parameter sanity checking
    if not isinstance(T,ete.Tree):
        err = "\nT must be an ete4.Tree instance\n\n"
        raise ValueError(err)

    if name_dict is not None:
        if not issubclass(type(name_dict),dict):
            err = "\nname_dict must be a dictionary\n\n"
            raise ValueError(err)
    
    if out_file is not None:
        if not issubclass(type(out_file),str):
            err = "\nout_file must be a string pointing to a file to write out\n\n"
            raise ValueError(err)
        
        overwrite = check_bool(overwrite,"overwrite")

        if os.path.exists(out_file):
            if os.path.isfile(out_file):
                if overwrite:
                    os.remove(out_file)
                else:
                    err = f"\nout_file '{out_file}' exists. Either delete or set overwrite to True\n\n"
                    raise FileExistsError(err)
            else:
                err = f"\nout_file '{out_file}' exists but is a directory. Cannot write output.\n\n"
                raise FileExistsError(err)

    anc_pp = check_bool(anc_pp,"anc_pp")
    anc_label = check_bool(anc_label,"anc_label")
    bs_support = check_bool(bs_support,"bs_support")
    event = check_bool(event,"event")

    # --------------------------------------------------------------------------
    # Set up output

    out_trees = []

    # Work on a copy
    T = T.copy()

    # If name dict is specified
    if name_dict is not None:
        for n in T.traverse():
            if n.is_leaf:
                n.name = name_dict[n.name]

    # --------------------------------------------------------------------------
    # Create bs_supports, events, anc_label, and anc_pp

    # bs_supports
    if bs_support:
        write_bs_supports = False
        for n in T.traverse():
            if not n.is_leaf:
                if n.get_prop("bs_support") is not None:
                    n.support = n.get_prop("bs_support")
                    write_bs_supports = True
                else:
                    n.support = 0.0

        # parser=0 --> all branches + leaf names + supports as internal names
        if write_bs_supports:
            out_trees.append(T.write(parser=0))

    # event
    if event:
        write_events = False
        for n in T.traverse():
            if not n.is_leaf:
                if n.get_prop("event") is not None:
                    n.name = n.get_prop("event")
                    write_events = True
                else:
                    n.name = ""

        # parser=1 --> all branches + leaf names + internal names
        if write_events:
            out_trees.append(T.write(parser=1))

    # anc_label
    if anc_label:
        write_anc_label = False
        for n in T.traverse():
            if not n.is_leaf:
                if n.get_prop("anc_label") is not None:
                    n.name = n.get_prop("anc_label")
                    write_anc_label = True
                else:
                    n.name = ""

        # parser=1 --> all branches + leaf names + internal names
        if write_anc_label:
            out_trees.append(T.write(parser=1))

    # anc_pp
    if anc_pp:
        write_anc_pp = False
        for n in T.traverse():
            if not n.is_leaf:
                if n.get_prop("anc_pp") is not None:
                    n.support = n.get_prop("anc_pp")
                    write_anc_pp = True
                else:
                    n.support = 0.0

        # parser=0 --> all branches + leaf names + supports as internal names
        if write_anc_pp:
            out_trees.append(T.write(parser=0))

    # --------------------------------------------------------------------------
    # Finalize output

    out_trees = "\n".join(out_trees)

    if out_file is not None:
        f = open(out_file,'w')
        f.write(out_trees)
        f.close()

    return out_trees 
    