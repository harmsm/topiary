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
    all_leaves = set(T.leaf_names())
    ref_leaf = min(all_leaves)

    node_dict = {}
    for node in T.traverse():
        d1 = set(node.leaf_names())
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
            node_dict[split] = []
        node_dict[split].append(node)

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
            node_dict[split] = []
        node_dict[split].append(node)

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

    # Load trees from the directory
    if directory is not None:

        tree_files = glob.glob(os.path.join(directory,"*.newick"))
        to_path = dict([(os.path.split(t)[-1],t) for t in tree_files])

        if prefix is None:

            num_rec = len([t for t in tree_files if t.startswith("reconciled")])
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
        return None

    # Make sure all trees have the same descendants
    ref_leaves = set(list(T_list[0].leaf_names()))
    for T in T_list[1:]:
        test_leaves = set(list(T.leaf_names()))
        if len(set(ref_leaves) - set(test_leaves)) != 0:
            err = "All trees must have the same leaves.\n"
            raise ValueError(err)

    # Strip root node support properties as they cause ETE4 assertion errors during rooting
    for T in T_list:
        try:
            T.root.del_prop("support")
        except (AttributeError, KeyError):
            pass
        try:
            T.root.del_prop("dist")
        except (AttributeError, KeyError):
            pass

    # If this is not a reconciled tree, strip the root node name. ETE4 requires 
    # a clean root to avoid issues during rooting/drawing, and most gene tree
    # formats don't have meaningful data at the root junction. For reconciled
    # trees, we keep it because it stores the root evolutionary event (D, S, etc.)
    if prefix != "reconciled":
        T.root.name = ""

    # If we have an event tree, root all trees on that rooted tree
    if prefix == "reconciled":

        if T_event is not None:
            root_tree = T_event
        else:
            root_tree = T_clean

        # Get left and right descendants of the root node
        root_on = []
        for n in root_tree.root.descendants():
            leaves = n.leaf_names()
            leaves = list(leaves)
            leaves.sort()
            root_on.append(tuple(leaves))
            if len(root_on) == 2:
                break

    # If not a reconciled tree, do midpoint rooting using first tree in list.
    else:
        # ete4 midpoint rooting changes tree in place and doesn't return outgroup easily
        # but we can get it from children if we root it.
        # We must do this on a COPY because set_midpoint_outgroup resolves multifurcations
        # permanently, causing topology mismatches down the line.
        tmp_T = T_list[0].copy()
        for n in tmp_T.traverse():
            if not n.is_leaf:
                try:
                    n.del_prop("support")
                except (AttributeError, KeyError):
                    pass
                n.name = ""
        
        # Root MUST NOT have distance or support for ETE4 rooting
        try:
            tmp_T.root.del_prop("dist")
        except (AttributeError, KeyError):
            pass

        tmp_T.set_midpoint_outgroup()
        
        root_on = []
        for n in tmp_T.children:
            leaves = list(n.leaf_names())
            leaves.sort()
            root_on.append(tuple(leaves))

    # Synchronize rooting for all trees in T_list based on root_on. Only do this 
    # for non-reconciled trees. Reconciled trees are already rooted correctly
    # on their event root. 
    if prefix != "reconciled":
        for T in T_list:
            if T is None:
                continue
            T.unroot()
            try:
                left = T.common_ancestor(root_on[0])
                T.set_outgroup(left)
            except Exception:
                try:
                    right = T.common_ancestor(root_on[1])
                    T.set_outgroup(right)
                except Exception:
                    pass

    # Make new tree from first tree in list. This will be our output tree.
    out_tree = T_list[0].copy()

    # Clean out_tree internal nodes of all properties that might interfere with 
    # rooting. Since we map properties back later, this is safe. 
    for n in out_tree.traverse():
        if not n.is_leaf:
            try:
                n.del_prop("support")
            except (AttributeError, KeyError):
                pass
            try:
                n.del_prop("dist")
            except (AttributeError, KeyError):
                if not n.is_root: # Root might not have dist
                    pass
            n.name = ""

    # Root out_tree!
    out_tree.unroot()
    left = out_tree.common_ancestor(root_on[0])
    right = out_tree.common_ancestor(root_on[1])
    if left is right:
        if len(root_on[0]) == 1:
            left = root_on[0][0]
        elif len(root_on[1]) == 1:
            right = root_on[1][0]
    
    try:
        if not left.is_root:
            out_tree.set_outgroup(left)
        elif not right.is_root:
            out_tree.set_outgroup(right)
    except (ete.TreeError,AssertionError) as e:
        try:
            if not right.is_root:
                out_tree.set_outgroup(right)
            elif not left.is_root:
                out_tree.set_outgroup(left)
        except (ete.TreeError,AssertionError) as e2:
            print(f"Failed to root out_tree: {e2}")
    for n in out_tree.traverse():

        # Create empty features
        if not n.is_leaf:
            n.add_prop("event",None)
            n.add_prop("anc_pp",None)
            n.add_prop("anc_label",None)
            n.add_prop("bs_support",None)


    # features_to_load maps trees with information to copy (keys) to what
    # feature we should extract from that tree. Values are (feature_in_tree,
    # name_of_feature_in_out_tree,whether_to_allow_root_value).
    features_to_load = {}
    if T_clean is not None:
        features_to_load[T_clean] = [("dist","dist",True),("support","support",False)]

    # Map other features
    root_ok_val = (prefix == "reconciled")
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

    # Only copy from trees that are not None.
    trees = [k for k in features_to_load.keys()]

    stash_values = {}
    for T in trees:

        # Since all have the same root and descendants, tree nodes should be
        # identical between trees and uniquely identified by their descendants
        shared, T_alone, out_alone = _map_tree_to_tree(T,out_tree,rooted=False)
        
        # Pull out_feature for error message if topology mismatch. This is just
        # picking the first one in the list for this tree. 
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

            # Node to copy from to
            in_node = s[0]
            out_node = s[1]

            for in_feature, out_feature, root_allowed in features_to_load[T]:

                # Get value from in
                value = in_node.get_prop(in_feature)
                if value is None:
                    value = in_node.get_prop(f"_{in_feature}")

                # If still None, check if this is a core attribute (e.g. name, support, dist)
                if value is None:
                    value = getattr(in_node,in_feature,None)

                # small hack --> anc to a
                if out_feature == "anc_label" and value is not None:
                    value = re.sub("anc","a",str(value))

                # Add value to out. We only add if not None because ETE4
                # properties like dist and support will crash on access if 
                # they are explicitly set to None (vs just being deleted or 
                # default). 
                if value is not None:
                    out_node.add_prop(out_feature,value)

                # If root node, pull out value if not allowed on root
                if out_node.is_root:
                    if not root_allowed:
                        stash_values[out_feature] = value
                        if out_feature in ["dist","support","name"]:
                            try:
                                out_node.del_prop(out_feature)
                            except (KeyError,AttributeError):
                                pass
                        else:
                            out_node.add_prop(out_feature,None)


    # Copy ancestor to correct node because displayed by rooting
    if len(stash_values) > 0 and "anc_label" in stash_values:

        if stash_values["anc_label"] != "":

            for n in out_tree.traverse():
                if n.is_leaf:
                    continue

                if n.get_prop("anc_label") is None or n.get_prop("anc_label") == "":
                    n.add_prop("anc_label",stash_values["anc_label"])
                    if "anc_pp" in stash_values:
                        n.add_prop("anc_pp",stash_values["anc_pp"])

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

        # parser=2 --> all branches + leaf names + internal supports
        if write_bs_supports:
            out_trees.append(T.write(parser=2))

    # event
    if event:
        write_events = False
        for n in T.traverse():
            if not n.is_leaf:
                if n.get_prop("event") is not None:
                    n.name = n.get_prop("event")
                    write_events = True

        # parser=3 --> all branches + all names
        if write_events:
            out_trees.append(T.write(parser=3))

    # anc_label
    if anc_label:
        write_anc_label = False
        for n in T.traverse():
            if not n.is_leaf:
                if n.get_prop("anc_label") is not None:
                    n.name = n.get_prop("anc_label")
                    write_anc_label = True

        # parser=3 --> all branches + all names
        if write_anc_label:
            out_trees.append(T.write(parser=3))

    # anc_pp
    if anc_pp:
        write_anc_pp = False
        for n in T.traverse():
            if not n.is_leaf:
                if n.get_prop("anc_pp") is not None:
                    n.support = n.get_prop("anc_pp")
                    write_anc_pp = True

        # parser=2 --> all branches + leaf names + internal supports
        if write_anc_pp:
            out_trees.append(T.write(parser=2))

    # --------------------------------------------------------------------------
    # Finalize output

    out_trees = "\n".join(out_trees)

    if out_file is not None:
        f = open(out_file,'w')
        f.write(out_trees)
        f.close()

    return out_trees 
    