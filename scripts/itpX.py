################################
dh1 = [9,10,11,12]		# What is the first dihedral in the trimer?
dh2 = [12,13,14,15]		# What is the second dihedral in the trimer?
hs_on_monomer = [27,26]  # What are the two terminal hydrogen's on the monomer? (order matters!!)
itp_3mer = 'trimer.itp'     # What is the name of the trimer file?
itp_1mer = 'monomer.itp'    # What is the name of the monomer file?
output_itp = "p3ht_76n.itp" # what do you want to name hte output file?
n = 76  # Length of the Polymer
################################

import matplotlib.pyplot as plt
from itp_parser import Itp_parser
from MolGraph import MolecularGraph
from MolGraphIso import MolGraphIsomorphism
import networkx as nx
from collections import Counter
import re
from collections import defaultdict
import pandas as pd

def split_G3(itp_file, dh1, dh2):
    """
    Identifies the three units of a trimer from an ITP file using graph theory.

    Args:
        itp_file (str): Path to the ITP file.
        dh1 (list): First dihedral defining the first cut (edge).
        dh2 (list): Second dihedral defining the second cut (edge).

    Returns:
        tuple: (Updated atoms dictionary with 'nmer' assignments, Molecular Graph, Subgraphs)
    """

    # Step 1: Initialize MolecularGraph Class (✅ FIXED)
    mol_graph = MolecularGraph(itp_file)
    mol_graph.assign_atoms()  # Ensure atom labels are assigned ✅
    G = mol_graph.get_graph()  # Get NetworkX graph ✅

    # Step 2: Extract the atoms dictionary
    parser = mol_graph.parser  # Get the parser from the class
    atoms_dict = parser.DF["atoms"]

    if not atoms_dict:
        raise ValueError("No atoms found in the ITP file.")

    # Step 3: Initialize "nmer" as a dictionary with default 0 values
    atoms_dict["nmer"] = {int(atom[0]): 0 for atom in atoms_dict["atoms"]}

    # Step 4: Define the cut edges (middle edges of dh1 and dh2)
    cut_edges = [(dh1[1], dh1[2]), (dh2[1], dh2[2])]

    # Step 5: Check if the edges exist before removing them
    valid_cut_edges = [edge for edge in cut_edges if G.has_edge(*edge)]

    if not valid_cut_edges:
        raise ValueError(f"Cut edges {cut_edges} not found in the graph!")

    print(f"Removing edges: {valid_cut_edges}")
    G.remove_edges_from(valid_cut_edges)

    # Step 6: Identify connected components (separate units)
    components = list(nx.connected_components(G))

    # Step 7: Assign `nmer` values based on specific rules
    node_to_nmer = {}

    for i, component in enumerate(components, start=1):
        component_atoms = set(component)  # Convert to set for easy checks

        contains_dh1 = {dh1[1], dh1[2]} & component_atoms
        contains_dh2 = {dh2[1], dh2[2]} & component_atoms

        print(f"\nChecking Component {i}: {component_atoms}")
        print(f"Contains dh1 atoms: {contains_dh1}")
        print(f"Contains dh2 atoms: {contains_dh2}")

        # Identify nmer based on the correct logic
        if len(contains_dh1) == 1 and not contains_dh2:
            nmer_value = 1  # First fragment
        elif len(contains_dh2) == 1 and not contains_dh1:
            nmer_value = 3  # Third fragment
        elif len(contains_dh1) == 1 and len(contains_dh2) == 1:
            nmer_value = 2  # Middle fragment
        else:
            nmer_value = 0  # Unexpected case

        print(f"Assigned nmer {nmer_value} to Component {i}")

        # Assign `nmer` to each node in the component
        for node in component:
            node_to_nmer[node] = nmer_value

    # Step 8: Update the 'nmer' field in atoms_dict
    for atom in atoms_dict["atoms"]:
        atom_id = int(atom[0])  # Ensure it's an integer
        atoms_dict["nmer"][atom_id] = node_to_nmer.get(atom_id, 0)  # Default to 0 if not found

    # Create subgraphs based on `nmer` values
    G3_n1 = G.subgraph([node for node, nmer in node_to_nmer.items() if nmer == 1]).copy()
    G3_n2 = G.subgraph([node for node, nmer in node_to_nmer.items() if nmer == 2]).copy()
    G3_n3 = G.subgraph([node for node, nmer in node_to_nmer.items() if nmer == 3]).copy()

    # Restore G3 by adding back the removed edges
    G3_restored = G.copy()
    G3_restored.add_edges_from(valid_cut_edges)

    # 🔹 FIX: Copy node attributes explicitly
    def copy_node_attributes(source_graph, target_graph):
        """Copies node attributes from source_graph to target_graph."""
        for node in target_graph.nodes():
            if node in source_graph.nodes():
                target_graph.nodes[node].update(source_graph.nodes[node])  # Copy all attributes

    # Apply fix to all subgraphs
    for subgraph in [G3_n1, G3_n2, G3_n3, G3_restored]:
        copy_node_attributes(G, subgraph)

    return atoms_dict, G, G3_n1, G3_n2, G3_n3, G3_restored
def plot_graphs(graphs, titles):
    """
    Plots multiple graphs in a grid layout with node labels displaying atom type and node number.

    Args:
        graphs (list): List of NetworkX graphs.
        titles (list): List of titles corresponding to each graph.
    """
    num_graphs = len(graphs)
    fig, axes = plt.subplots(1, num_graphs, figsize=(5 * num_graphs, 5))

    if num_graphs == 1:
        axes = [axes]  # Ensure iterable when only one graph exists

    for ax, graph, title in zip(axes, graphs, titles):
        ax.set_title(title)

        # Step 1: Create node labels in the format "node_number (atom_type)"
        labels = {node: f"{node} ({graph.nodes[node].get('atom_name', '?')})" for node in graph.nodes()}

        # Step 2: Draw the graph with custom labels
        nx.draw(
            graph,
            ax=ax,
            with_labels=True,
            labels=labels,
            node_color="lightblue",
            edge_color="gray",
            node_size=80,
            font_size=6  # Adjust for readability
        )

    plt.show()
def process_G1(G1, hs_on_monomer):
    """
    Processes G1 by renaming specific nodes and renumbering the remaining ones in ascending order.

    Args:
        G1 (networkx.Graph): The original graph.
        hs_on_monomer (list): A two-element list specifying nodes to be renamed.

    Returns:
        networkx.Graph: The processed graph with updated node names.
    """
    a, b = hs_on_monomer  # Extract values
    G1_processed = nx.relabel_nodes(G1, {a: "s", b: "e"}, copy=True)  # Rename nodes

    # Get remaining integer nodes (excluding "s" and "e"), sorted in ascending order
    remaining_nodes = sorted(node for node in G1_processed.nodes if isinstance(node, int))

    # Create a mapping to renumber these nodes sequentially
    renumber_map = {old: new for new, old in enumerate(remaining_nodes, start=1)}

    # Apply renumbering
    G1_processed = nx.relabel_nodes(G1_processed, renumber_map, copy=True)

    return G1_processed
def identify_binding_nodes(G1_processed, G3_n1, G3_n2, G3_n3, dh1, dh2):
    """
    Identifies key binding nodes for isomorphism preprocessing.

    - Finds `e_bind`: The node connected to `e` in G1_processed.
    - Finds `s_bind`: The node connected to `s` in G1_processed.
    - Finds `g3n1_bind`: The node in G3_n1 that is either dh1[1] or dh1[2].
    - Finds `g3n3_bind`: The node in G3_n3 that is either dh2[1] or dh2[2].
    - Finds `g3n2_bind1`: The node in G3_n2 that is either dh1[1] or dh1[2].
    - Finds `g3n2_bind2`: The node in G3_n2 that is either dh2[1] or dh2[2].

    Args:
        G1_processed (nx.Graph): The processed version of G1.
        G3_n1 (nx.Graph): First fragment of G3.
        G3_n2 (nx.Graph): Second fragment of G3.
        G3_n3 (nx.Graph): Third fragment of G3.
        dh1 (list): First dihedral defining the first cut (edge).
        dh2 (list): Second dihedral defining the second cut (edge).

    Prints:
        Identified key nodes with their labels.
    """

    # Identify e_bind (node in G1_proc that had an edge to e)
    e_bind = None
    if "e" in G1_processed:
        e_neighbors = list(G1_processed.neighbors("e"))
        if len(e_neighbors) == 1:
            e_bind = e_neighbors[0]
        else:
            print(f"⚠️ Unexpected number of neighbors for 'e': {e_neighbors}")

    # Identify s_bind (node in G1_proc that had an edge to s)
    s_bind = None
    if "s" in G1_processed:
        s_neighbors = list(G1_processed.neighbors("s"))
        if len(s_neighbors) == 1:
            s_bind = s_neighbors[0]
        else:
            print(f"⚠️ Unexpected number of neighbors for 's': {s_neighbors}")

    # Identify g3n1_bind (dh1[1] or dh1[2] in G3_n1)
    g3n1_bind = dh1[1] if dh1[1] in G3_n1 else dh1[2] if dh1[2] in G3_n1 else None

    # Identify g3n3_bind (dh2[1] or dh2[2] in G3_n3)
    g3n3_bind = dh2[1] if dh2[1] in G3_n3 else dh2[2] if dh2[2] in G3_n3 else None

    # Identify g3n2_bind1 (dh1[1] or dh1[2] in G3_n2)
    g3n2_bind1 = dh1[1] if dh1[1] in G3_n2 else dh1[2] if dh1[2] in G3_n2 else None

    # Identify g3n2_bind2 (dh2[1] or dh2[2] in G3_n2)
    g3n2_bind2 = dh2[1] if dh2[1] in G3_n2 else dh2[2] if dh2[2] in G3_n2 else None

    # Print the identified values
    print("\n🔍 Identified Binding Nodes:")
    print(f"e_bind (Connected to 'e' in G1_processed): {e_bind}")
    print(f"s_bind (Connected to 's' in G1_processed): {s_bind}")
    print(f"g3n1_bind (dh1[1] or dh1[2] in G3_n1): {g3n1_bind}")
    print(f"g3n3_bind (dh2[1] or dh2[2] in G3_n3): {g3n3_bind}")
    print(f"g3n2_bind1 (dh1[1] or dh1[2] in G3_n2): {g3n2_bind1}")
    print(f"g3n2_bind2 (dh2[1] or dh2[2] in G3_n2): {g3n2_bind2}")

    return {
        "e_bind": e_bind,
        "s_bind": s_bind,
        "g3n1_bind": g3n1_bind,
        "g3n3_bind": g3n3_bind,
        "g3n2_bind1": g3n2_bind1,
        "g3n2_bind2": g3n2_bind2
    }
def process_G3(G1_processed, G3_n1, G3_n2, G3_n3, dh1, dh2):
    """
    Runs graph isomorphism checks between modified versions of G1_processed and the three G3 subgraphs.

    - G3_n1: Uses a copy of G1_processed with `e` removed and initial mapping `{e_bind: g3n1_bind}`.
    - G3_n3: Uses a copy of G1_processed with `s` removed and initial mapping `{s_bind: g3n3_bind}`.
    - G3_n2: Uses a copy of G1_processed with both `e` and `s` removed and initial mapping `{e_bind: g3n2_bind1, s_bind: g3n2_bind2}`.

    Args:
        G1_processed (nx.Graph): The processed version of G1.
        G3_n1 (nx.Graph): First fragment of G3.
        G3_n2 (nx.Graph): Second fragment of G3.
        G3_n3 (nx.Graph): Third fragment of G3.
        dh1 (list): First dihedral defining the first cut (edge).
        dh2 (list): Second dihedral defining the second cut (edge).

    Returns:
        dict: A dictionary containing the isomorphism results and mappings.
    """

    # Step 1: Identify the binding nodes
    binding_nodes = identify_binding_nodes(G1_processed, G3_n1, G3_n2, G3_n3, dh1, dh2)

    e_bind = binding_nodes["e_bind"]
    s_bind = binding_nodes["s_bind"]
    g3n1_bind = binding_nodes["g3n1_bind"]
    g3n3_bind = binding_nodes["g3n3_bind"]
    g3n2_bind1 = binding_nodes["g3n2_bind1"]
    g3n2_bind2 = binding_nodes["g3n2_bind2"]

    # Step 2: Define initial mappings and modified graphs for each isomorphism test
    tests = {
        "G3_n1": {
            "graph": G3_n1,
            "G1_copy": G1_processed.copy(),
            "nodes_to_remove": ["e"],
            "initial_mapping": {e_bind: g3n1_bind} if e_bind and g3n1_bind else None
        },
        "G3_n3": {
            "graph": G3_n3,
            "G1_copy": G1_processed.copy(),
            "nodes_to_remove": ["s"],
            "initial_mapping": {s_bind: g3n3_bind} if s_bind and g3n3_bind else None
        },
        "G3_n2": {
            "graph": G3_n2,
            "G1_copy": G1_processed.copy(),
            "nodes_to_remove": ["e", "s"],
            "initial_mapping": {e_bind:g3n2_bind2 , s_bind: g3n2_bind1}
            if e_bind and s_bind and g3n2_bind1 and g3n2_bind2 else None
        }
    }

    results = {}

    for name, test in tests.items():
        print(f"\n🔍 Running Isomorphism Check: {name} vs G1_processed")

        # Remove the required nodes from the copied G1_processed graph
        G1_modified = test["G1_copy"]
        for node in test["nodes_to_remove"]:
            if node in G1_modified:
                G1_modified.remove_node(node)
                print(f"🛑 Removed node '{node}' from G1_processed copy for {name}")

        initial_mapping = test["initial_mapping"]
        if initial_mapping is None:
            print(f"⚠️ Skipping {name} due to missing required mappings.")
            results[name] = {"is_isomorphic": False, "mapping": {}}
            continue

        print(f"🔹 Using Initial Mapping: {initial_mapping}")

        # Instantiate the isomorphism checker
        mg_iso = MolGraphIsomorphism(G1_modified, test["graph"], initial_mapping)

        # Run isomorphism check
        is_iso, mapping_dict = mg_iso.check_isomorphism()

        # Store results
        results[name] = {
            "is_isomorphic": is_iso,
            "mapping": mapping_dict
        }

    return results
def finalize_G3_mappings(results):
    """
    Modifies the isomorphism mappings by appending 'A' for G3_n1, 'B' for G3_n2, and 'C' for G3_n3.
    Keeps non-numeric keys (like 's' and 'e') unchanged.

    Args:
        results (dict): The output dictionary from `process_G3` containing mappings.

    Returns:
        dict: A dictionary with modified mappings.
    """

    # Define suffixes for each mapping
    suffixes = {
        "G3_n1": "A",
        "G3_n2": "B",
        "G3_n3": "C"
    }

    # Create a new dictionary for modified results
    finalized_results = {}

    for name, result in results.items():
        suffix = suffixes.get(name, "")
        mapping = result["mapping"]

        # Modify the mapping
        modified_mapping = {
            f"{key}{suffix}" if isinstance(key, int) else key: value
            for key, value in mapping.items()
        }

        # Store the modified mapping
        finalized_results[name] = {
            "is_isomorphic": result["is_isomorphic"],
            "mapping": modified_mapping
        }

    return finalized_results
def finalize_G3(G3, finalized_results):
    """
    Finalizes G3 by renaming all nodes based on the reversed mappings from finalized_results.
    This converts numeric nodes in G3 back to their original names and creates a DataFrame
    with the renaming mappings.

    Args:
        G3 (nx.Graph): The full G3 graph.
        finalized_results (dict): The mappings from G3_n1, G3_n2, and G3_n3 obtained from `finalize_G3_mappings`.

    Returns:
        nx.Graph: The finalized G3 graph with renamed nodes.
        pd.DataFrame: A DataFrame (G3_remapped) containing the old and new node names.
    """

    # Create a reversed mapping dictionary to rename nodes in G3
    node_mapping = {}

    # Loop through finalized results to get renaming mappings
    for fragment_name, result in finalized_results.items():
        mapping = result["mapping"]

        for original_name, numeric_id in mapping.items():
            numeric_id_str = str(numeric_id)  # Ensure keys are strings for NetworkX compatibility

            # Store the reverse mapping {numeric_id: original_name}
            node_mapping[numeric_id_str] = original_name

    # Log the final renaming mapping
    print("\n✅ Final G3 Node Renaming Mapping:")
    for old_name, new_name in node_mapping.items():
        print(f"  {old_name} → {new_name}")

    # Apply the renaming to G3
    G3 = nx.relabel_nodes(G3, node_mapping, copy=True)

    # Create a DataFrame with the renaming mappings
    G3_remapped = pd.DataFrame(list(node_mapping.items()), columns=["Old", "New"])

    return G3, G3_remapped
def normalize_dict_length(data_dict):
    """
    Normalizes dictionary values to the same length by padding shorter lists with None.

    Args:
        data_dict (dict): Dictionary where values are lists.

    Returns:
        dict: Dictionary with equal-length lists.
    """
    max_length = max(len(v) if isinstance(v, list) else 0 for v in data_dict.values())

    for key in data_dict:
        if isinstance(data_dict[key], list):
            data_dict[key] += [None] * (max_length - len(data_dict[key]))  # Pad with None

    return data_dict
def process_itp_to_df(parser):
    """
    Processes an ITP file using Itp_parser, extracts all relevant sections, converts them to DataFrames,
    and prints metadata including the name, type, column names, and data types.

    Args:
        parser (Itp_parser): An instance of Itp_parser with parsed data.

    Returns:
        dict: A dictionary where keys are section names and values are Pandas DataFrames.
    """
    # Extract dictionaries from the ITP parser
    extracted_dicts = parser.DF  # Get the extracted dictionaries
    dataframes_dict = {}

    print("\n📂 Extracted ITP Sections:")

    # Convert extracted dictionaries to Pandas DataFrames
    for dict_name, data_dict in extracted_dicts.items():
        try:
            normalized_dict = normalize_dict_length(data_dict)  # Ensure uniform length
            df = pd.DataFrame(normalized_dict)  # Convert to DataFrame
            dataframes_dict[dict_name] = df  # Store in dictionary
        except Exception as e:
            print(f"⚠️ Could not convert {dict_name} to DataFrame: {e}")

    # Identify DataFrame types based on Itp_parser.TOTAL_DIRECTIVES
    directive_types = Itp_parser.TOTAL_DIRECTIVES  # Access known ITP section types

    print("\n📝 DataFrame Types & Metadata:")
    for df_name, df in dataframes_dict.items():
        df_type = "Unknown"  # Default type if not found in TOTAL_DIRECTIVES

        for category, section_dict in directive_types.items():
            if df_name in section_dict:
                df_type = category
                break  # Stop searching once we find a match

        print(f"\n📌 DataFrame Name: {df_name}")
        print(f"🔹 Type: {df_type}")
        print(f"🔸 Columns: {list(df.columns)}")
        print(f"🔹 Data Types:\n{df.dtypes}")

    return dataframes_dict
def remap_itp_data(dataframes_dict, G3_remapped, parser):
    """
    Remaps atom indices in the relevant DataFrames using the G3_remapped mapping.

    Args:
        dataframes_dict (dict): Dictionary where keys are section names and values are Pandas DataFrames.
        G3_remapped (pd.DataFrame): DataFrame containing "Old" and "New" mapping for atom renaming.
        parser (Itp_parser): An instance of Itp_parser to classify DataFrame types.

    Returns:
        dict: Updated dictionary of DataFrames with remapped "atoms" columns.
    """
    # Define the sections that require atom renaming
    atom_related_types = {
        "BONDS", "DIHEDRALS", "ANGLES", "CONSTRAINTS", "PAIRS", "PAIRS_NB",
        "EXCLUSIONS", "VIRTUAL_SITES1", "VIRTUAL_SITES2", "VIRTUAL_SITES3",
        "VIRTUAL_SITES4", "VIRTUAL_SITESn", "SETTLE", "POSITION_RESTRAINT",
        "DISTANCE_RESTRAINTS", "ORIENTATION_RESTRAINTS", "ANGLES_RESTRAINTS"
    }

    # Create a dictionary for quick mapping lookup
    rename_mapping = dict(zip(G3_remapped["Old"].astype(str), G3_remapped["New"]))

    # Access known ITP section types from the parser
    directive_types = Itp_parser.TOTAL_DIRECTIVES

    # Process each DataFrame in dataframes_dict
    for df_name, df in dataframes_dict.items():
        # Identify the type of DataFrame
        df_type = "Unknown"
        for category, section_dict in directive_types.items():
            if df_name in section_dict:
                df_type = category
                break

        # Check if this DataFrame needs atom remapping
        if df_name == "atoms" or df_type in atom_related_types:
            if "atoms" in df.columns:
                # Apply mapping to the "atoms" column, converting atom indices
                def map_atoms(atom_list):
                    if isinstance(atom_list, list):
                        return [rename_mapping.get(str(atom), atom) for atom in atom_list]
                    return atom_list  # Keep original value if not a list

                # Create new "atoms_new" column with remapped values
                dataframes_dict[df_name]["atoms_new"] = df["atoms"].apply(map_atoms)

                print(f"✅ Remapped 'atoms' column in '{df_name}' DataFrame.")

    return dataframes_dict
def assign_df_set(dataframes_dict):
    """
    Assigns a 'df_set' column to each DataFrame based on unique atom types in 'atoms_new',
    following a strict order of 'sCBAe'.

    Args:
        dataframes_dict (dict): Dictionary where keys are section names and values are Pandas DataFrames.

    Returns:
        dict: Updated dictionary with 'df_set' column added to relevant DataFrames.
    """

    # Function to determine df_set value for a row
    def determine_df_set(atom_list):
        if not isinstance(atom_list, list):
            return None  # Skip if atoms_new is missing or not a list

        # Define order precedence
        order = ["s", "A", "B", "C", "e"]

        # Extract unique label characters from atom_list
        unique_labels = set()
        for atom in atom_list:
            if isinstance(atom, str):  # Ensure it's a string before processing
                for char in atom:  # Extract individual characters
                    if char in order:
                        unique_labels.add(char)

        # Maintain strict order: sCBAe
        ordered_labels = "".join([char for char in order if char in unique_labels])

        return ordered_labels if ordered_labels else None

    # Process each DataFrame
    for df_name, df in dataframes_dict.items():
        if "atoms_new" in df.columns:  # Ensure 'atoms_new' exists
            df["df_set"] = df["atoms_new"].apply(determine_df_set)
            print(f"✅ Assigned 'df_set' column for '{df_name}' DataFrame.")

    return dataframes_dict
def generate_polymer_list(n): # This is where we plan out the sections required to build the polymer. each string here is listing the unit it is connecting. S = starting Hydrogen, E = Ending Hydrogen, A = unit 1, no Hydrogen, C = last unit, no Hydrogen. B1 = unit 2, B2 = unit 3, ..., Bn-2 = second to last unit, So if a string here is B1, that is all forces needed for unit 1 with no connecting parameters. If a string here is SAB1 then these are all parameters going from the first hydrogen, to unit 1, to unit 2.
    """
    Generates a list of polymer unit strings based on input n, ensuring the presence of required values.

    Args:
        n (int): The number of polymer units.

    Returns:
        list: A list containing the required values and additional B-related values.
    """
    if n < 3:
        return "Please choose a polymer with at least 3 units."

    # Required values
    required_values = ["s", "A", "C", "e", "sA", "Ce"]

    # Generate B values from B1 to B(n-2)
    b_values = [f"B{i}" for i in range(1, n - 1)]

    # Core additional required combinations (always present)
    additional_values = [
        "sB1", "sAB1",
        f"B{n-2}e", f"B{n-2}Ce",
        "AB1", f"B{n-2}C"
    ]

    # Special case for n == 3
    if n == 3:
        additional_values.extend(["se", "sCe", "sB1e", "sAe", "sAB1e", "sACe", "sB1Ce",
                                  "sAB1Ce", "AC", "AB1C", "sC", "sB1C", "sAC", "sAB1C",
                                  "Ae", "AB1e", "ACe", "AB1Ce"])

    # Extra values when n >= 4
    if n >= 4:
        extra_values = [
            "sAB2", "sAB1B2", "sB2",
            f"B{n-3}B{n-2}Ce",
            "AB1B2", "AB2",
            f"B{n-3}B{n-2}C", f"B{n-3}C",
            f"sB1B2", f"B1B2",f"B{n-3}e",
            f"B{n-3}Ce",f"B{n-3}B{n-2}e"
        ]
        additional_values += extra_values

    # **New B-Combination Rules for n >= 5**
    if n >= 4:
        for a in range(1, n - 3):  # Iterate from 1 to (n-3)
            b1 = f"B{a+1}B{a+2}"
            b2 = f"B{a}B{a+2}"
            b3 = f"B{a}B{a+1}B{a+2}"

            additional_values.extend([b1, b2, b3])
    # Combine all lists
    polymer_list = required_values + b_values + additional_values

    return polymer_list
def analyze_polymer_sections(dataframes_dict, polymer_sections, n):
    """
    Identifies:
    - Unique df_set values present in the DataFrames.
    - Unique polymer section labels from polymer_sections.
    - Which polymer_sections labels are missing from df_set values.
    - Generates copies of '_B', '_AB', and '_BC' DataFrames as needed.
    - Generates additional copies of '_BC' as 'B(a)B(a+1)' for a in range(1, n-1).
    - Ensures newly created sections are no longer missing.
    - Prints n.

    Args:
        dataframes_dict (dict): Dictionary where keys are section names and values are Pandas DataFrames.
        n (int): An integer value to determine the number of copies.

    Returns:
        dict: Contains 'df_set_values', 'polymer_section_labels', 'missing_polymer_labels', 'n', and updated 'dataframes_dict'.
    """

    # Extract unique df_set values from all DataFrames
    df_set_values = set()
    for df in dataframes_dict.values():
        if 'df_set' in df.columns:
            df_set_values.update(df['df_set'].dropna().unique())  # Collect full unique df_set values

    # Convert polymer_sections list to a set for comparison
    polymer_section_labels = set(polymer_sections)

    # Generate copies of '_B' DataFrames for B1 to B(n-2)
    b_copies = {}
    b_count = n - 2  # We want B1 to B(n-2)

    for df_name, df in dataframes_dict.items():
        if df_name.endswith("_B"):  # Identify '_B' DataFrames
            for i in range(1, b_count + 1):  # Generate B1 to B(n-2)
                new_name = f"{df_name}{i}"  # Example: atoms_B1, atoms_B2, ...
                b_copies[new_name] = df.copy()  # Create a copy of the DataFrame
                df_set_values.add(f"B{i}")  # Add new df_set values

    # Generate copies of '_AB' DataFrames as '_AB1'
    ab_copies = {}
    for df_name, df in dataframes_dict.items():
        if df_name.endswith("_AB"):  # Identify '_AB' DataFrames
            new_name = f"{df_name}1"  # Example: atoms_AB1
            ab_copies[new_name] = df.copy()  # Create a copy of the DataFrame
            df_set_values.add("AB1")  # Add 'AB1' to df_set values

    # Generate copies of '_BC' DataFrames as '_B(n-2)C'
    bc_copies = {}
    for df_name, df in dataframes_dict.items():
        if df_name.endswith("_BC"):  # Identify '_BC' DataFrames
            new_name = f"{df_name[:-2]}B{b_count}C"  # Example: atoms_B6C if n=8
            bc_copies[new_name] = df.copy()  # Create a copy of the DataFrame
            df_set_values.add(f"B{b_count}C")  # Add 'B(n-2)C' to df_set values

    # Generate additional copies of '_BC' as 'B(a)B(a+1)' for a in range(1, n-1)
    bc_transition_copies = {}
    for df_name, df in dataframes_dict.items():
        if df_name.endswith("_BC"):  # Identify '_BC' DataFrames
            for a in range(1, n-2):  # Generate B1B2 to B(n-2)B(n-1)
                new_name = f"{df_name[:-2]}B{a}B{a+1}"  # Example: atoms_B1B2, atoms_B2B3, ..., atoms_B(n-2)B(n-1)
                bc_transition_copies[new_name] = df.copy()  # Create a copy of the DataFrame
                df_set_values.add(f"B{a}B{a+1}")  # Add 'B(a)B(a+1)' to df_set values

    # Merge newly created copies into the original dictionary
    dataframes_dict.update(b_copies)
    dataframes_dict.update(ab_copies)
    dataframes_dict.update(bc_copies)
    dataframes_dict.update(bc_transition_copies)

    # Identify missing values (present in polymer_sections but missing from updated df_set)
    missing_polymer_labels = polymer_section_labels - df_set_values

    # Print results
    print(f"✅ Unique df_set values in DataFrames: {sorted(df_set_values)}")
    print(f"✅ Labels present in polymer_sections list: {sorted(polymer_section_labels)}")
    print(f"❌ Missing polymer_sections labels (not found in df_set): {sorted(missing_polymer_labels)}")
    print(f"✅ Value of n: {n}")
    print(f"✅ Created {b_count} copies for each '_B' DataFrame (B1 to B{b_count}).")
    print(f"✅ Created copies for '_AB' DataFrames as '_AB1'.")
    print(f"✅ Created copies for '_BC' DataFrames as '_B{b_count}C'.")
    print(f"✅ Created additional '_BC' copies as 'B(a)B(a+1)' for a in range(1, {n-2}).")

    # Return results as a dictionary for further use if needed
    return {
        "df_set_values": sorted(df_set_values),
        "polymer_section_labels": sorted(polymer_section_labels),
        "missing_polymer_labels": sorted(missing_polymer_labels),
        "n": n,
        "updated_dataframes_dict": dataframes_dict
    }
def clean_dfs_AB1(dataframes_dict):
    """
    Finds all DataFrames matching 'AB1', modifies them by:
    - Replacing 'B' with 'B1' in the 'atoms_new' column.
    - Setting all values in the 'df_set' column to 'AB1'.
    - Printing the modified DataFrames.

    Args:
        dataframes_dict (dict): Dictionary where keys are section names and values are Pandas DataFrames.

    Returns:
        None
    """
    section_type = "AB1"  # Remove the `$` since it's a regex anchor

    # Use re.fullmatch() to ensure an exact match at the end
    matching_keys = [key for key in dataframes_dict if re.fullmatch(f".*_{section_type}", key)]

    if matching_keys:
        print(f"✅ Found {len(matching_keys)} DataFrames matching '{section_type}':")
        for key in matching_keys:
            df = dataframes_dict[key]

            # Ensure 'atoms_new' and 'df_set' columns exist before modifying
            if 'atoms_new' in df.columns and 'df_set' in df.columns:
                # Replace 'B' with 'B1' in atoms_new
                df['atoms_new'] = df['atoms_new'].apply(lambda x: [atom.replace('B', 'B1') for atom in x] if isinstance(x, list) else x)

                # Change all values in df_set to 'AB1'
                df['df_set'] = 'AB1'

                print(f"🔹 Modified {key}:")
                print(df.head())  # Print first few rows to verify changes
                print("\n" + "-" * 40 + "\n")
    else:
        print(f"❌ No DataFrames found for '{section_type}'.")
def clean_dfs_Bn2C(dataframes_dict, n):
    """
    Finds all DataFrames matching 'B(n-2)C', modifies them by:
    - Replacing 'B' with 'B(n-2)' in the 'atoms_new' column.
    - Setting all values in the 'df_set' column to 'B(n-2)C'.
    - Printing the modified DataFrames.

    Args:
        dataframes_dict (dict): Dictionary where keys are section names and values are Pandas DataFrames.
        n (int): The predefined value used to compute 'B(n-2)'. Default is 8.

    Returns:
        None
    """
    b_n2 = f"B{n-2}"  # Compute B(n-2), e.g., B6 if n=8
    section_type = f"{b_n2}C"  # Example: _B6C if n=8

    # Use re.fullmatch() to ensure an exact match at the end
    matching_keys = [key for key in dataframes_dict if re.fullmatch(f".*_{section_type}", key)]

    if matching_keys:
        print(f"✅ Found {len(matching_keys)} DataFrames matching '{section_type}':")
        for key in matching_keys:
            df = dataframes_dict[key]

            # Ensure 'atoms_new' and 'df_set' columns exist before modifying
            if 'atoms_new' in df.columns and 'df_set' in df.columns:
                # Replace 'B' with 'B(n-2)' in atoms_new
                df['atoms_new'] = df['atoms_new'].apply(lambda x: [atom.replace('B', b_n2) for atom in x] if isinstance(x, list) else x)

                # Change all values in df_set to 'B(n-2)C'
                df['df_set'] = section_type

                print(f"🔹 Modified {key}:")
                print(df.head())  # Print first few rows to verify changes
                print("\n" + "-" * 40 + "\n")
    else:
        print(f"❌ No DataFrames found for '{section_type}'.")
def clean_dfs_Ba(dataframes_dict, n):
    """
    Finds all DataFrames matching 'Ba' where a is an integer from 1 to n-2.
    Modifies them by:
    - Replacing 'B' with 'Ba' in the 'atoms_new' column.
    - Setting all values in the 'df_set' column to 'Ba'.
    - Printing the modified DataFrames.

    Args:
        dataframes_dict (dict): Dictionary where keys are section names and values are Pandas DataFrames.
        n (int): The predefined value used to determine the range of 'Ba'. Default is 8.

    Returns:
        None
    """
    for a in range(1, n-1):  # a ranges from 1 to n-2
        ba_label = f"B{a}"  # Example: _B1, _B2, ..., _B(n-2)

        # Use re.fullmatch() to ensure the entire string matches exactly
        matching_keys = [key for key in dataframes_dict if re.fullmatch(f".*_{ba_label}", key)]

        if matching_keys:
            print(f"✅ Found {len(matching_keys)} DataFrames matching '{ba_label}':")
            for key in matching_keys:
                df = dataframes_dict[key]

                # Ensure 'atoms_new' and 'df_set' columns exist before modifying
                if 'atoms_new' in df.columns and 'df_set' in df.columns:
                    # Replace 'B' with 'Ba' (B1, B2, ..., B(n-2)) in atoms_new
                    df['atoms_new'] = df['atoms_new'].apply(lambda x: [atom.replace('B', ba_label) for atom in x] if isinstance(x, list) else x)

                    # Change all values in df_set to 'Ba'
                    df['df_set'] = ba_label

                    print(f"🔹 Modified {key}:")
                    print(df.head())  # Print first few rows to verify changes
                    print("\n" + "-" * 40 + "\n")
        else:
            print(f"❌ No DataFrames found for '{ba_label}'.")
def clean_dfs_Ba1Ba2(dataframes_dict, n):
    """
    Finds all DataFrames matching 'B(a)B(a+1)' where a ranges from 1 to n-3.
    Modifies them by:
    - Replacing 'B' with 'B(a)' and 'C' with 'B(a+1)' in the 'atoms_new' column.
    - Setting all values in the 'df_set' column to 'B(a)B(a+1)'.
    - Printing the modified DataFrames.

    Args:
        dataframes_dict (dict): Dictionary where keys are section names and values are Pandas DataFrames.
        n (int): The predefined value used to determine the range of 'B(a)B(a+1)'. Default is 8.

    Returns:
        None
    """
    for a in range(1, n-2):  # a ranges from 1 to n-3
        ba_ba1_label = f"_B{a}B{a+1}"  # Example: _B1B2, _B2B3, ..., _B(n-3)B(n-2)

        # Use re.fullmatch() to ensure an exact match at the end
        matching_keys = [key for key in dataframes_dict if re.fullmatch(f".*{ba_ba1_label}", key)]

        if matching_keys:
            print(f"✅ Found {len(matching_keys)} DataFrames matching '{ba_ba1_label}':")
            for key in matching_keys:
                df = dataframes_dict[key]

                # Ensure 'atoms_new' and 'df_set' columns exist before modifying
                if 'atoms_new' in df.columns and 'df_set' in df.columns:
                    # Replace 'B' with 'B{a}' and 'C' with 'B{a+1}' in atoms_new
                    df['atoms_new'] = df['atoms_new'].apply(
                        lambda x: [atom.replace('B', f'B{a}').replace('C', f'B{a+1}') for atom in x] if isinstance(x, list) else x
                    )

                    # Change all values in df_set to 'B(a)B(a+1)'
                    df['df_set'] = f"B{a}B{a+1}"

                    print(f"🔹 Modified {key}:")
                    print(df.head())  # Print first few rows to verify changes
                    print("\n" + "-" * 40 + "\n")
        else:
            print(f"❌ No DataFrames found for '{ba_ba1_label}'.")
def build_polymer_mapping(dataframes_dict, n):
    """
    Constructs and processes atomic mapping for a polymer by concatenating dataframes in a specific order,
    extracting atomic identifiers, assigning sequence positions, and ensuring proper sorting.

    Args:
        dataframes_dict (dict): Dictionary containing atomic DataFrames.
        n (int): Defines the number of B-sections in the polymer chain (default is 8).

    Returns:
        pd.DataFrame: Processed mapping DataFrame with final atom assignments.
    """

    # Step 1: Concatenate Atomic DataFrames
    ordered_names = ["atoms_s", "atoms_A"] + [f"atoms_B{i}" for i in range(1, n-1)] + ["atoms_C", "atoms_e"]
    dataframes_to_concat = [dataframes_dict[name] for name in ordered_names if name in dataframes_dict]

    if not dataframes_to_concat:
        print("❌ No matching DataFrames found to append.")
        return None

    atoms_df = pd.concat(dataframes_to_concat, ignore_index=True).dropna(subset=['atoms'])

    # Create Initial Mapping
    mapping = pd.DataFrame({
        "atoms": atoms_df["atoms"],
        "atoms_new": atoms_df["atoms_new"],
        "atoms_final": "TBD"  # Placeholder for final assignments
    })

    # Step 2: Extract Atom Identifiers
    mapping['atoms_new'] = mapping['atoms_new'].apply(lambda x: x[0] if isinstance(x, list) else x)

    def split_atom_identifier(value):
        """Splits atom identifier into integer and letter components."""
        if value in ['s', 'e']:  # Special cases
            return (1, value)
        match = re.match(r'(\d+)?([A-Z]+\d*)', value)
        if match:
            return (int(match.group(1)) if match.group(1) else None, match.group(2))
        return (None, None)

    mapping[['integer', 'letter']] = mapping['atoms_new'].apply(lambda x: pd.Series(split_atom_identifier(x)))

    # Step 3: Assign Polymer Sequence Positions (nmer)
    nmer_values = {'s': 0, 'A': 1, 'C': n, 'e': n + 1}

    for index, row in mapping.iterrows():
        if row['letter'].startswith('B'):
            match = re.search(r'B(\d+)', row['letter'])
            if match:
                nmer_values[row['letter']] = int(match.group(1)) + 1

    mapping['nmer'] = mapping['letter'].map(nmer_values)

    # Step 4: Compute Number of Atoms per Unit
    num_atoms = len(mapping)
    atoms_per_unit = (num_atoms - 2) / n

    # Step 5: Assign Final Atom Order
    def compute_final_atom_number(row):
        """Computes the final atom number assignment."""
        if row['letter'] == 's':
            return 1
        if row['letter'] == 'e':
            return num_atoms
        return (row['nmer'] - 1) * atoms_per_unit + row['integer'] + 1

    mapping.loc[mapping['atoms_final'] == 'TBD', 'atoms_final'] = mapping.apply(compute_final_atom_number, axis=1)
    mapping['atoms_final'] = mapping['atoms_final'].astype(int)

    # Step 6: Sort by Final Atom Order
    mapping = mapping.sort_values(by='atoms_final').reset_index(drop=True)

    # Step 7: Final Cleaning
    mapping.drop(columns=['integer', 'letter'], inplace=True)
    mapping['atoms'] = mapping['atoms'].apply(lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x)

    # Step 8: Adjust nmer Labels
    max_nmer = mapping['nmer'].max()
    mapping['nmer'] = mapping['nmer'].replace({0: 'H1', max_nmer: 'H2'})

    return mapping
def clear_irrelevant_dfs(dataframes_dict):
    """
    Categorizes and prints DataFrame names into three groups:
    1. Names **without** an underscore `_`.
    2. Names **with `_` and containing 'B'** **not immediately followed** by a number at some point after `_`.
    3. Remaining DataFrames.

    Args:
        dataframes_dict (dict): Dictionary of DataFrame names.
    """
    no_underscore = []
    b_not_followed_by_number = []
    other = []

    for name in dataframes_dict.keys():
        if "_" not in name:
            no_underscore.append(name)
        elif re.search(r"_.*B(?!\d)", name):  # Looks for 'B' after '_' not followed by a number
            b_not_followed_by_number.append(name)
        else:
            other.append(name)

    # Print results
    print("\n📂 **DataFrames without an underscore (_):**")
    print(no_underscore if no_underscore else "None found.")

    print("\n🔹 **DataFrames with 'B' after `_`, NOT followed by a number:**")
    print(b_not_followed_by_number if b_not_followed_by_number else "None found.")

    print("\n📁 **Other DataFrames:**")
    print(other if other else "None found.")
def cut_final_dfs(dataframes_dict):
    """
    Removes all DataFrames from the dictionary where:
    - The name contains an underscore `_`.
    - The name has 'B' appearing after `_`, and it is **not immediately followed by a number**.

    Args:
        dataframes_dict (dict): Dictionary containing DataFrames.

    Returns:
        dict: A new dictionary with the filtered DataFrames.
    """
    filtered_dict = {
        name: df for name, df in dataframes_dict.items()
        if not re.search(r"_.*B(?!\d)", name)  # Matches 'B' after '_' not followed by a number
    }

    print("✅ Removed DataFrames where 'B' appears after '_' and is NOT followed by a number.")
    return filtered_dict
def find_invalid_B_entries(split_dataframes):
    # Iterate through each DataFrame in the dictionary
    for df_name, df in split_dataframes.items():
        # Use regex to find rows where 'B' is present but not followed by a number
        mask = df['atoms_new'].astype(str).str.contains(r'B(?!\d)', regex=True)

        # Extract the rows that match the condition
        invalid_rows = df[mask]

        # Print results if any invalid rows exist
        if not invalid_rows.empty:
            print(f"🔍 DataFrame: {df_name}")
            print(invalid_rows)
            print("-" * 40)
def concatenate_dfs(dataframes_dict):
    """
    Groups and concatenates DataFrames that share the same prefix before the first underscore '_'.

    Args:
        dataframes_dict (dict): Dictionary containing DataFrames with names formatted as 'category_subcategory'.

    Returns:
        dict: A new dictionary with concatenated DataFrames grouped by their common prefix.
    """
    grouped_dataframes = defaultdict(list)

    # Group DataFrames based on the prefix before the first underscore
    for name, df in dataframes_dict.items():
        match = re.match(r"([^_]+)_.*", name)  # Extract prefix before first '_'
        if match:
            prefix = match.group(1)
            grouped_dataframes[prefix].append(df)

    # Concatenate DataFrames within each group
    final_dfs_dict = {}
    for prefix, dfs in grouped_dataframes.items():
        final_dfs_dict[prefix] = pd.concat(dfs, ignore_index=True)

    return final_dfs_dict
def process_dataframes(dataframes_dict, mapping):
    """
    Processes each DataFrame in the given dictionary by:
    1. Removing the 'atoms' column.
    2. Renaming 'atoms_new' to 'atoms'.
    3. Removing the 'df_set' column.
    4. Replacing values in the 'atoms' column based on the mapping DataFrame.

    Args:
        dataframes_dict (dict): Dictionary of DataFrames to process.
        mapping (pd.DataFrame): Mapping DataFrame with 'atoms_new' and 'atoms_final' columns.

    Returns:
        dict: A dictionary with updated DataFrames.
    """

    # Create a mapping dictionary from atoms_new to atoms_final
    mapping_dict = dict(zip(mapping['atoms_new'], mapping['atoms_final']))

    processed_dfs_dict = {}

    for name, df in dataframes_dict.items():
        print(f"🔹 Processing DataFrame: {name}")

        # Drop 'atoms' and 'df_set' columns if they exist
        df = df.drop(columns=['atoms'], errors='ignore')
        df = df.drop(columns=['df_set'], errors='ignore')

        # Rename 'atoms_new' to 'atoms' if it exists
        if 'atoms_new' in df.columns:
            df = df.rename(columns={'atoms_new': 'atoms'})

        # Replace values in the 'atoms' column using the mapping
        if 'atoms' in df.columns:
            df['atoms'] = df['atoms'].apply(lambda x: [mapping_dict.get(item, item) for item in x] if isinstance(x, list) else x)

        # Store the processed DataFrame
        processed_dfs_dict[name] = df

    return processed_dfs_dict
def reorder_processed_dfs(dataframes_dict):
    """
    Sorts and reindexes each DataFrame in processed_dfs_dict based on the first value of the 'atoms' column.
    Converts all values in the 'atoms' lists to integers.

    Args:
        dataframes_dict (dict): Dictionary containing DataFrames.

    Returns:
        dict: Updated dictionary with sorted and reindexed DataFrames.
    """
    updated_dfs_dict = {}

    for df_name, df in dataframes_dict.items():
        if 'atoms' in df.columns:
            # Ensure 'atoms' values are lists of integers
            df['atoms'] = df['atoms'].apply(lambda lst: [int(x) for x in lst] if isinstance(lst, list) else lst)

            # Create a temporary column for sorting based on the first value in the list
            df['_sort_key'] = df['atoms'].apply(lambda lst: lst[0] if isinstance(lst, list) and len(lst) > 0 else float('inf'))

            # Sort DataFrame based on the temporary column
            df = df.sort_values(by='_sort_key')

            # Drop the temporary column after sorting
            df = df.drop(columns=['_sort_key'])

            # Reset index after sorting
            df = df.reset_index(drop=True)

        updated_dfs_dict[df_name] = df  # Store the processed DataFrame

    print("✅ Successfully sorted and reindexed all DataFrames in processed_dfs_dict.")
    return updated_dfs_dict
def renumber_atom_names(dataframes_dict):
    """
    Processes the 'atoms' DataFrame in processed_dfs_dict:
    - Extracts the atomic symbol from the 'atom_name' column (a list with one value).
    - Generates a unique hex-style numbering sequence.
    - Updates the 'atom_name' column with the new numbering while maintaining the list format.
    - Keeps all other columns intact.

    Args:
        dataframes_dict (dict): Dictionary containing DataFrames, including 'atoms'.

    Returns:
        dict: Updated dictionary with the processed 'atoms' DataFrame.
    """
    if 'atoms' not in dataframes_dict:
        print("❌ 'atoms' DataFrame not found in processed_dfs_dict.")
        return dataframes_dict

    df = dataframes_dict['atoms'].copy()

    if 'atom_name' not in df.columns:
        print("❌ 'atom_name' column not found in 'atoms' DataFrame.")
        return dataframes_dict

    # List of atomic symbols
    atomic_symbols = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ]

    def split_atom_name(atom_list):
        """Extracts the atomic symbol and additional numbering."""
        if not isinstance(atom_list, list) or len(atom_list) != 1:
            raise ValueError(f"Invalid atom_name format: {atom_list}")

        atom_name = atom_list[0]  # Extract string from list
        for symbol in atomic_symbols:
            if atom_name.startswith(symbol):
                remaining = atom_name[len(symbol):]
                return symbol, remaining if remaining else None  # Return extracted parts

        raise ValueError(f"Atomic symbol not found in atom_name: {atom_name}")

    # Generate hex-style numbering
    def hex_style_atom_type(index):
        """Generate a three-character hex-style atom type based on index."""
        hex_digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        base = len(hex_digits)
        result = ''
        index -= 1
        while index > 0:
            index, remainder = divmod(index, base)
            result = hex_digits[remainder] + result
        return result.zfill(3)  # Ensure three-character output

    # Extract atomic symbol and numbering
    temp_split = df['atom_name'].apply(lambda x: pd.Series(split_atom_name(x)))
    df['atom_type_temp'] = temp_split[0]  # Atomic symbol
    df['old_number_temp'] = temp_split[1]  # Original numbering

    # Generate unique hex-style numbering
    df['new_number_temp'] = [hex_style_atom_type(i) for i in range(1, len(df) + 1)]

    # Update atom_name column while keeping it as a list
    df['atom_name'] = df.apply(lambda row: [row['atom_type_temp'] + row['new_number_temp']], axis=1)

    # Remove temporary columns
    df = df.drop(columns=['atom_type_temp', 'old_number_temp', 'new_number_temp'])

    # Update the dictionary
    dataframes_dict['atoms'] = df

    print("✅ Successfully renumbered atom_name in the 'atoms' DataFrame while keeping all other columns intact.")
    return dataframes_dict
def comment_nmer(dataframes_dict, mapping):
    """
    Updates the 'atoms' DataFrame with an 'nmer' column.
    - Extracts the integer from the 'atoms' column (list with one value).
    - Finds the corresponding value in 'atoms_final' in the 'mapping' DataFrame.
    - Retrieves and assigns the associated 'nmer' value from 'mapping' to 'atoms'.
    - Prefixes each 'nmer' value with '; '.

    Args:
        dataframes_dict (dict): Dictionary containing DataFrames, including 'atoms'.
        mapping (pd.DataFrame): DataFrame containing 'atoms_final' and 'nmer' columns.

    Returns:
        dict: Updated dictionary with the processed 'atoms' DataFrame.
    """
    if 'atoms' not in dataframes_dict:
        print("❌ 'atoms' DataFrame not found in processed_dfs_dict.")
        return dataframes_dict

    df_atoms = dataframes_dict['atoms'].copy()

    if 'atoms' not in df_atoms.columns:
        print("❌ 'atoms' column not found in 'atoms' DataFrame.")
        return dataframes_dict

    if 'atoms_final' not in mapping.columns or 'nmer' not in mapping.columns:
        print("❌ 'mapping' DataFrame is missing required columns: 'atoms_final' or 'nmer'.")
        return dataframes_dict

    # Extract integer from the atoms column (which is a list with one integer value)
    def extract_integer(atom_list):
        """Extracts the integer from the list."""
        if isinstance(atom_list, list) and len(atom_list) == 1 and isinstance(atom_list[0], int):
            return atom_list[0]
        return None

    # Apply extraction function to get the number from atoms column
    df_atoms['atom_number'] = df_atoms['atoms'].apply(extract_integer)

    # Merge with mapping DataFrame based on atom_number matching atoms_final
    df_atoms = df_atoms.merge(mapping[['atoms_final', 'nmer']], left_on='atom_number', right_on='atoms_final', how='left')

    # Prefix each nmer value with "; "
    df_atoms['nmer'] = df_atoms['nmer'].apply(lambda x: f"; nmer = {x}" if pd.notna(x) else "; nmer = N/A")

    # Drop temporary columns
    df_atoms = df_atoms.drop(columns=['atom_number', 'atoms_final'])

    # Update the dictionary
    dataframes_dict['atoms'] = df_atoms

    print("✅ Successfully updated 'atoms' with 'nmer' column based on 'mapping'.")
    return dataframes_dict
def display_all_dataframes(dataframes_dict):
    """
    Prints the names and full contents of all DataFrames in the given dictionary.

    Args:
        dataframes_dict (dict): Dictionary containing DataFrames.
    """
    if not dataframes_dict:
        print("❌ No DataFrames to display.")
        return

    print("📂 Displaying all DataFrames:\n")

    for name, df in dataframes_dict.items():
        print(f"🔹 DataFrame: {name} (Shape: {df.shape})\n")
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print(df)
            print("\n" + "="*80 + "\n")  # Separator for readability
def organize_dfs(dataframes_dict, parser_class):
    """
    Sorts DataFrame columns based on the expected structure from Itp_parser.
    - Moves any unexpected columns to the end.
    - Reports which DataFrames contain unexpected columns.

    Args:
        dataframes_dict (dict): Dictionary containing processed DataFrames.
        parser_class (Itp_parser): The Itp_parser class to extract expected columns (not an instance).

    Returns:
        dict: Updated dictionary with sorted DataFrames.
        dict: Report of extra columns found in each DataFrame.
    """
    sorted_dfs_dict = {}
    unexpected_columns_report = {}

    # Extract expected DataFrame structures from the parser class
    expected_df_structures = parser_class.TOTAL_DIRECTIVES

    for df_name, df in dataframes_dict.items():
        expected_columns = None

        # Check if the DataFrame name matches a known section
        for directive, section_dict in expected_df_structures.items():
            if df_name in section_dict:
                expected_columns = ["atoms", "function_type", "params", "comments"]  # Default structure from parser
                break

        if expected_columns is None:
            print(f"⚠️ Warning: {df_name} is not a recognized section in Itp_parser. Skipping column sorting.")
            sorted_dfs_dict[df_name] = df
            continue

        # Identify unexpected columns
        df_columns = list(df.columns)
        extra_columns = [col for col in df_columns if col not in expected_columns]

        if extra_columns:
            unexpected_columns_report[df_name] = extra_columns
            print(f"📌 Unexpected columns in {df_name}: {extra_columns}")

        # Create the new sorted column order (expected first, extras last)
        sorted_columns = [col for col in expected_columns if col in df_columns] + extra_columns

        # Reorder DataFrame columns
        sorted_dfs_dict[df_name] = df[sorted_columns]

    print("\n✅ Successfully sorted DataFrames based on Itp_parser.")
    return sorted_dfs_dict, unexpected_columns_report
def sort_and_identify_extra_columns(dataframes_dict, parser_class):
    """
    Sorts DataFrame columns based on the expected structure from Itp_parser.
    - Moves any unexpected columns to the end.
    - Reports which DataFrames contain unexpected columns.

    Args:
        dataframes_dict (dict): Dictionary containing processed DataFrames.
        parser_class (Itp_parser): The Itp_parser class to extract expected columns (not an instance).

    Returns:
        dict: Updated dictionary with sorted DataFrames.
        dict: Report of extra columns found in each DataFrame.
    """
    sorted_dfs_dict = {}
    unexpected_columns_report = {}

    # Extract expected DataFrame structures from the parser class
    expected_df_structures = parser_class.TOTAL_DIRECTIVES

    # Manually add atoms expected columns (since it isn't inside TOTAL_DIRECTIVES)
    expected_atoms_columns = [
        "atoms", "atom_types", "resodue#", "residue_name", "atom_name",
        "chargeGroups#", "charge", "mass", "comments", "Comments_top", "Coments_bottom"
    ]

    for df_name, df in dataframes_dict.items():
        expected_columns = None

        # Check if the DataFrame name matches a known section in TOTAL_DIRECTIVES
        for directive, section_dict in expected_df_structures.items():
            if df_name in section_dict:
                expected_columns = ["atoms", "function_type", "params", "comments"]  # Default structure for directives
                break

        # Special handling for the 'atoms' DataFrame
        if df_name == "atoms":
            expected_columns = expected_atoms_columns

        if expected_columns is None:
            print(f"⚠️ Warning: {df_name} is not a recognized section in Itp_parser. Skipping column sorting.")
            sorted_dfs_dict[df_name] = df
            continue

        # Identify unexpected columns
        df_columns = list(df.columns)
        extra_columns = [col for col in df_columns if col not in expected_columns]

        if extra_columns:
            unexpected_columns_report[df_name] = extra_columns
            print(f"📌 Unexpected columns in {df_name}: {extra_columns}")

        # Create the new sorted column order (expected first, extras last)
        sorted_columns = [col for col in expected_columns if col in df_columns] + extra_columns

        # Reorder DataFrame columns
        sorted_dfs_dict[df_name] = df[sorted_columns]

    print("\n✅ Successfully sorted DataFrames based on Itp_parser.")
    return sorted_dfs_dict, unexpected_columns_report
def remove_none_lines(input_itp: str, output_itp: str = None) -> None:
    """
    Removes any lines from an ITP file that contain the string 'None'.

    Args:
        input_itp (str): Path to the input ITP file.
        output_itp (str, optional): Path to save the cleaned ITP file.
                                    If None, it will overwrite the input file.
    """
    output_itp = output_itp or input_itp

    with open(input_itp, "r") as f:
        lines = f.readlines()

    cleaned_lines = []
    for line in lines:
        if "None" in line and not line.strip().startswith(";"):
            continue  # Skip lines with 'None', except comments
        cleaned_lines.append(line)

    with open(output_itp, "w") as f:
        f.writelines(cleaned_lines)

    print(f"✅ Cleaned ITP saved to {output_itp} (removed lines with 'None')")
def extract_and_paste_sections(original_itp: str, output_itp: str) -> None:
    """
    Extracts [ atomtypes ] and [ moleculetype ] sections from an ITP file and inserts
    them after the POLYX header in another ITP output file.

    Args:
        original_itp (str): Path to the source ITP file.
        output_itp (str): Path to the output ITP file (must already exist with POLYX header).
    """
    # Read original ITP to find and store sections
    with open(original_itp, "r") as f:
        original_lines = f.readlines()

    target_sections = {"[ atomtypes ]": [], "[ moleculetype ]": []}
    found = {"[ atomtypes ]": False, "[ moleculetype ]": False}
    current_section = None

    for line in original_lines:
        stripped = line.strip()

        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped.lower()
            continue

        if current_section in target_sections:
            if not stripped or stripped.startswith("["):
                current_section = None
            else:
                target_sections[current_section].append(line.rstrip())
                found[current_section] = True

    # Read output ITP (assumes it already exists with POLYX header)
    with open(output_itp, "r") as f:
        output_lines = f.readlines()

    # Insert the extracted sections after the POLYX header (assumes header ends with two blank lines)
    insert_index = 0
    blank_line_count = 0
    for i, line in enumerate(output_lines):
        if line.strip() == "":
            blank_line_count += 1
        else:
            blank_line_count = 0
        if blank_line_count == 2:
            insert_index = i + 1
            break

    # Build insertion block
    insertion_block = []
    for section_name in ["[ atomtypes ]", "[ moleculetype ]"]:
        if found[section_name]:
            insertion_block.append(f"{section_name}\n")
            insertion_block += [line + "\n" for line in target_sections[section_name]]
            insertion_block.append("\n")

    # Inject and rewrite
    new_lines = output_lines[:insert_index] + insertion_block + output_lines[insert_index:]
    with open(output_itp, "w") as f:
        f.writelines(new_lines)

    print(f"✅ Inserted found sections into {output_itp}")


# Run the function with your data
atoms_dict, G3, G3_n1, G3_n2, G3_n3, G3_restored = split_G3(itp_3mer, dh1, dh2)

# Print the dictionary with assigned `nmer` values
print("\nUpdated Atoms Dictionary:")
print(atoms_dict)


# Plot all graphs
G3 = G3_restored

mol_graph = MolecularGraph(itp_1mer)
mol_graph.assign_atoms()  # Ensure atom labels are assigned ✅
G1 = mol_graph.get_graph()  # Get NetworkX graph ✅

plot_graphs(
    [G3, G3_n1, G3_n2, G3_n3, G1],
    ["G3", "G3_n1 (nmer=1)", "G3_n2 (nmer=2)", "G3_n3 (nmer=3)", "G1"]
)

# Example usage
G1_processed = process_G1(G1, hs_on_monomer)

plot_graphs(
    [G1, G1_processed],
    ["G1", "G1_processed"]
)

# Run the function with dh1 and dh2
results = process_G3(G1_processed, G3_n1, G3_n2, G3_n3, dh1, dh2)

# Print results
for name, result in results.items():
    print(f"\n🔹 {name} Isomorphism Result:")
    print(f"Is Isomorphic? {result['is_isomorphic']}")
    print(f"Mapping: {result['mapping']}")

# Run finalize_G3_mappings to modify the mappings
finalized_results = finalize_G3_mappings(results)

# Print modified mappings
for name, result in finalized_results.items():
    print(f"\n🔹 {name} Finalized Mapping:")
    print(result["mapping"])

# Step 3: Apply Finalization Function
G3_final,G3_remapped  = finalize_G3(G3, finalized_results)

plot_graphs(
    [G3_final],
    ["G3"]
)

print(G3_remapped)

# Run the function with your ITP parser
parser = Itp_parser(itp_3mer)  # Load the ITP file
dataframes_dict = process_itp_to_df(parser)

# Now you can access individual DataFrames like this:
print(dataframes_dict["moleculetype"])  # Example of accessing a specific DataFrame
print(dataframes_dict["atoms"])         # Example of another DataFrame

# Apply the remapping function
dataframes_dict = remap_itp_data(dataframes_dict, G3_remapped, parser)

# Check the updated atoms DataFrame
print(dataframes_dict["atoms"])

# Apply the function to add 'df_set' to each DataFrame
dataframes_dict = assign_df_set(dataframes_dict)

# Example: Print the updated atoms DataFrame
print(dataframes_dict["atoms"])

# List DataFrames that have the 'df_set' column
df_with_df_set = [key for key, df in dataframes_dict.items() if 'df_set' in df.columns]

# Print the names of DataFrames that have 'df_set'
print(df_with_df_set)

# Initialize a counter to track occurrences of 'df_set' values
df_set_counter = Counter()

# Iterate through DataFrames that have 'df_set' and count occurrences
for key, df in dataframes_dict.items():
    if 'df_set' in df.columns:
        df_set_counter.update(df['df_set'])

# Convert the counter to a DataFrame for better visualization
df_set_counts = pd.DataFrame(df_set_counter.items(), columns=['df_set', 'count'])

# Display the result
print(df_set_counts)

# Dictionary to store the new split DataFrames
split_dataframes = {}

# Iterate through each DataFrame in the original dictionary
for df_name, df in dataframes_dict.items():
    if 'df_set' in df.columns:
        # Get unique values in 'df_set'
        unique_df_sets = df['df_set'].unique()

        # Split DataFrame by unique 'df_set' values
        for df_set_value in unique_df_sets:
            new_df_name = f"{df_name}_{df_set_value}"  # Create new name
            split_dataframes[new_df_name] = df[df['df_set'] == df_set_value].copy()  # Store split DataFrame

# Display the names of the new DataFrames
print("Newly created DataFrames:", list(split_dataframes.keys()))


# Example usage
polymer_sections = generate_polymer_list(n)
print(polymer_sections)

# Call the function
test_results = analyze_polymer_sections(split_dataframes, polymer_sections, n)

# Example usage
clean_dfs_AB1(split_dataframes)

# Example usage
clean_dfs_Bn2C(split_dataframes, n)

# Example usage
clean_dfs_Ba(split_dataframes, n)

# Example usage
clean_dfs_Ba1Ba2(split_dataframes, n)

# Example Usage
mapping = build_polymer_mapping(split_dataframes, n)

# Run the function on split_dataframes
clear_irrelevant_dfs(split_dataframes)

# Apply the function to split_dataframes
split_dataframes = cut_final_dfs(split_dataframes)

# Print remaining DataFrames
print("📂 Remaining DataFrames:", list(split_dataframes.keys()))

# Usage
find_invalid_B_entries(split_dataframes)

# Apply the function to split_dataframes
final_dfs_dict = concatenate_dfs(split_dataframes)

# Print result summary
print("✅ Successfully concatenated DataFrames into final_dfs_dict:")
for key in final_dfs_dict:
    print(f" - {key}: {final_dfs_dict[key].shape}")

# Run the function to process final_dfs_dict
processed_dfs_dict = process_dataframes(final_dfs_dict, mapping)

# Print summary
print("✅ Successfully processed all DataFrames in processed_dfs_dict.")

# Apply the function
processed_dfs_dict = reorder_processed_dfs(processed_dfs_dict)

# Apply the function
processed_dfs_dict = renumber_atom_names(processed_dfs_dict)

# Print the updated 'atoms' DataFrame
print(processed_dfs_dict['atoms'])

# Apply the function
processed_dfs_dict = comment_nmer(processed_dfs_dict, mapping)

# Print the updated 'atoms' DataFrame
print(processed_dfs_dict['atoms'])

# Run the function on final_dfs_dict
display_all_dataframes(processed_dfs_dict)

# Apply the function using the Itp_parser class (without needing an instance)
processed_dfs_dict, extra_columns_report = organize_dfs(processed_dfs_dict, Itp_parser)

# Print unexpected columns report
if extra_columns_report:
    print("\n🔍 Extra Columns Found in DataFrames:")
    for df_name, extra_cols in extra_columns_report.items():
        print(f"- {df_name}: {extra_cols}")
else:
    print("\n✅ No unexpected columns found.")

# Apply the function using the Itp_parser class (without needing an instance)
processed_dfs_dict, extra_columns_report = sort_and_identify_extra_columns(processed_dfs_dict, Itp_parser)

# Print unexpected columns report
if extra_columns_report:
    print("\n🔍 Extra Columns Found in DataFrames:")
    for df_name, extra_cols in extra_columns_report.items():
        print(f"- {df_name}: {extra_cols}")
else:
    print("\n✅ No unexpected columns found.")

parser.DF = processed_dfs_dict
parser.save_itp(output_itp)

# Example usage
remove_none_lines(output_itp)  # Overwrites original

extract_and_paste_sections(itp_3mer, output_itp)
