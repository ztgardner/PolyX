import matplotlib.pyplot as plt
import networkx as nx
from collections import Counter
import re
from collections import defaultdict
import pandas as pd
from itp_parser import Itp_parser
import copy
import numpy as np
from networkx.algorithms import isomorphism
from networkx.algorithms import isomorphism as iso
from functools import partial

class AnchoredGraphMatcher(iso.GraphMatcher):
    def __init__(self, G1, G2, node_match=None, edge_match=None, anchor_mapping=None):
        super().__init__(G1, G2, node_match=node_match, edge_match=edge_match)
        self.anchor_mapping = anchor_mapping or {}

    def semantic_feasibility(self, G1_node, G2_node):
        if not super().semantic_feasibility(G1_node, G2_node):
            return False

        # Strict anchor constraints
        g1_expected = self.anchor_mapping.get(G2_node)
       #print(f"Anchor constraint failed: G2[{G2_node}] must match G1[{g1_expected}], got G1[{G1_node}]")
        if g1_expected is not None and G1_node != g1_expected:
            return False

        if G1_node in self.anchor_mapping.values():
            expected_g2 = [k for k, v in self.anchor_mapping.items() if v == G1_node]
            if G2_node not in expected_g2:
                #rint(f"Anchor constraint failed: G1[{G1_node}] must match G2{expected_g2}, got G2[{G2_node}]")
                return False

        return True


class Swapper:

    def __init__(self, achor_mapping:list, sc_1: list, sc_rep: list, hs_on_monomer: list,hs_ontrimer:list, itp_file_trimer: str, gro_file_trimer: str, itp_file_1_mer: str, gro_file_1_mer: str,):
        self.ITP_TRI= Itp_parser(itp_file_trimer)
        self.ITP_SOLO= Itp_parser(itp_file_1_mer)
        self.TRI_GRO_LOAD=True
        self.SOLO_GRO_LOAD=True
        try:
            self.ITP_TRI.load_gro(gro_file_trimer)
        except:
            self.TRI_GRO_LOAD=False
            raise RuntimeError("The tri-mer gro did not load")
        try:
            self.ITP_SOLO.load_gro(gro_file_1_mer)
        except:
            self.SOLO_GRO_LOAD=False
            raise RuntimeError("The one-mer gro did not load")
        self.ITP_Graph_TRI = Swapper.Build_ITP_Nodes(self.ITP_TRI)
        self.ITP_Graph_SOLO = Swapper.Build_ITP_Nodes(self.ITP_SOLO)
        self.ITP_Graph_TRI = Swapper.Build_ITP_Vertices(self.ITP_Graph_TRI, self.ITP_TRI)
        self.ITP_Graph_SOLO = Swapper.Build_ITP_Vertices(self.ITP_Graph_SOLO, self.ITP_SOLO)
        self.sidechain_atoms_solo = Swapper.get_sidechain_atoms(self.ITP_Graph_SOLO, sc_1)

        self.side_chain_atoms_with_backbone, self.SIDE_CHAIN_ITP_SOLO = Swapper.get_side_ITP(self.ITP_SOLO, self.sidechain_atoms_solo)
        self.sidechain_atoms_tri = Swapper.get_sidechain_atoms(self.ITP_Graph_TRI, sc_rep)

        self.sidechain_atoms_tri.difference_update(
            set(atom for sublist in sc_rep for atom in sublist))  # removes all bridges from sidechain



        self.BACKBONE_ITP_TRI = Swapper.get_backbone_ITP(self.ITP_TRI, self.ITP_Graph_TRI, self.sidechain_atoms_tri)

        self.BACKBONE_GRAPH_TRI = self.ITP_Graph_TRI.copy()
        self.BACKBONE_GRAPH_TRI.remove_nodes_from(self.sidechain_atoms_tri)

        self.BACKBONE_GRAPH_TRI, self.renumber_map_BACKBONE_TRI = self.renumber_graph(self.BACKBONE_GRAPH_TRI, hs_ontrimer[0], hs_ontrimer[
            1])
        self.BACKBONE_ITP_TRI = self.renumber_itp(self.BACKBONE_ITP_TRI, self.renumber_map_BACKBONE_TRI)



        self.Side_chain_graph = self.ITP_Graph_SOLO.subgraph(self.side_chain_atoms_with_backbone).copy()

        self.Side_atoms_only = self.side_chain_atoms_with_backbone - self.sidechain_atoms_solo
        self.Side_atoms_only.add(sc_1[0][1])
        self.Side_atoms_only.add(sc_1[0][0])
        self.Side_atoms_only_Graph = self.Side_chain_graph.copy()
        self.Side_atoms_only_Graph.remove_nodes_from(self.Side_atoms_only)

        self.Backbone_Only = self.side_chain_atoms_with_backbone - self.Side_atoms_only
        self.Backbone_atoms_only_Graph = self.Side_chain_graph.copy()
        self.Backbone_atoms_only_Graph.remove_nodes_from(self.Backbone_Only)

    @staticmethod
    def Build_ITP_Nodes(ITP):
        ITP_Graph = nx.Graph()

        try:
            for atom in ITP.DF["atoms"]["atoms"]:
                Atom_Index = int(atom[0])
                # print(Atom_Index)

                attributes = {}
                for section in ITP.DF["atoms"].keys():
                    if section in ["Comments_top", "Coments_bottom",
                                   "atoms"]:  # Want to skip top and bottom comments and atoms section (it is the node index)
                        continue
                        # print(section,ITP.DF["atoms"][section][Atom_Index -1 ] )
                    attributes[section] = ITP.DF["atoms"][section][Atom_Index - 1]

                    try:
                        attributes['cord'] = ITP.coordinates[Atom_Index - 1]
                    except:
                        attributes['cord'] = np.nan
                    # print(ITP.DF["atoms"][section][Atom_Index - 1])

                ITP_Graph.add_node(Atom_Index, **attributes)


        except:
            raise RuntimeError("Failed to Turn Atom Numbers to Nodes")

        return ITP_Graph

    @staticmethod
    def Build_ITP_Vertices(ITP_Graph, ITP_OBJ):
        for bond in ITP_OBJ.DF["bond"]["atoms"]:
            i, j = int(bond[0]), int(bond[1])
            ITP_Graph.add_edge(i, j)
        if nx.is_connected(ITP_Graph):
            return ITP_Graph

        for angle in ITP_OBJ.DF["angle"]["atoms"]:
            i, j, k = int(angle[0]), int(angle[1]), int(angle[2])
            ITP_Graph.add_edge(i, j)
            ITP_Graph.add_edge(j, k)

        if nx.is_connected(ITP_Graph):
            return ITP_Graph

        for dih in ITP_OBJ.DF["Ryckaert-Bellemans dihedral"]["atoms"]:
            i, j, k, l = int(dih[0]), int(dih[1]), int(dih[2]), int(dih[3])
            ITP_Graph.add_edge(i, j)
            ITP_Graph.add_edge(j, k)
            ITP_Graph.add_edge(k, l)

        if nx.is_connected(ITP_Graph):
            return ITP_Graph

        raise RuntimeError(
            "ITP is not a connected Graph, Only Bonds/Angles/Ryckaert-Bellemans dihedral Are Check."
            "Please reach out for more Implementation")



    @staticmethod
    def renumber_itp( itp, renumber_map):
        """
        Returns a new ITP with *all* atom indices remapped according to renumber_map.
        Assumes renumber_map: { old_int: new_int }.
        """
        NEW = copy.deepcopy(itp)

        for section in NEW.sections:
            if section in ('defaults', 'atomtypes', 'moleculetype', 'dihedral_restraints'):
                continue
            # print(itp[section]["atoms"])
            for i, atoms in enumerate((NEW[section]["atoms"])):

                Entry = []
                for a in atoms:
                    Entry.append(str(renumber_map[int(a)]))
                NEW[section]["atoms"][i] = Entry
                # print(NEW[section]["atoms"])

        return NEW

    @staticmethod
    def renumber_graph(G, start_node, end_node):
        # 1) Compute the main path
        try:
            main_path = nx.shortest_path(G, source=start_node, target=end_node)
        except nx.NetworkXNoPath:
            raise ValueError(f"No path between {start_node} and {end_node}")

        main_set = set(main_path)
        mapping = {}
        next_label = 1

        # 2) Recursively map a side‐branch subtree in DFS order
        def traverse_subtree(u):
            nonlocal next_label
            mapping[u] = next_label
            next_label += 1
            # sort for deterministic order; drop any neighbors already mapped
            for v in sorted(G.neighbors(u)):
                if v not in mapping and v not in main_set:
                    traverse_subtree(v)

        # 3) Walk the main path, mapping each node, then its side branches
        def traverse_main(i):
            nonlocal next_label
            u = main_path[i]
            mapping[u] = next_label
            next_label += 1

            # for each neighbor off the main path, fully traverse its branch
            for v in sorted(G.neighbors(u)):
                if v not in mapping and v not in main_set:
                    traverse_subtree(v)

            # then move to the next node on the main path
            if i + 1 < len(main_path):
                traverse_main(i + 1)

        traverse_main(0)

        # 4) Any leftover (disconnected) nodes
        for n in sorted(G.nodes()):
            if n not in mapping:
                mapping[n] = next_label
                next_label += 1

        # 5) Relabel and return
        G_renumbered = nx.relabel_nodes(G, mapping, copy=True)
        return G_renumbered, mapping

    @staticmethod
    def get_sidechain_atoms(ITP_Graph, Bridge):
        """
        Given a graph of the molecule and a list of (backbone_atom, sidechain_atom) pairs,
        return a set of all atom indices belonging to the sidechain(s).
        """
        sidechain_atoms = set()
        for backbone_atom, sidechain_atom in Bridge:
            visited = set()
            queue = [sidechain_atom]

            # sidechain_atoms.add(Bridge[0][0])
            # sidechain_atoms.add(Bridge[0][1])

            while queue:
                current = queue.pop()
                if current in visited or current == backbone_atom:
                    continue
                visited.add(current)
                queue.extend(ITP_Graph.neighbors(current))

            sidechain_atoms.update(visited)

        # sidechain_atoms.discard(Bridge[0]) #Drops the starting bridge, it will get readded in get_SIDE_ITP but dont want to include it as a "relavent" sidechain

        return sidechain_atoms

    @staticmethod
    def get_side_ITP(ITP_OBJ, sidechain_atoms):
        """
        Given an ITP object, graph, and a list of bside_chain_Atoms, return a new ITP object
        containing only the atoms and interactions related to the sidechain(s).
        I apologize in advance for the number of loops, it makes me sad to see
        """

        ### SET UP NEW ITP
        sections = ITP_OBJ.sections
        sections = [entry for entry in sections if
                    entry not in ('defaults', 'atomtypes', "moleculetype")]  # Filter sections

        # Initialize NEW_ITP and clear non-atom data sections
        NEW_ITP = copy.deepcopy(ITP_OBJ)  # Use deepcopy to copy the structure without references to the original data
        for section in ITP_OBJ.sections:
            if section in ('defaults', 'atomtypes', "moleculetype"):  # This doesnt change
                continue
            NEW_ITP[section].clear()  # Clear the atom data

        ### Find Every Atom Relavent to SideChains
        relavent_atoms = set()
        for section in sections:
            n_entries = len(ITP_OBJ[section]["atoms"]) if "atoms" in ITP_OBJ[section] else len(
                next(iter(ITP_OBJ[section].values())))

            for i in range(n_entries):
                atoms = ITP_OBJ[section]["atoms"][i]
                if any(int(atom) in sidechain_atoms for atom in atoms):
                    relavent_atoms |= set(int(a) for a in atoms)

        for section in sections:
            n_entries = len(ITP_OBJ[section]["atoms"]) if "atoms" in ITP_OBJ[section] else len(
                next(iter(ITP_OBJ[section].values())))

            if section == "atoms":  # I need to look at "relavent atoms in the atom section"
                for i in range(n_entries):
                    atoms = ITP_OBJ[section]["atoms"][i]
                    if any(int(atom) in relavent_atoms for atom in atoms):
                        # Build new entry
                        entry = {key: ITP_OBJ[section][key][i]
                                 for key in ITP_OBJ[section]
                                 if isinstance(ITP_OBJ[section][key], list) and len(ITP_OBJ[section][key]) > i}

                        # Append to NEW_ITP
                        for key, value in entry.items():
                            if key not in NEW_ITP[section]:
                                NEW_ITP[section][key] = []
                            NEW_ITP[section][key].append(value)

            for i in range(n_entries):
                atoms = ITP_OBJ[section]["atoms"][i]
                if any(int(atom) in sidechain_atoms for atom in atoms):
                    # Build new entry
                    entry = {key: ITP_OBJ[section][key][i]
                             for key in ITP_OBJ[section]
                             if isinstance(ITP_OBJ[section][key], list) and len(ITP_OBJ[section][key]) > i}

                    # Append to NEW_ITP
                    for key, value in entry.items():
                        if key not in NEW_ITP[section]:
                            NEW_ITP[section][key] = []
                        NEW_ITP[section][key].append(value)

        return relavent_atoms, NEW_ITP


    @staticmethod
    def get_backbone_ITP(ITP_OBJ, ITP_Graph, sidechain_atoms):
        """
        Returns a copy of the ITP with all atoms and interactions related to sidechains removed.

        Parameters:
        - ITP_OBJ: original ITP object
        - ITP_Graph: graph object representing atom connectivity
        - sc_rep: list of (backbone_atom, sidechain_atom) tuples
        """

        # Get sidechain atom indices (excluding backbone atoms)

        # Create deep copy of ITP_OBJ
        NEW_ITP = copy.deepcopy(ITP_OBJ)
        SKIP = {'defaults', 'atomtypes', 'moleculetype'}

        for section in ITP_OBJ.sections:
            if section in SKIP:
                continue

            if "atoms" not in ITP_OBJ[section]:
                NEW_ITP[section] = copy.deepcopy(ITP_OBJ[section])
                continue

            NEW_ITP[section].clear()
            n_entries = len(ITP_OBJ[section]["atoms"])

            for i in range(n_entries):
                atom_list = ITP_OBJ[section]["atoms"][i]
                atom_indices = {int(a) for a in atom_list}

                # Keep entry only if it contains no sidechain atoms
                if not (atom_indices & sidechain_atoms):
                    for key, valuelist in ITP_OBJ[section].items():
                        if isinstance(valuelist, list) and len(valuelist) > i:
                            NEW_ITP[section].setdefault(key, []).append(valuelist[i])

        return NEW_ITP

    @staticmethod
    def node_match2(n1, n2, mass_tolerance=0.5):
        def unwrap(v):
            if isinstance(v, list) and len(v) == 1:
                return v[0]
            return v

        a1 = unwrap(n1.get("atom_name"))
        a2 = unwrap(n2.get("atom_name"))

        m1_str = unwrap(n1.get("mass"))
        m2_str = unwrap(n2.get("mass"))

        if a1 is None or a2 is None or m1_str is None or m2_str is None:
            return False

        # Check first letter or exact atom_name match
        atom_name_match = (a1 == a2) or (a1[0] == a2[0])

        # Compare mass as floats with tolerance
        try:
            m1 = float(m1_str)
            m2 = float(m2_str)
        except Exception:
            return False

        mass_match = abs(m1 - m2) <= mass_tolerance

        return atom_name_match and mass_match
