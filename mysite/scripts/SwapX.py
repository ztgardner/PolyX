import matplotlib.pyplot as plt
import networkx as nx
from collections import Counter
import re
from collections import defaultdict
import pandas as pd
from scripts.itp_parser import Itp_parser
import copy
import numpy as np

import networkx as nx
from networkx.algorithms import isomorphism as iso

from scipy.spatial.transform import Rotation as R

class RaiseError(Exception):
    pass




class POLY_Swap:
    def __init__(
        self,
        itp_file_path_tri: str,
        gro_file_path_tri: str,
        itp_file_path_1mer: str,
        gro_file_path_1mer: str,
        sc_1: list,
        sc_rep: list,
        hs_on_monomer: list,
        hs_ontrimer: list,
        Out_Name:str
    ):
        self.itp_file_path_tri = itp_file_path_tri
        self.gro_file_path_tri = gro_file_path_tri

        self.itp_file_path_1mer = itp_file_path_1mer
        self.gro_file_path_1mer = gro_file_path_1mer

        self.sc_1 = sc_1
        self.sc_rep = sc_rep
        self.hs_on_monomer = hs_on_monomer
        self.hs_ontrimer = hs_ontrimer

        # Load ITP and GRO files
        try:
            ITP_TRI = Itp_parser(itp_file_path_tri)
            ITP_SOLO = Itp_parser(itp_file_path_1mer)
        except Exception as e:
            print(e)
            raise RaiseError("ITP files did not load into ITP parser properly") from e

        # Attempt to load trimer GRO
        try:
            ITP_TRI.load_gro(gro_file_path_tri)
            self.tri_gro_loaded = True
        except Exception as e:
            self.tri_gro_loaded = False
            print(e)
            raise RaiseError(f"TRIMER: '{gro_file_path_tri}' GRO did not load properly") from e

        # Attempt to load monomer GRO
        try:
            ITP_SOLO.load_gro(gro_file_path_1mer)
            self.solo_gro_loaded = True
        except Exception as e:
            self.solo_gro_loaded = False
            print(e)
            raise RaiseError(f"MONOMER: '{gro_file_path_1mer}' GRO did not load properly") from e

        # Build Graphs
        try:
            ITP_Graph_TRI = self.Build_ITP_Nodes(ITP_TRI)
            ITP_Graph_SOLO = self.Build_ITP_Nodes(ITP_SOLO)
        except Exception as e:
            print(e)
            raise RaiseError("Failed to build ITP Graphs from ITP Parsers") from e
        
        #Build Edges
        try:         
            ITP_Graph_TRI = self.Build_ITP_Vertices(ITP_Graph_TRI,ITP_TRI)
            ITP_Graph_SOLO = self.Build_ITP_Vertices(ITP_Graph_SOLO,ITP_SOLO)
        except Exception as e:
            print(e)
            raise RaiseError("Failed to build ITP Edges from ITP Parsers") from e

        #Find Sidechian Atoms
        try:
            
            sidechain_atoms_solo = self.get_sidechain_atoms(ITP_Graph_SOLO, sc_1)
            print("Sidechain atoms in monomer:", sidechain_atoms_solo, "   Does this look correct? Compared to Input: ", sc_1)
        except Exception as e:
                    print(e)
                    raise RaiseError("Failed to get sidechain atoms from MONOMER") from e
        
        #Find Sidechain ITP
        try:
            side_chain_atoms_with_backbone, SIDE_CHAIN_ITP_SOLO = self.get_side_ITP(ITP_SOLO, sidechain_atoms_solo)
        except Exception as e:
            print(e)
            raise RaiseError("Failed to get sidechain ITP from MONOMER") from e
        

        #Find Sidechain Atoms on Trimer remove them get new ITP with just the Backbone
        try:
            sidechain_atoms_tri = self.get_sidechain_atoms(ITP_Graph_TRI, sc_rep)
            
            sidechain_atoms_tri.difference_update(set(atom for sublist in sc_rep for atom in sublist)) #removes all bridges from sidechain
            
            print("Sidechain atoms in trimer (excluding bridging input):", sidechain_atoms_tri, "   Does this look correct? Compared to Input: ", sc_rep)

            BACKBONE_ITP_TRI = self.get_backbone_ITP(ITP_TRI, ITP_Graph_TRI, sidechain_atoms_tri)
        
        except Exception as e:
            print(e)
            raise RaiseError("Failed to get backbone ITP from TRIMER") from e

        #drop sidechain nodes on the trimer, renumber trimer to go from left to right,
        try:
            BACKBONE_GRAPH_TRI = ITP_Graph_TRI.copy()
            BACKBONE_GRAPH_TRI.remove_nodes_from(sidechain_atoms_tri)

            BACKBONE_GRAPH_TRI,renumber_map_BACKBONE_TRI = self.renumber_graph(BACKBONE_GRAPH_TRI, hs_ontrimer[0], hs_ontrimer[1]) #Renumber TRIMER To go from left to right
            BACKBONE_ITP_TRI = self.renumber_itp(BACKBONE_ITP_TRI,renumber_map_BACKBONE_TRI)

            Side_chain_graph =  ITP_Graph_SOLO.subgraph(side_chain_atoms_with_backbone).copy() #graph of just sidechain


            Side_atoms_only = side_chain_atoms_with_backbone - sidechain_atoms_solo
      
            Side_atoms_only.add(sc_1[0][1])
            Side_atoms_only.add(sc_1[0][0])
          
            Side_atoms_only_Graph = Side_chain_graph.copy()
            Side_atoms_only_Graph.remove_nodes_from(Side_atoms_only)
        

            Backbone_Only = side_chain_atoms_with_backbone - Side_atoms_only
            Backbone_atoms_only_Graph = Side_chain_graph.copy()
            Backbone_atoms_only_Graph.remove_nodes_from(Backbone_Only)
         

        except Exception as e:
            print(e)
            raise RaiseError("Failed to renumber TRIMER backbone and isolate MONOMER sidechain graphs") from e
   
        try:
            ### I need to renumber the atoms in ONLY the side chain such that I can just add it to whatever the number of backbone atoms there are plus the number of times we have seen it.
            sc_rep = [(renumber_map_BACKBONE_TRI[backbone], renumber_map_BACKBONE_TRI[sidechain]) for backbone,sidechain in  sc_rep]


            node_match = iso.categorical_node_match(['type', 'label'], [None, None])
            Final_Graph = BACKBONE_GRAPH_TRI.copy()
            Final_ITP = BACKBONE_ITP_TRI
            processed_indexes = set()
            for backbone,sidechain in sc_rep:
                #print(backbone,sidechain)
                
                anchor_mapping = {sc_1[0][0] : backbone,
                                sc_1[0][1] : sidechain}


                matcher = self.AnchoredGraphMatcher(BACKBONE_GRAPH_TRI, Backbone_atoms_only_Graph, node_match=self.node_match2, anchor_mapping=anchor_mapping)
              
                if matcher.subgraph_is_isomorphic():
                    #print("Anchored subgraph isomorphism found!")
                    backbone_mapping = matcher.mapping
                    #print(backbone_mapping)
                    #all_maps = list(matcher.subgraph_isomorphisms_iter())
                    #print("Number of anchored isomorphisms:", len(all_maps))
                    
                else:
                    raise RaiseError("No match found that satisfies anchor constraints.")
                
                backbone_mapping = {v: k for k, v in backbone_mapping.items()} #reverse mapping
                New_Backbone = nx.relabel_nodes(Backbone_atoms_only_Graph, backbone_mapping, copy=True)

                
                side_atoms = sorted(Side_atoms_only_Graph.nodes())
                side_renumber_map = {}
                offset = max(Final_Graph.nodes()) + 1
                for i, node in enumerate(side_atoms):
                    side_renumber_map[node] = offset + i
                
                side_renumbered = nx.relabel_nodes(Side_atoms_only_Graph, side_renumber_map, copy=True)
                side_renumbered = self.align_side_chain(BACKBONE_GRAPH_TRI, Backbone_atoms_only_Graph, side_renumbered,anchor_mapping)

                dict_to_update = side_renumber_map | backbone_mapping
                

                #gafting sidechain back on backbone
                graph_to_append  = nx.compose(New_Backbone, side_renumbered)

                ##Adding any edge that was removed when splitting backbone and sidechain
                Side_chain_and_backbone = Side_chain_graph.copy()
                Side_chain_and_backbone = nx.relabel_nodes(Side_chain_and_backbone, dict_to_update, copy=True)

                E_bridge = E_bridge = set(Side_chain_and_backbone.edges()) - set(side_renumbered.edges()) - set(New_Backbone.edges())
                        
        
                #updates the ITP
                NEW_SIDE_ITP = self.renumber_itp(SIDE_CHAIN_ITP_SOLO,dict_to_update) #This just needs to be appended to the BACKBONE_ITP

                #Add new sidechian to the graph 
                Final_Graph = nx.compose(Final_Graph, side_renumbered)
                Final_Graph.add_edges_from(E_bridge)    
        



                #Appending to BACKBONE_ITP
                sections = NEW_SIDE_ITP.sections
                sections = [entry for entry in sections if entry not in ('defaults', 'atomtypes', "moleculetype")]  # Filter sections

                
                for section in sections:
                    for i, atom in enumerate(NEW_SIDE_ITP.DF[section]["atoms"]):
                            Entry = {key: NEW_SIDE_ITP[section][key][i] for key in NEW_SIDE_ITP[section] if isinstance(NEW_SIDE_ITP[section][key], list)}
                
                                    
                            if section == "atoms":
                                atom_index = int(Entry["atoms"][0])
                                if atom_index in processed_indexes or atom_index in New_Backbone.nodes():
                                    continue      
                                processed_indexes.add(atom_index)             
                            for key, value in Entry.items():
                                if key not in NEW_SIDE_ITP[section]:
                                    Final_ITP[section][key] = []

                                if isinstance(value, list):
                                    Final_ITP[section][key].append(value)
        except Exception as e:
            print(e)    
            raise RaiseError("Failed to graft sidechain onto backbone and update ITP") from e

        try:
            

            # 1) Get a numerically sorted list of all node‐labels
            sorted_nodes = sorted(Final_Graph.nodes())

            # 2) Build coords in that exact order
            coords = []
            for idx in sorted_nodes:
                xyz = Final_Graph.nodes[idx]['cord']   # ⬅ double‐check your attribute key
                coords.append([xyz[0], xyz[1], xyz[2]])

            # 3) Stack into an (N×3) array where row i is atom i (or atom i+1, if 1‐based)
            Final_ITP.coordinates = np.vstack(coords)

            cords =[ [data['cord'][0], data['cord'][1],data['cord'][2]] for node, data in Final_Graph.nodes(data=True)]
            Final_ITP.coordinates = np.vstack(cords)


            #Renumber final From left to right 
            hs_ontrimer = [renumber_map_BACKBONE_TRI[x] for x  in hs_ontrimer]

            Final_Graph,map = self.renumber_graph(Final_Graph, hs_ontrimer[0], hs_ontrimer[1])
            Final_ITP = self.renumber_itp(Final_ITP,map)

            # 1) Get a numerically sorted list of all node‐labels
            sorted_nodes = sorted(Final_Graph.nodes())

            # 2) Build coords in that exact order
            coords = []
            for idx in sorted_nodes:
                xyz = Final_Graph.nodes[idx]['cord']   # ⬅ double‐check your attribute key
                coords.append([xyz[0], xyz[1], xyz[2]])

            # 3) Stack into an (N×3) array where row i is atom i (or atom i+1, if 1‐based)
            Final_ITP.coordinates = np.vstack(coords)

        except Exception as e:
            print(e)    
            raise RaiseError("Failed to update final ITP coordinates") from e
        


        try:

            #Reorder_ITP in ascending order
            atom_list = Final_ITP['atoms']['atoms']

            atom_list = [atom[0] for atom in atom_list]
            Ascending_order_list = {int(value): index + 1 for index, value in enumerate(atom_list)}

            Final_ITP = self.renumber_itp(Final_ITP,Ascending_order_list)


            #Reorder Cords in Ascending order
            index_to_old_index = {new: old for old, new in Ascending_order_list.items()}

            # 2. Create reordered coordinates
            reordered_coordinates = [None] * len(Final_ITP.coordinates)
            for new_index in range(len(Final_ITP.coordinates)):
                old_index = index_to_old_index[new_index + 1]  # +1 if indices are 1-based
                reordered_coordinates[new_index] = Final_ITP.coordinates[old_index - 1]  # -1 to match Python 0-indexing

            # 3. Replace the original coordinates
            Final_ITP.coordinates = reordered_coordinates




            Final_ITP.save_itp(f"{Out_Name}.itp")

            coords = Final_ITP.coordinates
            box    = None
            atomnames = Final_ITP["atoms"]["atom_name"] # must match coords length
            atomnames = [atom[0] for atom in atomnames]
            Final_ITP.write_gro(
                f"{Out_Name}.gro",
                coords,
                box,
                title="Polymer Block",
                resname="BOBX",
                atomname="X   ",        # fallback if atomnames=None
                atomnames=atomnames,
                residue_number=1,
                start_index=1)
            
        except Exception as e:
            print(e)    
            raise RaiseError("Failed to print final ITP and GRO") from e



################################################### Build Graphs from ITPs
    def Build_ITP_Nodes(self,ITP):
        ITP_Graph = nx.Graph()
        try:
            for atom in ITP.DF["atoms"]["atoms"]:
                Atom_Index = int(atom[0])
                #print(Atom_Index)

                attributes = {}
                for section in ITP.DF["atoms"].keys():
                    if section in ["Comments_top", "Coments_bottom", "atoms"]: #Want to skip top and bottom comments and atoms section (it is the node index)
                        continue 
                    #print(section,ITP.DF["atoms"][section][Atom_Index -1 ] )
                    attributes[section] = ITP.DF["atoms"][section][Atom_Index -1 ]


                    try:
                        attributes['cord'] = ITP.coordinates[Atom_Index -1]
                    except:
                        attributes['cord'] = np.nan
                    #print(ITP.DF["atoms"][section][Atom_Index - 1])

                ITP_Graph.add_node(Atom_Index, **attributes)
                
                
        except:
            raise RaiseError("Failed to Turn Atom Numbers to Nodes")
        
        return ITP_Graph

    #Build edges based on bonds/angles/dihedrals
    def Build_ITP_Vertices(self,ITP_Graph, ITP_OBJ):
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
        
        RaiseError("ITP is not a connected Graph, Only Bonds/Angles/Ryckaert-Bellemans dihedral Are Check. Please reach out for more Implementation")


    #Renumbers ITP based on renumber_map        
    def renumber_itp(self,itp, renumber_map):
        """
        Returns a new ITP with *all* atom indices remapped according to renumber_map.
        Assumes renumber_map: { old_int: new_int }.
        """
        NEW = copy.deepcopy(itp)

        for section in NEW.sections:
            if section in ('defaults','atomtypes','moleculetype','dihedral_restraints'):
                continue
            #print(itp[section]["atoms"])
            for i, atoms in enumerate((NEW[section]["atoms"])):
                
                
                Entry = []
                for a in atoms:
                    Entry.append(str(renumber_map[int(a)]))
                NEW[section]["atoms"][i] = Entry 
            #print(NEW[section]["atoms"])        


        return NEW


    #Renumber Graph based on main path from start_node to end_node
    def renumber_graph(self,G, start_node, end_node):
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
    


    def get_sidechain_atoms(self,ITP_Graph, Bridge):
        """
        Given a graph of the molecule and a list of (backbone_atom, sidechain_atom) pairs,
        return a set of all atom indices belonging to the sidechain(s).
        """
        sidechain_atoms = set()
        for backbone_atom, sidechain_atom in Bridge:
            visited = set()
            queue = [sidechain_atom]
    
            #sidechain_atoms.add(Bridge[0][0])
            #sidechain_atoms.add(Bridge[0][1])

            while queue:
                current = queue.pop()
                if current in visited or current == backbone_atom:
                    continue
                visited.add(current)
                queue.extend(ITP_Graph.neighbors(current))
            
            sidechain_atoms.update(visited)

        #sidechain_atoms.discard(Bridge[0]) #Drops the starting bridge, it will get readded in get_SIDE_ITP but dont want to include it as a "relavent" sidechain
        
        return sidechain_atoms


    def get_side_ITP(self,ITP_OBJ, sidechain_atoms):
        """
        Given an ITP object, graph, and a list of side_chain_Atoms, return a new ITP object
        containing only the atoms and interactions related to the sidechain(s).
        I apologize in advance for the number of loops, it makes me sad to see
        """

        ### SET UP NEW ITP
        sections = ITP_OBJ.sections
        sections = [entry for entry in sections if entry not in ('defaults', 'atomtypes', "moleculetype")]  # Filter sections

        # Initialize NEW_ITP and clear non-atom data sections
        NEW_ITP = copy.deepcopy(ITP_OBJ)  # Use deepcopy to copy the structure without references to the original data
        for section in ITP_OBJ.sections:
            if section in ('defaults', 'atomtypes', "moleculetype"): #This doesnt change 
                continue
            NEW_ITP[section].clear()  # Clear the atom data
        

        ### Find Every Atom Relavent to SideChains
        relavent_atoms = set()
        for section in sections:
            n_entries = len(ITP_OBJ[section]["atoms"]) if "atoms" in ITP_OBJ[section] else len(next(iter(ITP_OBJ[section].values())))
            
            for i in range(n_entries):
                atoms = ITP_OBJ[section]["atoms"][i]
                if any(int(atom) in sidechain_atoms for atom in atoms):
                    relavent_atoms |= set(int(a) for a in atoms)

        
        for section in sections:
            n_entries = len(ITP_OBJ[section]["atoms"]) if "atoms" in ITP_OBJ[section] else len(next(iter(ITP_OBJ[section].values())))

            if section == "atoms": #I need to look at "relavent atoms in the atom section"
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
        
        return relavent_atoms,NEW_ITP

    def get_backbone_ITP(self,ITP_OBJ, ITP_Graph, sidechain_atoms):
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
  
    def node_match2(self,n1, n2, mass_tolerance=0.5):
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

    def get_coords(self,graph, nodes):
        """Extract coordinates as an Nx3 numpy array from nodes in the given order."""
        return np.array([graph.nodes[n]['cord'] for n in nodes])

    def align_side_chain(self,New_Backbone, Backbone_atoms_only_Graph, side_renumbered,anchor_mapping):
        """
        Repositions `side_renumbered` so that it has the same spatial relationship
        to `New_Backbone` as `Side_atoms_only_Graph` had to `Backbone_atoms_only_Graph`.
        """
        old_list = []
        new_list = []
        for old_bb_node, new_bb_node in anchor_mapping.items():
            old_list.append(Backbone_atoms_only_Graph.nodes[old_bb_node]['cord'])
            new_list.append(New_Backbone.nodes[new_bb_node]['cord'])

        old_coords = np.array(old_list)  # shape (K, 3)
        new_coords = np.array(new_list)
        
        r, _ = R.align_vectors(new_coords, old_coords)
        R_mat = r.as_matrix()

        for node, data in side_renumbered.nodes(data=True):
            original_cord = np.array(side_renumbered.nodes[node]['cord'])  # use un-renumbered source
            centered = original_cord - old_coords.mean(axis=0)
            rotated = R_mat.dot(centered)
            transformed = rotated + r.apply(centered) + new_coords.mean(axis=0)
            data['cord'] = transformed.tolist()
        
        
        return side_renumbered
        

################################################ EXAMPLE
'''
sc_1 = [(14,16)]    # What bond binds the sidechain to the backbone on the MONOMER. [(A,B)] A would be the atom on the backbone, B would be the first atom on the sidechain
sc_rep = [(89,95),(53,59),(83,86),(47,50),(17,23),(11,14)]  # What bonds bind the sidechains, that you wish to replace, to the backbone on the POLYMER [(A,B),(C,D)] A and C would be the atoms on the backbone, B and D would be the atom on the sidechains
hs_on_monomer = [1,54]  # These are the hydrogens the polymer will propagate from if it were to extend
hs_ontrimer = [109,110]
itp_file_path_tri = 'Tirmer.itp'
gro_file_path_tri = "Trimer.gro"

itp_file_path_1mer = 'Monomer_Sidechain.itp'
gro_file_path_1mer = 'Monomer_Sidechain.gro'



POLY_Swap_Instance = POLY_Swap(
    itp_file_path_tri,
    gro_file_path_tri,
    itp_file_path_1mer,
    gro_file_path_1mer,
    sc_1,
    sc_rep,
    hs_on_monomer,
    hs_ontrimer
)
'''