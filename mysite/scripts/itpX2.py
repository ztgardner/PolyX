
import networkx as nx
from scripts.itp_parser import Itp_parser
import copy
import numpy as np


class RaiseError(Exception):
    pass

class POLY_Extend:
    def __init__(
        self,
        ITP_FILE: str,
        GRO_FILE: str,
        Repeat: int,
        H_Start: int,
        H_End: int,
        Bridge_1: list,
        Bridge_2: list,
        Out_Name: str
    ):
        """
        Initializes the object with input and output files, repeat count, hydrogen indices,
        and bridge atom labels.

        Parameters:
            ITP_FILE (str): Path to the .itp file.
            GRO_FILE (str): Path to the .gro file.
            Repeat (int): Number of repeat units.
            H_Start (int): Index hydrogen atoms starting polymer.
            H_End (int): Index hydrogen atoms ending polymer.
            Bridge_1 (str): Name of the first bridge atom.
            Bridge_2 (str): Name of the second bridge atom.
            Out_Name (str): Output file base name.
        """
        self.ITP_FILE = ITP_FILE
        self.GRO_FILE = GRO_FILE
        self.Repeat = Repeat
        self.H_Start = H_Start
        self.H_end = H_End
        self.Bridge_1 = Bridge_1
        self.Bridge_2 = Bridge_2
        self.Out_Name = Out_Name


        #LoadITP
        ITP = Itp_parser(self.ITP_FILE)
        ITP.load_gro(self.GRO_FILE)
        #Build Nodes
        ITP_Graph = self.Build_ITP_Nodes(ITP)
        #Build Vertices (to make graph complete)
        ITP_Graph = self.Build_ITP_Vertices(ITP_Graph,ITP)

        #Reumber Graph and ITP
        ITP_Graph,renumber_map = self.renumber_graph(ITP_Graph,self.H_Start,self.H_end)
        B1 = [renumber_map[Bridge_1[0]], renumber_map[Bridge_1[1]]]
        B2 = [renumber_map[Bridge_2[0]], renumber_map[Bridge_2[1]]]
        start_subgraph, middle_subgraph, end_subgraph ,Division_map = self.Divide_Graph(ITP_Graph,B1,B2)
        ITP = self.renumber_itp(ITP,renumber_map)

        #Extend
        NEW_ITP = self.repeat_middle_section_ITP(ITP, Division_map, repeat=self.Repeat)
        NEW_ITP.coordinates = self.repeat_middle_section_cords(NEW_ITP,ITP_Graph,Division_map,repeat=self.Repeat)

        #Reorder_ITP in ascending order
        atom_list = NEW_ITP['atoms']['atoms']
        atom_list = [atom[0] for atom in atom_list]
        Ascending_order_list = {int(value): index + 1 for index, value in enumerate(atom_list)}
        NEW_ITP = self.renumber_itp(NEW_ITP,Ascending_order_list)

        #Reorder Cords in Ascending order
        index_to_old_index = {new: old for old, new in Ascending_order_list.items()}
        reordered_coordinates = [None] * len(NEW_ITP.coordinates)
        for new_index in range(len(NEW_ITP.coordinates)):
            old_index = index_to_old_index[new_index + 1]  # +1 if indices are 1-based
            reordered_coordinates[new_index] = NEW_ITP.coordinates[old_index - 1]  # -1 to match Python 0-indexing
        NEW_ITP.coordinates = reordered_coordinates


        # Calculate box size
        coords = np.array(NEW_ITP.coordinates)

        min_coords = np.min(coords, axis=0)
        max_coords = np.max(coords, axis=0)

        # Calculate box size


          ### WRITE
        box = max_coords - min_coords
        atomnames = NEW_ITP["atoms"]["atom_name"] # must match coords length
        atomnames = [atom[0] for atom in atomnames]

        NEW_ITP.write_gro(
            f"{self.Out_Name}.gro",
            coords,
            box,
            title="GRO made through BOBCAT",
            resname="BOBX",
            atomname="X   ",        # fallback if atomnames=None
            atomnames=atomnames,
            residue_number=1,
            start_index=1
        )


        NEW_ITP.save_itp(f"{self.Out_Name}.itp")

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
                    attributes[section] = ITP.DF["atoms"][section][Atom_Index -1 ]
                    attributes['cord'] = ITP.coordinates[Atom_Index -1]
                    #print(ITP.DF["atoms"][section][Atom_Index - 1])
                ITP_Graph.add_node(Atom_Index, **attributes)


        except:
            raise RaiseError("Failed to Turn Atom Numbers to Nodes")

        return ITP_Graph

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


    def Divide_Graph(self,ITP_Graph, Bond1,Bond2):
        ITP_Graph_temp = ITP_Graph.copy()


        ITP_Graph_temp.remove_edge(*Bond1)
        ITP_Graph_temp.remove_edge(*Bond2)

        components = list(nx.connected_components(ITP_Graph_temp))


        if len(components) != 3:
          raise ValueError(f"Could not break ITP into 3 Section, Broke into: {len(components)}.")


        section_map = {}
        comp_labels = {}

        #Plug in Node: Get its section (RN its 0 1 2)
        node_to_comp = {}
        for idx, comp in enumerate(components):
            for node in comp:
                node_to_comp[node] = idx


        b1_comp = {node_to_comp[Bond1[0]], node_to_comp[Bond1[1]]}
        b2_comp = {node_to_comp[Bond2[0]], node_to_comp[Bond2[1]]}

        # Ensure each bond bridges two *different* components
        if len(b1_comp) != 2 or len(b2_comp) != 2:
            raise ValueError("Either Bond1 or Bond2 connects nodes in the same component. They should bridge different components.")

        all_comps = {0, 1, 2}
        start_comp = list(b1_comp - b2_comp)[0] #Starting Block cannot have both
        end_comp = list(b2_comp - b1_comp)[0] #Ending Block Cannot have both
        middle_comp = list(all_comps - {start_comp, end_comp})[0] #Middle Block needs to span Both Start and End Block

        #Now I need to make sure the starting Block has the lowest node in it
        min_start = min(components[start_comp])
        min_end = min(components[end_comp])
        if min_end < min_start:
          start_comp, end_comp = end_comp, start_comp
          #Just switching the sections if the end block has a larger min value than the start

        section_map = {}
        for node, comp_idx in node_to_comp.items():
            if comp_idx == start_comp:
                section_map[node] = "Start"
            elif comp_idx == middle_comp:
                section_map[node] = "Middle"
            elif comp_idx == end_comp:
                section_map[node] = "End"


        #Divides everything up into sugraphs
        sections = {"Start": [], "Middle": [], "End": []}

        # Group nodes by their section
        for node, section in section_map.items():
            sections[section].append(node)

        # Create subgraphs for each section
        subgraphs = {
            section: ITP_Graph.subgraph(nodes).copy()
            for section, nodes in sections.items()
        }

        return subgraphs["Start"], subgraphs["Middle"], subgraphs["End"], section_map


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


    def repeat_middle_section_ITP(self,ITP_obj, section_map, repeat):

          Offset = sum(1 for section in section_map.values() if section == "Middle")  # This is the number of atoms in the Middle Section

          sections = ITP_obj.sections
          sections = [entry for entry in sections if entry not in ('defaults', 'atomtypes', "moleculetype", 'dihedral_restraints')]  # Filter sections

          # Initialize NEW_ITP and clear non-atom data sections
          NEW_ITP = copy.deepcopy(ITP_obj)  # Use deepcopy to copy the structure without references to the original data
          for section in ITP_obj.sections:
              if section in ('defaults', 'atomtypes', "moleculetype", 'dihedral_restraints'):
                  continue
              NEW_ITP[section].clear()  # Clear the atom data

          # Process sections and repeat atoms
          for section in sections:
              for n in range(repeat):
                  for i, atom in enumerate(ITP_obj.DF[section]["atoms"]):
                      # Map atom nodes to their respective section types (e.g., Start, Middle, End)
                      local_map = {node: section_map.get(int(node), "Unknown") for node in atom}

                      # Create entry from the current atom data
                      Entry = {key: ITP_obj[section][key][i] for key in ITP_obj[section] if isinstance(ITP_obj[section][key], list) and len(ITP_obj[section][key]) > i}


                      # Check if there is any "Unknown" section
                      if any(value == "Unknown" for value in local_map.values()):
                          raise ValueError("Unknown Section Found")

                      # Modify indices based on section and iteration
                      if any(value == "Start" for value in local_map.values()):
                          if n == 0:
                              #write atoms as is if we are at the start and ready for the start
                              pass
                          else:
                              continue
                      elif any(value == "Middle" for value in local_map.values()):
                          # Case 2: Entries in "STart to Middle" or "Middle"
                          Entry['atoms'] = [str(int(atom) + Offset * n) for atom in Entry['atoms']]

                      elif any(value == "End" for value in local_map.values()):
                          # Case 3: Entries in "End" – only modify at the last iteration (transition from Middle to End) or just the end block
                          if  n == repeat - 1:
                              Entry['atoms'] = [str(int(atom) + Offset * n) for atom in Entry['atoms']]
                          else:
                              #if we have an atom at the end but arent ready for the end just skip
                              continue


                      # Append modified Entry data to NEW_ITP
                      for key, value in Entry.items():
                          if key not in NEW_ITP[section]:
                              NEW_ITP[section][key] = []

                          if isinstance(value, list):
                              NEW_ITP[section][key].append(value)

          return NEW_ITP


    def repeat_middle_section_cords(self,ITP_obj,ITP_map, section_map, repeat):
          # 1) Grab your template coords as an (N,3) array
          coords = np.array(ITP_obj.coordinates)

          coords = {node: (data['cord'][0], data['cord'][1],data['cord'][2]) for node, data in ITP_map.nodes(data=True)}
          coords = np.array([coords[i] for i in sorted(coords.keys())])

          # 2) Extract the 1-based indices of the middle section
          middle_ids = sorted(i for i, m in section_map.items() if m == "Middle")
          if not middle_ids:
              raise ValueError("No atoms labeled 'Middle' in section_map")

          first_mid_i = middle_ids[0] - 1
          last_mid_i  = middle_ids[-1] - 1

          # 3) Compute the step‐vector of just the middle block

          offset = coords[last_mid_i + 1] - coords[first_mid_i]  # shape (3,)

          #offset = offset  = Extra_Shift
          # 4) Build your new coordinate list
          new_coords = []
          for n in range(repeat):
              for atom_idx, coord in enumerate(coords):
                  role = section_map[atom_idx+1]

                  if role == "Start":
                      # only include your template Start once, at the very beginning
                      if n == 0:
                          new_coords.append(coord.copy())

                  elif role == "Middle":
                      # stamp it out each time, shifted by n⋅offset
                      new_coords.append(coord + n*offset )

                  elif role == "End":
                      # include the End only on the last repeat
                      if n == repeat - 1:
                          new_coords.append(coord + n*offset )

          new_coords = np.vstack(new_coords)
          return (new_coords.tolist())




