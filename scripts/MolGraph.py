import networkx as nx
from itp_parser import Itp_parser

class MolecularGraph:
    """
    A class to parse an ITP file, construct a molecular graph, and assign atomic symbols to nodes.
    """

    # Define atomic symbols inside the class
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

    def __init__(self, itp_file):
        """
        Initializes the class by parsing the ITP file.

        Args:
            itp_file (str): Path to the ITP file.
        """
        self.parser = Itp_parser(itp_file)  # Parse the ITP file
        self.G = nx.Graph()  # Initialize the graph
        self.build_graph()  # Construct the molecular graph

    def build_graph(self):
        """
        Constructs the molecular graph using bonds and constraints from the ITP file.
        """
        relevant_sections = {"BONDS", "CONSTRAINTS"}
        edges_set = set()

        for section in self.parser.DF:
            section_type = None
            for category, sub_sections in self.parser.TOTAL_DIRECTIVES.items():
                if section in sub_sections and category in relevant_sections:
                    section_type = category
                    break

            if section_type is not None:
                # Extract atom connections (edges)
                atoms_data = self.parser.DF[section]["atoms"]

                for atom_list in atoms_data:
                    if isinstance(atom_list, str):
                        atom_list = atom_list.split()

                    atom_list = [int(atom) for atom in atom_list]

                    for i in range(len(atom_list) - 1):
                        edge = tuple(sorted([atom_list[i], atom_list[i + 1]]))
                        edges_set.add(edge)

        # Add edges to graph
        self.G.add_edges_from(edges_set)

    def extract_atomic_symbol(self, atom_name):
        atomic_set = set(self.atomic_symbols)  # Convert list to set for fast lookup
        for length in range(1, 3):  # Atomic symbols are 1 or 2 letters long
            symbol_candidate = atom_name[:length]
            if symbol_candidate in atomic_set:
                return symbol_candidate
        raise ValueError(f"⚠️ Atomic symbol not found in atom_name: {atom_name}")

    def assign_atoms(self):
        """
        Assigns cleaned atom types (atomic symbols) to each node in the graph.
        Optimized to reduce redundant lookups and unnecessary printing.
        """
        if "atoms" in self.parser.DF:
            atom_data = self.parser.DF["atoms"]

            if "atom_name" in atom_data and "atoms" in atom_data:
                # Convert atomic symbol list to a set for fast lookup
                atomic_set = set(self.atomic_symbols)

                atom_id_symbol_map = {}

                for atom, name in zip(atom_data["atoms"], atom_data["atom_name"]):
                    # Unwrap if list
                    if isinstance(atom, list):
                        atom = atom[0]
                    if isinstance(name, list):
                        name = name[0]

                    try:
                        atom_idx = int(atom)
                    except Exception:
                        print(f"⚠️ Skipping invalid atom index: {atom}")
                        continue

                    for length in range(1, 3):
                        symbol_candidate = name[:length]
                        if symbol_candidate in atomic_set:
                            atom_id_symbol_map[atom_idx] = symbol_candidate
                            break
                    else:
                        atom_id_symbol_map[atom_idx] = "Unknown"
                        print(f"❓ Could not extract atomic symbol from: {name}")

                # Set attributes
                nx.set_node_attributes(self.G, atom_id_symbol_map, "atom_name")

                print("✅ Assigned Atomic Symbols to Graph Nodes (First 10):", list(self.G.nodes(data=True))[:10])

            else:
                print("⚠️ 'atom_name' or 'atoms' key not found inside 'atoms'.")
        else:
            print("⚠️ 'atoms' section not found in ITP file.")

    def get_graph(self):
        """
        Returns the NetworkX graph with assigned atomic symbols.

        Returns:
            nx.Graph: The molecular graph.
        """
        return self.G
