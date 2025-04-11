import networkx as nx
from networkx.algorithms import isomorphism
from itp_parser import Itp_parser
from MolGraph import MolecularGraph

class MolGraphIsomorphism:
    """
    A class to check if two molecular graphs are isomorphic, considering both
    atomic structure (connectivity) and atom types.
    """

    def __init__(self, G1, G2, initial_mapping=None):
        """ Initializes the molecular graph isomorphism checker. """
        self.G1 = G1
        self.G2 = G2
        self.initial_mapping = initial_mapping or {}

    def node_match(self, n1, n2):
        """ Ensures that atom types are the same for two nodes. """
        return n1.get("atom_name", "?") == n2.get("atom_name", "?")

    def iso_node_check(self, g1_node, g2_node):
        """ Ensures that g1_node and g2_node have identical atom types and bonded neighbors. """
        if not self.node_match(self.G1.nodes[g1_node], self.G2.nodes[g2_node]):
            return False  # Atom types must match

        # Get sorted lists of neighbor atom types
        g1_neighbors = sorted(self.G1.nodes[n]["atom_name"] for n in self.G1.neighbors(g1_node))
        g2_neighbors = sorted(self.G2.nodes[n]["atom_name"] for n in self.G2.neighbors(g2_node))

        return g1_neighbors == g2_neighbors  # Ensure neighbors match in type & count

    def check_isomorphism(self):
        """ Checks if G1 and G2 are isomorphic under the given constraints. """

        # Initialize the GraphMatcher
        GM = isomorphism.GraphMatcher(self.G1, self.G2, node_match=self.node_match)

        # 🔍 Step 1: Validate Initial Mapping
        if self.initial_mapping:
            print("\n🔍 Using Initial Mapping:", self.initial_mapping)

            # Ensure the given nodes exist in both graphs
            for g1_node, g2_node in self.initial_mapping.items():
                if g1_node not in self.G1.nodes or g2_node not in self.G2.nodes:
                    print(f"❌ ERROR: Invalid node in initial mapping! {g1_node} ↔ {g2_node} not found.")
                    return False, {}

                # Ensure pre-mapped nodes are isomorphic in **both atom type and connectivity**
                if not self.iso_node_check(g1_node, g2_node):
                    print(f"❌ ERROR: Nodes {g1_node} ↔ {g2_node} do not have matching connectivity!")
                    return False, {}

            # **Step 2: Enforce Initial Mapping as a Constraint**
            # Set GM.mapping directly before running is_isomorphic()
            GM.mapping = self.initial_mapping.copy()

            # Sort nodes so initial mapping nodes are processed first
            GM.order = sorted(self.initial_mapping.keys())

            print("\n✅ Initial Mapping Verified: Proceeding with Constrained Isomorphism Search.")

        # 🔍 Step 3: Compute Full Isomorphism with Forced Mapping
        if GM.is_isomorphic():
            print("\n✅ Graphs are isomorphic!")

            # **Enforce initial mapping on final result**
            node_mapping = self.initial_mapping.copy()

            # Add remaining mappings from GraphMatcher
            for g1_node, g2_node in GM.mapping.items():
                if g1_node not in node_mapping:
                    node_mapping[g1_node] = g2_node  # Only add if not already fixed

            print("\n🔹 Node Mapping:")
            for g1_node, g2_node in node_mapping.items():
                atom_type = self.G1.nodes[g1_node]["atom_name"]
                print(f"G1 Node {g1_node} ({atom_type}) ↔ G2 Node {g2_node} ({atom_type})")

            return True, node_mapping
        else:
            print("\n❌ Graphs are NOT isomorphic, even with the given initial mapping!")
            return False, {}
