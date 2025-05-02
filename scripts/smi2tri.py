from rdkit import Chem
import networkx as nx
import matplotlib.pyplot as plt
from rdkit.Chem import Draw

class Smi2Tri:
    def __init__(self, smiles, n=3):
        self.smiles = smiles
        self.n = n
        self.base_graph, self.atom_labels = self.smiles_to_nx(smiles)
        self.base_named_nodes = {}
        self.G_poly = None
        self.repeated_named_nodes = {}

    def smiles_to_nx(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        mol = Chem.AddHs(mol)
        G = nx.Graph()
        atom_labels = {}
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            G.add_node(idx, label=symbol)
            atom_labels[idx] = symbol
        for bond in mol.GetBonds():
            G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=str(bond.GetBondType()))
        return G, atom_labels

    def draw_monomer(self):
        mol = Chem.AddHs(Chem.MolFromSmiles(self.smiles))
        for atom in mol.GetAtoms():
            atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
        img = Draw.MolToImage(mol, size=(400, 400), kekulize=True)
        plt.imshow(img)
        plt.axis('off')
        plt.title("Initial Molecule with Atom Indices")
        plt.show()

    def assign_nodes(self, user_input):
        label_str, index_str = user_input.split("=")
        labels = label_str.strip(" {}\n").split(",")
        indices = list(map(int, index_str.strip(" {}\n").split(",")))
        if len(labels) != len(indices):
            raise ValueError("Mismatch between labels and indices")
        selected_nodes = {}
        for lbl, idx in zip(labels, indices):
            atom = self.atom_labels.get(idx, "Unknown")
            selected_nodes[lbl] = (idx, atom)
            globals()[lbl] = f"Node {lbl}: index={idx}, atom={atom}"
        self.base_named_nodes = selected_nodes
        self.G_poly, self.repeated_named_nodes = self.repeat_graph()

    def repeat_graph(self):
        G_poly = nx.Graph()
        expanded_labels = {}
        nodes_per_copy = self.base_graph.number_of_nodes()
        for i in range(self.n):
            offset = i * nodes_per_copy
            mapping = {node: node + offset for node in self.base_graph.nodes()}
            G_copy = nx.relabel_nodes(self.base_graph, mapping)
            for node, data in G_copy.nodes(data=True):
                G_poly.add_node(node, **data)
            for u, v, data in G_copy.edges(data=True):
                G_poly.add_edge(u, v, **data)
            for label, (base_idx, atom) in self.base_named_nodes.items():
                new_idx = base_idx + offset
                new_label = f"{label}_{i+1}"
                expanded_labels[new_label] = (new_idx, atom)
                globals()[new_label] = f"Node {new_label}: index={new_idx}, atom={atom}"
        return G_poly, expanded_labels

    def polymerize(self):
        for a in range(1, self.n):
            C_a_key = f"C_{a}"
            B_next_key = f"B_{a+1}"
            D_a_key = f"D_{a}"
            A_next_key = f"A_{a+1}"
            try:
                C_idx = self.repeated_named_nodes[C_a_key][0]
                B_next_idx = self.repeated_named_nodes[B_next_key][0]
                if self.G_poly.has_node(C_idx) and self.G_poly.has_node(B_next_idx):
                    self.G_poly.add_edge(C_idx, B_next_idx, bond_type="custom")
                if D_a_key in self.repeated_named_nodes:
                    d_node = self.repeated_named_nodes[D_a_key][0]
                    if self.G_poly.has_node(d_node):
                        self.G_poly.remove_node(d_node)
                if A_next_key in self.repeated_named_nodes:
                    a_node = self.repeated_named_nodes[A_next_key][0]
                    if self.G_poly.has_node(a_node):
                        self.G_poly.remove_node(a_node)
            except KeyError as e:
                print(f"⚠️ Missing symbolic label: {e}")
            except nx.NetworkXError as e:
                print(f"⚠️ Graph operation failed: {e}")

    def draw_graph(self):
        existing_nodes = set(self.G_poly.nodes)
        safe_highlights = [v[0] for k, v in self.repeated_named_nodes.items() if v[0] in existing_nodes]
        safe_label_map = {k: v for k, v in self.repeated_named_nodes.items() if v[0] in existing_nodes}
        pos = nx.kamada_kawai_layout(self.G_poly)
        default_nodes = list(set(self.G_poly.nodes) - set(safe_highlights))
        labels = {i: str(i) for i in self.G_poly.nodes}
        edge_labels = nx.get_edge_attributes(self.G_poly, 'bond_type')
        nx.draw_networkx_nodes(self.G_poly, pos, nodelist=default_nodes, node_color='lightblue', node_size=500)
        nx.draw_networkx_nodes(self.G_poly, pos, nodelist=safe_highlights, node_color='red', node_size=500)
        nx.draw_networkx_labels(self.G_poly, pos, labels, font_color='black')
        overlay = {v[0]: k for k, v in safe_label_map.items()}
        nx.draw_networkx_labels(self.G_poly, pos, overlay, font_color='white', font_weight='bold')
        nx.draw_networkx_edges(self.G_poly, pos)
        nx.draw_networkx_edge_labels(self.G_poly, pos, edge_labels=edge_labels)
        plt.axis('off')
        plt.title("Molecular Graph (Kamada-Kawai Layout)")
        plt.show()

    def to_smiles(self):
        mol = Chem.RWMol()
        index_map = {}
        for node in sorted(self.G_poly.nodes()):
            atom_symbol = self.G_poly.nodes[node].get("label", "C")
            atom = Chem.Atom(atom_symbol)
            new_idx = mol.AddAtom(atom)
            index_map[node] = new_idx
        for u, v, data in self.G_poly.edges(data=True):
            bond_type_str = data.get("bond_type", "SINGLE").upper()
            try:
                bond_type = {
                    "SINGLE": Chem.rdchem.BondType.SINGLE,
                    "DOUBLE": Chem.rdchem.BondType.DOUBLE,
                    "TRIPLE": Chem.rdchem.BondType.TRIPLE,
                    "AROMATIC": Chem.rdchem.BondType.AROMATIC,
                    "CUSTOM": Chem.rdchem.BondType.SINGLE,
                }.get(bond_type_str, Chem.rdchem.BondType.SINGLE)
                mol.AddBond(index_map[u], index_map[v], bond_type)
            except Exception as e:
                print(f"⚠️ Failed to add bond ({u}, {v}): {e}")
        try:
            Chem.SanitizeMol(mol)
            final_mol = mol.GetMol()
            smiles = Chem.MolToSmiles(final_mol, canonical=True)
            mol_no_h = Chem.RemoveHs(final_mol)
            img = Draw.MolToImage(mol_no_h, size=(400, 400))
            plt.imshow(img)
            plt.axis('off')
            plt.title("Final Molecule (No Hydrogens)")
            plt.show()
            return smiles
        except Exception as e:
            print(f"⚠️ RDKit sanitization failed: {e}")
            return None

# Step 1: Initialize and visualize the original molecule
g = Smi2Tri("c1ccccc1")  # n=3 by default
g.draw_monomer()  # shows Hs and indices

# Step 2: Assign symbolic labels
g.assign_nodes("{A,B,C,D} = {10,4,0,6}")
g.polymerize()
final_smiles = g.to_smiles()
print("Final SMILES:", final_smiles)