import pandas as pd
import networkx as nx
import numpy as np
import pandas as pd
from itp_parser import *
dh1 = [5,6,7,10]		# What is the first dihedral?
dh2 = [10,11,12,15]		# What is the second dihedral?
itp_file_path = 'polymer.itp'
n = 5
#data= Itp_parser(File=itp_file_path)
def parse_itp_file(file_path):
    data = {}
    current_section = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            elif line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                data[current_section] = []
                if current_section in ['angles', 'moleculetype', 'atoms']:
                    next(file)  # Skip one line for certain sections
            else:
                line_data = line.split(';')[0].strip()
                if line_data:
                    if current_section == 'dihedrals':
                        if 'IMPROPER DIHEDRAL' in line_data:
                            next(file)  # Skip the next line for improper dihedrals
                            continue
                        elif 'PROPER DIHEDRAL' in line_data:
                            next(file)  # Skip the next line for proper dihedrals
                            continue
                    data[current_section].append(line_data.split())
    return data

def process_dihedrals_sections(file_path):
    improper_data = []
    proper_data = []
    with open(file_path, 'r') as file:
        section_flag = None  # Tracks the current section
        dihedrals_flag = False  # Flag to indicate dihedrals section
        section_type = None  # Tracks the type of dihedrals section
        for line in file:
            line = line.strip()
            if line.startswith('[ dihedrals ]'):
                dihedrals_flag = True
            elif dihedrals_flag and line.startswith('; IMPROPER DIHEDRAL ANGLES'):
                section_flag = 'improper'
                section_type = 'Improper Dihedrals'
            elif dihedrals_flag and line.startswith('; PROPER DIHEDRAL ANGLES'):
                section_flag = 'proper'
                section_type = 'Proper Dihedrals'
            elif '[ pairs ]' in line:
                break  # Stop collecting data after encountering [ pairs ]
            elif section_flag == 'improper' and line.strip() and not line.startswith(';'):
                improper_data.append(line.split())
            elif section_flag == 'proper' and line.strip() and not line.startswith(';'):
                proper_data.append(line.split())

    # Convert data to dataframes
    df_improper = pd.DataFrame(improper_data, columns=['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2'])
    df_proper = pd.DataFrame(proper_data, columns=['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5'])

    return df_improper, df_proper

def update_nmer_using_graph(df_atoms, df_bonds, dihedral_centers):
    G = nx.Graph()

    # Add edges from bonds
    for _, row in df_bonds.iterrows():
        G.add_edge(int(row['ai']), int(row['aj']))

    # Identify components
    for dihedral_center in dihedral_centers:
        G.remove_edge(dihedral_center[0], dihedral_center[1])

    components = list(nx.connected_components(G))
    atom_to_nmer = {}
    for nmer, component in enumerate(components, start=1):
        for atom in component:
            atom_to_nmer[atom] = nmer

    # Update nmer in df_atoms
    df_atoms['nmer'] = df_atoms['nr'].astype(int).map(atom_to_nmer)

    return df_atoms

# Example usage

dihedral_centers = [(dh1[1], dh1[2]), (dh2[1], dh2[2])]
parsed_data = parse_itp_file(itp_file_path)
df_improper, df_proper = process_dihedrals_sections(itp_file_path)

# Convert each section to a DataFrame
df_atomtypes = pd.DataFrame(parsed_data.get('atomtypes', []), columns=['type', 'type_name', 'mass', 'charge', 'element', 'sigma', 'epsilon'])
df_moleculetype = pd.DataFrame(parsed_data.get('moleculetype', []), columns=['name', 'nrexcl'])
df_atoms = pd.DataFrame(parsed_data.get('atoms', []), columns=['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass'])
df_bonds = pd.DataFrame(parsed_data.get('bonds', []), columns=['ai', 'aj', 'funct', 'c0', 'c1'])
df_angles = pd.DataFrame(parsed_data.get('angles', []), columns=['ai', 'aj', 'ak', 'funct', 'c0', 'c1'])
df_pairs = pd.DataFrame(parsed_data.get('pairs', []), columns=['ai', 'aj', 'funct'])

df_atoms = update_nmer_using_graph(df_atoms, df_bonds, dihedral_centers)

# Print the DataFrames
print("DataFrame df_atomtypes:")
print(df_atomtypes)

print("\nDataFrame df_moleculetype:")
print(df_moleculetype)

print("\nDataFrame df_atoms:")
print(df_atoms)

print("\nDataFrame df_bonds:")
print(df_bonds)

print("\nDataFrame df_angles:")
print(df_angles)

print("\nDataFrame df_pairs:")
print(df_pairs)

print("\nDataFrame df_improper:")
print(df_improper)

print("\nDataFrame df_proper:")
print(df_proper)

 
def assign_nmer_columns(df_atoms, df_bonds, df_angles, df_improper, df_proper):
    # Ensure the columns in df_atoms are integers
    df_atoms['nr'] = df_atoms['nr'].astype(int)
    df_atoms['nmer'] = df_atoms['nmer'].astype(int)

    # Create a mapping of atom number to nmer value
    atom_to_nmer = df_atoms.set_index('nr')['nmer'].to_dict()

    # Ensure the columns in other DataFrames are integers
    df_bonds[['ai', 'aj']] = df_bonds[['ai', 'aj']].astype(int)
    df_angles[['ai', 'aj', 'ak']] = df_angles[['ai', 'aj', 'ak']].astype(int)
    df_improper[['ai', 'aj', 'ak', 'al']] = df_improper[['ai', 'aj', 'ak', 'al']].astype(int)
    df_proper[['ai', 'aj', 'ak', 'al']] = df_proper[['ai', 'aj', 'ak', 'al']].astype(int)

    # Helper function to map nmer values to the respective columns
    def map_nmer(df, columns):
        for col in columns:
            nmer_col = 'nmer_' + col
            df[nmer_col] = df[col].map(atom_to_nmer)
        return df

    # Assign nmer values to df_bonds
    df_bonds = map_nmer(df_bonds, ['ai', 'aj'])

    # Assign nmer values to df_angles
    df_angles = map_nmer(df_angles, ['ai', 'aj', 'ak'])

    # Assign nmer values to df_improper
    df_improper = map_nmer(df_improper, ['ai', 'aj', 'ak', 'al'])

    # Assign nmer values to df_proper
    df_proper = map_nmer(df_proper, ['ai', 'aj', 'ak', 'al'])

    return df_bonds, df_angles, df_improper, df_proper

# Debugging step to print mappings and verify data
print("\nSample df_atoms:")
print(df_atoms)

print("\nSample df_bonds before mapping:")
print(df_bonds)

df_bonds, df_angles, df_improper, df_proper = assign_nmer_columns(df_atoms, df_bonds, df_angles, df_improper, df_proper)

print("\nSample df_bonds after mapping:")
print(df_bonds)

print("\nSample df_angles after mapping:")
print(df_angles)

print("\nSample df_improper after mapping:")
print(df_improper)

print("\nSample df_proper after mapping:")
print(df_proper)
 
def filter_nmer(df_atoms, df_bonds, df_angles, df_improper, df_proper, nmer_value=2):
    df_atoms_mm = df_atoms[df_atoms['nmer'] == nmer_value].copy()

    df_bonds_mm = df_bonds[(df_bonds['nmer_ai'] == nmer_value) | (df_bonds['nmer_aj'] == nmer_value)].copy()

    df_angles_mm = df_angles[(df_angles['nmer_ai'] == nmer_value) |
                             (df_angles['nmer_aj'] == nmer_value) |
                             (df_angles['nmer_ak'] == nmer_value)].copy()

    df_improper_mm = df_improper[(df_improper['nmer_ai'] == nmer_value) |
                                 (df_improper['nmer_aj'] == nmer_value) |
                                 (df_improper['nmer_ak'] == nmer_value) |
                                 (df_improper['nmer_al'] == nmer_value)].copy()

    df_proper_mm = df_proper[(df_proper['nmer_ai'] == nmer_value) |
                             (df_proper['nmer_aj'] == nmer_value) |
                             (df_proper['nmer_ak'] == nmer_value) |
                             (df_proper['nmer_al'] == nmer_value)].copy()

    return df_atoms_mm, df_bonds_mm, df_angles_mm, df_improper_mm, df_proper_mm

def adjust_mm_dataframes(df_atoms_mm, df_bonds_mm, df_angles_mm, df_improper_mm, df_proper_mm):
    # Calculate the minimum value in the nr column of df_atoms_mm
    min_nr_atoms = df_atoms_mm['nr'].min()

    # Define the adjustment value
    adjustment_value = min_nr_atoms - 1

    # Adjust the nr column in df_atoms_mm
    df_atoms_mm['nr'] = df_atoms_mm['nr'] - adjustment_value

    # Adjust the ai and aj columns in df_bonds_mm
    df_bonds_mm['ai'] = df_bonds_mm['ai'] - adjustment_value
    df_bonds_mm['aj'] = df_bonds_mm['aj'] - adjustment_value

    # Adjust the ai, aj, and ak columns in df_angles_mm
    df_angles_mm['ai'] = df_angles_mm['ai'] - adjustment_value
    df_angles_mm['aj'] = df_angles_mm['aj'] - adjustment_value
    df_angles_mm['ak'] = df_angles_mm['ak'] - adjustment_value

    # Adjust the ai, aj, ak, and al columns in df_improper_mm
    df_improper_mm['ai'] = df_improper_mm['ai'] - adjustment_value
    df_improper_mm['aj'] = df_improper_mm['aj'] - adjustment_value
    df_improper_mm['ak'] = df_improper_mm['ak'] - adjustment_value
    df_improper_mm['al'] = df_improper_mm['al'] - adjustment_value

    # Adjust the ai, aj, ak, and al columns in df_proper_mm
    df_proper_mm['ai'] = df_proper_mm['ai'] - adjustment_value
    df_proper_mm['aj'] = df_proper_mm['aj'] - adjustment_value
    df_proper_mm['ak'] = df_proper_mm['ak'] - adjustment_value
    df_proper_mm['al'] = df_proper_mm['al'] - adjustment_value

    return df_atoms_mm, df_bonds_mm, df_angles_mm, df_improper_mm, df_proper_mm

# Example usage:
df_atoms_mm, df_bonds_mm, df_angles_mm, df_improper_mm, df_proper_mm = filter_nmer(df_atoms, df_bonds, df_angles, df_improper, df_proper)
df_atoms_mm, df_bonds_mm, df_angles_mm, df_improper_mm, df_proper_mm = adjust_mm_dataframes(
    df_atoms_mm, df_bonds_mm, df_angles_mm, df_improper_mm, df_proper_mm
)

print(df_atoms_mm)
print(df_bonds_mm)
print(df_angles_mm)
print(df_improper_mm)
print(df_proper_mm)
 
def filter_nmer(df_atoms, df_bonds, df_angles, df_improper, df_proper, nmer_value=1):
    df_atoms_start = df_atoms[df_atoms['nmer'] == nmer_value].copy()

    df_bonds_start = df_bonds[(df_bonds['nmer_ai'] == nmer_value) & (df_bonds['nmer_aj'] == nmer_value)].copy()

    df_angles_start = df_angles[(df_angles['nmer_ai'] == nmer_value) &
                                (df_angles['nmer_aj'] == nmer_value) &
                                (df_angles['nmer_ak'] == nmer_value)].copy()

    df_improper_start = df_improper[(df_improper['nmer_ai'] == nmer_value) &
                                    (df_improper['nmer_aj'] == nmer_value) &
                                    (df_improper['nmer_ak'] == nmer_value) &
                                    (df_improper['nmer_al'] == nmer_value)].copy()

    df_proper_start = df_proper[(df_proper['nmer_ai'] == nmer_value) &
                                (df_proper['nmer_aj'] == nmer_value) &
                                (df_proper['nmer_ak'] == nmer_value) &
                                (df_proper['nmer_al'] == nmer_value)].copy()

    return df_atoms_start, df_bonds_start, df_angles_start, df_improper_start, df_proper_start

# Example usage:
df_atoms_start, df_bonds_start, df_angles_start, df_improper_start, df_proper_start = filter_nmer(df_atoms, df_bonds, df_angles, df_improper, df_proper)
df_atoms_end, df_bonds_end, df_angles_end, df_improper_end, df_proper_end = filter_nmer(df_atoms, df_bonds, df_angles, df_improper, df_proper, nmer_value=3)

print(df_atoms_start)
print(df_bonds_start)
print(df_angles_start)
print(df_improper_start)
print(df_proper_start)

print(df_atoms_end)
print(df_bonds_end)
print(df_angles_end)
print(df_improper_end)
print(df_proper_end )

def initialize_final_dataframes(df_atoms_start, df_bonds_start, df_angles_start, df_improper_start, df_proper_start):
    df_atoms_fin = df_atoms_start.copy()
    df_bonds_fin = df_bonds_start.copy()
    df_angles_fin = df_angles_start.copy()
    df_improper_fin = df_improper_start.copy()
    df_proper_fin = df_proper_start.copy()

    max_nr_atoms = df_atoms_fin['nr'].max()
    print(f"Max nr value in df_atoms_fin: {max_nr_atoms}")

    return df_atoms_fin, df_bonds_fin, df_angles_fin, df_improper_fin, df_proper_fin

# Example usage:
df_atoms_fin, df_bonds_fin, df_angles_fin, df_improper_fin, df_proper_fin = initialize_final_dataframes(
    df_atoms_start, df_bonds_start, df_angles_start, df_improper_start, df_proper_start
)

print(df_atoms_fin)
print(df_bonds_fin)
print(df_angles_fin)
print(df_improper_fin)
print(df_proper_fin)
 
def repeat_process(df_atoms_start, df_bonds_start, df_angles_start, df_improper_start, df_proper_start,
                   df_atoms_mm, df_bonds_mm, df_angles_mm, df_improper_mm, df_proper_mm, n):
    # Initialize the final dataframes
    df_atoms_fin, df_bonds_fin, df_angles_fin, df_improper_fin, df_proper_fin = initialize_final_dataframes(
        df_atoms_start, df_bonds_start, df_angles_start, df_improper_start, df_proper_start
    )

    for i in range(1, n - 1):
        # 1. Calculate the max value in the nr column of df_atoms_fin
        max_nr_atoms = df_atoms_fin['nr'].max()

        # 2. Copy all data from each df_mm
        df_atoms_mm_copy = df_atoms_mm.copy()
        df_bonds_mm_copy = df_bonds_mm.copy()
        df_angles_mm_copy = df_angles_mm.copy()
        df_improper_mm_copy = df_improper_mm.copy()
        df_proper_mm_copy = df_proper_mm.copy()

        # 3. Adjust the copied data by adding the max value to nr, ai, aj, ak, al
        df_atoms_mm_copy['nr'] += max_nr_atoms
        df_bonds_mm_copy['ai'] += max_nr_atoms
        df_bonds_mm_copy['aj'] += max_nr_atoms
        df_angles_mm_copy['ai'] += max_nr_atoms
        df_angles_mm_copy['aj'] += max_nr_atoms
        df_angles_mm_copy['ak'] += max_nr_atoms
        df_improper_mm_copy['ai'] += max_nr_atoms
        df_improper_mm_copy['aj'] += max_nr_atoms
        df_improper_mm_copy['ak'] += max_nr_atoms
        df_improper_mm_copy['al'] += max_nr_atoms
        df_proper_mm_copy['ai'] += max_nr_atoms
        df_proper_mm_copy['aj'] += max_nr_atoms
        df_proper_mm_copy['ak'] += max_nr_atoms
        df_proper_mm_copy['al'] += max_nr_atoms

        # 4. Update nmer value
        nmer_value = 1 + i
        df_atoms_mm_copy['nmer'] = nmer_value
        df_bonds_mm_copy['nmer_ai'] = nmer_value
        df_bonds_mm_copy['nmer_aj'] = nmer_value
        df_angles_mm_copy['nmer_ai'] = nmer_value
        df_angles_mm_copy['nmer_aj'] = nmer_value
        df_angles_mm_copy['nmer_ak'] = nmer_value
        df_improper_mm_copy['nmer_ai'] = nmer_value
        df_improper_mm_copy['nmer_aj'] = nmer_value
        df_improper_mm_copy['nmer_ak'] = nmer_value
        df_improper_mm_copy['nmer_al'] = nmer_value
        df_proper_mm_copy['nmer_ai'] = nmer_value
        df_proper_mm_copy['nmer_aj'] = nmer_value
        df_proper_mm_copy['nmer_ak'] = nmer_value
        df_proper_mm_copy['nmer_al'] = nmer_value

        # 5. Append the copied and adjusted data to df_fin
        df_atoms_fin = pd.concat([df_atoms_fin, df_atoms_mm_copy])
        df_bonds_fin = pd.concat([df_bonds_fin, df_bonds_mm_copy])
        df_angles_fin = pd.concat([df_angles_fin, df_angles_mm_copy])
        df_improper_fin = pd.concat([df_improper_fin, df_improper_mm_copy])
        df_proper_fin = pd.concat([df_proper_fin, df_proper_mm_copy])

    # Return the final dataframes
    return df_atoms_fin, df_bonds_fin, df_angles_fin, df_improper_fin, df_proper_fin

# Example usage:
df_atoms_fin, df_bonds_fin, df_angles_fin, df_improper_fin, df_proper_fin = repeat_process(
    df_atoms_start, df_bonds_start, df_angles_start, df_improper_start, df_proper_start,
    df_atoms_mm, df_bonds_mm, df_angles_mm, df_improper_mm, df_proper_mm, n
)

print(df_atoms_fin)
print(df_bonds_fin)
print(df_angles_fin)
print(df_improper_fin)
print(df_proper_fin)
 
def adjust_end_dataframes(df_atoms_end, df_bonds_end, df_angles_end, df_improper_end, df_proper_end):
    # Calculate the minimum value in the nr column of df_atoms_end
    min_nr_atoms = df_atoms_end['nr'].min()

    # Define the adjustment value
    adjustment_value = min_nr_atoms - 1

    # Adjust the nr column in df_atoms_end
    df_atoms_end['nr'] = df_atoms_end['nr'] - adjustment_value

    # Adjust the ai and aj columns in df_bonds_end
    df_bonds_end['ai'] = df_bonds_end['ai'] - adjustment_value
    df_bonds_end['aj'] = df_bonds_end['aj'] - adjustment_value

    # Adjust the ai, aj, and ak columns in df_angles_end
    df_angles_end['ai'] = df_angles_end['ai'] - adjustment_value
    df_angles_end['aj'] = df_angles_end['aj'] - adjustment_value
    df_angles_end['ak'] = df_angles_end['ak'] - adjustment_value

    # Adjust the ai, aj, ak, and al columns in df_improper_end
    df_improper_end['ai'] = df_improper_end['ai'] - adjustment_value
    df_improper_end['aj'] = df_improper_end['aj'] - adjustment_value
    df_improper_end['ak'] = df_improper_end['ak'] - adjustment_value
    df_improper_end['al'] = df_improper_end['al'] - adjustment_value

    # Adjust the ai, aj, ak, and al columns in df_proper_end
    df_proper_end['ai'] = df_proper_end['ai'] - adjustment_value
    df_proper_end['aj'] = df_proper_end['aj'] - adjustment_value
    df_proper_end['ak'] = df_proper_end['ak'] - adjustment_value
    df_proper_end['al'] = df_proper_end['al'] - adjustment_value

    return df_atoms_end, df_bonds_end, df_angles_end, df_improper_end, df_proper_end

# Example usage:
df_atoms_end_adjusted, df_bonds_end_adjusted, df_angles_end_adjusted, df_improper_end_adjusted, df_proper_end_adjusted = adjust_end_dataframes(
    df_atoms_end, df_bonds_end, df_angles_end, df_improper_end, df_proper_end
)

print(df_atoms_end_adjusted)
print(df_bonds_end_adjusted)
print(df_angles_end_adjusted)
print(df_improper_end_adjusted)
print(df_proper_end_adjusted)
 
def add_end_dataframes(df_atoms_fin, df_bonds_fin, df_angles_fin, df_improper_fin, df_proper_fin,
                                  df_atoms_end, df_bonds_end, df_angles_end, df_improper_end, df_proper_end):
      # Calculate the max value in the nr column of df_atoms_fin
    max_nr_atoms = df_atoms_fin['nr'].max()

    # Adjust the copied data by adding the max value to nr, ai, aj, ak, al
    df_atoms_end['nr'] += max_nr_atoms
    df_bonds_end['ai'] += max_nr_atoms
    df_bonds_end['aj'] += max_nr_atoms
    df_angles_end['ai'] += max_nr_atoms
    df_angles_end['aj'] += max_nr_atoms
    df_angles_end['ak'] += max_nr_atoms
    df_improper_end['ai'] += max_nr_atoms
    df_improper_end['aj'] += max_nr_atoms
    df_improper_end['ak'] += max_nr_atoms
    df_improper_end['al'] += max_nr_atoms
    df_proper_end['ai'] += max_nr_atoms
    df_proper_end['aj'] += max_nr_atoms
    df_proper_end['ak'] += max_nr_atoms
    df_proper_end['al'] += max_nr_atoms

    # Update nmer value
    iteration_value = n - 1  # This is the first run-through
    nmer_value = 1 + iteration_value
    df_atoms_end['nmer'] = nmer_value
    df_bonds_end['nmer_ai'] = nmer_value
    df_bonds_end['nmer_aj'] = nmer_value
    df_angles_end['nmer_ai'] = nmer_value
    df_angles_end['nmer_aj'] = nmer_value
    df_angles_end['nmer_ak'] = nmer_value
    df_improper_end['nmer_ai'] = nmer_value
    df_improper_end['nmer_aj'] = nmer_value
    df_improper_end['nmer_ak'] = nmer_value
    df_improper_end['nmer_al'] = nmer_value
    df_proper_end['nmer_ai'] = nmer_value
    df_proper_end['nmer_aj'] = nmer_value
    df_proper_end['nmer_ak'] = nmer_value
    df_proper_end['nmer_al'] = nmer_value

    # Append the copied and adjusted data to df_fin
    df_atoms_fin = pd.concat([df_atoms_fin, df_atoms_end])
    df_bonds_fin = pd.concat([df_bonds_fin, df_bonds_end])
    df_angles_fin = pd.concat([df_angles_fin, df_angles_end])
    df_improper_fin = pd.concat([df_improper_fin, df_improper_end])
    df_proper_fin = pd.concat([df_proper_fin, df_proper_end])

    # Return the final dataframes
    return df_atoms_fin, df_bonds_fin, df_angles_fin, df_improper_fin, df_proper_fin

# Example usage:
df_atoms_fin, df_bonds_fin, df_angles_fin, df_improper_fin, df_proper_fin = add_end_dataframes(
    df_atoms_fin, df_bonds_fin, df_angles_fin, df_improper_fin, df_proper_fin,
    df_atoms_end, df_bonds_end, df_angles_end, df_improper_end, df_proper_end
)

print(df_atoms_fin)
print(df_bonds_fin)
print(df_angles_fin)
print(df_improper_fin)
print(df_proper_fin)
 
# Edit the df_atoms_fin dataframe to make the cgnr column the same as the nmer column
df_atoms_fin['cgnr'] = df_atoms_fin['nmer']

# Print the updated dataframe to verify changes
print(df_atoms_fin)
 
df_atoms_fin.reset_index(drop=True, inplace=True)

 
print(df_atoms_fin)
 
def add_charges_to_units(df_atoms_fin):
    while True:
        # Ask the user if they would like to add charges to a unit
        add_charges = input("Would you like to add charges to a unit? (y/n): ").strip().lower()

        if add_charges == 'n':
            break

        if add_charges == 'y':
            # Ask the user for the units they would like to add charge to
            units_input = input("Which unit(s) would you like to add charge to? (provide numbers separated by commas): ").strip()
            units = [int(unit) for unit in units_input.split(',')]

            # Ask the user for the name of the charges file
            charges_file_name = input("What is the name of the charges file? (don't include .txt): ").strip()
            charges_file_path = f"{charges_file_name}.txt"

            try:
                # Read the charges from the text file and convert to floats rounded to 4 decimal places
                with open(charges_file_path, 'r') as file:
                    charges = []
                    for line in file:
                        line = line.strip()
                        # Ignore lines with unexpected characters
                        try:
                            charge = round(float(line), 4)
                            charges.append(charge)
                        except ValueError:
                            print(f"Warning: Ignoring invalid line: '{line}'")

                for unit in units:
                    # Get the index for the specified unit
                    unit_indices = df_atoms_fin[df_atoms_fin['nmer'] == unit].index
                    num_atoms_in_unit = len(unit_indices)

                    # Check if the number of charges matches the number of atoms in the unit
                    if num_atoms_in_unit != len(charges):
                        print(f"Error: The number of charges provided ({len(charges)}) does not match the number of atoms ({num_atoms_in_unit}) in unit {unit}.")
                        continue

                    # Update the charges for the specified unit
                    for idx, charge in zip(unit_indices, charges):
                        df_atoms_fin.at[idx, 'charge'] = f"{charge:.4f}"

                print(f"Charges updated successfully for units: {', '.join(map(str, units))}")

            except FileNotFoundError:
                print(f"Error: The file {charges_file_path} was not found.")
            except Exception as e:
                print(f"An unexpected error occurred: {e}")


# Call the function to add charges to units
add_charges_to_units(df_atoms_fin)
 
print(df_atoms_fin)
 
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

def split_atom_name(df):
    atom_col = []
    atom_namnum_col = []

    for atom_name in df['atom']:
        atom = None
        atom_namnum = None

        # Check if the atom name starts with an atomic symbol
        for symbol in atomic_symbols:
            if atom_name.startswith(symbol):
                # Check if there are additional characters after the atomic symbol
                if len(atom_name) > len(symbol):
                    next_char = atom_name[len(symbol)]
                    if next_char.isdigit() or (next_char.isalpha() and next_char.isupper()):
                        atom = symbol
                        atom_namnum = atom_name[len(symbol):]
                        break
                    elif next_char.isalpha() and next_char.islower():
                        # Check if the full atom name matches any atomic symbol
                        if atom_name in atomic_symbols:
                            atom = atom_name
                            break

        if atom is None:
            raise ValueError(f"Atomic symbol not found in atom_name: {atom_name}")

        atom_col.append(atom)
        atom_namnum_col.append(atom_namnum)

    df['atom'] = atom_col
    df['atom_namnum'] = atom_namnum_col

    return df

# Update the atom_nam and atom_namnum columns using a custom hex style
def hex_style_atom_type(index):
    """Generate hex-style atom type based on index."""
    hex_digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    base = len(hex_digits)
    result = ''
    index -= 1
    while index > 0:
        index, remainder = divmod(index, base)
        result = hex_digits[remainder] + result
    return result.zfill(3)  # Pad with zeros to ensure three characters

def update_dataframe_with_hex_style(df):
    # Update atom_num to be sequential starting from 1
    df['atom_num'] = range(1, len(df) + 1)

    # Create a custom sequence for atom_namnum using hex_style_atom_type function
    df['atom_namnum'] = [hex_style_atom_type(i) for i in range(1, len(df) + 1)]

    return df

# Call the function to split the atom_name column
df_atoms_fin = split_atom_name(df_atoms_fin)

# Call the function to update the dataframe with the custom sequence
df_atoms_fin = update_dataframe_with_hex_style(df_atoms_fin)
df_atoms_fin['atom'] = df_atoms_fin['atom'] + df_atoms_fin['atom_namnum']

# Print the updated dataframe
print(df_atoms_fin)
 
print(df_atoms_fin)
print(df_bonds_fin)
print(df_angles_fin)
print(df_improper_fin)
print(df_proper_fin)
 
import pandas as pd

# Assuming the dataframes are already defined in your Colab notebook

# Helper function to remove duplicates based on sorted tuples of relevant columns
def remove_duplicates_by_columns(df, columns):
    # Create a sorted tuple for the specified columns
    df['sorted_tuple'] = df.apply(lambda row: tuple(sorted([row[col] for col in columns])), axis=1)
    # Drop duplicate rows based on the sorted tuple
    df_cleaned = df.drop_duplicates(subset='sorted_tuple')
    # Drop the temporary sorted tuple column
    df_cleaned = df_cleaned.drop(columns=['sorted_tuple'])
    return df_cleaned

# Remove duplicates from df_bonds_fin based on ai, aj
df_bonds_fin_cleaned = remove_duplicates_by_columns(df_bonds_fin, ['ai', 'aj'])

# Remove duplicates from df_angles_fin based on ai, aj, ak
df_angles_fin_cleaned = remove_duplicates_by_columns(df_angles_fin, ['ai', 'aj', 'ak'])

# Remove duplicates from df_proper_fin based on ai, aj, ak, al
df_proper_fin_cleaned = remove_duplicates_by_columns(df_proper_fin, ['ai', 'aj', 'ak', 'al'])

# Remove duplicates from df_improper_fin based on ai, aj, ak, al
df_improper_fin_cleaned = remove_duplicates_by_columns(df_improper_fin, ['ai', 'aj', 'ak', 'al'])

# Now, cleaned dataframes contain the cleaned versions of each dataframe
# Assign back to original variables if needed
df_bonds_fin = df_bonds_fin_cleaned
df_angles_fin = df_angles_fin_cleaned
df_proper_fin = df_proper_fin_cleaned
df_improper_fin = df_improper_fin_cleaned

print(df_bonds_fin)
print(df_angles_fin)
print(df_proper_fin)
print(df_improper_fin)
 
def remove_columns(df, columns_to_remove):
    df.drop(columns=columns_to_remove, inplace=True, errors='ignore')

def format_atoms_df(df):
    formatted_lines = []
    for index, row in df.iterrows():
        formatted_line = f"{row['nr']:>6}{row['type']:>11}{row['resnr']:>7}{row['residue']:>7}{row['atom']:>6}{row['cgnr']:>7}{row['charge']:>11}{row['mass']:>11}"
        formatted_lines.append(formatted_line)
    return '\n'.join(formatted_lines)

def format_bonds_str(df):
    formatted_lines = []
    for index, row in df.iterrows():
        formatted_line = f"{row['ai']:>5}{row['aj']:>6}{row['funct']:>6}{row['c0']:>12}{row['c1']:>11}"
        formatted_lines.append(formatted_line)
    return '\n'.join(formatted_lines)

def format_angles_str(df):
    formatted_lines = []
    for index, row in df.iterrows():
        formatted_line = f"{row['ai']:>5}{row['aj']:>6}{row['ak']:>6}{row['funct']:>6}{row['c0']:>11}{row['c1']:>11}"
        formatted_lines.append(formatted_line)
    return '\n'.join(formatted_lines)

def format_improper_str(df):
    formatted_lines = []
    for index, row in df.iterrows():
        formatted_line = f"{row['ai']:>6}{row['aj']:>6}{row['ak']:>6}{row['al']:>6}{row['funct']:>5}{row['c0']:>15}{row['c1']:>11}{row['c2']:>6}"
        formatted_lines.append(formatted_line)
    return '\n'.join(formatted_lines)

def format_proper_str(df):
    formatted_lines = []
    for index, row in df.iterrows():
        formatted_line = f"{row['ai']:>5}{row['aj']:>5}{row['ak']:>5}{row['al']:>5}{row['funct']:>9}{row['c0']:>12}{row['c1']:>8}{row['c2']:>8}{row['c3']:>8}{row['c4']:>8}{row['c5']:>8}"
        formatted_lines.append(formatted_line)
    return '\n'.join(formatted_lines)

###############################################
# AFTER EXTENSION
###############################################
# list of columns to remove
columns_to_remove = [ "nmer_ai",  "nmer_aj",  "nmer_ak",  "nmer_al" ]
# Call the function to remove the specified columns from the DataFrame
remove_columns(df_proper_fin, columns_to_remove)
# Print the updated DataFrame
print(df_proper_fin)

# list of columns to remove
columns_to_remove = [ "nmer_ai",  "nmer_aj",  "nmer_ak",  "nmer_al" ]
# Call the function to remove the specified columns from the DataFrame
remove_columns(df_improper_fin, columns_to_remove)
# Print the updated DataFrame
print(df_improper_fin)

# list of columns to remove
columns_to_remove = [ "nmer_ai",  "nmer_aj",  "nmer_ak" ]
# Call the function to remove the specified columns from the DataFrame
remove_columns(df_angles_fin, columns_to_remove)
# Print the updated DataFrame
print(df_angles_fin)

# list of columns to remove
columns_to_remove = [ "nmer_ai",  "nmer_aj"]
# Call the function to remove the specified columns from the DataFrame
remove_columns(df_bonds_fin, columns_to_remove)
# Print the updated DataFrame
print(df_bonds_fin)

# list of columns to remove
columns_to_remove = ['nmer']
# Call the function to remove the specified columns from the DataFrame
remove_columns(df_atoms_fin, columns_to_remove)
# Define the desired column order
desired_columns = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']
# Reorder the columns of df_atoms using reindex
df_atoms_fin = df_atoms_fin.reindex(columns=desired_columns)
# Print the rearranged DataFrame
print(df_atoms_fin)

# Convert the specified columns to strings with three decimal places
df_atoms_fin[['charge','mass']] = df_atoms_fin[['charge','mass']].astype(float)
df_bonds_fin[['c0', 'c1']] = df_bonds_fin[['c0', 'c1']].astype(float)
df_angles_fin[['c0', 'c1']] = df_angles_fin[['c0', 'c1']].astype(float)
df_improper_fin[['c0', 'c1']] = df_improper_fin[['c0', 'c1']].astype(float)
df_proper_fin[['c0', 'c1', 'c2', 'c3', 'c4', 'c5']] = df_proper_fin[['c0', 'c1', 'c2', 'c3', 'c4', 'c5']].astype(float)


df_atoms_fin[['charge','mass']] = df_atoms_fin[['charge','mass']].applymap(lambda x: f"{x:.4f}")
df_bonds_fin[['c0']] = df_bonds_fin[['c0']].applymap(lambda x: f"{x:.4f}")
df_bonds_fin[['c1']] = df_bonds_fin[['c1']].applymap(lambda x: f"{x:.3f}")
df_angles_fin[['c0', 'c1']] = df_angles_fin[['c0', 'c1']].applymap(lambda x: f"{x:.3f}")
df_improper_fin[['c0', 'c1']] = df_improper_fin[['c0', 'c1']].applymap(lambda x: f"{x:.3f}")
df_proper_fin[['c0', 'c1', 'c2', 'c3', 'c4', 'c5']] = df_proper_fin[['c0', 'c1', 'c2', 'c3', 'c4', 'c5']].applymap(lambda x: f"{x:.3f}")

# Define the input and output file paths
itp_file_path = 'test.itp'
output_file_path = 'output.itp'

# Initialize an empty list to store the modified lines
modified_lines = []

# Flag to start clearing data after encountering the "[ atoms ]" section
clear_data = False

# Read the file line by line
with open(itp_file_path, 'r') as file:
    for line in file:
        # Check if the line is a header or comment
        if line.strip() == '' or line.strip().startswith(';') or line.strip().startswith('['):
            # Retain headers and comments
            modified_lines.append(line)

            # Check if the current line is the "[ atoms ]" section
            if '[ atoms ]' in line:
                clear_data = True
        else:
            # If clear_data is True, ignore data lines
            if clear_data:
                continue
            else:
                modified_lines.append(line)

# Remove lines containing the string "pair"
modified_lines = [line for line in modified_lines if 'pair' not in line]

# Write the modified content back to the file
with open(output_file_path, 'w') as file:
    file.writelines(modified_lines)

# Print a message indicating that the process is complete
print(f"Data has been cleared after '[ atoms ]', headers are retained, and lines containing 'pair' have been removed in {output_file_path}")

# Read the original file
with open('output.itp', 'r') as file:
    lines = file.readlines()

# Identify the line containing "[ bonds ]" and add a blank line before it
new_lines = []
for line in lines:
    if '[ bonds ]' in line:
        new_lines.append('\n')  # Add a blank line
        new_lines.append(line)  # Add the original line
    else:
        new_lines.append(line)  # Add the original line as is

# Write the modified content back to the file
with open('output.itp', 'w') as file:
    file.writelines(new_lines)

# Convert DataFrames to string format with correct spacing
atoms_str = format_atoms_df(df_atoms_fin)
bonds_str = format_bonds_str(df_bonds_fin)
angles_str = format_angles_str(df_angles_fin)
improper_str = format_improper_str(df_improper_fin)
proper_str = format_proper_str(df_proper_fin)

# Read the original file
with open('output.itp', 'r') as file:
    lines = file.readlines()

# Identify sections and insert DataFrame strings
new_lines = []
i = 0
while i < len(lines):
    line = lines[i]
    if '[ atoms ]' in line:
        print("Found [ atoms ] section")
        new_lines.append(line)
        i += 1
        # Skip header line
        new_lines.append(lines[i])
        i += 1
        # Insert formatted DataFrame
        new_lines.append(atoms_str + '\n')
    elif '[ bonds ]' in line:
        print("Found [ bonds ] section")
        new_lines.append(line)
        i += 1
        # Insert formatted DataFrame
        new_lines.append(bonds_str + '\n')
    elif '[ angles ]' in line:
        print("Found [ angles ] section")
        new_lines.append(line)
        i += 1
        # Skip header line
        new_lines.append(lines[i])
        i += 1
        # Insert formatted DataFrame
        new_lines.append(angles_str + '\n')
    elif 'IMPROPER' in line:
        print("Found IMPROPER section")
        new_lines.append(line)
        i += 1
        # Skip headers until a blank line or next section
        while i < len(lines) and lines[i].strip() != '' and not any(keyword in lines[i] for keyword in ['[', 'IMPROPER', 'PROPER']):
            new_lines.append(lines[i])
            i += 1
        # Insert formatted DataFrame
        new_lines.append(improper_str + '\n')
    elif 'PROPER' in line:
        print("Found PROPER section")
        new_lines.append(line)
        i += 1
        # Skip headers until a blank line or next section
        while i < len(lines) and lines[i].strip() != '' and not any(keyword in lines[i] for keyword in ['[', 'IMPROPER', 'PROPER']):
            new_lines.append(lines[i])
            i += 1
        # Insert formatted DataFrame
        new_lines.append(proper_str + '\n')
    else:
        new_lines.append(line)
    i += 1

# Write the modified content back to the file
with open('output.itp', 'w') as file:
    file.writelines(new_lines)
