import pandas as pd
import numpy as np
from scipy.spatial.transform import Rotation as R
import re



#####################################################################################################################
# Load and prepare the DataFrame
file_path = 'input.gro'
dihedral_file_path = 'dihedrals.txt'
start_atoms = [2,1,5]  # Going into the polymer H-C-S, this is the start of each polymer unit.
end_atoms = [5,6,7]    # Coming out the polymer S-C-H, this is where the polymer will propagate from.
n = 5
bond_length = 0.1380
#####################################################################################################################


def parse_gro(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract the number of atoms from the second line
    num_atoms = int(lines[1].strip())

    # Initialize an empty list to store atom data
    atom_data = []

    # Loop through the lines that contain atom information
    for i in range(2, 2 + num_atoms):
        line = lines[i].strip()

        # Use regular expression to split the line by whitespace
        parts = re.split(r'\s+', line)

        # Extract relevant parts
        res = parts[0]
        atom_name = parts[1]
        atom_num = int(parts[2])
        x = float(parts[3])
        y = float(parts[4])
        z = float(parts[5])

        # Append the extracted data to the list
        atom_data.append([res, atom_name, atom_num, x, y, z])

    # Create a DataFrame with the extracted data
    df = pd.DataFrame(atom_data, columns=['res', 'atom_name', 'atom_num', 'x', 'y', 'z'])

    return df

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

def extract_atomic_symbols(df):
    atom_col = []
    atom_namnum_col = []

    for atom_name in df['atom_name']:
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
    df.drop(columns=['atom_name'], inplace=True)

    # Rearrange the columns
    df = df[['res', 'atom', 'atom_namnum', 'atom_num', 'x', 'y', 'z']]

    return df

def monomer(df, start_atoms, end_atoms):
    # Create 'start' column based on specified start_atoms list
    df['start'] = 0  # Initialize all values to 0
    for i, atom_num in enumerate(start_atoms, 1):
        df.loc[df['atom_num'] == atom_num, 'start'] = i

    # Create 'end' column based on specified end_atoms list
    df['end'] = 0  # Initialize all values to 0
    for i, atom_num in enumerate(end_atoms, 1):
        df.loc[df['atom_num'] == atom_num, 'end'] = i

    return df

def prep_nmer(df, n):
    # Initialize an empty list to store DataFrame copies
    df_list = []

    # Loop through n times to create the copies with different nmer values
    for i in range(1, n + 1):
        df_copy = df.copy()  # Create a copy of the DataFrame
        df_copy['nmer'] = i  # Set the nmer column to the current value of i
        df_list.append(df_copy)  # Append the copy to the list

    # Concatenate all copies into a single DataFrame
    result_df = pd.concat(df_list, ignore_index=True)

    return result_df

def translate_coordinates(df):
    # Find the coordinates of nmer=1, end=2
    nmer_1_end_2_coords = df[(df['nmer'] == 1) & (df['end'] == 2)][['x', 'y', 'z']].values[0]

    # Calculate the translation values
    dx = df.loc[(df['nmer'] == 2) & (df['start'] == 1), 'x'].values[0] - nmer_1_end_2_coords[0]
    dy = df.loc[(df['nmer'] == 2) & (df['start'] == 1), 'y'].values[0] - nmer_1_end_2_coords[1]
    dz = df.loc[(df['nmer'] == 2) & (df['start'] == 1), 'z'].values[0] - nmer_1_end_2_coords[2]

    # Translate coordinates for all nmer=2 rows
    df.loc[df['nmer'] == 2, ['x', 'y', 'z']] -= [dx, dy, dz]

    return df

def rotate_points(df):
    # Extract coordinates for nmer=1, end=2, nmer=1, end=3, and nmer=2, start=2
    nmer_1_end_2_coords = df[(df['nmer'] == 1) & (df['end'] == 2)][['x', 'y', 'z']].values
    nmer_1_end_3_coords = df[(df['nmer'] == 1) & (df['end'] == 3)][['x', 'y', 'z']].values
    nmer_2_start_2_coords = df[(df['nmer'] == 2) & (df['start'] == 2)][['x', 'y', 'z']].values

    # Calculate the midpoint between nmer=1, end=3 and nmer=2, start=2
    midpoint = (nmer_1_end_3_coords + nmer_2_start_2_coords) / 2

    # Calculate the axis of rotation from nmer=1, end=2 to the midpoint
    axis = midpoint - nmer_1_end_2_coords

    # Rotate nmer=2 points 180 degrees around the calculated axis
    r = R.from_rotvec(axis * np.pi/2)
    nmer_2_points = df[df['nmer'] == 2][['x', 'y', 'z']].values
    rotated_points = r.apply(nmer_2_points - nmer_1_end_2_coords) + nmer_1_end_2_coords

    # Update the DataFrame with rotated points
    df.loc[df['nmer'] == 2, ['x', 'y', 'z']] = rotated_points

    return df

def align_monomer_orientation(df):
    # Extract coordinates for nmer=1, end=3, nmer=1, end=2, and nmer=2, start=2
    nmer_1_end_3_coords = df[(df['nmer'] == 1) & (df['end'] == 3)][['x', 'y', 'z']].values
    nmer_1_end_2_coords = df[(df['nmer'] == 1) & (df['end'] == 2)][['x', 'y', 'z']].values
    nmer_2_start_2_coords = df[(df['nmer'] == 2) & (df['start'] == 2)][['x', 'y', 'z']].values

    # Ensure that the shapes are correct
    if nmer_1_end_3_coords.shape != (1, 3) or nmer_1_end_2_coords.shape != (1, 3) or nmer_2_start_2_coords.shape != (1, 3):
        raise ValueError("Vector dimensions are not as expected.")

    # Calculate vectors between the points
    vec_1_3 = nmer_1_end_3_coords - nmer_1_end_2_coords
    vec_1_2 = nmer_2_start_2_coords - nmer_1_end_2_coords

    # Reshape the vectors to ensure correct dimensions for the dot product
    vec_1_3 = vec_1_3.reshape(3,)  # Reshape to a 1D array
    vec_1_2 = vec_1_2.reshape(3,)  # Reshape to a 1D array

    # Calculate the rotation angle as the angle between vec_1_3 and vec_1_2
    cos_angle = np.dot(vec_1_3, vec_1_2) / (np.linalg.norm(vec_1_3) * np.linalg.norm(vec_1_2))
    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))

    # Calculate the normal vector to the plane defined by vec_1_3 and vec_1_2
    normal_vec = np.cross(vec_1_3, vec_1_2)
    norm = np.linalg.norm(normal_vec)
    if norm != 0:
        normal_vec /= norm  # Normalize the normal vector

    # Rotate nmer=2 points around the normal vector at nmer=1, end=2 by the calculated angle
    r = R.from_rotvec(normal_vec * -angle)
    nmer_2_points = df[df['nmer'] == 2][['x', 'y', 'z']].values
    nmer_1_end_2_coords_repeated = np.repeat(nmer_1_end_2_coords, len(nmer_2_points), axis=0)
    rotated_points = r.apply(nmer_2_points - nmer_1_end_2_coords_repeated) + nmer_1_end_2_coords_repeated

    # Update the DataFrame with rotated points for nmer=2
    df.loc[df['nmer'] == 2, ['x', 'y', 'z']] = rotated_points

    return df

def adjust_bond_length(df, target_distance):
    # Extract coordinates for nmer=1, end=2 and nmer=2, start=2
    nmer_1_end_2_coords = df[(df['nmer'] == 1) & (df['end'] == 2)][['x', 'y', 'z']].values
    nmer_2_start_2_coords = df[(df['nmer'] == 2) & (df['start'] == 2)][['x', 'y', 'z']].values

    # Ensure that the shapes are correct
    if nmer_1_end_2_coords.shape != (1, 3) or nmer_2_start_2_coords.shape != (1, 3):
        raise ValueError("Vector dimensions are not as expected.")

    # Calculate the current distance between the two points
    current_distance = np.linalg.norm(nmer_2_start_2_coords - nmer_1_end_2_coords)

    # Calculate the direction vector from nmer=1, end=2 to nmer=2, start=2
    direction_vector = (nmer_2_start_2_coords - nmer_1_end_2_coords).reshape(3,)

    # Calculate the required translation vector to achieve the target distance
    scale_factor = (target_distance - current_distance) / current_distance
    translation_vector = scale_factor * direction_vector

    # Translate all points of nmer=2
    df.loc[df['nmer'] == 2, ['x', 'y', 'z']] += translation_vector

    return df

def calculate_dihedral_angle(A, B, C, D):
    """Calculate the dihedral angle between four points A, B, C, and D."""
    AB = B - A
    BC = C - B
    CD = D - C

    n1 = np.cross(AB, BC)
    n2 = np.cross(BC, CD)

    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)

    m1 = np.cross(n1, BC / np.linalg.norm(BC))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return np.degrees(np.arctan2(y, x))

def rotate_to_target_dihedral(df, target_dihedral_angle):
    # Extract coordinates for the specified points
    nmer_1_end_1_coords = df[(df['nmer'] == 1) & (df['end'] == 1)][['x', 'y', 'z']].values
    nmer_1_end_2_coords = df[(df['nmer'] == 1) & (df['end'] == 2)][['x', 'y', 'z']].values
    nmer_2_start_2_coords = df[(df['nmer'] == 2) & (df['start'] == 2)][['x', 'y', 'z']].values
    nmer_2_start_3_coords = df[(df['nmer'] == 2) & (df['start'] == 3)][['x', 'y', 'z']].values

    # Ensure that the shapes are correct
    if (nmer_1_end_1_coords.shape != (1, 3) or nmer_1_end_2_coords.shape != (1, 3) or
        nmer_2_start_2_coords.shape != (1, 3) or nmer_2_start_3_coords.shape != (1, 3)):
        raise ValueError("Vector dimensions are not as expected.")

    # Reshape the points to 1D arrays
    A = nmer_1_end_1_coords.reshape(3,)
    B = nmer_1_end_2_coords.reshape(3,)
    C = nmer_2_start_2_coords.reshape(3,)
    D = nmer_2_start_3_coords.reshape(3,)

    # Calculate the current dihedral angle
    current_dihedral_angle = calculate_dihedral_angle(A, B, C, D)

    # Calculate the required rotation angle to achieve the target dihedral angle
    rotation_angle = target_dihedral_angle - current_dihedral_angle

    # Normalize rotation angle to be within -180 to 180 degrees
    rotation_angle = (rotation_angle + 180) % 360 - 180

    # Calculate the rotation vector
    rotation_vector = (C - B) / np.linalg.norm(C - B)
    r = R.from_rotvec(np.radians(-rotation_angle) * rotation_vector)

    # Rotate nmer=2 points around the axis defined by B and C
    nmer_2_points = df[df['nmer'] == 2][['x', 'y', 'z']].values
    B_coords_repeated = np.repeat(B.reshape(1, 3), len(nmer_2_points), axis=0)
    rotated_points = r.apply(nmer_2_points - B_coords_repeated) + B_coords_repeated

    # Update the DataFrame with rotated points for nmer=2
    df.loc[df['nmer'] == 2, ['x', 'y', 'z']] = rotated_points

    # Recalculate the dihedral angle after rotation
    new_A = A
    new_B = B
    new_C = df[(df['nmer'] == 2) & (df['start'] == 2)][['x', 'y', 'z']].values.reshape(3,)
    new_D = df[(df['nmer'] == 2) & (df['start'] == 3)][['x', 'y', 'z']].values.reshape(3,)

    new_dihedral_angle = calculate_dihedral_angle(new_A, new_B, new_C, new_D)

    return df

def remove_specific_rows(df):
    # Remove the row where nmer=1 and start=3
    condition_1 = (df['nmer'] == 1) & (df['end'] == 3)
    df = df[~condition_1]

    # Remove the row where nmer=2 and end=2
    condition_2 = (df['nmer'] == 2) & (df['start'] == 1)
    df = df[~condition_2]

    return df

def subtract_one_from_nmer(df):
    # Subtract 1 from the 'nmer' column
    df['nmer'] = df['nmer'] - 1

    return df

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

def read_dihedral_angles(file_path):
    with open(file_path, 'r') as file:
        angles = [float(line.strip()) for line in file.readlines()]
    return angles

def perform_polymerization(df, n, dihedral_angles):
    num_angles = len(dihedral_angles)

    # Run the process n times
    for i in range(n-1):
        dihedral_angle = dihedral_angles[i % num_angles]  # Select angle in order

        # Perform operations
        df = translate_coordinates(df)
        df = rotate_points(df)
        df = align_monomer_orientation(df)
        df = adjust_bond_length(df, bond_length)
        df = rotate_to_target_dihedral(df, dihedral_angle)  # Use the selected angle
        df = remove_specific_rows(df)
        df = subtract_one_from_nmer(df)

    return df

def process_dataframe(df):
    # Concatenate "atom" and "atom_namnum" into a new "atom_names" column
    df['atom_names'] = df['atom'] + df['atom_namnum']

    # Remove the columns "nmer", "end", and "start"
    df = df.drop(columns=['nmer', 'end', 'start', 'atom', 'atom_namnum'])

    return df

def write_gro_file(df, input_file_path, output_file_path):
    # Read the first line from the input.gro file
    with open(input_file_path, 'r') as file:
        first_line = file.readline().strip()

    # Calculate the total number of atoms
    total_atoms = df['atom_num'].max()

    # Check and convert x, y, z columns to numeric if possible
    numeric_columns = ['x', 'y', 'z']
    for col in numeric_columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except ValueError:
            return

    # Prepare the DataFrame content for writing
    df['x'] = df['x'].map('{:.3f}'.format)
    df['y'] = df['y'].map('{:.3f}'.format)
    df['z'] = df['z'].map('{:.3f}'.format)

    # Write to the output.gro file
    with open(output_file_path, 'w') as file:
        file.write(f"{first_line}\n")
        file.write(f"  {total_atoms}\n")

        for _, row in df.iterrows():
            res = row['res']
            atom_names = row['atom_names']
            atom_num = row['atom_num']
            x = row['x']
            y = row['y']
            z = row['z']

            file.write(f"{res:>8}{atom_names:>7}{atom_num:>5}{x:>8}{y:>8}{z:>8}\n")
        file.write(f"   1.00000   1.00000   1.00000\n")

# Read dihedral angles from file
dihedral_angles = read_dihedral_angles(dihedral_file_path)

# Load the DataFrame and prepare the system
df = parse_gro(file_path)
df = extract_atomic_symbols(df)
df = monomer(df, start_atoms, end_atoms)
df = prep_nmer(df, n)

# Run the polymerization process using the dihedral angles from the file
df = perform_polymerization(df, n, dihedral_angles)

# Final update after the loop
df = update_dataframe_with_hex_style(df)

# Process and write the final DataFrame to a GRO file
processed_df = process_dataframe(df)
write_gro_file(processed_df, file_path, 'output.gro')