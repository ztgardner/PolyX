import matplotlib
matplotlib.use('Agg') #Turns of GUI sup


import pandas as pd
import math
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import sys
# Redirect print statements to the system output
sys.stdout = sys.__stdout__


def get_old_paramters():
    try:
        index = 1
        directory_name = "run_1"
        run_check = False
        while os.path.exists("../"+directory_name): #checks what directory to look
            index += 1
            run_check = True
            directory_name = f"run_{index}"
            
            

        if run_check is False:
            print("No Old Run Directories Found Containing Parameters.csv")
            Parameters = [0,0,0,0,0,0]
            print("Using Paramters: " + str(Parameters))
        else:
            directory_name = "../"+ f"run_{index -2}"
            print("Looking for Parameters in: " + directory_name)
            P_path = os.path.join(os.getcwd(), directory_name + "/Parameters.csv")
            Parameters = []
            with open(P_path, "r") as file:
                
                for line in file:
                    line = line.strip()  # Remove leading and trailing whitespaces
                    if line:  # Check if the line is not empty
                        Parameters.append(float(line))
                print("Old Paramters Found: " + str(Parameters))
    except Exception as e:
        print(e)
        Parameters = [0,0,0,0,0,0]
        print("Failed to extract Old Parameters from Parameters.csv. Using Paramters: " + str(Parameters))
    return(Parameters)


def extract_potential_energy(log_file):
    print(log_file)
    with open(log_file, 'r') as file:
        content = file.read()
        match = re.search(r'Potential Energy\s+=\s+([0-9.+-e]+)', content)
        if match:
            potential_energy = match.group(1)
            print(f"Extracted potential energy from {log_file}: {potential_energy}")  # Debugging output
            return potential_energy
    print(f"Potential energy not found in {log_file}")  # Debugging output
    return None

def extract_atom_coordinates(gro_file, atom_number):
    with open(gro_file, 'r') as file:
        lines = file.readlines()
        for line in lines[2:]:  # Start from line 3 to skip header lines
            if line.strip():  # Check if the line is not empty
                fields = line.split()
                if len(fields) >= 6 and fields[2] == str(atom_number):
                    # Extract x, y, z coordinates
                    x = float(fields[3])
                    y = float(fields[4])
                    z = float(fields[5])
                    return [x, y, z]
    return None

def calculate_dihedral_angle(atom1_xyz, atom2_xyz, atom3_xyz, atom4_xyz):
    # Convert coordinates to numpy arrays
    atom1 = np.array(atom1_xyz, dtype=float)
    atom2 = np.array(atom2_xyz, dtype=float)
    atom3 = np.array(atom3_xyz, dtype=float)
    atom4 = np.array(atom4_xyz, dtype=float)

    # Vectors representing the bonds
    b1 = atom2 - atom1
    b2 = atom3 - atom2
    b3 = atom4 - atom3

    # Normal vectors to the planes defined by the three bonds
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    # Angle between the normal vectors (dihedral angle)
    angle_rad = np.arccos(np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)))

    # Convert angle to degrees
    angle_deg = np.degrees(angle_rad)

    return angle_deg

def process_dft_results(file_path):
    # Read the DFT results file into a DataFrame.
    df = pd.read_csv(file_path, sep=' ')

    # Calculate relative DFT energy in kJ/mol.
    df['dft_rel_E(kJ/mol)'] = (df.iloc[:, 0] - df.iloc[:, 0].min()) * 2625.5

    # Create a new column called "angle(degrees)".
    num_rows = len(df)
    angles = [i * 10 for i in range(num_rows)]
    df['angle(degrees)'] = angles

    # Load MD results and calculate relative MD energy.
    md_results = pd.read_csv('md_results.txt', sep="	")
    df['md_results'] = md_results.iloc[:, 1]
    df['md_rel_E(kJ/mol)'] = df['md_results'] - df['md_results'].min()

    # Calculate the difference between DFT and MD energies.
    df['dft-md(kJ/mol)'] = df['dft_rel_E(kJ/mol)'] - df['md_rel_E(kJ/mol)']

    # Calculate psi in radians.
    df['psi'] = (df['angle(degrees)'] + 180) / 180 * math.pi

    print(df)

    return df

def calculate_fit(params, psi):
    c0, c1, c2, c3, c4, c5 = params
    fit = c0 + c1 * math.cos(psi)**1 + c2 * math.cos(psi)**2 + c3 * math.cos(psi)**3 + c4 * math.cos(psi)**4 + c5 * math.cos(psi)**5
    return fit

def optimize_parameters(df):
    initial_guess = [0, 0, 0, 0, 0, 0]

    def objective_function(params, psi_values, dft_md):
        fit_values = [calculate_fit(params, psi) for psi in psi_values]
        sq_err_values = (dft_md - fit_values) ** 2
        return sum(sq_err_values)

    psi_values = df['psi']  # Convert Series to numpy array         CHECKING TO SEE IF the first element18 fit better
    dft_md = df['dft-md(kJ/mol)']  # Convert Series to numpy array

    result = minimize(objective_function, initial_guess, args=(psi_values, dft_md), method='L-BFGS-B')
    return result

def get_dihedral_index():
    itp = "000.itp" #Checking the first one
    with open(itp) as fin:
        lines = fin.readlines()
    
    for line_num,line in enumerate(lines, start=1):
        
        match = re.search(r'\[ dihedral_restraints \]', line)
        if match:
            index = line_num +1

    
    index  = lines[index].split()
    index = index[0:4]
    index = [int(i) for i in index]
    return index
    





def main():
    folder_path = os.getcwd()
    [id_1, id_2, id_3, id_4] = get_dihedral_index()
    with open('md_results.txt', 'w') as output_file:
        output_file.write("file\tE(kJ/mol)\t\tDihedral\n")
        log_files = [file for file in os.listdir(os.getcwd()) if file.endswith('.log')]
        log_files = sorted(log_files, key=lambda x: int(re.search(r'\d+', x).group()))


        for file_name in log_files:
            if file_name.endswith('.log'):
                log_file = os.path.join(folder_path, file_name)
                gro_file = os.path.join(folder_path, file_name.replace('.log', '.gro'))

                file_number = re.search(r'em_(\d+)\.log', file_name).group(1)
                potential_energy = extract_potential_energy(log_file)
                
                # Extract atom coordinates from the .gro file for the dihedral calculation
                atom1_xyz = extract_atom_coordinates(gro_file, id_1)
                atom2_xyz = extract_atom_coordinates(gro_file, id_2)
                atom3_xyz = extract_atom_coordinates(gro_file, id_3)
                atom4_xyz = extract_atom_coordinates(gro_file, id_4)

                # Calculate dihedral angle
                dihedral_angle = calculate_dihedral_angle(atom1_xyz, atom2_xyz, atom3_xyz, atom4_xyz)

                if potential_energy is not None and dihedral_angle is not None:
                    output_file.write(f"{file_number}\t{potential_energy}\t{dihedral_angle}\n")

    print("Output file 'md_results.txt' created with file number, potential energy, and dihedral angle.")
    # Call the function with the provided file path.
    df = process_dft_results('../DFT_Energy.txt')
 
    # Optimize the parameters.
    result = optimize_parameters(df)

    #Writing to File
    old_params = get_old_paramters() 
    params_list = result.x.tolist()
    params_list = [a + b for a, b in zip(params_list, old_params)]
    df_params = pd.DataFrame(params_list)
    df_params.to_csv("Parameters.csv", index=False,header=False)

    # Print the optimized parameters and final sum of squared errors.
    print("Final Optimized Parameters:", result.x)
    print("Final Sum of Squared Errors:", result.fun)

    # Apply the optimized parameters to calculate the fit and squared error.
    df['fit'] = [calculate_fit(result.x, psi) for psi in df['psi']]
    df['sq_error'] = (df['dft-md(kJ/mol)'] - df['fit']) ** 2

    # Sum the squared errors.
    total_sq_error = df['sq_error'].sum()
    print(f"Total squared error: {total_sq_error}")

    # Save the updated DataFrame to a new text file.
    df.to_csv('Summary.txt', sep='\t', index=False)

    # Plot the first graph: Dihedral Angle vs Relative Energy (kJ/mol)

    plt.figure(figsize=(10, 6))
    plt.plot(df['angle(degrees)'], df['dft_rel_E(kJ/mol)'], color='blue', label='DFT Relative Energy')
    plt.plot(df['angle(degrees)'], df['md_rel_E(kJ/mol)'], color='orange', label='MD Relative Energy')
    plt.xlabel('Dihedral Angle')
    plt.ylabel('Relative Energy (kJ/mol)')
    plt.legend()
    plt.savefig('DFT_and_MD.png')  # Save the figure as 'dihedral_vs_energy_1.png'

    # Plot the second graph: Dihedral Angle vs Relative Energy (kJ/mol)
    plt.figure(figsize=(10, 6))
    plt.plot(df['angle(degrees)'], df['dft-md(kJ/mol)'], color='blue', label='DFT-MD')
    plt.plot(df['angle(degrees)'], df['fit'], color='orange', label='Fit')
    plt.xlabel('Dihedral Angle')
    plt.ylabel('Relative Energy (kJ/mol)')
    plt.legend()
    plt.savefig('Fit.png')  # Save the figure as 'dihedral_vs_energy_2.png'

if __name__ == "__main__":
    main()
