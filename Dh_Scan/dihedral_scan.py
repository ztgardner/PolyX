import re
import os
import pandas as pd
import sys
import platform


# Redirect print statements to the system output
sys.stdout = sys.__stdout__

class Dh_Info: #This class builds all the information needed to create the dihedral scan either by scrapping or user prompting
    def __init__(self, backbone_indexes=None, scan_index=None, charges=None , DFT_Energy=None, MD_Energy=None, Parameters=None):
        
        if backbone_indexes or scan_index is None:
            try:
                print("\n----------------------------------------------------------------\n")
                print("            Found Dihedral Scan Info From Scan.log\n")
                print("----------------------------------------------------------------\n")
                [backbone_indexes,scan_index] = self.get_indexes()
            except Exception as e:
                print(e)
                print("Failed to get backbone indexes and scan index From log file. Please provide them.")
                backbone_indexes_input = input("Enter backbone indexes (comma-separated, each sublist in square brackets): ")
                scan_index_input = input("Enter scan index (comma-separated): ")

                # Convert input strings to appropriate data types
                backbone_indexes = [list(map(int, sublist.strip().split(','))) for sublist in backbone_indexes_input.strip('[]').split('],[')]
                scan_index = list(map(int, scan_index_input.split(',')))   

        if charges is None:
            #Try to read charges  from Charges.txt
            try:
                print("\nAttempting to read charges from 'charges.txt' file")
                with open("charges.txt", 'r') as f:
                    charges = f.readlines()
                    # Remove the first line (header) and any leading or trailing whitespace
                    charges = [charge.strip() for charge in charges[1:]]
                print("Charges Found in 'charges.txt'")
            except:
                #Trys to extract them from .log
                print("Failed to read charges from 'charges.txt' attempting to extract from Opt.log File...")
                try:
                    charges = self.read_charges_from_log()
                    print("Charges Found in Optimized Log File")
                    print("Writing charges to 'charges.txt'")
                    charges.to_csv("charges.txt", index=False, header=True, sep=" ")
                    
                    #Once Written to charges.txt reads from file for consitency 
                    with open("charges.txt", 'r') as f:
                        charges = f.readlines()
                        charges = [charge.strip() for charge in charges[1:]]
                except Exception as e:
                    print("Failed to extract charges from Opt.log File. Please provide charges manually.")
                    print(e)
                    charges_input = input("Enter charges (comma-separated): ")
                    charges = list(map(float, charges_input.split(',')))
        
            if DFT_Energy is None:
                try:
                    print("\nAttempting to read DFT Energy from 'DFT_Energy.txt")
                    with open("DFT_Energy.txt", 'r') as f:
                        DFT_Energy = f.readlines()
                        # Remove the first line (header) and any leading or trailing whitespace
                        DFT_Energy = [DFT_Energy.strip() for DFT_Energy in DFT_Energy[1:]]
                    print("Energies Found in 'DFT_Energy.txt'")
                except:
                    print("Failed to read DFT Energy from 'DFT_Energy.txt' attempting to extract from Scan.log File...")
                    try:
                        E = self.read_energies_from_log()
                        print("DFT_Energy Found in Scan File")
                        print("Writing charges to 'DFT_Energy.txt'")
                        E.to_csv("DFT_Energy.txt", index=False, header=True, sep=" ")
                        
                        #Once Written to charges.txt reads from file for consitency 
                        with open("DFT_Energy.txt", 'r') as f:
                            DFT_Energy = f.readlines()
                            DFT_Energy = [DFT_Energy.strip() for DFT_Energy in DFT_Energy[1:]]
                    except Exception as e:
                        print("Failed to extract DFT_Energy from Scan.log File. Please provide Energy manually.")
                        print(e)
                        charges_input = input("Enter charges (comma-separated): ")
                        charges = list(map(float, charges_input.split(',')))
        if Parameters is None:
            Parameters = []
            print('----------------------------------------------------------------')
            print("\nAttempting to read Parameters from 'Parameters.csv'")
            try:
                index = 1
                directory_name = "run_1"
                run_check = False
                while os.path.exists(directory_name): #checks what directory to look
                    index += 1
                    run_check = True
                    directory_name = f"run_{index}"
                    
                    

                if run_check is False:
                    print("No Run Directories Found Containing Parameters.csv")
                    Parameters = [0,0,0,0,0,0]
                    print("Using Paramters: " + str(Parameters))
                else:
                    directory_name = f"run_{index -1}"
                    print("Looking for Parameters in: " + directory_name)
                    P_path = os.path.join(os.getcwd(), directory_name + "/Parameters.csv")
                    with open(P_path, "r") as file: 
                        for line in file:
                            line = line.strip()  # Remove leading and trailing whitespaces
                            if line:  # Check if the line is not empty
                                Parameters.append(float(line))
                        print("Using Paramters: " + str(Parameters))
            except Exception as e:
                print(e)
                Parameters = [0,0,0,0,0,0]
                print("Failed to extract Parameters from Parameters.csv. Using Paramters: " + str(Parameters))


        temp_dh = Scan.check_dihedral(self,dihedral_atoms=scan_index,base_itp="base.itp") #This makes sure that the scan paramters are in the correct order of the itp
        if temp_dh == []:
            temp_dh = Scan.check_dihedral(self,dihedral_atoms=scan_index.reverse(),base_itp="base.itp")

            if temp_dh != []:
                scan_index = temp_dh
                temp = temp[0] #not sure why this is need but it keeps the output clean :---/

        
        print('----------------------------------------------------------------')
        print("       Attempting to get Geometries from " + "Scan.log")
        print('----------------------------------------------------------------')
        self.extract_geometries()

        
        self.backbone_indexes = backbone_indexes
        self.scan_indexes = scan_index
        self.charges = charges
        self.DFT_Energy = DFT_Energy
        self.MD_Energy = MD_Energy
        self.Parameters = Parameters
    ################################
    "Reads Scan.log file and extracts the geometry of each scan, then writes as .gjf files"
    ################################
    def extract_geometries(self):
        log_file = "Scan.log"  # The name of any Gaussian Log File Scan Output
        def atomic_number_to_symbol(atomic_number):
            element_symbols = {
                1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B',
                6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
                11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
                16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
                21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn',
                26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
                31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br',
                36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
                41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh',
                46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
                51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs',
                56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
                61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb',
                66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
                71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re',
                76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
                81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
                86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
                91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am',
                96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
                101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db',
                106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
                111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc',
                116: 'Lv', 117: 'Ts', 118: 'Og'
            }
            return element_symbols.get(atomic_number, "Unknown Symbol")



        def write_geometries(file_out, start,stop): #writes gjf
            
            if os.path.isfile(f'GJF/{file_out}.gjf'):
                print("Warning: " + f'GJF/{file_out}.gjf' + " already exists. Skipping")
                return
            
            with open(f'GJF/{file_out}.gjf', "w") as fout:
                fout.write("#File auto generated by dihedral_scan.py -Lok 9/7/2024\n")
                with open(log_file, "r") as fin:
                    lines = fin.readlines()
                    for line_number, line in enumerate(lines, start=0):
                        if line_number >= start and line_number <= stop:
                            line = line.split()
                            #print(line)
                            symbol = atomic_number_to_symbol(int(line[1]))
                            x = float(line[3])
                            y = float(line[4])
                            z = float(line[5])
                            line = f"{symbol:2}  {x:15.8f}  {y:15.8f}  {z:15.8f}\n"                
                            fout.write(line)
                            
                                

        def transform_number(n): #Truns integers into the correct 000.gro format
            if n % 10 == 0:
                return f'{n:03}'
            else:
                return f'{n}'
            
        
        
        new_opt_ptr = " Optimization completed."#This is where geometries start 
        start_new_geom_ptr = "                         Standard orientation:                         " #This is where cordinate data starts
        end_new_geom_ptr = " ---------------------------------------------------------------------"
        new_opt = False
        in_geom = False

        scan_rate = 10 #Assuming the scan rate is 10 degrees 
        angle = 0
        print("Searching for geometries in " + log_file +  "...")
        print("Found: \n")
        print("Angle|Start Line|End Line")
        with open(log_file, 'r') as input_file:
            lines = input_file.readlines()
        
        
        for line_number, line in enumerate(lines, start=0):
            match_new_opt = re.search(new_opt_ptr, line)
            match_new_geom = re.search(start_new_geom_ptr, line)
            match_end_geom = re.search(end_new_geom_ptr, line)
            if match_new_opt: #Double checks we have seen a new optimzation
                new_opt = True
            
            if new_opt and match_new_geom: #Finds where the geometry data from the new optimization section starts       
                start_geom = line_number + 5 #+4 to remove heading
                in_geom = True
            
            if new_opt and in_geom and match_end_geom and line_number > start_geom + 5: #Finds where the geometry data from the new optimization section ends (+6 to skip heading)
                end_geom = line_number
                new_opt = False
                in_geom = False
                file_name = transform_number(angle)
                print(file_name,start_geom,end_geom) 
                write_geometries(file_name,start_geom,end_geom - 1)
                angle += 10
           
        
        

    ################################
    "Finds the Dihedral Angel in Scan.log"
    ################################
    def find_dihedral(self,atom_list):
        log_file = "Scan.log"  # The name of any Gaussian Log File Scan line
        [a1, a2, a3, a4] = atom_list
        pattern = r'D\((\d+),(\d+),(\d+),(\d+)\)' #This is the Line Pattern
        pattern_atom = f"D({a1},{a2},{a3},{a4})" #this is the atom index

        with open(log_file, 'r') as input_file:
            lines = input_file.readlines()
        
        for line_number, line in enumerate(lines, start=1):
            match = re.search(pattern, line)
            if match:
                temp = line.split()
                if temp[2] == pattern_atom:
                    return(float(temp[3]))
                
        
    ################################
    "opens the Scan.log file to get required info: Energy"
    ################################
    def read_energies_from_log(self):
        file_path = "Scan.log"
        energies = []
        energy = False
        # Regular expression pattern to find energy in log file (SCF Done energy for example)
        energy_pattern = re.compile(r'SCF Done:')
        opt_pattern = re.compile(r" Optimization completed." )     
        with open(file_path, 'r') as file:
            for line in file:
                # Search for energy in each line using the pattern
                match_e = energy_pattern.search(line)
                match_grad = opt_pattern.search(line)
                if match_e:
                    energy = line.split()[4]
                elif match_grad:
                    energies.append(float(energy))
                    energy = False
        print("Found Energy in Scan.log")
        df = pd.DataFrame(energies)
        print(df)
        return df
################################
    "opens the opt.log file to get required info: charges"
    ################################
    def read_charges_from_log(self):
        charges = []
        log_file_opt = "Opt.log" #The name of any Gaussian Log File Opt Output
        pattern_start = r"Hirshfeld charges, spin densities, dipoles, and CM5 charges using IRadAn=\s*5"
        pattern_end = r"Hirshfeld charges with hydrogens summed into heavy atoms:"
        print("Looking for Pattern in " + log_file_opt + ":")
        print("Pattern Start: " + pattern_start)
        print("Pattern End: " + pattern_end)
        print('----------------------------------------------------------------')
        with open(log_file_opt, 'r') as input_file:
            lines = input_file.readlines()

        Read = False
        for line_number, line in enumerate(lines, start=1):
            match_start = re.search(pattern_start, line)
            match_end = re.search(pattern_end, line)
            
            if match_start:
                Read = True
                continue 
            if Read and match_end:
                break
            if Read:
                temp = line.split()
                #print(f"Line {line_number}: {temp}")
                charges.append([temp[2],temp[1]])
        df = pd.DataFrame(charges, columns=["charges", "atom"])
        df = df.drop(df.index[0]) #Drops first line which labels the columns
        df = df.drop(df.index[-1]) #Drops last line which is the totals
        print(df)
        print('----------------------------------------------------------------')
        charges = df.drop(columns=["atom"])
        return charges
    
    ################################
    "opens the scan.log file to get required info: scanned cords, backbone indexes"
    ################################
    def get_indexes(self):
        
        log_file_scan = "Scan.log" #The name of any Gaussian Log File Scan Output
        print("Atempting to read indexes from "+ log_file_scan )
        freeze_cords_data = []
        scanned_cords = []
        pattern = r"D\s+\d+\s+\d+\s+\d+\s+\d+\s+[FS]\s*(?:\s+\d+\s+\d+\s*)?"
        #Read .log
        print("Looking for Pattern in " + log_file_scan + ":")
        print(pattern)
        print('----------------------------------------------------------------')
        with open(log_file_scan, 'r') as input_file:
            lines = input_file.readlines()
        
        
        for line_number, line in enumerate(lines, start=1):
            match = re.search(pattern, line)
            if match:
                temp = line.split()  
                freeze_cords = [temp[1], temp[2], temp[3], temp[4]]
                freeze_cords.append(self.find_dihedral(freeze_cords)) #Adds the Dihedral
                freeze_cords.append(len(temp)>6) #Checks if we are scanning the cordinate
                freeze_cords_data.append(freeze_cords)
                #print(f"Line {line_number}: {temp}")
                if len(temp)>6:
                    scanned_cords = [temp[1], temp[2], temp[3], temp[4]]
        
        # Create DataFrame from freeze_cords_data
        df = pd.DataFrame(freeze_cords_data, columns=['Atom1', 'Atom2', 'Atom3', 'Atom4', "Dihedral Angel", "Scan_Cords?"])
        print("BackBone Index To Freeze Extracted From Log File")
        print(df)
        print("\n")
        print("BackBone Index To Scan Extracted From Log File")
        print(scanned_cords)
        print('----------------------------------------------------------------')
        return [df,scanned_cords]
    
    def display_info(self):
        print("Backbone Indexes To Freeze\n")
        for sublist in self.backbone_indexes:
            print(sublist)
        
        print("Indexes to Scan")
        print(self.scan_index)


############################################################################################
"This class deals with the scan.log files to building the md simulations"
"At this Point Scan.log, Opt.log, base.itp, and 000-360.gjf files must be present in the same directory"
############################################################################################
class Scan:
    def __init__(self,dh_info_class):
        self.log_file_scan = "Scan.itp"
        self.itp_file = "base.itp"
        self.charges_file = "charges.txt"
        self.dh_info =dh_info_class

    
    #Run this to Set up ITP Files (ie 000.itp-360.itp)
    def make_itp_files(self):
            input_file = self.itp_file
            self.build_gro()
            
            print("\n----------------------------------------------------------------\n")
            print("                 EDITING " +  input_file + "\n")
            print("----------------------------------------------------------------\n")
            
            self.clean_itp(input_file) #Cleans Base.itp file (unless changed in dh_INfo)
            self.add_charges()
            self.find_dihedral_to_scan_and_update()
            self.fix_dihedral_spaces()
            self.write_dihedral_restraints()
            

            print("\n----------------------------------------------------------------\n")
            print("                 Finished EDITING " +  input_file)
            print("  Using " +  input_file + " as a template, building 000-360.itp Files\n")
            print("----------------------------------------------------------------\n")

            self.propagate_itp()
    ################################
    "Helper Functions"
    ################################
    def create_directory(self,directory_path):
        if os.path.exists(directory_path):
            print("Directory already exists: " + directory_path)
            print("Clearing Directory")
            for file in os.listdir(directory_path):
                file_path = os.path.join(directory_path, file)
                if os.path.isfile(file_path):
                    os.remove(file_path)
        else:
            print("Creating directory " + directory_path)
            os.makedirs(directory_path)
    def clean_itp(self,file_path):
        print("Cleaning ITP file: " + file_path)
        with open(file_path, 'r') as file:
            lines = file.readlines()
        # Remove lines that are only whitespace
        lines = [line for line in lines if line.strip() != '']



        with open(file_path, 'w') as file:
            file.writelines(lines)
    def hex_style_atom_type(self,index):
        """Generate hex-style atom type based on index."""
        hex_digits = '0123456789ABCDEFGHIJKMNOPQRSTUVWXYZ'
        base = len(hex_digits)
        result = ''
        index -= 1
        while index > 0:
            index, remainder = divmod(index, base)
            result = hex_digits[remainder] + result
        return result.zfill(2)  # Pad with zeros to ensure three characters
    
    def read_charges_from_file(self,filename):
        with open(filename, 'r') as f:
            charges = f.readlines()
        # Remove the first line (header) and any leading or trailing whitespace
        charges = [charge.strip() for charge in charges[1:]]
        return charges
    
    def update_charge_column(self,input_file, charges):
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Find the index range for the section to update
        start_idx = lines.index("[ atoms ]\n")
        end_idx = lines.index("[ bonds ]\n")

        # Update the charge column
        for idx, charge in zip(range(start_idx + 2, end_idx), charges):
            line_parts = lines[idx].split()
            line_parts[-2] = f"{float(charge):.4f}"  # Replace the charge value with 4 decimal places
            lines[idx] = " ".join(line_parts) + '\n'

        # Write the updated content back to the file
        with open(input_file, 'w') as f:
            f.writelines(lines)
    
    def fix_spacing_atoms(self,input_file):
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Find the index range for the "atoms" section
        start_idx = lines.index("[ atoms ]\n") + 1
        end_idx = lines.index("[ bonds ]\n")

        # Update the spacing for each line in the "atoms" section
        for idx in range(start_idx, end_idx):
            line_parts = lines[idx].split()
            line_parts[0] = line_parts[0].rjust(4)  # Right-align nr column
            line_parts[1] = line_parts[1].rjust(9)  # Right-align type column
            line_parts[2] = line_parts[2].rjust(6)  # Right-align resnr column
            line_parts[3] = line_parts[3].rjust(8)  # Right-align residue column
            line_parts[4] = line_parts[4].rjust(7)  # Right-align atom column
            line_parts[5] = line_parts[5].rjust(6)  # Right-align cgnr column
            line_parts[6] = line_parts[6].rjust(10)  # Right-align charge column
            line_parts[7] = line_parts[7].rjust(10)  # Right-align mass column
            lines[idx] = " ".join(line_parts) + '\n'
        # Write the updated content back to the file
        with open(input_file, 'w') as f:
            f.writelines(lines)


    def remove_semicolon_after_atoms(self,input_file):
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Find the index range for the "[ atoms ]" section
        start_idx = lines.index("[ atoms ]\n")
        end_idx = lines.index("[ bonds ]\n")

        # Check if the line after "[ atoms ]" has a semicolon and remove it if present
        if lines[start_idx + 1].strip().startswith(";"):
            lines[start_idx + 1] = lines[start_idx + 1].lstrip(";").lstrip()

        # Write the updated content back to the file
        with open(input_file, 'w') as f:
            f.writelines(lines)

    def add_semicolon_after_atoms(self,input_file):
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Find the index range for the "[ atoms ]" section
        start_idx = lines.index("[ atoms ]\n")
        end_idx = lines.index("[ bonds ]\n")

        # Find the index of the line containing "  nr"
        nr_line_idx = lines.index("  nr      type  resnr  residue    atom   cgnr     charge       mass\n", start_idx, end_idx)

        # Replace the first character of the line with a semicolon
        lines[nr_line_idx] = ";" + lines[nr_line_idx][1:].replace("resnr", "resnr ")

        # Write the updated content back to the file
        with open(input_file, 'w') as f:
            f.writelines(lines)
    ################################
    "Takes the generated .gjf files and converts them to .gro"
    ################################
    def gjf_to_gro(self,gjf_file, gro_file):
        atom_pattern = re.compile(r'^\s*([A-Z][a-z]?)(?:\s+[A-Z][a-z]?)*\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')

        # Initialize lists to store atom information
        atom_data = []

        # Open the .gjf file for reading
        with open(gjf_file, 'r') as input_file:
            lines = input_file.readlines()

        # Parse the .gjf file to extract atom information using regular expression
        for line in lines:
            match = atom_pattern.match(line)
            if match:
                atom_symbol = match.group(1)
                x_angstrom = float(match.group(2))
                y_angstrom = float(match.group(3))
                z_angstrom = float(match.group(4))
                atom_data.append((atom_symbol, x_angstrom, y_angstrom, z_angstrom))

        # Convert atomic coordinates from angstroms to nanometers and write to .gro file
        with open(gro_file, 'w') as output_file:
            output_file.write(f'GRO FILE GENERATED BY GJF_TO_GRO.PY BY GARDNER, Z. T.\n')
            output_file.write(f'   {len(atom_data)}\n')
            for index, atom in enumerate(atom_data, start=1):
                atom_type, x_angstrom, y_angstrom, z_angstrom = atom
                x_nm = x_angstrom / 10.0
                y_nm = y_angstrom / 10.0
                z_nm = z_angstrom / 10.0
                atom_type_hex = self.hex_style_atom_type(index)
                output_file.write(f'    1UNK{atom_type:>5}{atom_type_hex}{index:>5}'
                                f'{x_nm:>8.3f}{y_nm:>8.3f}{z_nm:>8.3f}\n')

            output_file.write('   0.00000   0.00000   0.00000\n')  # Write box dimensions
    


    ################################
    "Converts all the .gjf files in the GJF directroy to .gro files in the GRO directory"
    ################################
    def build_gro(self):
        print("----------------------------------------------------------------\n")
        print("Converting all the .gjf files in the GJF directroy to .gro files in the GRO directory")
        print("----------------------------------------------------------------\n")
        cwd = os.getcwd()
        GJF_path =  os.path.join(cwd, 'GJF')
        GRO_path =  os.path.join(cwd, 'GRO')
        gjf_files = sorted(os.listdir(GJF_path), key=lambda x: int(re.search(r'\d+', x).group()))
        gjf_files = [file for file in gjf_files if file.endswith('.gjf')]
        gjf_files = [os.path.join(GJF_path, file) for file in gjf_files] #Getting the Path to all GJF files in the GJF directory

        
        for gjf in gjf_files:
            print("Using " + os.path.basename(gjf) + " To Build " + os.path.basename(gjf.replace(".gjf", ".gro")))
            self.gjf_to_gro(gjf ,gjf.replace(".gjf", ".gro"))
    
        gro_files = [file for file in os.listdir(GJF_path) if file.endswith(".gro")]
        for gro_file in gro_files:
            src_file = os.path.join(GJF_path, gro_file)
            dest_file = os.path.join(GRO_path, gro_file)
            
            if os.path.isfile(dest_file):
                print("Warning: " + dest_file + " already exists. Overwriting...")
                os.remove(dest_file)


            os.rename(src_file, dest_file)
        print("All .gro files have been moved to the GRO directory.")

    ################################
    "Replaces Base.itp with the correct charge"
    ################################  
    def add_charges(self):
            input_file = self.itp_file  
            charges_file = self.charges_file 
            print("Reading charges from: " +  charges_file)
            charges = self.read_charges_from_file(charges_file)
            print("Successfully read charges.")
            self.update_charge_column(input_file, charges)
            print("Successfully updated charges in " + input_file)
            self.remove_semicolon_after_atoms(input_file)
            self.fix_spacing_atoms(input_file)
            self.add_semicolon_after_atoms(input_file)
            print("Successfully reformated " + input_file)
    ################################
    "Clears Params of Scanned Dihedral (does not reformat)"
    ################################
    def find_dihedral_to_scan_and_update(self,dihedral_atoms=None):
        itp_file = "base.itp"
        if dihedral_atoms is None: #if not scan indexes provided will update dihedral of dh_Info
            dihedral_atoms = self.dh_info.scan_indexes
             
        

        with open(itp_file, 'r') as f:
            lines = f.readlines()
        # Find the start and end indices of the dihedrals section
        start_idx = -1
        end_idx = -1

        for i, line in enumerate(lines):
            line = line.split()
            if len(line) > 1 and line[1] == "PROPER":
    
                start_idx = i + 1  # Skip the comment line
            elif start_idx != -1 and line[0] == ";IMPROPER":
                end_idx = i - 1 #removing the [dihedral] line
                break

        if start_idx == -1 or end_idx == -1:
            raise ValueError("Dihedrals section not found in the correct format.")

        # Convert integers to strings
        dihedral_atoms = [str(atom) for atom in dihedral_atoms]

        # Find dihedral lines with the specified atoms in the middle two positions
        dihedral_lines = [
            line.strip() for line in lines[start_idx:end_idx]
            if (line.split()[1:3] == [dihedral_atoms[1], dihedral_atoms[2]]) or
            (line.split()[1:3] == [dihedral_atoms[2], dihedral_atoms[1]])
        ]

        if dihedral_lines:
            print("Dihedral(s) found. Updating the parameters and applying changes to the .itp file...")
            updated_lines = []
            for line in lines:
                if line.strip() in dihedral_lines:
                    # Update dihedral line with zeros for parameters
                    line_parts = line.split()
                    updated_params = ["{:>12}".format("0.000")] * 6  # Right-align with the same width
                    updated_line = " ".join(line_parts[:-6] + updated_params) + "\n"  # Add newline character
                    updated_lines.append(updated_line)
                else:
                    updated_lines.append(line)

            # Write the updated lines back to the .itp file
            with open(itp_file, 'w') as f:
                f.writelines(updated_lines)
            
            print("Changes applied successfully.")
        else:
            raise ValueError("Something went wrong.      Dihedral(s) not found")

    ################################
    "Fix Dihedral Spaces"
    ################################
    def fix_dihedral_spaces(self):
        input_file = self.itp_file
        output_file = self.itp_file
        print("\nAttempting to reformat " + input_file + "...")
        # Open the input file to read data
        with open(input_file, 'r') as fin:
            lines = fin.readlines()

        # Open the output file in write mode
        with open(output_file, 'w') as fout:
            # Flag to identify the [ dihedrals ] section
            in_dihedrals = False

            # Counter for tracking the number of lines skipped in the [ dihedrals ] section
            lines_skipped = 0

            # Loop through each line in the input file
            for line in lines:
                # Check if the line starts with "[ dihedrals ]"
                if line.startswith('[ dihedrals ]'):
                    in_dihedrals = True
                    lines_skipped = 0  # Reset the lines skipped counter
                    fout.write(line)  # Write the header line as-is
                elif in_dihedrals and lines_skipped < 3:
                    lines_skipped += 1  # Increment lines skipped counter
                    fout.write(line)  # Write the line as-is (skip the first three lines after the header)
                elif in_dihedrals:
                    # Split the line into columns based on whitespace
                    columns = line.split()

                    # Check if the line has at least 11 columns (assuming the format you provided)
                    if len(columns) >= 11:
                        # Adjust the spacing and right-align each column
                        adjusted_line = ' {:>4}{:>5}{:>5}{:>5}{:>9}{:>12}{:>8}{:>8}{:>8}{:>8}{:>8}\n'.format(
                            columns[0], columns[1], columns[2], columns[3],
                            columns[4], columns[5], columns[6], columns[7],
                            columns[8], columns[9], columns[10]
                        )
                        fout.write(adjusted_line)
                    else:
                        fout.write(line)  # Write the line as-is with a leading space if it doesn't match the expected format
                else:
                    fout.write(line)  # Write the line as-is with a leading space if not in the [ dihedrals ] section

        print(f'Refromatted Successfully. Data has been written to {output_file}')
        
    ################################
    "Adds Dihedral Parameters found in the Dh_Info class to the end of each itp"
    ################################
    def check_dihedral(self,dihedral_atoms, base_itp=None):
        if base_itp == None:
            base_itp = self.itp_file

        with open(base_itp, 'r') as f:
            itp_lines = f.readlines()

        dihedrals_found = []

        for line in itp_lines:
            if line.strip().startswith(";") or "ai" in line.lower():
                continue  # Skip comment lines and header line
            dihedral_info = line.split()
            if len(dihedral_info) >= 8:
                try:
                    atoms = [int(atom) for atom in dihedral_info[0:4]]
                    if atoms == dihedral_atoms:
                        dihedrals_found.append(dihedral_info)
                except ValueError:
                    continue  # Skip lines that cannot be converted to integers

        return dihedrals_found
    def write_dihedral_restraints(self):
        itp_file = self.itp_file

        with open(itp_file) as fin:
            lines = fin.readlines()
        
        for line in lines:
            
            match = re.search(r'\[ dihedral_restraints \]', line)
            if match:
                print("Dihedral Restraints Already Found in " + itp_file + "not adding...")
                return
        
        print("\nAttempting to Add Dihedral Restraints...\n")
        dihedrals = self.dh_info.backbone_indexes
        dihedrals = dihedrals.sort_values(by="Scan_Cords?" ,ascending=False) #Puts the scan indexes first

        itp_dh = []
        print("Dihedrals in ITP Found:")
        print("[Atom1, Atom2, Atom3, Atom4, funct,c0,c1,c2,c3,c4,c5,Dihedral Angel]")
        for index, r in dihedrals.iterrows():
            dh = [int(r["Atom1"]), int(r["Atom2"]), int(r["Atom3"]), int(r["Atom4"])]
            temp = self.check_dihedral(dh)
            if temp == []:
                dh.reverse() 
                temp = self.check_dihedral(dh) #checks to see if the other direction of dihedrals can be found
            temp = temp[0] #not sure why this is need but it keeps the output clean :---/
            temp.append(r["Dihedral Angel"]) #Add dihedral angle from scan.log to end of the list
            print(temp)
            itp_dh.append(temp)
        
        print("\nWriting Found Dihedrals to " + self.itp_file + "...")
        
        with open(itp_file, 'a') as f:
            f.write(f"[ dihedral_restraints ]\n")
            f.write(";   ai    aj    ak    al  type   phi    dphi kfac\n")
            print(f"[ dihedral_restraints ]") # I like a pretty log file dont judge
            print(";   ai    aj    ak    al  type   phi    dphi kfac")
            for dh in itp_dh:
                [Atom1,Atom2,Atom3,Atom4,angle] = [dh[0], dh[1], dh[2], dh[3],dh[-1]]
                f.write(f"{Atom1:>5} {Atom2:>5} {Atom3:>5} {Atom4:>5}     1 {angle:>6.1f}     0  9000\n")
                print(f"{Atom1:>5} {Atom2:>5} {Atom3:>5} {Atom4:>5}     1 {angle:>6.1f}     0  9000")
            
            f.write("; Include Position restraint file\n")
            f.write("#ifdef POSRES\n")
            f.write("#include \"posre.itp\"\n")
            f.write("#endif\n")

    def propagate_itp(self):
        input_itp = self.itp_file
        cwd = os.getcwd()
        GRO_path =  os.path.join(cwd, 'GRO')
        gro_files = sorted(os.listdir(GRO_path), key=lambda x: int(re.search(r'\d+', x).group()))
        gro_files = [file for file in gro_files if file.endswith('.gro')]
        itp_files = [file.replace('.gro', '.itp') for file in gro_files]

        dihedral_restarint_pattern = re.compile((r'\[ dihedral_restraints \]'))
        #Finding what line the scanned dihedral is on (always two under the dihedral restraints)
        with open(input_itp, 'r') as f:
            lines = f.readlines()
            for line_number, line in enumerate(lines, start=0):
                if re.search(dihedral_restarint_pattern, line):
                    dh_line_number = line_number + 2 # Skip the header line 
                    print(f"Found the Dihedral Restraint to Update in {input_itp} @ line: {dh_line_number}")
        dh_line = lines[dh_line_number]
        print(f"Line {dh_line_number} looks like: {dh_line}")
        dh_line = dh_line.strip().split()
        [Atom1,Atom2,Atom3,Atom4,angle] = [dh_line[0], dh_line[1], dh_line[2], dh_line[3],float(dh_line[6])] 
        
        for itp in itp_files:
            #An important note is made here. We are not taking the name of each itp file and writing the angle from the name, but instead adding 10 to whatever the first inputed angle is in the base.itp
            #Hopefully this allows us to deal with other naming conventions.
            new_dh_line = f"{Atom1:>5} {Atom2:>5} {Atom3:>5} {Atom4:>5}     1 {angle:>6.1f}     0  9000\n"
            lines[dh_line_number] = new_dh_line
            with open(itp, 'w') as fout:
                fout.writelines(lines)
            angle += 10
        
        #Moving all the new ITP files
        self.create_directory("ITP") #Making a new directory for the new ITP files.
        ITP_PATH = os.path.join(cwd, 'ITP')
        for itp_file in itp_files:
            dest_file = os.path.join(ITP_PATH, itp_file)
            if os.path.isfile(dest_file):
                print("Warning: " + dest_file + " already exists. Overwriting...")
                os.remove(dest_file)
            os.rename(itp_file, dest_file)
        print("All .itp files have been moved to the ITP directory.")
    def set_up_run(self):
        dihedral_atoms = self.dh_info.scan_indexes
        paramters = self.dh_info.Parameters
        paramters = list(map(lambda x: round(float(x), 3), paramters))
        print("----------------------------------------------------------------")
        print("                  Setting up MD Run with parameters")
        print("                  " + str(paramters))
        print("                  On Dihedrals:")
        print("                  " + str(dihedral_atoms))
        print("----------------------------------------------------------------")
        index = 1
        directory_name = "run_1"
        while os.path.exists(directory_name): #checks what directory to make
            directory_name = f"run_{index}"
            index += 1
        os.mkdir(directory_name) #Making a new directory for the run
        print("Running in:" + directory_name)

        ITP_PATH = os.path.join(os.getcwd(), 'ITP')
        GRO_PATH = os.path.join(os.getcwd(),"GRO")
        new_PATH = os.path.join(os.getcwd(), directory_name)
        itp_files = sorted(os.listdir(ITP_PATH), key=lambda x: int(re.search(r'\d+', x).group()))
        gro_files = sorted(os.listdir(GRO_PATH), key=lambda x: int(re.search(r'\d+', x).group()))
        itp_files = [file for file in itp_files if file.endswith('.itp')]
        gro_files = [file for file in gro_files if file.endswith('.gro')]
        
        

        for itp_file in itp_files:
            src_file = os.path.join(ITP_PATH, itp_file)
            dest_file = os.path.join(new_PATH, itp_file)

            if os.path.isfile(dest_file):
                print(f"Warning: {dest_file} already exists. Overwriting...")
                os.remove(dest_file)

            if os.name == 'nt':  # Windows
                cmd = f'copy "{src_file}" "{dest_file}"'
            else:  # Unix/Linux
                cmd = f'cp "{src_file}" "{dest_file}"'
            
            os.system(cmd)

        for gro_file in gro_files:
            src_file = os.path.join(GRO_PATH, gro_file)
            dest_file = os.path.join(new_PATH, gro_file)

            if os.path.isfile(dest_file):
                print(f"Warning: {dest_file} already exists. Overwriting...")
                os.remove(dest_file)

            if os.name == 'nt':  # Windows
                cmd = f'copy "{src_file}" "{dest_file}"'
            else:  # Unix/Linux
                cmd = f'cp "{src_file}" "{dest_file}"'
            
            os.system(cmd)
        
        #Now that everything is moved into the new directory will update to the paramters from Dh_info
        
        for itp_file in itp_files:
            itp_file_path = os.path.join(new_PATH, itp_file)
            print(f"Processing file: {itp_file_path}")
            self.find_dihedral_line_and_replace(itp_file_path, dihedral_atoms, paramters)
            self.adjust_dihedral_section(itp_file_path)
            #Building .top
            self.build_top(itp_file.split('.')[0],new_PATH)
            
        
    def find_dihedral_line_and_replace(self,itp_file, dihedral_atoms, new_parameters):
        dihedral_atoms = [int(dihedral_atom) for dihedral_atom in dihedral_atoms]
        new_parameters = [str(new_parameter) for new_parameter in new_parameters]
       
        try:
            with open(itp_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"Error: File {itp_file} not found.")
            return

        updated_lines = []  # List to store updated lines

        # Search for a line that starts with the specified dihedral_atoms
        found_line = False
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    current_atoms = [int(part) for part in parts[:4]]
                except ValueError:
                    updated_lines.append(line)
                    continue
                
                if current_atoms == dihedral_atoms:
                    # Ensure we have enough parts to replace the last six numbers
                    if len(parts) >= 10:
                        parts[-6:] = new_parameters
                        
                        updated_line = " ".join(parts) + "\n"
                        updated_lines.append(updated_line)
                        found_line = True
                    else:
                        updated_lines.append(line)
                else:
                    updated_lines.append(line)
            else:
                updated_lines.append(line)

        if not found_line:
            print(f"Error: No line with specified dihedral atoms found in {itp_file}.")
            return

        # Write the updated lines back to the input file
        with open(itp_file, 'w') as f:
            f.writelines(updated_lines)

    def adjust_dihedral_section(self,itp_file):
        try:
            with open(itp_file, 'r') as fin:
                lines = fin.readlines()
        except FileNotFoundError:
            print(f"Error: File {itp_file} not found.")
            return

        output_file = itp_file  # Since we want to overwrite the same file

        with open(output_file, 'w') as fout:
            in_dihedrals = False
            lines_skipped = 0

            for line in lines:
                if line.startswith('[ dihedrals ]'):
                    in_dihedrals = True
                    lines_skipped = 0
                    fout.write(line)
                elif in_dihedrals and lines_skipped < 3:
                    lines_skipped += 1
                    fout.write(line)
                elif in_dihedrals:
                    columns = line.split()
                    if len(columns) >= 11:
                        adjusted_line = ' {:>4}{:>5}{:>5}{:>5}{:>9}{:>12}{:>8}{:>8}{:>8}{:>8}{:>8}\n'.format(
                            columns[0], columns[1], columns[2], columns[3],
                            columns[4], columns[5], columns[6], columns[7],
                            columns[8], columns[9], columns[10]
                        )
                        fout.write(adjusted_line)
                    else:
                        fout.write(line)
                else:
                    fout.write(line)

    def build_top(self,name,out_path): 
        
        path = os.path.join(out_path, f"{name}.top")
        with open(path, "w") as top_file:
            top_file.write(
                "#define _FF_OPLS\n"
                "#define _FF_OPLSAA\n"
                "\n"
                ";[ defaults ]\n"
                "; nbfunc    comb-rule   gen-pairs   fudgeLJ    fudgeQQ\n"
                ";1            3         yes       0.5        0.5\n"
                "\n"
                "; Include force field parameters\n"
                f"#include \"{name}.itp\"\n"
                "\n"
                "[system]\n"
                "system\n"
                "\n"
                "[ molecules ]\n"
                "; Name      Number\n"
                "UNK   1\n"
            )
        print("Built topology: " +f"{name}.top" )


def copy_script_to_run_dir():
    src_file = os.path.join(os.getcwd(), "Codes/paramaterization.py")
    if not os.path.isfile(src_file):
        print(f"Error: Source file '{src_file}' does not exist.")
        return

    # Get a list of all directories starting with "run_"
    run_directories = [d for d in os.listdir() if d.startswith("run_")]

    for directory in run_directories:
        dest_file = os.path.join(directory, os.path.basename(src_file))

        # Check if the destination file already exists
        if os.path.isfile(dest_file):
            continue

        if os.name == 'nt':  # Windows
            cmd = f'copy "{src_file}" "{dest_file}"'
        else:  # Unix/Linux
            cmd = f'cp "{src_file}" "{dest_file}"'
        os.system(cmd)



Info = Dh_Info() 
S = Scan(Info)
S.make_itp_files() 
S.set_up_run()
copy_script_to_run_dir()

print("\n----------------------------------------------------------------\n")
print("                  dihdral_scan.py Completed Successfully")
print("\n----------------------------------------------------------------\n")