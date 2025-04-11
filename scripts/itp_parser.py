import os
from random import randint




class Itp_parser:
    """
    A class to parse and manipulate gromacs itp files.

    Attributes:
        load_path (str): Loads the data to the data frame.
        sections (tuple): A tuple of section headers extracted from the file.
        sections_to_data_dic (dict): A dictionary mapping section headers to raw data lines.
        sections_to_comments (dict): A dictionary mapping section headers to comments.
        number_of_sections (int): The number of sections in the ITP file.
        sections_to_pure_data_dic (dict): A dictionary mapping section headers to pure data lines (excluding comments).

    """
    # The information is organized as following: (number of atoms, fuc type, number of parameters, total entires-- sum of num atom, number of paramters +1)
    BONDS = {
        "bond": (2, 1, 2, 5),
        "G96 bond": (2, 2, 2, 5),
        "Morse": (2, 3, 3, 6),
        "cubic bond": (2, 2, 2, 5),
        "connection": (2, 5, 0, 3),
        "harmonic potentia": (2, 6, 2, 5),
        "FENE bond": (2, 7, 2, 5),
        "tabulated bond1": (2, 8, 2, 5),
        "tabulated bond2": (2, 9, 2, 5),
        "restraint potential": (2, 10, 2, 5),
    }
    DIHEDRALS = {
        "proper dihedral": (4, 1, 2, 5),
        "Ryckaert-Bellemans dihedral": (4, 3, 6, 11),
        "improper_dihedral": (4, 2, 2, 5),
        "periodic improper dihedral": (4, 4, 6, 11),
        "Fourier dihedral": (4, 5, 2, 5),
        "tabulated dihedral": (4, 8, 6, 11),
        "proper dihedral (multiple)": (4, 9, 2, 5),
        "restricted dihedral": (4, 10, 6, 11),
        "combined bending-torsion potential": (4, 11, 2, 5),
    }
    PAIRS = {
        "extra LJ or Coulomb1": (2, 1, 2, 5),
        "extra LJ or Coulomb2": (2, 2, 4, 7),
    }
    PAIRS_NB = {
        "extra LJ or CoulombNB": (2, 1, 3, 6),
    }
    ANGLES = {
        "angle": (3, 1, 2, 6),
        "G96 angle": (3, 2, 2, 5),
        "cross bond-bond": (3, 3, 2, 5),
        "cross bond-angle": (3, 4, 2, 5),
        "Urey-Bradley": (3, 5, 2, 5),
        "quartic angle": (3, 6, 2, 5),
        "tabulated angle": (3, 8, 2, 5),
        "linear angle": (3, 9, 2, 5),
        "restricted bending potential": (3, 10, 2, 5),
    }
    EXCLUSIONS = {
        "exclusions": (1, 0, 1, 0),
    }
    CONSTRAINTS = {
        "constraints1": (2, 1, 1, 4),
        "constraints2": (2, 2, 1, 4),
    }
    SETTLES = {"SETTLE": (1, 1, 0, 0)}
    VIRTUAL_SITES1 = {
        "1-body virtual site": (2, 1, 0, 3),
    }
    VIRTUAL_SITES2 = {
        "2-body virtual site": (3, 1, 2, 5),
        "2-body virtual site (fd)": (3, 2, 2, 5),
    }
    VIRTUAL_SITES3 = {
        "3-body virtual site": (4, 1, 2, 7),
        "3-body virtual siteFD": (4, 2, 2, 7),
        "3-body virtual siteFAD": (4, 3, 2, 7),
        "3-body virtual siteOUT": (4, 4, 2, 7),
    }
    VIRTUAL_SITES4 = {
        "4-body virtual site (fdn)": (5, 2, 2, 8),
    }
    VIRTUAL_SITESn = {
        "N-body virtual site (COG)": (1, 1, 0, 2),
        "N-body virtual site (COM)": (1, 2, 0, 2),
        "N-body virtual site (COW)": (1, 3, 0, 2),
    }
    POSITION_RESTRAINT = {
        "position restraint": (1, 1, 1, 3),
        "flat-bottomed position restraint": (1, 2, 1, 3),
    }
    DISTANCE_RESTRAINTS = {
        "distance restraint1": (2, 1, 4, 7),
        "distance restraint2": (4, 1, 3, 8),
    }
    ORIENTATION_RESTRAINTS = {
        "orientation restraint": (2, 1, 6, 9),
    }
    ANGLES_RESTRAINTS = {
        "angle restraint": (4, 1, 3, 8),
        "angle restraint (z)": (2, 1, 3, 6),
    }
    TOTAL_DIRECTIVES = {
        "BONDS": BONDS,
        "DIHEDRALS": DIHEDRALS,
        "PAIRS": PAIRS,
        "PAIRS_NB": PAIRS_NB,
        "ANGLES": ANGLES,
        "EXCLUSIONS": EXCLUSIONS,
        "CONSTRAINTS": CONSTRAINTS,
        "SETTLE": SETTLES,
        "VIRTUAL_SITES1": VIRTUAL_SITES1,
        "VIRTUAL_SITES2": VIRTUAL_SITES2,
        "VIRTUAL_SITES3": VIRTUAL_SITES3,
        "VIRTUAL_SITES4": VIRTUAL_SITES4,
        "VIRTUAL_SITESn": VIRTUAL_SITESn,
        "POSITION_RESTRAINT": POSITION_RESTRAINT,
        "DISTANCE_RESTRAINTS": DISTANCE_RESTRAINTS,
        "ORIENTATION_RESTRAINTS": ORIENTATION_RESTRAINTS,
        "ANGLES_RESTRAINTS": ANGLES_RESTRAINTS,
    }

    def __init__(self, File: str = "Base.itp"):
        """
        Initializes the Itp_parser instance with a file path and parses the file.

        Args:
            File (str, optional): The path to the ITP file to load. Defaults to Base.itp.
        """
        if not os.path.isfile(File):
            raise FileNotFoundError(File)
        self.load_path = File
        file_lines, sections = self._load()
        self.sections = tuple(sections)
        self.sections_to_data_dic = {}
        self.sections_to_comments = {}
        self.number_of_sections = len(self.sections)
        self.sections_to_pure_data_dic = {}
        self.converte_internal_sections_to_gromacs_sections={"atoms":"atoms","moleculetype":"moleculetype"}
        self._extract_data(sections, file_lines)
        self._clean_data()
        sections = self._match_tempsection_to_proper_section()
        self.sections = tuple(sections)
        self.DF = self._make_sections_data_to_df()

    def __repr__(self) -> str:
        """
        Returns a string representation of the Itp_parser instance.

        Returns:
            str: A string representation of the sections.
        """
        return str(self.DF)

    def __getitem__(self, i: str) -> list or None:
        """
        Retrieves the pure data for a given section.

        Args:
            i (str): The name of the section to retrieve.

        Returns:
            list: The list of pure data lines for the specified section.
        """
        key = i
        if key in self.sections:
            return self.DF[key]
        else:
            return None

    def __len__(self) -> int:
        """

        Returns:
            int: The number of sections.
        """
        return len(self.sections)

    def __setitem__(self, i: str, j: list)-> None:
        """
        Sets the pure data for a given section.

        Args:
            i (str): The name of the section to update.
            j (list): The list of pure data lines to set for the section.
        """
        internal_section=i
        info_aboutsection=self.converte_internal_sections_to_gromacs_sections[i][1]
        splited=[]
        a = int(info_aboutsection[0])
        d = int(info_aboutsection[3])
        for lines in j:
            splited.append(lines.split())


        atoms_to_extend=[slip[:a] for slip in splited]
        params_to_extend=[

                (
                    slip[a + 1:]
                    if not self.check_if_comment(slip)
                    else slip[a + 1: d - 1]
                )
                for slip in splited

            ]
        func_to_extend=[[str(info_aboutsection[1])] for _ in range(len(splited))]

        self.DF[internal_section]["atoms"].extend(atoms_to_extend)
        self.DF[internal_section]["params"].extend(params_to_extend)
        self.DF[internal_section]["function_type"].extend(func_to_extend)


    def __iter__(self)-> None:
        for key,val in self.DF.items():
            yield key, val

    def raw(self, section) -> list:
        """
        Retrieves the raw data (including comments) for a given section.

        Args:
            section (str): The name of the section to retrieve.

        Returns:
            list: The list of raw data lines (including comments) for the specified section.
        """
        key = section
        return self.sections_to_data_dic[key]

    def save_itp(self, name: str)-> None:
        """
        Saves the current state of the ITP data to a file.

        Args:
            name (str): The name of the file to save the data to.
        """
        with open(name, "w") as f:
            f.write("; \n")
            f.write("; \n")
            f.write("; \n")
            f.write(";              An extender, a parser, a website -- are you getting it?????"+"\n")
            f.write(";                                  ITP made through POLYX"+"\n")
            f.write("; \n")
            f.write("; \n")
            f.write("; \n")
            f.write(" \n")
            f.write(" \n")

            for key, data in self.DF.items():
                section_as_a_string=str(key).replace("[ ", "").replace("]","").strip()
                match section_as_a_string:
                    case "moleculetype" |"atoms":
                        f.writelines(f"[ {key} ]")
                        f.write("\n")
                    case _:
                        actual_section=self.find_topology_section(key)
                        f.writelines(f"{actual_section}")
                        f.write("\n")

                strings_to_write = None
                bottom_coments=[]
                for inside_key,inside_data in data.items():
                    if len(inside_data) > 0:
                        match inside_key:
                            case "Comments_top":
                                for top_coment in inside_data:
                                    f.writelines(str(top_coment)+"\n")
                            case "Coments_bottom":
                                for bottom_coment in inside_data:
                                    bottom_coments.append(str(bottom_coment))
                            case _:
                                if not strings_to_write:
                                    strings_to_write=['' for i in range(len(inside_data))]
                                for param_data, index in zip(inside_data, range(len(inside_data))):
                                    current_string=strings_to_write[index]
                                    current_string+=str(param_data).replace(","," ").strip("[").strip("]").replace("'"," ").center(10)
                                    current_string+=" "
                                    if "POLYX" in current_string:
                                        current_string=current_string.replace("POLYX","   ")
                                    strings_to_write[index]=current_string
                for i in strings_to_write:
                    f.write(str(i)+"\n")
                for bottom_coment in bottom_coments:
                    f.write(str(bottom_coment)+"\n")
                f.write(" \n")


    def fill_in_data(self) -> dict:
        """
        Fills in the data for each section, combining comments and pure data.

        Returns:
            dict: A dictionary where keys are section headers and values are lists of lines combining comments and data.
        """
        current_itp = {}
        for i, j in self.sections_to_pure_data_dic.items():
            current_itp[i] = (
                self.sections_to_comments[i]
                + self.sections_to_pure_data_dic[i]
                + ["\n"]
            )
        return current_itp


    def _load(self)-> tuple:
        """
        Loads the ITP file and extracts lines and section headers. PRIVATE FUNCTION.

        Returns:
            tuple: A tuple containing a list of file lines and a list of section headers.
        """
        file_lines = []
        sections = []
        with open(self.load_path, "r") as f:
            file_lines.extend(f.readlines())
        for lines in file_lines:
            if lines.strip():
                match lines.strip()[0]:
                    case "[":
                        sections.append(lines.strip())
                    case _:
                        pass
        return file_lines, sections

    def _extract_data(self, sect, file_lines) -> None:
        """
        Extracts data for each section from the file lines. PRIVATE FUNCTION.

        Args:
            sect (list): A list of section headers.
            file_lines (list): A list of lines from the ITP file.
        """
        sections = sect
        keys = []
        self.index_reached=0

        for i in range(self.number_of_sections):
            file_left_to_parse = file_lines[self.index_reached:]
            section = sections[0]
            sections.pop(0)
            Found = False
            First = True
            Temp_key = section + f"&{randint(1,1000000)}"
            keys.append(Temp_key)
            self.sections_to_data_dic[Temp_key] = []
            times_iterated=0
            for index,line in enumerate(file_left_to_parse):
                times_iterated += 1
                if line.strip() in sections and not First:
                    if not self.sections_to_data_dic[Temp_key]:
                        raise(ValueError(f"section is empty{Temp_key}, line: {times_iterated}"))

                    First = True
                    self.index_reached+=times_iterated-1
                    break
                if Found:
                    self.sections_to_data_dic[Temp_key].append(line)
                if section in line:

                    Found = True
                    First = False
        self.sections = tuple(keys)


    def _clean_data(self) -> None:
        """
        Cleans up the data by separating comments from pure data. PRIVATE FUNCTION.
        """
        for i, j in self.sections_to_data_dic.items():
            pure_data = []
            comments = []
            for lines in j:
                if lines.strip():
                    match lines.strip()[0]:
                        case ";":
                            comments.append(lines)
                        case "[":
                            pass
                        case _:
                            pure_data.append(lines)
            self.sections_to_pure_data_dic[i] = pure_data
            self.sections_to_comments[i] = comments

    @staticmethod
    def find_repeating_sections(sections: list) -> list:
        append = []
        repeat = []
        for section in sections:
            if section not in append:
                sections.append(section)
            else:
                repeat.append(section)
        return repeat


    def _match_tempsection_to_proper_section(self) -> list:
        sections = []
        for i, j in self.sections_to_pure_data_dic.items():

            rud_directive = i.split("&")[0]
            directive = rud_directive[1 : len(rud_directive) - 1].strip().upper()
            if directive in self.TOTAL_DIRECTIVES.keys():
                directive_info = self.TOTAL_DIRECTIVES[directive]
                if len(directive_info) != 1 and len(j)!=0:

                    for k, l in directive_info.items():
                        direc_type = k
                        number_atoms_in_param = l[0]
                        fun_type = j[0].split()[number_atoms_in_param] if len(j[0].split())>=number_atoms_in_param else j[1].split()[number_atoms_in_param]
                        if int(fun_type) == l[1]:
                            sections.append(direc_type)
                            break
                else:
                    sections.append(directive.lower())
            else:
                sections.append(directive.lower())
        for i, j in zip(self.sections, sections):
            self.sections_to_data_dic[j] = self.sections_to_data_dic.pop(i)
            self.sections_to_pure_data_dic[j] = self.sections_to_pure_data_dic.pop(i)
            self.sections_to_comments[j] = self.sections_to_comments.pop(i)
        return sections

    def _make_sections_data_to_df(self)-> dict:

        self.check_if_comment = lambda s: any(i.startswith("#") for i in s)
        Dataframe_dic = {}

        for sections in self.sections:
            match sections:
                case "moleculetype":
                    spliter = []
                    sections_data = self.sections_to_data_dic[sections]
                    top_comment = []
                    bottom_comment = []
                    data_found = False
                    for lines in sections_data:
                        if lines.startswith(";"):
                            if data_found:
                                bottom_comment.append(lines)
                            else:
                                top_comment.append(lines)
                        else:
                            data_found = True
                            if len(lines.split()) != 0:
                                spliter.append(lines.split())
                    data_for_moletype = {
                        "Comments_top": top_comment,
                        "ResidueName": [[spliter[0][0]]],
                        "nrexcl": [[spliter[0][1]]],
                        "comments": [spliter[0][2:]],
                        "Coments_bottom": bottom_comment,
                    }
                    Dataframe_dic[sections] = data_for_moletype

                case "atoms":
                    slipters = []
                    found_data = False
                    bottom_comment = []
                    top_comment = []
                    for lines in self.sections_to_data_dic[sections]:
                        if lines.startswith(";"):
                            if found_data:
                                bottom_comment.append(lines)
                            else:
                                top_comment.append(lines)
                        else:
                            found_data = True
                            if len(lines.split()) != 0:
                                slipters.append(lines.split())
                    to_add = {
                        "Comments_top": top_comment,
                        "atoms": [[slip[0]] if len(slip)>0 else ["POLYX"] for slip in slipters],
                        "atom_types": [[slip[1]] if len(slip)>1 else ["POLYX"] for slip in slipters],
                        "resodue#": [[slip[2]] if len(slip)>2 else ["POLYX"] for slip in slipters],
                        "residue_name": [[slip[3]] if len(slip)>3 else ["POLYX"] for slip in slipters],
                        "atom_name": [[slip[4]] if len(slip)>4 else ["POLYX"] for slip in slipters],
                        "chargeGroups#": [[slip[5]] if len(slip)>5 else ["POLYX"] for slip in slipters],
                        "charge": [[slip[6]] if len(slip)>6 else ["POLYX"] for slip in slipters],
                        "mass": [[slip[7]] if len(slip)>7 else ["POLYX"] for slip in slipters],
                        "comments": [["".join(slip[8:])] if len(slip)>8 else ""for slip in slipters ],
                        "Coments_bottom": bottom_comment,
                    }
                    Dataframe_dic[sections] = to_add

                case _:
                    info_aboutsection = None

                    for key, val in self.TOTAL_DIRECTIVES.items():
                        for sect, stuff in val.items():
                            if sections == sect:
                                info_aboutsection = stuff
                                self.converte_internal_sections_to_gromacs_sections[str(sect)] = [str(key).lower(),stuff]
                                break
                            if info_aboutsection:
                                break
                        if info_aboutsection:
                            break

                    found_data = False
                    top_comment = []
                    bottom_comment = []
                    current_data = self.sections_to_data_dic[sections]
                    slipters = []

                    for lines in current_data:
                        if info_aboutsection:
                            a = int(info_aboutsection[0])
                            d = int(info_aboutsection[3])
                            if lines.startswith(";"):
                                if found_data:
                                    bottom_comment.append(lines)
                                else:
                                    top_comment.append(lines)
                            else:
                                found_data = True
                                slipter = lines.split()
                                if len(slipter)>0:
                                    slipters.append(slipter)
                            to_add = {
                                "Comments_top": top_comment,
                                "atoms": [slip[:a] for slip in slipters],
                                "function_type": [[str(info_aboutsection[1])] for _ in slipters ],
                                "params": [

                                        (
                                            slip[a + 1 :]
                                            if not self.check_if_comment(slip)
                                            else slip[a + 1 : d - 1]
                                        )
                                        for slip in slipters

                                ],
                                "comments": [
                                    "".join(slip[d - 1 :]) if self.check_if_comment(slip) else ""
                                    for slip in slipters
                                ],
                                "Coments_bottom": bottom_comment,
                            }
                            Dataframe_dic[sections] = to_add
        return Dataframe_dic


    def load_gro(self, gro_file:str)-> None:
        if not isinstance(gro_file, str) or not gro_file.endswith(".gro"):
            raise ValueError("Cordinate File must be a .gro file'")
        with open(gro_file, "r") as f:
            gro_lines = f.readlines()
            cords = []
            for line_number, line in enumerate(gro_lines):
                line = line.split()
                if line_number <= 1 or len(line) < 4:
                    continue
                x, y, z = [float(line[3]), float(line[4]), float(line[5])]
                cords.append([x, y, z])
        self.coordinates = cords
        print(f"Loaded Coordinates From {gro_file}")

    def load_pdb(self, pdb_file:str) -> None:
        if not isinstance(pdb_file, str) or not pdb_file.endswith(".pdb"):
            raise ValueError("Cordinate File must be a .pdb file'")
        with open(pdb_file, "r") as f:
            pdb_lines = f.readlines()
            cords = []
            for line_number, line in enumerate(pdb_lines):
                line = line.split()
                if line_number <= 4 or len(line) < 4:
                    print(f"Skipping Line {line} Containing: {line}")
                    continue
                x, y, z = [float(line[5]), float(line[6]), float(line[7])]
                cords.append([x, y, z])
        self.coordinates = cords
        print(f"Loaded Coordinates From {pdb_file}")



    def  find_topology_section(self,name:str)-> str:
        name=name.replace("["," ").replace("]"," ").strip()
        for key, inner_dic in self.TOTAL_DIRECTIVES.items():
            if name in inner_dic.keys():
                return f"[ {str(key).lower()} ]"
        raise ValueError(f"Topology section:{name}  does not correspond to an gromacs topology section")


    def add_to_comments(self, section: str, extend: list,bottom: bool=False) -> None:
        internal_section=section
        sect="Coments_bottom" if bottom else "Comments_top"
        self.DF[internal_section][sect].extend(extend)