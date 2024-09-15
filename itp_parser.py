import os


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
        self.sections_to_data_dic = {i: [] for i in self.sections}
        self.sections_to_comments = {i: [] for i in self.sections}
        self.number_of_sections = len(self.sections)
        self.sections_to_pure_data_dic = {}
        self._extract_data(sections, file_lines)
        self._clean_data()

    def __repr__(self) -> str:
        """
        Returns a string representation of the Itp_parser instance.

        Returns:
            str: A string representation of the sections.
        """
        return str(self.sections)

    def __getitem__(self, i: str) -> list:
        """
        Retrieves the pure data for a given section.

        Args:
            i (str): The name of the section to retrieve.

        Returns:
            list: The list of pure data lines for the specified section.
        """
        key = "[" + i + "]" if "[" + i + "]" in self.sections else "[ " + i + " ]"
        return self.sections_to_pure_data_dic[key]

    def __len__(self) -> int:
        """

        Returns:
            int: The number of sections.
        """
        return len(self.sections)

    def __setitem__(self, i: str, j: list):
        """
        Sets the pure data for a given section.

        Args:
            i (str): The name of the section to update.
            j (list): The list of pure data lines to set for the section.
        """
        key = "[" + i + "]" if "[" + i + "]" in self.sections else "[ " + i + " ]"
        self.sections_to_pure_data_dic[key] = j

    def raw(self, section) -> list:
        """
        Retrieves the raw data (including comments) for a given section.

        Args:
            section (str): The name of the section to retrieve.

        Returns:
            list: The list of raw data lines (including comments) for the specified section.
        """
        key = (
            "[" + section + "]"
            if "[" + section + "]" in self.sections
            else "[ " + section + " ]"
        )
        return self.sections_to_data_dic[key]

    def save_itp(self, name: str):
        """
        Saves the current state of the ITP data to a file.

        Args:
            name (str): The name of the file to save the data to.
        """
        current_itp = self.fill_in_data()
        with open(name, "w") as f:
            for i, j in current_itp.items():
                f.write(i)
                f.write("\n")
                f.writelines(j)

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

    def set_charge(self, charge_list: list):
        """
        Updates the charge values in the 'atoms' section with the provided charge list.

        Args:
            charge_list (list): A list of new charge values to update.
        """
        key = (
            "[" + "atoms" + "]"
            if "[" + "atoms" + "]" in self.sections
            else "[ " + "atoms" + " ]"
        )
        index_to_change = len(self.sections_to_pure_data_dic[key][2].split()) - 2
        new_lines = self.set_value(
            self.sections_to_pure_data_dic[key], index_to_change, charge_list
        )
        self.sections_to_pure_data_dic[key] = new_lines

    @staticmethod
    def set_value(to_change: list, index_to_change: int, value_to_set: list) -> list:
        """
        Updates specific values in a list of strings based on the provided index and new values.

        Args:
            to_change (list): The list of lines to update.
            index_to_change (int): The index where the values should be updated.
            value_to_set (list): The list of new values to set.

        Returns:
            list: The updated list of lines with new values applied.
        """
        new_lines = []
        for i, j in zip(to_change, value_to_set):
            words = i.split()
            len_word = []
            index = []

            for word in words:
                index.append(i.index(word))
                len_word.append(len(word))
            new_line = (
                str(i[: index[index_to_change]])
                + str(j).rjust(len_word[index_to_change])
                + str(i[index[index_to_change] + len_word[index_to_change] :])
            )
            new_lines.append(new_line)
        return new_lines

    def _load(self):
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

    def _extract_data(self, sect, file_lines):
        """
        Extracts data for each section from the file lines. PRIVATE FUNCTION.

        Args:
            sect (list): A list of section headers.
            file_lines (list): A list of lines from the ITP file.
        """
        sections = sect
        for i in range(self.number_of_sections):
            section = sections[0]
            sections.pop(0)
            Found = False
            Done = False
            for line in file_lines:
                if line.strip() in sections:
                    Done = True
                if Found and not Done:
                    self.sections_to_data_dic[section].append(line)
                if section in line:
                    Found = True

    def _clean_data(self):
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
                        case _:
                            pure_data.append(lines)
            self.sections_to_pure_data_dic[i] = pure_data
            self.sections_to_comments[i] = comments
