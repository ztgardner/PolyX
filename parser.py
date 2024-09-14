

class Itp_parser:
    def __init__(self, File: str = None):
        self.load_path = File
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
        self.sections = tuple(sections)
        self.sections_to_data_dic = {i: [] for i in self.sections}
        self.sections_to_comments = {i: [] for i in self.sections}

        self.number_of_sections = len(self.sections)
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
        self.sections_to_pure_data_dic = {}
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

    def __repr__(self) -> str:
        return str(self.sections)

    def __getitem__(self, i:str) -> list:
        key = "[" + i + "]" if "[" + i + "]" in self.sections else "[ " + i + " ]"
        return self.sections_to_pure_data_dic[key]

    def __len__(self) -> int:
        return len(self.sections_to_pure_data_dic)

    def __setitem__(self, i:str, j:list):
        key = "[" + i + "]" if "[" + i + "]" in self.sections else "[ " + i + " ]"
        self.sections_to_pure_data_dic[key] = j

    def raw(self, section) -> list:
        key = (
            "[" + section + "]"
            if "[" + section + "]" in self.sections
            else "[ " + section + " ]"
        )
        return self.sections_to_data_dic[key]

    def save_itp(self, name:str):
        current_itp = self.fill_in_data()
        with open(name, "w") as f:
            for i, j in current_itp.items():
                f.write(i)
                f.write("\n")
                f.writelines(j)

    def fill_in_data(self) -> dict:
        current_itp = {}
        for i, j in self.sections_to_pure_data_dic.items():
            current_itp[i] = (
                self.sections_to_comments[i]
                + self.sections_to_pure_data_dic[i]
                + ["\n"]
            )
        return current_itp
