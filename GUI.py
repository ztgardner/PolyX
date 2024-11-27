import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from mpl_toolkits.mplot3d import proj3d
from itp_parser import Itp_parser 


class GUI():
    def __init__(self, root, itp_loader_object = None):
        self.root = root
       

        if itp_loader_object is None:
            #Example coordinates
            self.root.title("Example Data")
            self.coordinates = np.array([[1, 2, 3], [3, 5, 6], [6, 8, 7], [8, 3, 9]])  # 3D coordinates
            self.vectors = np.array([[2, 3, 1], [3, 1, -2], [0, -5, 4], [1, 2, -3]])   # 3D vectors
            print("Example ITP Loaded into GUI...")
        else:
            itp = itp_loader_object
            self.itp = itp_loader_object
            self.root.title(itp.load_path)
            self.sections = itp.sections
            self.itp_currently_shown = itp_loader_object
            print("Loaded ITP into GUI...")

            
            self.coordinates = np.array(itp_obj.coordinates)  # 3D coordinates
            self.atom_types = self.itp_currently_shown.DF["atoms"]['atom_name'][0] #The atom types in the same order as the plotted cords
            self.atom_colors = [self.atom_type_to_color(type) for type in self.atom_types]
            
            
            self.vectors = self.section_to_vectors() # 3D vectors



        self.setup_layout()
        # Draw initial plot
        self.plot_atoms_vectors()
        #Setting up clicks
        self.prev_clicked_atom = None
    def setup_layout(self):
        """Set up the GUI layout, including the 3D plot, toolbar, and controls."""


        # === Create the Display Area for the 3D Plot ===
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')  # Create a 3D axis
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=0, column=0, columnspan=2, sticky="nsew")
        self.canvas.draw()

        

        # === Add Toolbar for Interactive Features ===
        toolbar_frame = tk.Frame(self.root)
        toolbar_frame.grid(row=1, column=0, columnspan=2, sticky="ew")
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()

        # === Create the Textbox for Information ===
        self.info_box = tk.Text(self.root, height=5, width=50, wrap=tk.WORD)
        self.info_box.grid(row=2, column=0, sticky="nsew", padx=5, pady=5)

        # === Create Buttons and Dropdown Controls ===
        button_frame = tk.Frame(self.root)
        button_frame.grid(row=2, column=1, sticky="nsew", padx=5, pady=5)

        # Add a Dropdown menu
        self.dropdown_var = tk.StringVar()
        self.dropdown = ttk.Combobox(button_frame, textvariable=self.dropdown_var,style="TCombobox")

        # Populate dropdown with available sections, excluding "moleculetype"
        section = [entry for entry in self.sections if entry != "moleculetype"]
        self.dropdown["values"] = section

        # Set default selection if the list is not empty
        if section:
            self.dropdown.current(0)  # Default to the first available option
        self.dropdown.pack(fill=tk.X)
        self.dropdown.bind("<<ComboboxSelected>>", lambda event: self.handle_dropdown_selection())

        # Add a Back button
        self.back_button = tk.Button(button_frame, text="Back", command=self.go_back)
        self.back_button.pack(fill=tk.X, padx=2.5, pady=5)

       

        # === Configure Grid Layout for Resizing ===
        self.root.grid_rowconfigure(0, weight=3)  # 3D plot
        self.root.grid_rowconfigure(2, weight=1)  # Info box and controls
        self.root.grid_columnconfigure(0, weight=3)  # Info box
        self.root.grid_columnconfigure(1, weight=1)  # Buttons

        # === Connect Plot Click Event ===
        self.canvas.mpl_connect("button_press_event", self.on_click)
        self.canvas.mpl_connect("scroll_event", self.on_scroll)

    def on_scroll(self, event):
        """Handle mouse scroll event for zooming."""
        base_scale = 1.1  # Scale factor for zooming
        if event.button == 'up':
            # Zoom in
            scale_factor = 1 / base_scale
        elif event.button == 'down':
            # Zoom out
            scale_factor = base_scale
        else:
            # No action for other scroll events
            return

        # Adjust the 3D view limits
        self.ax.set_xlim([scale_factor * x for x in self.ax.get_xlim()])
        self.ax.set_ylim([scale_factor * y for y in self.ax.get_ylim()])
        self.ax.set_zlim([scale_factor * z for z in self.ax.get_zlim()])

        # Redraw the canvas with updated limits
        self.canvas.draw()       
    def plot_atoms_vectors(self):
        """Plot 3D atoms and vectors on the canvas."""
        self.ax.clear()
        self.ax.grid(False)
        self.ax.set_axis_off()
        # Plot the 3D atoms
        x, y, z = self.coordinates[:, 0], self.coordinates[:, 1], self.coordinates[:, 2]
        self.ax.scatter(x, y, z, color=self.atom_colors, s=50, zorder=5, marker='o')
        self.drawn_atom = range(1, len(x) + 1) #Length of x is how many atoms are plotted(this only works since I am plotting all atoms)
        # Plot the 3D vectors
        for i in range(len(self.vectors)):
            
            data = self.vectors[i]
            start = data[0]
            vec = data[1]
            
            self.ax.quiver(start[0], start[1], start[2], vec[0], vec[1], vec[2], color="red", zorder=4)
        self.canvas.draw()

    def on_click(self, event):
        if event.inaxes == self.ax:
            closest_atom = self.get_closest_atom(event)
            if closest_atom is not None:
                #print(f"Clicked atom index: {closest_atom}")
                #print(self.itp_currently_shown.DF["atoms"]['atom_name'][0][closest_atom]) #Shows the pressed atom type
                #print(self.atom_type_to_color(self.itp_currently_shown.DF["atoms"]['atom_name'][0][closest_atom])) #Shows the pressed color
                
                # Reset the color of the previously clicked atom
                if self.prev_clicked_atom is not None:
                    self.ax.scatter(*self.coordinates[self.prev_clicked_atom], color="blue", s=50, zorder=6, marker='o')  # Reset to original color
                
                
                # Update the textbox with information related to the clicked atom
                self.update_atom_info(closest_atom)
                
                # Store the clicked atom as the previous clicked atom for next time
                self.prev_clicked_atom = closest_atom
                
                # Redraw the figure to reflect changes
                self.ax.grid(False)
                self.ax.set_axis_off()
                self.fig.canvas.draw()

    def get_closest_atom(self, event):
        """Returns the index of the closest atom to the click position in 3D space."""
        screen_coords = []
        
        # Loop over each 3D point and project it into 2D
        for x, y, z in self.coordinates:
            x2, y2, _ = proj3d.proj_transform(x, y, z, self.ax.get_proj())  # Use proj3d module
            x2_screen, y2_screen = self.ax.transData.transform((x2, y2))
            screen_coords.append((x2_screen, y2_screen))

        # Calculate the distances between the click and each projected point
        distances = [np.sqrt((x2_screen - event.x) ** 2 + (y2_screen - event.y) ** 2) for x2_screen, y2_screen in screen_coords]
        closest_index = np.argmin(distances)

        # Set a tolerance to ensure the click is close enough to the atom
        if distances[closest_index] < 20:  # Adjust threshold as needed (in screen pixels)
            return closest_index
        return None


    def update_atom_info(self, atom_index):
        """Update the information displayed when a atom is clicked."""
        atom_index = atom_index #Starts at 0
        atom_info = f"Atom {atom_index + 1} selected at coordinates {self.coordinates[atom_index]}\n"
        self.info_box.delete("1.0", tk.END)  # Clear previous info
        self.info_box.insert(tk.END, atom_info)  # Display new info
        section = self.dropdown_var.get()

        # Highlight the clicked atom in green
      
        self.ax.clear()
        self.ax.scatter(*self.coordinates[atom_index], color="green", s=60, zorder=6, marker='o')
        self.drawn_atom = atom_index
        self.canvas.draw()

        """Draws atoms of intrest"""
        atom_info = self.index_info(atom_index)

        
        
        relavent_atoms = self.itp_currently_shown.DF[section]["atoms"]
        unique_atoms = []
        for l in relavent_atoms:
            if str(atom_index + 1) in l: #Adding one to be on the same  index as itp
                unique_atoms.extend(l)
            
        unique_atoms = set(unique_atoms) #itp index 
        #relavent_atoms = set(relavent_atoms)

        for atom in  unique_atoms:
            atom = int(atom) - 1 #puts in plotting index
            if atom == atom_index:
                continue 
                
            self.ax.scatter(*self.coordinates[atom], color=self.atom_colors[atom], s=50, zorder=5, marker='o')


        """Writing Atom Info Based on Dropdown"""
        
        atom_info = self.index_info(atom_index + 1)
        info = self.itp_currently_shown
        text = f"\nShowing {section} from .itp\n\n"
        self.info_box.insert(tk.END, text)
        
        try: #If a top comment exists it will write it
            self.info_box.insert(tk.END, info.sections_to_comments[section][0])
        except:
            None

        for i in atom_info[section]:
            line = f"   {info.sections_to_pure_data_dic[section][i]}"
            self.info_box.insert(tk.END, line)
        
        
        """Plotting Vectors"""
        if section != "atoms":
            vectors = []
            for i in atom_info[section]:
                #print(info.DF[section]["atoms"][i]) #These are the vectors to make
                pair = info.DF[section]["atoms"][i]
                
                vectors.extend((self.mk_vector(pair)))
            self.vectors = vectors
    
            Entire_itp = True #Set to false if only want to set atoms clicked
            if Entire_itp:
                self.plot_atoms_vectors()
            else:
                for i in range(len(self.vectors)):
                    data = self.vectors[i]
                    start = data[0]
                    vec = data[1]
                    self.ax.quiver(start[0], start[1], start[2], vec[0], vec[1], vec[2], color="red", zorder=4)
                self.canvas.draw()
            
            

      

    def reset_zoom(self):
        """Reset the zoom/translate features."""
        self.toolbar.home()  # Reset view to original zoom

    def go_back(self):
        
        self.info_box.delete("1.0", tk.END)
        self.info_box.insert(tk.END, "Returning...")
        self.ax.clear()
        self.handle_dropdown_selection()
        self.plot_atoms_vectors()

    def handle_dropdown_selection(self):
        selected_option = self.dropdown_var.get()
    
        self.info_box.delete("1.0", tk.END)
        self.info_box.insert(tk.END, f"ITP Section: {selected_option}\n")


        info = self.itp_currently_shown
        if isinstance(self.drawn_atom, np.int64): #if the atom is drawn not all of them drawn_atom will be an integer
            self.info_box.delete("1.0", tk.END)
            atom_info = f"Atom {self.drawn_atom + 1} selected at coordinates {self.coordinates[self.drawn_atom]}\n"
            self.info_box.insert(tk.END, atom_info)  # Display new info

            atom_info = self.index_info(self.drawn_atom + 1) #getting the info about the drawn atom
            
            text = f"\nShowing {selected_option} from .itp\n\n"
            self.info_box.insert(tk.END, text)
            try: #If a top comment exists it will write it
                self.info_box.insert(tk.END, info.sections_to_comments[selected_option][0])
            except:
                None

            for i in atom_info[selected_option]:
                line = f"   {info.sections_to_pure_data_dic[selected_option][i]}"
                self.info_box.insert(tk.END, line)


        else:
            try: #If a top comment exists it will write it
                self.info_box.insert(tk.END, info.sections_to_comments[selected_option][0])
            except:
                None
            
            for i in info.sections_to_pure_data_dic[selected_option]:
                self.info_box.insert(tk.END, i)


        ##Plots the vectors of the section
        self.vectors= self.section_to_vectors()
        self.plot_atoms_vectors()
        


    ################################################################
    "This section is for all functions that manipulate the ITP object"
    ################################################################
    def index_info(self, index: int):
        """
        Finds all the infromation in the ITP regarding the given index (itp index so the atom number)

        Args:
            index (interger): The atom index to look for 

        Returns:
            dict: returns a dictornry of all the locations of the specified index
        """
        indexed_itp = self.itp_currently_shown

        found_dict = {}  # Will store the lines to remove from each section
        for section in indexed_itp.sections:
            if section == "moleculetype":
                continue

            found= []
            for i, data in enumerate(indexed_itp.DF[section]["atoms"]):  
                if str(index) in data:
                    found.append(i) #all the lines to find
            found_dict[section] = found                     
        
    
        return(found_dict)


    def mk_vector(self,index):
        """
        Takes the indexes in the itp file and returns vectors. ie [1,2,3,4] will return vectors for atom indexes [1,2],[2,3]and,[3,4]

        Args:
            index (list of integers): The atom index to build vectors

        Returns:
            np.array: all the vectors ready to draw
        """
        try: #Need to check the section to draw Virtual Sites
            selected_option = self.dropdown_var.get()
        except:
            selected_option = "bond"
        
       
        
        if "virtual site" in selected_option: #Virtual Site vectors are plotted diffrent
            pairs = []
            for i in range(len(index) - 1):
                pair = [int(index[i+1]), int(index[0])]
                pairs.append(pair)
            vectors = []
            for pair in pairs:
                cord_1 = self.coordinates[pair[0]-1]
                cord_2= self.coordinates[pair[1]-1]

                #print(cord_1,cord_2,cord_2 - cord_1 )
                vectors.append([cord_1,cord_2 - cord_1])
            return vectors
    
        pairs = []
        for i in range(len(index) - 1):
            pair = [int(index[i]), int(index[i+1])]
            pairs.append(pair)
        vectors = []
        for pair in pairs:
            cord_1 = self.coordinates[pair[0]-1]
            cord_2= self.coordinates[pair[1]-1]

            #print(cord_1,cord_2,cord_2 - cord_1 )
            vectors.append([cord_1,cord_2 - cord_1])
        return vectors

    
    def section_to_vectors(self):
        try:
            selected_option = self.dropdown_var.get()
        except:
            selected_option = "bond"
        indexed_itp = self.itp_currently_shown
        
        if selected_option == "atoms":
            vectors = []
            return vectors

        vectors = []
        for pair in indexed_itp.DF[selected_option]["atoms"]:
            if pair == []:
                continue
            
            vectors.extend((self.mk_vector(pair)))
        return vectors
    ################################################################
    "This section is for all helper functions not directly interacting with the GUI"
    ################################################################
        
    def atom_type_to_color(self,atom_type):
        """
        Converts an atom type into a RGB color.
        """
        if len(atom_type) > 1:
            atom_type = atom_type[0] #Going to take the first letter of every atom type to get color


        atom_colors_rgba = {
    "H": (1.0, 0, 1.0, 1.0), "C": (0.0, 0.0, 0.0, 1.0), "N": (0.0, 0.0, 1.0, 1.0), "O": (1.0, 0.0, 0.0, 1.0),
    "S": (1.0, 1.0, 0.0, 1.0), "P": (1.0, 0.647, 0.0, 1.0), "F": (0.565, 0.933, 0.565, 1.0), "Cl": (0.0, 1.0, 0.0, 1.0),
    "Br": (0.647, 0.165, 0.165, 1.0), "I": (0.580, 0.0, 0.827, 1.0), "Mg": (0.678, 0.847, 0.902, 1.0), "Fe": (0.545, 0.271, 0.075, 1.0),
    "Zn": (0.663, 0.663, 0.663, 1.0), "Na": (0.0, 0.0, 0.545, 1.0), "K": (0.541, 0.169, 0.886, 1.0), "Ca": (0.133, 0.545, 0.133, 1.0),
    "Al": (0.753, 0.753, 0.753, 1.0), "Si": (0.824, 0.706, 0.549, 1.0), "Ti": (0.412, 0.412, 0.412, 1.0), "Cr": (1.0, 0.549, 0.0, 1.0),
    "Mn": (0.545, 0.0, 0.545, 1.0), "Co": (1.0, 0.078, 0.576, 1.0), "Ni": (0.251, 0.878, 0.816, 1.0), "Cu": (0.722, 0.451, 0.200, 1.0),
    "Pb": (0.467, 0.533, 0.600, 1.0)}
        
        return(atom_colors_rgba.get(atom_type, (0, 0, 0)))


        
                


        
        


    

            


# Main program

if __name__ == "__main__":
    root = tk.Tk()
    #itp_obj = Itp_parser("acn.itp")
    #itp_obj.load_gro("acn.gro")
    #itp_obj = Itp_parser("pf6.itp")
    #itp_obj.load_gro("pf6.gro")
    #itp_obj = Itp_parser("tba.itp")
    #itp_obj.load_gro("tba.gro")
    itp_obj = Itp_parser("Trimer_topology.itp")
    itp_obj.load_gro("Trimer.gro")   
    #itp_obj = Itp_parser("7mer_n.itp")
    #itp_obj.load_gro("7mer_n.gro")   
    gui = GUI(root = root,itp_loader_object = itp_obj)
    root.mainloop()
