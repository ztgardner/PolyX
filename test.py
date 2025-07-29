

from itp_parser import *
#
# Loading in molecule with the path
itp = Itp_parser(File="tba.itp")
#add comments
itp.add_to_comments("bond",[' new ','comments bottom'])
itp.add_to_comments("bond",[' new ','comments bottom'],bottom=True)
#extends the section wit the new line
itp["bond"]=['  18     2     1      0.1090 284512.000 ;aasdfas']

# # Outputs the data of that sections
print(itp["moleculetype"])

# # Setting the variable to the data from that section
moleculetype = itp["moleculetype"]
#
# # Adding to the data in the section, you can add whatever, i just chose to duplicate the data
print(itp["moleculetype"])
print(itp["atoms"]['charge'])
itp.set_charge("charges")
print(itp["atoms"]['charge'])
# # Saving itp with new charge
itp.save_itp("new_tba.itp")





