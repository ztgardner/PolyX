from jsonschema.validators import extend

from itp_parser import *
#
# Loading in molecule with the path
itp = Itp_parser(File="polymer.itp")
# Will print out what sections are in the itp
print(itp)
# Outputs the data of that sections minus comments
print(itp["moleculetype"])
# Outputs the data of that section plus comments
print(itp.raw("moleculetype"))
# Setting the variable to the data from that section
moleculetype = itp["moleculetype"]
# Adding to the data in the section, you can add whatever, i just chose to duplicate the data
itp.add_to_section("moleculetype",moleculetype)
print(itp["moleculetype"])
# Saving the itp, it saves the comments in each section too, but all of the comments are put at the top of the section
itp.save_itp("demo.itp")
# Changing the charges
new_charge = [
    *range(332),
]
itp.set_charge(new_charge)
# Saving itp with new charge
itp.save_itp("demo2.itp")




