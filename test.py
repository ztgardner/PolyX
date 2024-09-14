from parser import *
# Loading in molecule with the path
itp=Itp_parser(File="polymer.itp")
# Will print out what sections were in the itp
print(itp)
#Outputs the data of that sections minus comments
print(itp["moleculetype"])
#Ouputs the data of that section plus comments
print(itp.raw("moleculetype"))
# Setting the variable to the data from that section
moleculetype=itp["moleculetype"]
# Adding to the data in the section, you can add whatever, i just chose to duplicate the data
itp["moleculetype"].extend(moleculetype)
# Saving the itp, it saves the comments in each section too, but all of the comments are put at the top of the section
itp.save_itp("demo.itp")