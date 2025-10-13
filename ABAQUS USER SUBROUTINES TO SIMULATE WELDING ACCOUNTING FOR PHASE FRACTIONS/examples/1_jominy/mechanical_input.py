"""
abaqus python mechanical_input.py filename.inp

Creates a file filename_mesh.inp were the state variables
are inserted
"""

import sys

filename = sys.argv[1]

with open(filename, "r") as file:
	content = [line for line in file]

for iline, line in enumerate(content):
	if line[:7] == "*Depvar" and content[iline+2][0] == "*":

		content[iline+1] = "{:7},\n".format(30)

		content.insert(iline+2, "1, THE, THE\n")
		content.insert(iline+3, "2, PEEQ, PEEQ\n")
		content.insert(iline+4, "3, PPEEQ, PPEEQ\n")
		content.insert(iline+5, "4, TPEEQ, TPEEQ\n")
		content.insert(iline+6, "5, PE11, PE11\n")
		content.insert(iline+7, "6, PE22, PE22\n")
		content.insert(iline+8, "7, PE33, PE33\n")
		content.insert(iline+9, "8, PE12, PE12\n")
		content.insert(iline+10, "9, PPE11, PE11\n")
		content.insert(iline+11, "10, PPE22, PE22\n")
		content.insert(iline+12, "11, PPE33, PE33\n")
		content.insert(iline+13, "12, PPE12, PE12\n")
		content.insert(iline+14, "13, TPE11, TPE11\n")
		content.insert(iline+15, "14, TPE22, TPE22\n")
		content.insert(iline+16, "15, TPE33, TPE33\n")
		content.insert(iline+17, "16, TPE12, TPE12\n")
		content.insert(iline+18, "17, EE11, EE11\n")
		content.insert(iline+19, "18, EE22, EE22\n")
		content.insert(iline+20, "19, EE33, EE33\n")
		content.insert(iline+21, "20, EE12, EE12\n")
		content.insert(iline+22, "21, SPD, SPD\n")
		content.insert(iline+23, "22, Y, Y\n")
		content.insert(iline+24, "23, Y0, Y0\n")
		content.insert(iline+25, "24, K, K\n")
		content.insert(iline+26, "25, G, G\n")
		content.insert(iline+27, "26, RHO, RHO\n")
		content.insert(iline+28, "27, KTP, KTP\n")
		content.insert(iline+29, "28, RHOR, RHOR\n")
		content.insert(iline+30, "29, SYF, SYF\n")
		content.insert(iline+31, "30, PEEQT, PEEQT\n")


with open(filename[:-4] + "_mesh.inp", "w") as file:
	file.writelines(content)
