"""
abaqus python thermal_input.py filename.inp

Creates a file filename_mesh.inp were the state variables
are inserted
"""

import sys

filename = sys.argv[1]

with open(filename, "r") as file:
	content = [line for line in file]

for iline, line in enumerate(content):
	if line[:7] == "*Depvar" and content[iline+2][0] == "*":

		content[iline+1] = "{:7},\n".format(12)

		content.insert(iline+2, "1, Xf, XF\n")
		content.insert(iline+3, "2, XP, XP\n")
		content.insert(iline+4, "3, XB, XB\n")
		content.insert(iline+5, "4, XM, XM\n")
		content.insert(iline+6, "5, XA, XA\n")
		content.insert(iline+7, "6, GSIZE, GSIZE\n")
		content.insert(iline+8, "7, GGROW, GGROW\n")
		content.insert(iline+9, "8, NUCF, NUCF\n")
		content.insert(iline+10, "9, NUCP, NUCP\n")
		content.insert(iline+11, "10, NUCB, NUCB\n")
		content.insert(iline+12, "11, DT700, DT700\n")
		content.insert(iline+13, "12, TMAX, TMAX\n")
		
with open(filename[:-4] + "_mesh.inp", "w") as file:
	file.writelines(content)

