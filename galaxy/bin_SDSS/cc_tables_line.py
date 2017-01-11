import glob
import numpy as n
from os.path import join
import os

all_tabs = n.array([ n.array(glob.glob(str(jj)+"*_line.data")) for jj in n.arange(10, 74, 1) ])
all_names = n.array([str(jj)+"xx_line.data" for jj in n.arange(10, 74, 1) ])

bad_Ncol = []
bad_Ncol_shape = []

for out_name, tab_1 in zip(all_names, all_tabs):
	if os.path.isfile(join("summary",out_name)):
		pass
	else :
		print out_name, len(tab_1)
		if len(tab_1)==1:
			data = n.loadtxt(tab_1[0])
			n.savetxt(join("summary",out_name), data)	
		if len(tab_1)>1:
			data = n.loadtxt(tab_1[0])
			print data.shape
			for tab in tab_1[1:]:
				inter = n.loadtxt(tab)
				print inter.shape
				if len(inter.shape)==2 and inter.shape[1]==229:
					data = n.vstack((data, inter ))

				elif len(inter.shape)==1 and inter.shape[0]==229:
					data = n.vstack((data, inter ))

				else:
					print tab, inter.shape
					bad_Ncol.append(tab)
					bad_Ncol_shape.append(inter.shape[1])

			n.savetxt(join("summary",out_name), data)

n.savetxt(join("summary","wrong_col_N"),n.transpose(bad_Ncol_shape))
