import os

# define params
def runPMP(LBOX = "650.000", NROW = "1612", NGRID = "1612", DTHR = "120.000", NORM = "0.0001", nrm= "N1em4", SLOPE= "-1.000", slp = "am100", DMAX = "240.000", RR=1):
	# create outdir
	topdir = "/data2/DATA/eBOSS/PM_runs/PM_C"
	os.chdir(topdir)
	outdir = "PM_Nr"+NROW+"_L"+LBOX+"_g"+NGRID+"_Dgt"+DTHR[:4]+"_"+slp+"_Dlt"+DMAX[:4] +"_"+nrm
	os.mkdir(outdir)

	os.remove("Init.dat")

	f=open("Init.dat", 'w')
	f.write("""Box      =  """ + LBOX  + "\n")
	f.write("""Nrow     =      """ + NROW  + "\n")
	f.write("""Ngrid    =      """ + NGRID  + "\n")
	f.write("""sig8     =     0.828"""  + "\n")
	f.write("""z_init   =   100.000"""  + "\n")
	f.write("""step da  =  4.0000E-04"""  + "\n")
	f.write("""z_final  =     0.700"""  + "\n")
	f.write("""#outputs =         8"""  + "\n")
	f.write("""2.00 1.70 1.40 1.20 1.10 1.00 0.9 0.8 """  + "\n")
	f.write("""dens_thr =   """ + DTHR + "\n")
	f.write("""Vrms     =     0.000"""  + "\n")
	f.write("""#Params  =        10"""  + "\n")
	f.write("""Parametr =   """ + DTHR  + "\n")
	f.write("""Parametr =     """ +NORM + "\n")
	f.write("""Parametr =     """ + SLOPE + "\n")
	f.write("""Parametr =   """ +DMAX + "\n")
	f.write("""Parametr =      0.000"""  + "\n")
	f.write("""Parametr =      0.000"""  + "\n")
	f.write("""Parametr =      0.000"""  + "\n")
	f.write("""Parametr =      0.000"""  + "\n")
	f.write("""Parametr =      0.000"""  + "\n")
	f.write("""Parametr =      0.000"""  )
	f.close()

	os.system("make clean")
	os.system("make PMP2main")
	os.system("make PMP2start ")
	os.system("make PMP2init ")
	os.system("make PMPanalysis")
	os.system("./PMP2init.exe")

	os.chdir(outdir)
	for ii in range(RR):
                os.system("touch begin"+str(ii+1))
                os.system("echo "+str(ii+1)+" | ../PMP2start.exe")
                os.system("echo 1000000 | ../PMP2main.exe")
                os.system("rm -rf PMcr*.DAT")
                os.system("touch fin"+str(ii+1))


runPMP(LBOX="1040.000", NROW="1300", NGRID="2600", DTHR="10.000", NORM="0.0003", nrm= "N3em4", SLOPE= "0.110", slp="ap011", DMAX="400.000", RR = 100)


