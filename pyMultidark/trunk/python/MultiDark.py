
"""
.. class:: MultiDark

.. moduleauthor:: Johan Comparat <johan.comparat__at__gmail.com>

The class MultiDark is a wrapper to handle Multidark simulations results / outputs.

"""
import cPickle
import fileinput
import astropy.cosmology as co
c2 = co.Planck13
from scipy.interpolate import interp1d
from os.path import join

import astropy.units as uu
import numpy as n
import glob

class MultiDarkSimulation :
    """
    Loads the environement proper to the Multidark simulations. All these parameters can only be changed with the setters and obtained with the getters as we do not wnat them to be changed easily. This is the fixed framework of the simulation.
            
    :param __Lbox: length of the box in Mpc/h 
    :param __wdir: Path to the multidark lightcone directory
    :param __boxDir: box directory name
    :param __snl: list of snapshots available
    :param __zsl: list of redshift corresponding to the snapshots   
    :param __zArray: redshift array to be considered to interpolate the redshift -- distance conversion
    :param __Hbox: Hubble constant at redshift 0 of the box
    :param __Melement: Mass of the resolution element in solar masses.   
    """
	__Lbox = 2500.0 * uu.Mpc # Mpc/h
	__wdir = "/data2/DATA/eBOSS/Multidark-lightcones/"
	__boxDir = "MD_2.5Gpc" # can be any of MD_0.4Gpc  MD_1Gpc_new_BDM  MD_1Gpc_new_rockS  MD_1Gpc_old_BDM	MD_2.5Gpc  MD_4Gpc
	__snl =  n.array(glob.glob(join(__wdir , __boxDir , "snapshots/hlist_?.?????.list")))
	__zsl =  None
	__zArray = n.arange(0.2,2.4,1e-1)
	__Hbox = 67.77 * uu.km / (uu.s * uu.Mpc)
	__Melement = 23593750000.0
    __columnDict = {'scale': 0, 'id': 1, 'desc_scale': 2, 'desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Rs_Klypin': 35, 'Mmvir_all': 36, 'M200b': 37, 'M200c': 38, 'M500c': 39, 'M2500c': 40, 'Xoff': 41, 'Voff': 42, 'Spin_Bullock': 43, 'b_to_a': 44, 'c_to_a': 45, 'Ax': 46, 'Ay': 47, 'Az': 48, 'b_to_a_500c': 49, 'c_to_a_500c': 50, 'Ax_500c': 51, 'Ay_500c': 52, 'Az_500c': 53, 'TU': 54, 'M_pe_Behroozi': 55, 'M_pe_Diemer': 56, 'Halfmass_Radius': 57, 'Macc': 58, 'Mpeak': 59, 'Vacc': 60, 'Vpeak': 61, 'Halfmass_Scale': 62, 'Acc_Rate_Inst': 63, 'Acc_Rate_100Myr': 64, 'Acc_Rate_1Tdyn': 65, 'Acc_Rate_2Tdyn': 66, 'Acc_Rate_Mpeak': 67, 'Mpeak_Scale': 68, 'Acc_Scale': 69, 'First_Acc_Scale': 70, 'First_Acc_Mvir': 71, 'First_Acc_Vmax': 72, 'VmaxatMpeak': 73}

	def __init__(self,Lbox,wdir,boxDir,snl,zsl,zArray,Hbox,Melement ):
		self.__Lbox = Lbox # box length
		self.__Hbox = Hbox # Hubble constant at redshift 0 in the box
		self.__wdir = wdir # working directory
		self.__boxDir = boxDir # directory of the box where the snapshots a stored
		self.__snl = snl # snapshot list
		self.__zsl = zsl # corresponding redshift list
		self.__zArray = zArray # redshift for the dC - z conversion
		self.__Melement = Melement # mass of one particle in the box
        self.__columnDict = __columnDict

	def set_Melement(self,Melement):
		self.__Melement = Melement

	def get_Melement(self):
		return self.__Melement 

	def set_zArray(self,zArray):
		self.__zArray = zArray

	def get_zArray(self):
		return self.__zArray 

	def set_Lbox(self,Lbox):
		self.__Lbox = Lbox

	def get_Lbox(self):
		return self.__Lbox 

	def set_Hbox(self,Hbox):
		self.__Hbox = Hbox

	def get_Hbox(self):
		return self.__Hbox 

	def set_wdir(self,wdir):
		self.__wdir = wdir

	def get_wdir(self):
		return self.__wdir 

	def set_boxDir(self,boxDir):
		self.__boxDir = boxDir

	def get_boxDir(self):
		return self.__boxDir 

	def set_snl(self,snl):
		self.__snl = snl

	def get_snl(self):
		return self.__snl 

	def set_zsl(self,zsl):
		self.__zsl = zsl

	def get_zsl(self):
		return self.__zsl 
        
    def set_columnDict(self,columnDict):
        self.__columnDict = columnDict

    def get_columnDict(self):
        return self.__columnDict

	def computeSingleDistributionFunction(self, ii, name, bins ) :
		"""
		Computes the mass, velocity and concentration histograms for a rockstar snapshot.
        
        I oopens a snapshot and reads line by line extractin the quantity of interest.
        
		:param ii: index of the snapshot in the list self.get_snl()
        :param name: name of the quantity of interest, mass, velocity.
        :param index: of the quantity of interest in the snapshots.
        :param bins: binning scheme to compute the historgram.
		"""		
        index = get_columnDict(name)
        output_dir = join(self.get_wdir(),self.get_boxDir(),"properties",name)
        os.system('mkdir '+ output_dir)
        NperBatch = 10000000
		qtyCentral = n.empty(NperBatch)  # 10M array
		qtySat = n.empty(NperBatch)  # 10M array

		fl = fileinput.input(self.get_snl()[ii])
		nameSnapshot = self.get_snl()[ii].split('/')[-1][:-5]

		countCen,countSat,countFileCen,countFileSat = 0,0,0,0
        
		for line in fl:
			if line[0] == "#" :
				continue

			line = line.split()
			sat_or_cen = float(line[5])
			if sat_or_cen ! = -1 :
				countSat+ = 1					
				qtySat[countSat] = float(line[index])
                
			if sat_or_cen = = -1 :
				countCen+ = 1					
				qtyCentral[countCen] = float(line[index])
                
			if countCen == NperBatch-1 :
				nnM,bb = n.histogram(n.log10(qtyCentral),bins = bins)
				print "countCen",countCen
				f = open(join(output_dir, nameSnapshot + "_" + name + "_Central_" + str(countFileCen)+ ".pkl"),'w')
				cPickle.dump(nnM,f)
				f.close()
				countFileCen+ = 1
				countCen = 0

			if countSat == NperBatch-1 :
				nnM,bb = n.histogram(n.log10(qtySat),bins = bins)
				print "countSat", countSat
				f = open(join(output_dir, nameSnapshot + "_" + name+ "_Satellite_" + str(countFileSat)+ ".pkl"),'w')
				cPickle.dump(nnM,f)
				f.close()
				countFileSat+ = 1
				countSat = 0

        # and for the last batch :
        nnM,bb = n.histogram(n.log10(qtyCentral),bins = bins)
        f = open(join(output_dir, nameSnapshot + "_" + name +"_Central_" + str(countFileCen)+ ".pkl"),'w')
        cPickle.dump(nnM,f)
        f.close()

        nnM,bb = n.histogram(n.log10(qtySat),bins = bins)
        f = open(join(output_dir, nameSnapshot + "_" + name + "_Satellite_" + str(countFileSat)+ ".pkl"),'w')
        cPickle.dump(nnM,f)
        f.close()


	def combinesSingleDistributionFunction(self, name ) :
		"""
		Coombines the outputs of computeSingleDistributionFunction. TO BE WRITTEN.
        """
        return 0.


	def computeMassVelocityConcentrationFunction(self,ii) :
		"""
		computes the mass, velocity and concentration histograms for a rockstar snapshot.
		:param ii: index of the snapshot in the list self.get_snl()
		"""
		massB = n.arange(8,16,0.01)
		vcirB = n.arange(0,4.5,0.01)
		concB = n.arange(1,3,0.1)

		NperBatch = 10000000
		mvcCentralMatrix = n.empty((NperBatch,3))  # 1M matrixes
		mvcSatMatrix = n.empty((NperBatch,3))  # 1 M matrixes

		fl = fileinput.input(self.get_snl()[ii])
		name = self.get_snl()[ii].split('/')[-1][:-5]
		countCen,countSat,countFileCen,countFileSat = 0,0,0,0
		for line in fl:
			if line[0] == "#" :
				continue

			line = line.split()
			sat_or_cen = float(line[5])
			if sat_or_cen ! = -1 :
				countSat+ = 1					
				mvcSatMatrix[countSat] = float(line[10]), float(line[16]), float(line[11]) 
                
			if sat_or_cen = = -1 :
				countCen+ = 1					
				mvcCentralMatrix[countCen] = float(line[10]), float(line[16]), float(line[11])
                
			if countCen == NperBatch-1 :
				nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
				nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
				nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
				dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
				dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
				print "countCen",countCen
				f = open(join(self.get_wdir(),self.get_boxDir(),"properties", name+"_MVRmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
				cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
				f.close()
				countFileCen+ = 1
				countCen = 0

			if countSat == NperBatch-1 :
				nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
				nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
				nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
				dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
				dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
				print "countSat", countSat
				f = open(join(self.get_wdir(),self.get_boxDir() ,"properties" , 
name+"_MVRmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
				cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
				f.close()
				countFileSat+ = 1
				countSat = 0

		# and for the last batch :
		nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
		nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
		nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
		dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
		dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
		f = open(join(self.get_wdir(),self.get_boxDir(),"properties",name+ "_MVRmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
		cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
		f.close()

		nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
		nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
		nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
		dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
		dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
		f = open(join(self.get_wdir(),self.get_boxDir(),"properties",name+ "_MVRmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
		cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
		f.close()


	def computeMassVelocityPeakAccRateFunctions(self,ii) :
		"""
		computes the mass, velocity and concentration histograms for a rockstar snapshot.
		:param ii: index of the snapshot in the list self.get_snl()
		"""
		massB = n.arange(8,16,0.01)
		vcirB = n.arange(0,4.5,0.01)
		concB = n.arange(-5e4,5e4+1,1e3)

		NperBatch = 10000000
		mvcCentralMatrix = n.empty((NperBatch,3))  # 1M matrixes
		mvcSatMatrix = n.empty((NperBatch,3))  # 1 M matrixes

		fl = fileinput.input(self.get_snl()[ii])
		name = self.get_snl()[ii].split('/')[-1][:-5]
		countCen,countSat,countFileCen,countFileSat = 0,0,0,0
		for line in fl:
			if line[0] == "#" :
				continue

			line = line.split()
			sat_or_cen = float(line[5])
			if sat_or_cen ! = -1 :
				countSat+ = 1					
				#print mvcSatMatrix[countSat]
				#print line[59], line[61], line[67]
				mvcSatMatrix[countSat] = float(line[59]), float(line[61]), float(line[67]) # check the right indices ... MASS velocity concentration

			if sat_or_cen = = -1 :
				countCen+ = 1					
				#print mvcCentralMatrix[countCen]
				#print line[59], line[61], line[67]
				mvcCentralMatrix[countCen] = float(line[59]), float(line[61]), float(line[67]) # check the right indices ... MASS velocity concentration

			if countCen == NperBatch-1 :
				nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
				nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
				nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
				dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
				dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
				print "countCen",countCen
				f = open(join(self.get_wdir(),self.get_boxDir(),"properties", name+"_MVAmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
				cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
				f.close()
				countFileCen+ = 1
				countCen = 0

			if countSat == NperBatch-1 :
				nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
				nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
				nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
				dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
				dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
				print "countSat", countSat
				f = open(join(self.get_wdir(),self.get_boxDir() ,"properties" , 
name+"_MVAmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
				cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
				f.close()
				countFileSat+ = 1
				countSat = 0

		# and for the last batch :
		nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
		nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
		nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
		dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
		dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
		f = open(join(self.get_wdir(),self.get_boxDir(),"properties",name+ "_MVAmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
		cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
		f.close()

		nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
		nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
		nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
		dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
		dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
		f = open(join(self.get_wdir(),self.get_boxDir(),"properties",name+ "_MVAmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
		cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
		f.close()


