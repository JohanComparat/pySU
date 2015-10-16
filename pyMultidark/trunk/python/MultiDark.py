import cPickle
import fileinput
import astropy.cosmology as co
c2=co.Planck13
from scipy.interpolate import interp1d
from os.path import join

import astropy.units as uu
import numpy as n
import glob

class MultiDark :
	__Lbox = 2500.0 * uu.Mpc # Mpc/h
	__wdir = "/home2/jcomparat/database/Multidark"
	__boxDir = "MD_2.5Gpc"
	__snl =  n.array(glob.glob(join(__wdir , __boxDir , "hlist_?.?????.list")))
	__zsl =  None
	__zArray = n.arange(0.2,2.4,1e-1)
	__Hbox = 67.77 * uu.km / (uu.s * uu.Mpc)
	__Melement = 23593750000.0

	def __init__(self,Lbox,wdir,boxDir,snl,zsl,zArray,Hbox,Melement ):
		self.__Lbox = Lbox # box length
		self.__Hbox = Hbox # Hubble constant at redshift 0 in the box
		self.__wdir = wdir # working directory
		self.__boxDir = boxDir # directory of the box where the snapshots a stored
		self.__snl = snl # snapshot list
		self.__zsl = zsl # corresponding redshift list
		self.__zArray = zArray # redshift for the dC - z conversion
		self.__Melement = Melement # mass of one particle in the box

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

	def computeMassVelocityConcentrationFunction(self,ii) :
		"""
		computes the mass, velocity and concentration histograms for a rockstar snapshot.
		:param ii: index of the snapshot in the list self.get_snl()
		"""
		massB=n.arange(8,16,0.01)
		vcirB=n.arange(0,4.5,0.01)
		concB=n.arange(-5e4,5e4+1,1e3)

		NperBatch = 10000000
		mvcCentralMatrix = n.empty((NperBatch,3))  # 1M matrixes
		mvcSatMatrix = n.empty((NperBatch,3))  # 1 M matrixes

		fl=fileinput.input(self.get_snl()[ii])
		name = self.get_snl()[ii].split('/')[-1][:-5]
		countCen,countSat,countFileCen,countFileSat=0,0,0,0
		for line in fl:
			if line[0]=="#" :
				continue

			line = line.split()
			sat_or_cen = float(line[5])
			if sat_or_cen !=-1 :
				countSat+=1					
				#print mvcSatMatrix[countSat]
				#print line[59], line[61], line[67]
				mvcSatMatrix[countSat] = float(line[59]), float(line[61]), float(line[67]) # check the right indices ... MASS velocity concentration

			if sat_or_cen ==-1 :
				countCen+=1					
				#print mvcCentralMatrix[countCen]
				#print line[59], line[61], line[67]
				mvcCentralMatrix[countCen] = float(line[59]), float(line[61]), float(line[67]) # check the right indices ... MASS velocity concentration

			if countCen==NperBatch-1 :
				nnM,bb=n.histogram(n.log10(mvcCentralMatrix.T[0]),bins=massB)
				nnV,bb=n.histogram(n.log10(mvcCentralMatrix.T[1]),bins= vcirB)
				nnC,bb=n.histogram(n.log10(mvcCentralMatrix.T[2]),bins= vcirB)
				dataMC=n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins=[massB,concB])
				dataVC=n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins=[vcirB,concB])
				print "countCen",countCen
				f=open(join(self.get_wdir(),self.get_boxDir(),"properties", name+"_MVAmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
				cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
				f.close()
				countFileCen+=1
				countCen=0

			if countSat==NperBatch-1 :
				nnM,bb=n.histogram(n.log10(mvcSatMatrix.T[0]),bins=massB)
				nnV,bb=n.histogram(n.log10(mvcSatMatrix.T[1]),bins= vcirB)
				nnC,bb=n.histogram(n.log10(mvcSatMatrix.T[2]),bins= vcirB)
				dataMC=n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins=[massB,concB])
				dataVC=n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins=[vcirB,concB])
				print "countSat", countSat
				f=open(join(self.get_wdir(),self.get_boxDir() ,"properties" , 
name+"_MVAmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
				cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
				f.close()
				countFileSat+=1
				countSat=0

		# and for the last batch :
		nnM,bb=n.histogram(n.log10(mvcCentralMatrix.T[0]),bins=massB)
		nnV,bb=n.histogram(n.log10(mvcCentralMatrix.T[1]),bins= vcirB)
		nnC,bb=n.histogram(n.log10(mvcCentralMatrix.T[2]),bins= vcirB)
		dataMC=n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins=[massB,concB])
		dataVC=n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins=[vcirB,concB])
		f=open(join(self.get_wdir(),self.get_boxDir(),"properties",name+ "_MVAmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
		cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
		f.close()

		nnM,bb=n.histogram(n.log10(mvcSatMatrix.T[0]),bins=massB)
		nnV,bb=n.histogram(n.log10(mvcSatMatrix.T[1]),bins= vcirB)
		nnC,bb=n.histogram(n.log10(mvcSatMatrix.T[2]),bins= vcirB)
		dataMC=n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins=[massB,concB])
		dataVC=n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins=[vcirB,concB])
		f=open(join(self.get_wdir(),self.get_boxDir(),"properties",name+ "_MVAmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
		cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
		f.close()


