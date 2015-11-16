
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
import os
import astropy.units as uu
import numpy as n
import glob

class MultiDarkSimulation :
    """
    Loads the environement proper to the Multidark simulations. This is the fixed framework of the simulation.
            
    :param Lbox: length of the box in Mpc/h 
    :param wdir: Path to the multidark lightcone directory
    :param boxDir: box directory name
    :param snl: list of snapshots available
    :param zsl: list of redshift corresponding to the snapshots   
    :param zArray: redshift array to be considered to interpolate the redshift -- distance conversion
    :param Hbox: Hubble constant at redshift 0 of the box
    :param Melement: Mass of the resolution element in solar masses.   
    :param columnDict: dictionnary to convert column name into the index to find it in the snapshots
    """

    def __init__(self,Lbox=2500.0 * uu.Mpc, wdir="/data2/DATA/eBOSS/Multidark-lightcones/", boxDir="MD_2.5Gpc", snl=n.array(glob.glob("/data2/DATA/eBOSS/Multidark-lightcones/MD_2.5Gpc/snapshots/hlist_?.?????.list")), zsl=None, zArray=n.arange(0.2,2.4,1e-1), Hbox = 67.77 * uu.km / (uu.s * uu.Mpc), Melement = 23593750000.0, columnDict = {'scale': 0, 'id': 1, 'desc_scale': 2, 'desc_id': 3, 'num_prog': 4, 'pid': 5, 'upid': 6, 'desc_pid': 7, 'phantom': 8, 'sam_mvir': 9, 'mvir': 10, 'rvir': 11, 'rs': 12, 'vrms': 13, 'mmp?': 14, 'scale_of_last_MM': 15, 'vmax': 16, 'x': 17, 'y': 18, 'z': 19, 'vx': 20, 'vy': 21, 'vz': 22, 'Jx': 23, 'Jy': 24, 'Jz': 25, 'Spin': 26, 'Breadth_first_ID': 27, 'Depth_first_ID': 28, 'Tree_root_ID': 29, 'Orig_halo_ID': 30, 'Snap_num': 31, 'Next_coprogenitor_depthfirst_ID': 32, 'Last_progenitor_depthfirst_ID': 33, 'Last_mainleaf_depthfirst_ID': 34, 'Rs_Klypin': 35, 'Mmvir_all': 36, 'M200b': 37, 'M200c': 38, 'M500c': 39, 'M2500c': 40, 'Xoff': 41, 'Voff': 42, 'Spin_Bullock': 43, 'b_to_a': 44, 'c_to_a': 45, 'Ax': 46, 'Ay': 47, 'Az': 48, 'b_to_a_500c': 49, 'c_to_a_500c': 50, 'Ax_500c': 51, 'Ay_500c': 52, 'Az_500c': 53, 'TU': 54, 'M_pe_Behroozi': 55, 'M_pe_Diemer': 56, 'Halfmass_Radius': 57, 'Macc': 58, 'Mpeak': 59, 'Vacc': 60, 'Vpeak': 61, 'Halfmass_Scale': 62, 'Acc_Rate_Inst': 63, 'Acc_Rate_100Myr': 64, 'Acc_Rate_1Tdyn': 65, 'Acc_Rate_2Tdyn': 66, 'Acc_Rate_Mpeak': 67, 'Mpeak_Scale': 68, 'Acc_Scale': 69, 'First_Acc_Scale': 70, 'First_Acc_Mvir': 71, 'First_Acc_Vmax': 72, 'VmaxatMpeak': 73} ):
        self.Lbox = Lbox # box length
        self.Hbox = Hbox # Hubble constant at redshift 0 in the box
        self.wdir = wdir # working directory
        self.boxDir = boxDir # directory of the box where the snapshots a stored
        self.snl = snl # snapshot list
        self.zsl = zsl # corresponding redshift list
        self.zArray = zArray # redshift for the dC - z conversion
        self.Melement = Melement # mass of one particle in the box
        self.columnDict = columnDict

    def computeSingleDistributionFunction(self, ii, name, bins ) :
        """
        Extracts the distribution of quantity 'name' out of all snapshots of the Multidark simulation.        
        :param ii: index of the snapshot in the list self.snl
        :param name: name of the quantity of interest, mass, velocity.
        :param index: of the quantity of interest in the snapshots.
        :param bins: binning scheme to compute the historgram.
        """		
        index = self.columnDict[name]
        output_dir = join(self.wdir,self.boxDir,"properties",name)
        os.system('mkdir '+ output_dir)
        NperBatch = 10000000
        qtyCentral = n.empty(NperBatch)  # 10M array
        qtySat = n.empty(NperBatch)  # 10M array
        print name, index, output_dir

        fl = fileinput.input(self.snl[ii])
        nameSnapshot = self.snl[ii].split('/')[-1][:-5]

        countCen,countSat,countFileCen,countFileSat = 0,0,0,0
        
        for line in fl:
            if line[0] == "#" :
                continue

            line = line.split()
            sat_or_cen = float(line[5])
            if sat_or_cen != -1 :
                countSat+= 1					
                qtySat[countSat] = float(line[index])
                
            if sat_or_cen == -1 :
                countCen+= 1					
                qtyCentral[countCen] = float(line[index])
                
            if countCen == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(qtyCentral),bins = bins)
                print "countCen",countCen
                f = open(join(output_dir, nameSnapshot + "_" + name + "_Central_" + str(countFileCen)+ ".pkl"),'w')
                cPickle.dump(nnM,f)
                f.close()
                countFileCen+= 1
                countCen = 0

            if countSat == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(qtySat),bins = bins)
                print "countSat", countSat
                f = open(join(output_dir, nameSnapshot + "_" + name+ "_Satellite_" + str(countFileSat)+ ".pkl"),'w')
                cPickle.dump(nnM,f)
                f.close()
                countFileSat+= 1
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


    def combinesSingleDistributionFunction(self, ii, name='Vpeak', bins=10**n.arange(0,3.5,0.01), type = "Central" ) :
        """
        Coombines the outputs of computeSingleDistributionFunction.
        :param ii: index of the snapshot
        :param name: name of the quantity studies
        :param bins: bins the histogram was done with
        :param type: "Central" or "Satellite"
        """
        output_dir = join(self.wdir,self.boxDir,"properties",name)
        nameSnapshot = self.snl[ii].split('/')[-1][:-5]
        pklList = n.array(glob.glob(join(output_dir, nameSnapshot + "_" + name +"_"+type+"_*.pkl")))

        nnM = n.empty( [len(pklList),len(bins)-1] ) 
        for jj in range(len(pklList)):
            f=open(pklList[jj],'r')
            nnMinter = cPickle.load(f)
            nnM[jj] = nnMinter
            f.close()

        n.savetxt(join(output_dir,"hist-"+type+"-"+name+"-"+nameSnapshot[6:]+".dat"),n.transpose([bins[:-1], bins[1:], nnM.sum(axis=0)]))


    def computeMassVelocityConcentrationFunction(self,ii) :
        """
        computes the mass, velocity and concentration histograms for a rockstar snapshot.
        :param ii: index of the snapshot in the list self.snl
        """
        massB = n.arange(8,16,0.01)
        vcirB = n.arange(0,4.5,0.01)
        concB = n.arange(1,3,0.1)

        NperBatch = 10000000
        mvcCentralMatrix = n.empty((NperBatch,3))  # 1M matrixes
        mvcSatMatrix = n.empty((NperBatch,3))  # 1 M matrixes

        fl = fileinput.input(self.snl[ii])
        name = self.snl[ii].split('/')[-1][:-5]
        countCen,countSat,countFileCen,countFileSat = 0,0,0,0
        for line in fl:
            if line[0] == "#" :
                continue

            line = line.split()
            sat_or_cen = float(line[5])
            if sat_or_cen != -1 :
                countSat+= 1					
                mvcSatMatrix[countSat] = float(line[10]), float(line[16]), float(line[11]) 
                
            if sat_or_cen == -1 :
                countCen+= 1					
                mvcCentralMatrix[countCen] = float(line[10]), float(line[16]), float(line[11])
                
            if countCen == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
                nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
                nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
                dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
                dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
                print "countCen",countCen
                f = open(join(self.wdir,self.boxDir,"properties", name+"_MVRmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
                cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
                f.close()
                countFileCen+= 1
                countCen = 0

            if countSat == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
                nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
                nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
                dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
                dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
                print "countSat", countSat
                f = open(join(self.wdir,self.boxDir ,"properties" , 
    name+"_MVRmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
                cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
                f.close()
                countFileSat+= 1
                countSat = 0

        # and for the last batch :
        nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
        nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
        nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
        dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
        dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
        f = open(join(self.wdir,self.boxDir,"properties",name+ "_MVRmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
        cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
        f.close()

        nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
        nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
        nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
        dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
        dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
        f = open(join(self.wdir,self.boxDir,"properties",name+ "_MVRmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
        cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
        f.close()


    def computeMassVelocityPeakAccRateFunctions(self,ii) :
        """
        computes the mass, velocity and concentration histograms for a rockstar snapshot.
        :param ii: index of the snapshot in the list self.snl()
        """
        massB = n.arange(8,16,0.01)
        vcirB = n.arange(0,4.5,0.01)
        concB = n.arange(-5e4,5e4+1,1e3)

        NperBatch = 10000000
        mvcCentralMatrix = n.empty((NperBatch,3))  # 1M matrixes
        mvcSatMatrix = n.empty((NperBatch,3))  # 1 M matrixes

        fl = fileinput.input(self.snl[ii])
        name = self.snl[ii].split('/')[-1][:-5]
        countCen,countSat,countFileCen,countFileSat = 0,0,0,0
        for line in fl:
            if line[0] == "#" :
                continue

            line = line.split()
            sat_or_cen = float(line[5])
            if sat_or_cen != -1 :
                countSat+= 1					
                #print mvcSatMatrix[countSat]
                #print line[59], line[61], line[67]
                mvcSatMatrix[countSat] = float(line[59]), float(line[61]), float(line[67]) # check the right indices ... MASS velocity concentration

            if sat_or_cen == -1 :
                countCen+= 1					
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
                f = open(join(self.wdir,self.boxDir,"properties", name+"_MVAmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
                cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
                f.close()
                countFileCen+= 1
                countCen = 0

            if countSat == NperBatch-1 :
                nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
                nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
                nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
                dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
                dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
                print "countSat", countSat
                f = open(join(self.wdir,self.boxDir ,"properties" , 
    name+"_MVAmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
                cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
                f.close()
                countFileSat+= 1
                countSat = 0

        # and for the last batch :
        nnM,bb = n.histogram(n.log10(mvcCentralMatrix.T[0]),bins = massB)
        nnV,bb = n.histogram(n.log10(mvcCentralMatrix.T[1]),bins =  vcirB)
        nnC,bb = n.histogram(n.log10(mvcCentralMatrix.T[2]),bins =  concB)
        dataMC = n.histogram2d(n.log10(mvcCentralMatrix.T[0]), mvcCentralMatrix.T[2] ,bins = [massB,concB])
        dataVC = n.histogram2d(n.log10(mvcCentralMatrix.T[1]), mvcCentralMatrix.T[2] , bins = [vcirB,concB])
        f = open(join(self.wdir,self.boxDir,"properties",name+ "_MVAmatrixCentral_" +str(countFileCen)+ ".pkl"),'w')
        cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
        f.close()

        nnM,bb = n.histogram(n.log10(mvcSatMatrix.T[0]),bins = massB)
        nnV,bb = n.histogram(n.log10(mvcSatMatrix.T[1]),bins =  vcirB)
        nnC,bb = n.histogram(n.log10(mvcSatMatrix.T[2]),bins =  concB)
        dataMC = n.histogram2d(n.log10(mvcSatMatrix.T[0]), mvcSatMatrix.T[2] ,bins = [massB,concB])
        dataVC = n.histogram2d(n.log10(mvcSatMatrix.T[1]), mvcSatMatrix.T[2] , bins = [vcirB,concB])
        f = open(join(self.wdir,self.boxDir,"properties",name+ "_MVAmatrixSatellite_" +str(countFileSat)+ ".pkl"),'w')
        cPickle.dump([nnM,nnV,nnC,dataMC,dataVC],f)
        f.close()


