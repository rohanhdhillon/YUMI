#This file is for post processing once all the runs are done.


#import the libraries
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.io import FortranFile

def ReadFortranFile(Location):
    '''
    Location is the path of the fortran file.

    The size of the file is expected to be part of the filename.
    '''
    if not(os.path.exists(Location)):
        print("File not found")
        return 0
    Size = int(Location.split("_")[-1][:-4])
    return FortranFile(Location, 'r').read_reals(dtype=np.float64).reshape(Size,Size).T
    

def findRunStatus(RunFolderLocation):
    #This is the list of the directories
    os.chdir("%s" %RunFolderLocation)

    AllDirectories = np.array(glob.glob("CO2*"))
    AllEnergies = np.array([float(Value.split("_")[-1]) for Value in AllDirectories])
    ArrangeIndex = np.argsort(AllEnergies)


    AllEnergies = AllEnergies[ArrangeIndex]
    AllDirectories = AllDirectories[ArrangeIndex]


    TrFilesNumber = []
    TiFilesNumber = []

    for Directory, Energy in zip(AllDirectories, AllEnergies):
        os.chdir(Directory)
        TrFileList = glob.glob("Tr*.dat")
        TiFileList = glob.glob("Tr*.dat")
        TrFilesNumber.append(len(TrFileList))
        TiFilesNumber.append(len(TiFileList))


        os.chdir("..")


    TrFilesNumber = np.array(TrFilesNumber)
    TiFilesNumber = np.array(TiFilesNumber)
    MaxTrFileNumber = np.max(TrFilesNumber)
    print(MaxTrFileNumber)

    GoodRunIndex1 = TrFilesNumber==MaxTrFileNumber
    GoodRunIndex2 = TiFilesNumber==MaxTrFileNumber
    GoodRunIndex = np.logical_and(GoodRunIndex1, GoodRunIndex2)

    print(np.sum(GoodRunIndex1), np.sum(GoodRunIndex2), np.sum(GoodRunIndex2))
    SuccessEnergies = AllEnergies[GoodRunIndex]
    FailureEnergies = AllEnergies[~GoodRunIndex]
    SuccessDirectories = AllDirectories[GoodRunIndex]

    print("The failure energies are given by:", FailureEnergies)
    return FailureEnergies, SuccessEnergies


def plotCrossSection(RunFolderLocation):
    '''

    '''
    os.chdir("%s" %RunFolderLocation)

    #Find and arrange directories

    AllDirectories = np.array(glob.glob("CO2*"))
    print("AllDirec")
    AllEnergies = np.array([float(Value.split("_")[-1]) for Value in AllDirectories])
    ArrangeIndex = np.argsort(AllEnergies)


    AllEnergies = AllEnergies[ArrangeIndex]
    AllDirectories = AllDirectories[ArrangeIndex]

    import matplotlib as mpl
    mpl.rc('font',**{'sans-serif':['Helvetica'], 'size':25,'weight':'bold'})
    mpl.rc('axes',**{'labelweight':'bold', 'linewidth':1.5, 'labelsize':20})
    mpl.rc('ytick',**{'major.pad':22, 'color':'k', 'major.width':1.5,'major.size':12.5})
    mpl.rc('xtick',**{'major.pad':10,'color':'k', 'major.width':1.5,'major.size':12.5})
    mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
    mpl.rc('text', **{'usetex':True})
    mpl.rc('contour', **{'negative_linestyle':'solid'})
    
    
    Data2PlotX = []
    Data2PlotY = []
    for Energy, Directory in zip(AllEnergies, AllDirectories):

        print("The Directory is:", Directory)
        
        #Go to the directory
        os.chdir(Directory)

        SMatrixFile = glob.glob("SMatrix_*.dat")

        Items2Plot = [[0,1], [1,0], [0,2], [2,0]]    

        if len(SMatrixFile)>0:
            SMatrixFile = SMatrixFile[0]
            print("The Sigma  is given by:", SMatrixFile)
            Data2Plot = ReadFortranFile(SMatrixFile)
            Data2PlotX.append(Energy)
            List2Append = []
            for Item in Items2Plot:
                PosX, PosY = Item
                List2Append.append(Data2Plot[PosX,PosY])
            Data2PlotY.append(List2Append)
           
            
        os.chdir("..")

    Data2PlotX = np.array(Data2PlotX)
    Data2PlotY = np.array(Data2PlotY).T
    print(Items2Plot)
    print("The shape of X", np.shape(Data2PlotX))
    print("The shape of Y", np.shape(Data2PlotY))
    input("Wait here...")
    plt.figure(figsize=(12,8))

    colorList = ["maroon", "blue", "green", "cyan", "red", "black", "orange"]
    for Counter, DataY in enumerate(Data2PlotY):
        LabelText = str(Items2Plot[Counter][0])+"$\\rightarrow$"+str(Items2Plot[Counter][1])
        plt.plot(Data2PlotX, DataY, color=colorList[Counter], lw=2, label=LabelText)
    #plt.plot(Data2PlotX, Data2PlotY[:,1], color="blue", lw=2, label="2 $\\rightarrow$ 2")
    #plt.plot(Data2PlotX, Data2PlotY[:,2], color="green", lw=2, label="3 $\\rightarrow$ 1")
    #plt.plot(Data2PlotX, Data2PlotY[:,3], color="cyan", lw=2, label="4 $\\rightarrow$ 2")
    plt.xlim(0,max(Data2PlotX))
    plt.xlabel("Energy [cm$^{\mathrm{-1}}$]", fontsize=30)
    plt.ylabel("Cross-Section [$\mathrm{\AA^2}$]", fontsize=30)
    plt.legend(loc=1)
    plt.tight_layout()
    plt.savefig("../CO2_H2_para_CrossSection.png")
    plt.savefig("../CO2_H2_para_CrossSection.pdf")
    plt.show()
    plt.close()

    pass
    #Check if the 




#Now create the diagnostic plots..
RunFolderLocation = "CO2_H2_para_2023-12-12_16_18_00"
plotCrossSection(RunFolderLocation)