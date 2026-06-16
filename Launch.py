import argparse
import os
import socket
import datetime
import numpy as np
import glob
import subprocess


MITSupercloudBaseText="#!/bin/bash\n#SBATCH -c NUMCORES -n 1\n#SBATCH -o Run_%j.log\nsource /etc/profile\nmodule load intel-oneapi/2023.1\nexport OMP_NUM_THREADS=NUMCORES\nexport MKL_NUM_THREADS=NUMCORES\n./yumi.x INPUTFILENAME"
AekalavyaBaseText="#!/bin/bash\nexport OMP_NUM_THREADS=NUMCORES\nexport MKL_NUM_THREADS=NUMCORES\n./yumi.x INPUTFILENAME"
LLMapBaseText="#!/bin/bash\nsource /etc/profile\nmodule load intel-oneapi/2023.1\nexport OMP_NUM_THREADS=NUMCORES\nexport MKL_NUM_THREADS=NUMCORES\ncd $1\n./yumi.x $2"


def Compile(Platform=''):
    #This is the step for compilation and linking it to the library.
    os.chdir("src")
    
    print("Removing the make file")
    os.system("rm Makefile")
    
    if "AEKALAVYA" == Platform:    
        os.system("cp Makefile.Desktop Makefile")
    elif "SUPERCLOUD.MIT" == Platform:
        os.system("cp Makefile.Supercloud Makefile")
    else:
        print("Please add the name of the platform here to run", Platform)
        assert 1==2
    
    #os.system("chmod u+x Makefile")
    os.system("make clean")
    print("Compiling the fortran code.")
    if "AEKALAVYA" == Platform:    
        os.system("make")
    elif "SUPERCLOUD.MIT" == Platform:
        os.system("chmod u+x Compile.sh")
        subprocess.run("./Compile.sh")

    MakeStatus = os.path.exists("yumi.x")

    if MakeStatus != 1:
        print("Could not complete the compilation. Look for compilation errors in fortran.")
        print("Exit here...")
        assert 1 == 2
    else:
        print("*"*20, "Successfully compiled","*"*20,"\n"*3)
    os.chdir("..")



def Run(Platform="", NumCores=8, BaseText="", CompileFlag=True, RunCase="CO2_H2_para", RELAUNCH=False):
    '''
    This code runs the code for the production run.

    Platform is for MIT Supercloud. 
    '''

    if not(RELAUNCH):
        #Create a folder for a run case. 
        print("Launching a fresh run")
        CurrentDatetime = str(datetime.datetime.now()).replace(" ", "_").replace(":","_").split(".")[0]
        SaveFolder = RunCase+"_"+CurrentDatetime
        RunFileLocation = "RunCase/"+RunCase
    else:
      
        RunFileLocation = "RunCase/"+RunCase

        FolderName = RunCase.split("_")[0]
        SaveFolder = glob.glob(RunCase+"*")
        assert len(SaveFolder)==1, "The folder should be unique."
        SaveFolder = SaveFolder[0]

    #Now find the right folder
    try:
        assert os.path.exists("RunCase/"+RunCase)
        print("Located the folder")
    except:
        print(RunFileLocation, " not found...")
        return 0
    
    InputFile = glob.glob(RunFileLocation+"/*.inp")[0].split("/")[-1][:-4]
    BaseText = BaseText.replace("INPUTFILENAME", InputFile).replace("NUMCORES",str(NumCores))

    #Use the following type of command
    if not(os.path.exists(SaveFolder)):
        os.system("mkdir %s" %(SaveFolder))
    

    print("*************************************")
    print("*       Compiling the code.         *")
    print("*************************************")

    if CompileFlag:
        CompilationError = Compile(Platform=Platform)
    else:
        print("Using previously compiled code.")

    #The energy levels for performing the run.   
    AllEnergyLevels = list(np.arange(30,300,2))+list(np.arange(300,1200,5))
    
    #Change the location of the folder.
    os.chdir("%s" %SaveFolder)
    
    NJobs = 0
    
    print("\nLaunching the code in the triple mode.\n")
    for EnergyL in AllEnergyLevels:
        #Launch in the triple mode setting

        #Create subfolder
        CurrentSubFolder = RunCase+"_E_"+str(EnergyL)
            
        #Make the subfolder
        if not(os.path.exists(CurrentSubFolder)):
            os.system("mkdir %s" %CurrentSubFolder)

        if RELAUNCH:
            #Check if the output file is present.
            OutputFiles = glob.glob(CurrentSubFolder+"/*.out")
            if len(OutputFiles)==1:
                OutputFileLocation = OutputFiles[0]
                LastLines = open(OutputFileLocation,'r').readlines()[-6:]
                LastLines = "".join(LastLines).upper()
                if "WALL TIME" in LastLines:  
                    print(" Skipping. ", CurrentSubFolder, " already finished. ")
                    continue
                else:
                    print("The folder", CurrentSubFolder, " is already present.")
                    print("But the run was not completed...")
        print("Launching: ", CurrentSubFolder)
        with open("input.txt", "a") as f:
            f.write(CurrentSubFolder+" "+InputFile+"\n")
        os.system("cp ../%s/* %s" %(RunFileLocation, CurrentSubFolder))

        #Modify the inputfile
        InputFileSubFolder = glob.glob("%s/*.inp" %CurrentSubFolder)[0]
        
        #Compile the compiled file here.    
        os.system("cp ../src/yumi.x %s" %CurrentSubFolder)
              
        with open(InputFileSubFolder,'r') as f:
            InputFileContent = f.readlines()
        
        InputFileContent = "".join(InputFileContent)
        InputFileContent = InputFileContent.replace("ENERGYVAL",str(EnergyL))

        #Now replace the J1Max value
        if EnergyL<100:
            J1MaxVal = 20
        elif EnergyL>100 and EnergyL<200:
            J1MaxVal = 25    
        elif EnergyL>200 and EnergyL<400:
            J1MaxVal = 50
        elif EnergyL>400 and EnergyL<600:
            J1MaxVal = 65
        else:
            J1MaxVal = 75   
        InputFileContent = InputFileContent.replace("J1MAXVAL",str(J1MaxVal))

        with open(InputFileSubFolder,'w') as f:
            f.write(InputFileContent)      
        NJobs+=1


    #Submit the job
    if "SUPERCLOUD.MIT" == Platform:
        os.system("LLsub ./SubmitJob.sh")
    else:
        print("Did not find the platform.")
    
    #write to mapper.sh
    print("Jobs Submitted", NJobs)
    
    
    with open("mapper.sh", 'w') as f:
        f.write(BaseText)
    os.system("chmod u+x mapper.sh")    
    
    NumProcess = 48//NumCores
    #Change this to number of nodes.
    NumNodes = 5
    print("NumNodes is:", NumNodes)
    print("Number of Process is is:", NumProcess)
    Command = "LLMapReduce --mapper mapper.sh --input input.txt --np=[%s,%s,%s] --keep=true " %(NumNodes,NumProcess,NumCores)
    print("The Command is given by:", Command)
    os.system(Command)


def UnitTest(Platform="", NumCores="", CompileFlag=True, BaseText="", UnitTestCase=5):
    '''
    This function performs the unit test for the code.
    Unit Test case 1-6 have been implemented.
    '''
    
    #Run a smallest case possible in 4 cores
    print("\n\n*************************************")
    print("           Performing Unit Test")
    print("*************************************\n\n")
    BaseText=BaseText.replace("NUMCORES",str(NumCores))

    if CompileFlag:
        CompilationError = Compile(Platform=Platform)
    else:
        print("Using previously compiled code.")

    # Get the current date and time
    CurrentDatetime = str(datetime.datetime.now()).replace(" ", "_").replace(":","_").split(".")[0]
    print("The current date time is given by:", CurrentDatetime)
    SaveFolder = "Test_"+CurrentDatetime 

    #Change the directory to Un
    #Command = "cp -r UnitTests/UnitTest2 %s" %SaveFolder
    #os.system(Command)

    #Change the directory to Un
    Command = "cp -r UnitTests/UnitTest%d %s" %(UnitTestCase, SaveFolder)

    if UnitTestCase == 1:
        StartLine, StopLine = 80, 88
    elif UnitTestCase == 2:
        StartLine, StopLine = 126, 134
    elif UnitTestCase == 3:
        StartLine, StopLine = 46, 49
    elif UnitTestCase == 4:
        StartLine, StopLine = 10, 10
    elif UnitTestCase == 7:
        print("Using the test case 7")
        StartLine, StopLine = 80, 88
    elif UnitTestCase == 8:
        print("Using the test case 8")
        StartLine, StopLine = 80, 88
    elif UnitTestCase == 9 or UnitTestCase == 10 or UnitTestCase == 11:
        print("Understanding the effect of NPAS on the values obtained")
        StartLine, StopLine = 92, 109
    else:
        print("Only four test cases has been implemented.")

    
    os.system(Command)

    #Copy the properly compiled YUMI
    os.system("cp -r src/yumi.x %s" %SaveFolder)
    
    #print("Why does not the linking work")
    #os.system("ln -s src/yumi.x %s/yumi.x" %SaveFolder)

  
    #Change the directory to 
    os.chdir("%s" %SaveFolder)
    InputFile = glob.glob("*.inp")[0][:-4]
    print("Running UnitTest case for:", InputFile)
    BaseText = BaseText.replace("INPUTFILENAME", InputFile)

    with open("SubmitJob.sh",'w') as f:
       f.write(BaseText)
       
    #Submit the job
    if "AEKALAVYA" == Platform:    
        os.system("./SubmitJob.sh")
    elif "SUPERCLOUD.MIT" == Platform:
        os.system("LLsub SubmitJob.sh")

    
    OutputFiles = glob.glob("*.out")
    assert len(OutputFiles)==2, "The run should produce an output file which has an extension .out" 

    if "ORIGINAL" in OutputFiles[0].upper():
        File1 = OutputFiles[0]
        File2 = OutputFiles[1]
    else:
        File1 = OutputFiles[1]
        File2 = OutputFiles[0]
   

    #Perform detailed analysis of the output file only in the desktop mode.
    if "AEKALAVYA" == Platform:                
        OriginalFileContent = open(File1, 'r').readlines()[StartLine:StopLine]
        MatrixContent1 = np.genfromtxt(OriginalFileContent)
        
        #Now compare the data. 
        NewFileContent = open(File2, 'r').readlines()[StartLine:StopLine]
        MatrixContent2 = np.genfromtxt(NewFileContent)

        Tolerance = 1e-9
        RelativeDifference =   np.abs(MatrixContent1- MatrixContent2)/(MatrixContent1+Tolerance)*100.0
        MaxDeviation = round(np.max(RelativeDifference),2)
        RunStatus= str(MaxDeviation<1)
        MaxDeviation = str(MaxDeviation)
        MeanDeviation = str(round(np.mean(RelativeDifference),2))
        STDDeviation = str(round(np.std(RelativeDifference),2))

        print("The maximum deviation is:", MaxDeviation)

        #Given the maximum deviation from the expected value and if this should raise the sta
        #Get the distribution of the deviation. Mention what tolerance was used.

        #Run Summary
        SummaryFile = open('SummaryFile.txt','w')
        StatusText = ("Test Passed: "+RunStatus+"\n")
        SummaryFile.write(StatusText)
        SummaryFile.write("Max Deviation: "+ MaxDeviation+"\n")
        SummaryFile.write("Mean Deviation: "+ MeanDeviation+"\n")
        SummaryFile.write("STD Deviation: "+ STDDeviation+"\n")
        SummaryFile.write("Num Cores Used: "+ str(NumCores)+"\n")
        SummaryFile.close()
    else:
        print("Please check the two output files for the performance.")


    
        
def SpeedTest(Platform=None, BaseText=None, CompileFlag=True, NumCores=[1,2,4,6,8,12,24,48], UnitTestCase=2):
    '''
    This functions performs speed test with different number of cores.
    '''

    print("Running the SpeedTest mode")
    #for numCores in NumCores:

    if CompileFlag:
        CompilationError = Compile(Platform=Platform)
    else:
        print("Using previously compiled code.")


    for nCores in NumCores:
        #Copy the files to the 
        CurrentDatetime = str(datetime.datetime.now()).replace(" ", "_").replace(":","_").split(".")[0]
        SaveFolder = "TestNumCores_"+str(nCores)+"_"+CurrentDatetime

        os.system("cp -r UnitTests/UnitTest%s %s" %(str(UnitTestCase), SaveFolder))

        #Copy the properly compiled YUMI
        os.system("cp -r src/yumi.x %s" %SaveFolder)
        
        os.chdir(SaveFolder)
        InputFile = glob.glob("*.inp")[0][:-4]
        print("Running UnitTest case for:", InputFile)
        BaseText = BaseText.replace("INPUTFILENAME", InputFile)
        
        with open("SubmitJob.sh", 'w') as f:
            ModifiedText = BaseText.replace("NUMCORES", str(nCores))
            f.write(ModifiedText)

        if "AEKALAVYA" == Platform:    
            print("Submitted job in Aekalayva...")
            os.system("./SubmitJob.sh")
        elif "SUPERCLOUD.MIT" == Platform:
            os.system("LLsub SubmitJob.sh")  
            print("Submitted job in MIT Supercloud...")  
        
        os.chdir("..")



def SpeedTestLLTripleMode(Platform=None,  CompileFlag=True, BaseText=None, UnitTestCase=1):

    '''
    This functions performs speed test with different number of cores.

    Mode: Can be small or large matrix.

    '''
    

    print("\n\n************************************************************")
    print("           Performing Speed benchmark Test")
    print("************************************************************\n\n")
    
    #Making sure the speed test exists. 1, 2 and 3 are for different size of the matrices.
    assert os.path.exists("UnitTests/SpeedTest"+str(UnitTestCase)), "The speed test %s should be present." %str(UnitTestCase)
    SourceFolder = "UnitTests/SpeedTest"+str(UnitTestCase)
     

    print("Running the SpeedTest mode")
    #for numCores in NumCores:

    print("First Compiling the code.")
    if CompileFlag:
        CompilationError = Compile(Platform=Platform)
    else:
        print("Using previously compiled code.")

    #Different combination of  Cores per Process, Number of Processe Node is always going to be 1.
    PossibleCombinations = [[1,48],[2,24],[4,12],[6,8],[8,6],[12,4],[24,2],[48,1],
                            [1,96],[2,48],[4,24],[6,16],[8,12],[12,8],[24,4],[48,2]]
    
    #Look at the possible combinations of values here.
    #PossibleCombinations = [[4,12],[6,8],[8,6],[12,4]]
    #Create a subfolder within folder to do this test.
    CurrentDatetime = str(datetime.datetime.now()).replace(" ", "_").replace(":","_").split(".")[0]
    TopFolder = "speedbenchMark"+"_Case"+str(UnitTestCase)+"_"+CurrentDatetime
    print("Running speed benchmarks under:", TopFolder)
        
    os.system("mkdir %s" %TopFolder)
    os.chdir(TopFolder)

    
    for CaseCounter, Cases in enumerate(PossibleCombinations): 
        if CaseCounter==0:
            print("Launching job for :", Cases[0], "x", Cases[1],"\n")
        
        print("Job Number:", CaseCounter+1)
     

        #Now create subfolder for each case
        SubFolderName = "Case_"+str(Cases[0]).zfill(3)+"_"+str(Cases[1]).zfill(3)
        os.system("mkdir %s" %(SubFolderName))
             
        NumCores=str(Cases[0])
        TempBaseText = BaseText
        TempBaseText = TempBaseText.replace("NUMCORES",NumCores)
        
        #Change directory to the case running
        os.chdir(SubFolderName)

        #Make the twice the number of processes
        for NumProcess in range(1,Cases[1]+1):
            print("Submitting Job:", NumProcess)
            CaseFolder = "Case_"+str(NumProcess).zfill(3)
            os.system("mkdir %s" %CaseFolder)
            os.system("cp -r ../../%s/* %s" %(SourceFolder, CaseFolder))
            os.system("cp -r ../../src/yumi.x %s" %CaseFolder)    
            #Do you need to modify the code.

            #For each subfolder create a mapper.sh and input.txt file

            #get the name of the input file
            os.chdir(CaseFolder)
            InputFileName = glob.glob("*.inp")[0][:-4].strip()
            os.chdir("..")
            with open("input.txt", "a") as f:
                Text2Write = CaseFolder+" "+"%s\n" %InputFileName
                f.write(Text2Write)
            
        with open("mapper.sh", "a") as f:
            f.write(TempBaseText)
            
            
        os.system("chmod u+x mapper.sh")        
        #This is the command for running the code
        Command = "LLMapReduce --mapper mapper.sh --input input.txt --np=[%s,%s,%s] --keep=true " %(1,NumProcess,NumCores)
        os.system(Command)
        os.chdir("..")


        
        
        




def main():
    parser = argparse.ArgumentParser(description="Perform different operations.")

    # Define the mode as a choice
    parser.add_argument(
        "--mode",
        choices=["run", "unittest", "speedtest", "relaunch"],
        required=True,
        help="Choose a mode: Run, Test, or SpeedTest"
    )

    # Allow the choice for the number of cores.
    parser.add_argument(
        "--NumCores",
        type=int,
        default=None,
        required=False,
        help="Choose a mode: run, unittest, or speedTest"
    )

    #Define the test case to run. The unit tests are provided within the folder UnitTests
    parser.add_argument(
        "--testcase", 
        type=int, 
        default=1,
        required=False,
        help="Select which unit case to run.")

    #Define the test case to run. The unit tests are provided within the folder UnitTests
    parser.add_argument(
        "--compile", 
        type=int, 
        default=1,
        required=False,
        help="Select which unit case to run.")
    
    args = parser.parse_args()
    RunMode = args.mode
    TestCase = args.testcase
    NUMCORES = args.NumCores

    if args.compile == 0: 
      CompileFlag = False
    elif args.compile == 1: 
      CompileFlag = True
    else:
      assert 1==2, "The compile flag should be either 0 (False) or 1 (True). 1 is default."
   
    
    if not(NUMCORES):
        NUMCORES = 8

    
    Directory = socket.gethostname()+ socket.getfqdn()
    Platform = "NOT_RECOGNIZED"

    #Find which platform you are running in
    if "AEKALAVYA" in Directory.upper():
        print("Running on Aekalavya")
        Platform = "AEKALAVYA"
        BaseText = AekalavyaBaseText
        
    elif "SUPERCLOUD.MIT" in Directory.upper() or "TX-GREEN" in Directory.upper():
        print("Running on MIT Supercloud or Tx-Green")
        Platform = "SUPERCLOUD.MIT"
        BaseText =  MITSupercloudBaseText
    
    else:
        print("Please add the name of the platform here to run")
        assert 1==2


    #Run the code in the proper mode.
    if RunMode.upper() == "RUN":
        Run(Platform=Platform,  BaseText=LLMapBaseText, NumCores=NUMCORES, CompileFlag=CompileFlag, RELAUNCH=False)
    elif RunMode.upper() == "UNITTEST":
        #Use 12 cores regardless

        if "AEKALAVYA" == Platform:
            if not(NUMCORES):
                NUMCORES=12
        elif "SUPERCLOUD.MIT" == Platform:
            if not(NUMCORES):
                NUMCORES = 48
        UnitTest(Platform=Platform, NumCores=NUMCORES, BaseText=BaseText, CompileFlag=CompileFlag, UnitTestCase=TestCase)
    elif RunMode.upper() == "SPEEDTEST":
        print("Launching in the speed test case.")
        #SpeedTest(Platform=Platform,  BaseText=BaseText, UnitTestCase=TestCase)
        SpeedTestLLTripleMode(Platform=Platform,  BaseText=LLMapBaseText, CompileFlag=CompileFlag, UnitTestCase=TestCase)
    elif RunMode.upper()=="RELAUNCH":
        print("Run only the cases which are not completed.")
        Run(Platform=Platform,  BaseText=LLMapBaseText, NumCores=NUMCORES, CompileFlag=CompileFlag,  RELAUNCH=True)
    else:
        print("Please add the name of the mode to run")
        assert 1==2

if __name__ == "__main__":
    main()
