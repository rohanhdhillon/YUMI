#Created by: Prajwal Niraula
#Date: Nov 6, 2023
#Potentials for CO2-H2 were provided by Prof. Laurent Weisenfeld


#importing libraries
import glob
import numpy as np


#different types of CO2+H2 potential considered.
AllPotentials = ["avcbs_u","avtz_u","avqz_u","avqzf12_u"]


def FormatEntry(Array, NumEntryPerLine=10, format="Potential"):
    if format=="Potential":
        Array = np.array([str("%16.15e" %Item) for Item in Array])
    else:
        Array = np.array([str("%3.2f" %Item) for Item in Array])
    Text2Return = ""
    while len(Array)>0:
        LineEntry = " ".join(Array[:NumEntryPerLine])
      
        Array = Array[NumEntryPerLine:]
        Text2Return += LineEntry+"\n"
    return Text2Return

for PotentialType in AllPotentials:

    AllFiles = np.array(glob.glob("all_fits/*%s*.fit" %PotentialType))

    Distance = [Item.split("/")[1].split("_")[0] for Item in AllFiles]
    Distance = np.array([float(Item.replace("p", ".").replace("d","")) for Item in Distance])

    ArrangeIndex = np.argsort(Distance)

    AllFiles = AllFiles[ArrangeIndex]
    Distance = Distance[ArrangeIndex]
   
    SmallestDistance = np.min(Distance)
    MedianDistanceDiff = "%3.2f" %np.median(np.diff(Distance[:10]))
    SaveName = "Potential_%s.pot" %PotentialType
    print("Saving under:", SaveName)
    

    JNumbers = []
    

    #77 is the number of lines 
    for FileCounter, FileName in enumerate(AllFiles):
        FileContent = open(FileName, 'r').readlines()
        
        if FileCounter==0:
            NumLines = len(FileContent)-3
          
            #Initialize the Matrix
            SaveMatrix = np.zeros((len(Distance), NumLines), dtype=np.float64)
            for EntryCounter, Entry in enumerate(FileContent[3:]):
                JNumbers.append(Entry[:12]+"\n")

        #Start Populating the matrix
        for EntryCounter, Entry in enumerate(FileContent[3:]):
            SaveMatrix[FileCounter, EntryCounter] = float(Entry[13:])

  

    DistanceText = FormatEntry(Distance, format="Distance")

    with open(SaveName, 'w') as f:
        FirstLine=PotentialType+"\n"
        f.write(FirstLine)
        SecondLine = "77 27 4 %s\n" %(MedianDistanceDiff)
        f.write(FirstLine)

        #Save the distance
        f.write(DistanceText)

       
        for JCounter, JValues in enumerate(JNumbers):
            AllValues = SaveMatrix[:,JCounter]
            f.write(JValues) 

            AllValues = SaveMatrix[:,JCounter]
            Values2Write = FormatEntry(AllValues, NumEntryPerLine=5, format="Potential")
            f.write(Values2Write) 
