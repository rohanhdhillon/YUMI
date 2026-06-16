import numpy as np
import matplotlib.pyplot as plt
import os
import glob


#Now add the formatter of matplotlib
import matplotlib as mpl
mpl.rc('font', family='sans-serif', size=25)
mpl.rc('font', serif='Helvetica Neue')
mpl.rc('font', serif='Skia')
mpl.rc('text', usetex='True')
mpl.rc('ytick',**{'major.pad':5, 'color':'black', 'major.size':11,'major.width':1.5, 'minor.size':5,'minor.width':0.75})
mpl.rc('xtick',**{'major.pad':5, 'color':'black',  'major.size':11,'major.width':1.5, 'minor.size':5,'minor.width':0.75})
#mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('axes',**{'linewidth':1.0,'edgecolor':'black'})


XValues = np.array([1,2,3,4,5,6,7,8])
XLabelsText = np.array(['1$\\times$48', '2$\\times$24', '4$\\times$12', '6$\\times$8', '8$\\times$6', '12$\\times$4', '24$\\times$2', '48$\\times$1'])
NumJobs = np.array([48, 24, 12, 8, 6, 4, 2, 1])

NumCores = np.array([1,2,4,6,8,12,24,48])
print(NumCores)
input("Wait here...")
##Matrix Size 200
#001_048 Number of cores used: 1 Time: 6.501 pm 0.175
#002_024 Number of cores used: 2 Time: 3.813 pm 0.152
#004_012 Number of cores used: 4 Time: 2.199 pm 0.069
#006_008 Number of cores used: 6 Time: 1.815 pm 0.046
#008_006 Number of cores used: 8 Time: 1.689 pm 0.077
#012_004 Number of cores used: 12 Time: 1.432 pm 0.015
#024_002 Number of cores used: 24 Time: 1.202 pm 0.006
#048_001 Number of cores used: 48 Time: 1.293 pm 0.0
Case1 = np.array([[6.501, 0.175],
                  [3.813, 0.152],
                  [2.199, 0.069],
                  [1.815, 0.046],
                  [1.689, 0.077],
                  [1.432, 0.015],
                  [1.202, 0.006],
                  [1.293, 0.0]
                ])

#Matrix Size 500
#001_048 Number of cores used: 1 Time: 65.383 pm 0.857
#002_024 Number of cores used: 2 Time: 37.887 pm 1.013
#004_012 Number of cores used: 4 Time: 23.255 pm 0.76
#006_008 Number of cores used: 6 Time: 19.495 pm 0.425
#008_006 Number of cores used: 8 Time: 17.596 pm 0.298
#012_004 Number of cores used: 12 Time: 15.968 pm 0.06
#024_002 Number of cores used: 24 Time: 13.943 pm 0.026
#048_001 Number of cores used: 48 Time: 14.118 pm 0.0
Case2 = np.array([[65.383, 0.857],
                  [37.887, 1.013],
                  [23.255, 0.76],
                  [19.495, 0.425],
                  [17.596, 0.298],
                  [15.968, 0.06],
                  [13.943, 0.026],
                  [14.118, 0.0]
                 ])


#Matrix Size 1200
#001_048 Number of cores used: 1 Time: 297.812 pm 3.599
#002_024 Number of cores used: 2 Time: 174.535 pm 2.403
#004_012 Number of cores used: 4 Time: 113.243 pm 1.362
#006_008 Number of cores used: 6 Time: 98.618 pm 1.259
#008_006 Number of cores used: 8 Time: 89.374 pm 0.351
#012_004 Number of cores used: 12 Time: 80.195 pm 0.433
#024_002 Number of cores used: 24 Time: 73.246 pm 0.034
#048_001 Number of cores used: 48 Time: 74.119 pm 0.0
Case3 = np.array([[297.812, 3.599],
                  [174.535, 2.403],
                  [113.243, 1.362],
                  [98.618, 1.259],
                  [89.374, 0.351],
                  [80.195, 0.433],
                  [73.246, 0.034],
                  [74.119, 0.0]
                 ])


#Matrix Size 4000
#001_048 Number of cores used: 1 Time: 3516.744 pm 101.846
#002_024 Number of cores used: 2 Time: 2406.197 pm 84.067
#004_012 Number of cores used: 4 Time: 1321.561 pm 119.977
#006_008 Number of cores used: 6 Time: 832.298 pm 10.97
#008_006 Number of cores used: 8 Time: 631.838 pm 1.19
#012_004 Number of cores used: 12 Time: 436.649 pm 7.753
#024_002 Number of cores used: 24 Time: 245.46 pm 2.322
#048_001 Number of cores used: 48 Time: 179.728 pm 0.0
Case4 = np.array([[3516.744, 101.846],
                  [2406.197, 84.067],
                  [1321.561, 119.977],
                  [832.298, 10.97],
                  [631.838, 1.19],
                  [436.649, 7.753],
                  [245.46, 2.322],
                  [172.81, 3.506]
                ])


#Matrix Size 10000
#002_024 Number of cores used: 2 Time: 462.526 pm 0.0 ##Cannot 
#004_012 Number of cores used: 4 Time: nan pm nan
#006_008 Number of cores used: 6 Time: 7989.823 pm 0.0
#008_006 Number of cores used: 8 Time: 5984.354 pm 0.0
#012_004 Number of cores used: 12 Time: 4294.288 pm 0.0
#024_002 Number of cores used: 24 Time: 2848.946 pm 0.0
#048_001 Number of cores used: 48 Time: 1527.838 pm 0.0

Case5 = np.array([[35142.64,1241.027],
                 [20340.9585,288.314],   
                 [10525.543,129.22],
                 [7305.302, 83.89],
                 [5478.124, 104.456],
                 [3946.70925, 18.33], #This is for 012_004
                 [2667.0052, 138.90],
                 [1436.56, 5.95],
                 ])


####Using Amdahl's Law to calculate the value of P
#S = 1/((1-P)+P/N)
#S*((1-P)+P/N) = 1
#S*(1-P)+SP/N = 1
#S-SP+SP/N = 1
#S-1 = SP(1-1/N)
#P = (S-1)/(S(1-1/N)) 


N = NumCores[1:]
T0 = Case1[0,0]
T = Case1[1:,0]
S = T0/T
P = (S-1)/(S*(1-1/N))
print("Case 1")
print("Speed-Up: ", S)
print(np.mean(P), np.std(P))

T0 = Case2[0,0]
T = Case2[1:,0]
S = T0/T
P = (S-1)/(S*(1-1/N))
print("Case 2")
print("Speed-Up: ", S)
print(np.mean(P), np.std(P))

T0 = Case3[0,0]
T = Case3[1:,0]
S = T0/T
P = (S-1)/(S*(1-1/N))
print("Case 3")
print("Speed-Up: ", S)
print(np.mean(P), np.std(P))


T0 = Case4[0,0]
T = Case4[1:,0]
S = T0/T
P = (S-1)/(S*(1-1/N))
print("Case 4")
print("Speed-Up: ", S)
print(np.mean(P), np.std(P))

T0 = Case5[0,0]
T = Case5[1:,0]
S = T0/T
P = (S-1)/(S*(1-1/N))
print("Case 5")
print("Speed-Up: ", S)
print(np.mean(P), np.std(P))


#Now we can plot the performance of the code
plt.figure(figsize=(12,8))
plt.errorbar(XValues,Case1[:,0], yerr=Case1[:,1], linestyle='-', marker='o', markersize=10, color='maroon', label='Matrix Size 200')
plt.errorbar(XValues,Case2[:,0], yerr=Case2[:,1], linestyle=':', marker='d', markersize=10, color='blue', label='Matrix Size 500')
plt.errorbar(XValues,Case3[:,0], yerr=Case3[:,1], linestyle='-.', marker='p', markersize=10, color='cyan', label='Matrix Size 1200')
plt.errorbar(XValues,Case4[:,0], yerr=Case4[:,1], linestyle='--', marker='x', markersize=10, color='goldenrod', label='Matrix Size 4000')
plt.errorbar(XValues,Case5[:,0], yerr=Case5[:,1], linestyle=(0, (3, 5, 1, 5, 1, 5)), marker=4, markersize=10, color='green', label='Matrix Size 10000')
plt.xlabel("Performance Case [" + "Cores"+" $\\times$ "+"\#Jobs]", labelpad=30, fontsize=30)
plt.ylabel("Total Time Taken [s]", fontsize=30)
plt.yscale('log')
YLIM = plt.gca().get_ylim()
plt.ylim(YLIM[0], YLIM[1]*20)
plt.xticks(range(1,9), XLabelsText, rotation=90, fontsize=30)
plt.legend(loc='upper right', ncols=2, fontsize=20)
plt.tight_layout()
plt.savefig("PerformanceGraph.pdf")
plt.show()


print("The shape of Case1 is given by: ", Case1.shape)
plt.figure(figsize=(12,8))
plt.errorbar(XValues,Case1[:,0]*NumCores, yerr=Case1[:,1], linestyle='-', marker='o', markersize=10, color='maroon', label='Matrix Size 200')
plt.errorbar(XValues,Case2[:,0]*NumCores, yerr=Case2[:,1], linestyle=':', marker='d', markersize=10, color='blue', label='Matrix Size 500')
plt.errorbar(XValues,Case3[:,0]*NumCores, yerr=Case3[:,1], linestyle='-.', marker='p', markersize=10, color='cyan', label='Matrix Size 1200')
plt.errorbar(XValues,Case4[:,0]*NumCores, yerr=Case4[:,1], linestyle='--', marker='x', markersize=10, color='goldenrod', label='Matrix Size 4000')
plt.errorbar(XValues,Case5[:,0]*NumCores, yerr=Case5[:,1], linestyle=(0, (3, 5, 1, 5, 1, 5)), marker=4, markersize=10, color='green', label='Matrix Size 10000')
plt.xlabel("Performance Case [" + "Cores"+" $\\times$ "+"\#Jobs]", labelpad=30, fontsize=30)
plt.ylabel("Time x \# Cores [s]", fontsize=30)
plt.xticks(range(1,9), XLabelsText, rotation=90, fontsize=30)
plt.yscale('log')
YLIM = plt.gca().get_ylim()
plt.ylim(YLIM[0], YLIM[1]*50)
plt.legend(loc='upper right', ncols=2, fontsize=20)
plt.tight_layout()
plt.savefig("PerformanceNormalized.pdf")
plt.show()
