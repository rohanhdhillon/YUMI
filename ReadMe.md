# YUMI v2.0.0

## Version

2.0.0 - Current working stable version.

1.0.0 - SauveGage version which was not optimized, but was used for most of the research works prior to this work.


## Introduction
YUMI code was developed Prof. Nejmeddine Jaidane (At Tunisia) and Prof. Laurent Weisenfeld (CNRS), while being further optimized by Prajwal Niraula (MIT). This codes uses scattering matrices approach to calculating the pressure broadening parameters. The first use case particularly inspired by the exoplanetary needs will be CO2+H2. A number of exoplanetary atmospheres including WASP-39 b will be used as a testcase. And this computation will be largely carried out in MIT-Supercloud/MIT Lincoln supercomputers.


## Contact Us
If you have any questions in regards to the performance of YUMI Code, please contact:
- Prof. Laurent Weisenfeld: Laurent.wiesenfeld@universite-paris-saclay.fr 
- Prajwal Niraula: pniraula@mit.edu


## Output Files
There are four types of files that are produced by YUMI. Some of these allows us to access, others are used to save all the data that down the road potentially may be used on 
- T-Matrix is the transfer matrix, which can be used to calculate the Scattering matrix. 
- output.dat contains information on the convergence test.
- *.out file contains time information it took to run the code for different matrix sizes.
- SigmaMatrix is the S matrix that are calculated by summing the T-Matrix.

Post processing will be part of Paper 2, and will be accordingly added here which is necessary for getting cross-sections from T-Matrices.


## Things to perform
  1. Save SMatrix. [Check if this works since the convergence test.]
  2. Test if the useCholesky is working properly, particularly once it fails if still works.
  3. Do further convergence test at lower energy [J larger than 5]. start value of J to be larger than 5. Look at the potential effects of NPAS value. Make JTotal smaller for J Matrix... [Testing the microreversibility] Symmetry of TMatrix, NPas, and JTotal, J1Max. One of this could be wrong, too small? 
  4. Add plotting to the post-processing, which produces the reports for the runs and graphs for the speed tests.
  5. Once the impact of NPAS is understood, maybe modify the run version such that when NPAS and J are accordingly modified for different energies. [Implemented for J. Need to check for the impact of NPAS on the result]
  6. If the cholesky version does not work add  turn it off, and use alternative. Use it whenever it can be.  [Should be easy to implement in c0003.f90]
  7. Look at a faster inverse algorithm. Also look into GPU based code perhaps using CUDA or MAGMA for matrix, when large matrices are involved.
  8. Generate the performance graph like: https://arxiv.org/pdf/2209.05725.pdf


## Compiling the code.
  1. Compile the code using proper flags and library. (Required OpenMP, Intel MKL Library, LAPACK & BLAS). The flags are very important for optimal performance of the code.
  2. Make sure you have the potential. The ab-initio potentials are pre-calculated.
  3. Make sure the YUMI is properly compiled.
  4. Make sure the input parameter which feeds total angular momentum, J1, J2, and the name of the potential files is properly fed to the YUMI. More details on calculating the potential and its formatting will be provided later.


## Running the code.
  There are multiple modes to run the code. You need numpy, python library to run the python wrapper built around YUMI. The make files are named as Makefile.Desktop and Makefile.Supercloud. The python wrapper recognizes the platform, and compiles the files accordingly. The flags are extremely important for compiling the code, which makes big difference. 

  1. run - This is run different cases such as CO2+H2 
  ```
  >>> python Launch.py --mode run
  ```
  2. unittest - These are cases that have been previously tested and validated. Make sure the OriginalOutput.out and the newoutputFile are the same. 
  ```
  >>> python Launch.py --mode unittest
  ```
  3. speedtest - This runs with different number of cores and see which is the potentially the best combination for number of cores. This will be modified to create the performance graph that will go into the Paper II.
  ```
  >>> python Launch.py --mode speedtest
  ```
  4. relaunch - This runs with running the code in relaunch mode. This is particularly helpful to relaunch the run to complete the cases which have not been completed.   
  ```
  >>> python Launch.py --mode relaunch
  ```
  5. If you do not want to recompile the code use --compile flag as following:
  ```
  >>> python Launch.py --mode unittest --testcase 1 --compile 0
  ```

## Tests Performed:
  1. Benchmarked the YUMI with HCO+ that has been fully tested. This means potential fit, input, output agree with the molscat computation within a few percent. This test is used for unittest. Unittest1 is the small case which can be run by command:
  ```
  >>> python Launch.py --mode unittest --testcase 1
  ```
  2. Unittest2 is the small case which can be run using command 
  ```
  >>> python Launch.py --mode unittest --testcase 2
  ```
  3. Run the unittest case 3 to test more than single channel are open.
  ```
  >>> python Launch.py --mode unittest --testcase 3
  ```
  
  
## Potential files included
  The potential for CO2+H2 was performed by Prof. Laurent, and it is now included here. This is a tedious task, and different versions of the potential have been calculated some of which are expected to have better precision over others. These values are converted to the form that YUMI uses using ConvertPotential.py program included in the folder. 

