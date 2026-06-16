#!/bin/bash 

echo "Now compiling the files"
module load intel-oneapi/2023.1
echo "Removing previously compiled files."
rm *.o *.x
make
