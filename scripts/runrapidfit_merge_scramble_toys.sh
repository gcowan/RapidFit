#! /bin/bash
export ROOTSYS=/Disk/speyside4/ROOT/root_v5.28.00.Linux-slc5-gcc3.4/root
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
export PATH=${PATH}:${ROOTSYS}/bin
echo  'Running rapidfit with args: ' $@
echo 'merging and scrambling toys afterwards'
/phys/linux/s0127440/rapidfit/trunk/bin/fitting $@
hadd toydataset_merged.root toydataset*.root
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/phys/linux/s0127440/public/simpletools/simpletools_2.0c/lib
/phys/linux/s0127440/public/simpletools/simpletools_2.0c/bin/tuplescrambler toydataset_merged.root dataNTuple 0 toydataset_scrambled.root
