#! /bin/bash
export ROOTSYS=/Disk/speyside4/ROOT/root_v5.28.00.Linux-slc5-gcc3.4/root
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
export PATH=${PATH}:${ROOTSYS}/bin
echo  'Running rapidfit with args: ' $@
/phys/linux/s0127440/rapidfit/trunk/bin/fitting $@
