#! /bin/bash
export ROOTSYS=/Disk/lochnagar0/general/root/5.24.00_slc5_gcc34_bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}
export PATH=${PATH}:${ROOTSYS}/bin
~/rapidfit/trunk/bin/fitting -f ~/rapidfit/trunk/config/RaPDF_Bs2JpsiPhi_sWave_ConfigFile.xml --doPulls PullPlots.root

