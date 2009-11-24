#! /bin/bash
export ROOTSYS=/Disk/lochnagar0/general/root/5.24.00_slc5_gcc34_bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}
export PATH=${PATH}:${ROOTSYS}/bin
#~/lhcb/beta_s/fitting/Fitting/RapidFit_new/bin/fitting -f ~/lhcb/beta_s/fitting/Fitting/RapidFit_new/config/RaPDFToy_noBkg.xml
~/lhcb/beta_s/fitting/Fitting/RapidFit_new/bin/fitting
