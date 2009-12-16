#! /bin/bash
export ROOTSYS=/Disk/lochnagar0/general/root/5.25.04_gcc41/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib
export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}
export PATH=${PATH}:${ROOTSYS}/bin

#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0.1fixed_ds1.57fixed_Phis0.0368.xml  --doPulls PullPlots.root
#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0fixed_ds0fixed_Phis0.0368.xml
#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0.1fixed_ds0fixed_Phis0.0368.xml
~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0.1fixed_ds1.57fixed_Phis0.0368.xml
#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0.05fixed_ds0fixed_Phis0.0368.xml
#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0.05fixed_ds1.57fixed_Phis0.0368.xml
#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0fixed_ds0fixed_Phis0.5.xml
#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0.1fixed_ds0fixed_Phis0.5.xml
#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0.1fixed_ds1.57fixed_Phis0.5.xml
#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0.05fixed_ds0fixed_Phis0.5.xml
#~/rapidfit/trunk/bin/fitting -f /phys/linux/s0127440/rapidfit/trunk/config/sWave_note/sWave_Rs0.05fixed_ds1.57fixed_Phis0.5.xml

