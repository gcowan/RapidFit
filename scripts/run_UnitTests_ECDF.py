#!/usr/bin/env python
import os

os.system( 'cd .. && make clean && make -j50 lib && cd scripts' )
os.system( 'source /exports/work/physics_ifp_ppe/ganga/install/etc/setup-lhcb.sh 5.4.3' )
os.system( 'ganga ./Unit_Test_Scripts/Ganga_UnitTest_ECDF_cFit.py' )
os.system( 'ganga ./Unit_Test_Scripts/Ganga_UnitTest_ECDF_sFit.py' )
os.system( 'ganga Unit_Test_Scripts/Ganga_UnitTest_ECDF_cFit_TOYS.py' )
