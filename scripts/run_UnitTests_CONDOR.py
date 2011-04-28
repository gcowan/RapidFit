#!/usr/bin/env python
import os

os.system( 'cd .. && make -j3 lib && cd scripts' )
os.system( 'SetupProject Ganga latest' )
os.system( 'ganga ./Unit_Test_Scripts/Ganga_UnitTest_CONDOR_cFit.py' )
os.system( 'ganga ./Unit_Test_Scripts/Ganga_UnitTest_CONDOR_sFit.py' )
os.system( 'ganga ./Unit_Test_Scripts/Ganga_UnitTest_CONDOR_cFit_TOYS.py' )
