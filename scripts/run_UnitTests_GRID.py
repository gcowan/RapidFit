#!/usr/bin/env python
import os

os.system( 'cd .. && make -j3 lib && cd scripts' )
os.system( 'SetupProject ganga latest' )
os.system( 'ganga ./Unit_Test_Scripts/Ganga_UnitTest_GRID_cFit.py' )
os.system( 'ganga ./Unit_Test_Scripts/Ganga_UnitTest_GRID_sFit.py' )
os.system( 'ganga ./Unit_Test_Scripts/Ganga_UnitTest_GRID_cFit_TOYS.py' )
