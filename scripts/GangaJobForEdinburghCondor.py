#Ganga# File created by Ganga - Tue Oct  6 22:34:51 2009
#Ganga#
#Ganga# Object properties may be freely edited before reloading into Ganga
#Ganga#
#Ganga# Lines beginning #Ganga# are used to divide object definitions,
#Ganga# and must not be deleted

#Ganga# Job object (category: jobs)

################################################
# Submission to Condor of multiple toy studies #
################################################
# Note: For submission to Condor to work, you must be on a PPE machine 
# that has access to Condor client. jura.ph.ed.ac.uk is a good example.
# Usage:
# > setupLHCb
# > GangaEnv
# > ganga
# > load( 'path/to/GangaJobForEdinburghCondor.py' )
# > jobs[-1].submit()
#

#numberOfSubjobs = 10
#subjobList = [['%s'% i] for i in range( numberOfSubjobs )]
#pathToExecutable = '/Home/gcowan1/lhcb/beta_s/fitting/Fitting/RapidFit_new/scripts/runFit.sh'

Job (
 name = 'RapidFit Toy Study - mass' ,
 outputsandbox = ['PullPlots.root'] ,
 info = JobInfo (
    ) ,
 inputdata = None ,
 merger = RootMerger (
    files = ['PullPlots.root'] ,
    args = None ,
    ignorefailed = False ,
    overwrite = False 
    ) ,
 inputsandbox = [ ] ,
 application = Executable (
    exe = '/Home/gcowan/lhcb/beta_s/fitting/Fitting/RapidFit/scripts/runFit.sh' , #change to your home dir
    env = {} ,
    args = [] 
    ) ,
 outputdata = None ,
 splitter = ArgSplitter (
    args = [['%s'% i] for i in range( 10 )]
    ) ,
 backend = Condor (
    submit_options = [] ,
    rank = 'Memory' ,
    getenv = 'False' ,
    globusscheduler = '' ,
    env = {'DPM_HOST': 'srm.glite.ecdf.ed.ac.uk', 'RFIO_PORT': '5001', 'X509_CERT_DIR': '/Disk/lochnagar0/grid/glite/external/etc/grid-security/certificates', 'DPNS_HOST': 'srm.glite.ecdf.ed.ac.uk', 'BASH_ENV': '/Disk/lochnagar0/lhcb/lhcb-soft/scripts/lhcb-condorsetup.sh', 'LD_LIBRARY_PATH': '/Disk/lochnagar0/lhcb/lhcb-soft/cern/usr/lib:/Disk/lochnagar0/grid/glite/glite/lib:/Disk/lochnagar0/grid/glite/lcg/lib:/Disk/lochnagar0/grid/glite/globus/lib:/Disk/lochnagar0/lhcb/lhcb-soft/lcg/external/castor/2.1.1-9/slc4_ia32_gcc34/usr/lib'} ,
    universe = 'vanilla' ,
    globus_rsl = '' ,
    shared_filesystem = True ,
    requirements = CondorRequirements (
       opsys = 'LINUX' ,
       machine = '' ,
       other = [] ,
       excluded_machine = 'moyles.ph.ed.ac.uk brambles.ph.ed.ac.uk vance.ph.ed.ac.uk goodier.ph.ed.ac.uk fabio.ph.ed.ac.uk hobbs.ph.ed.ac.uk lowe.ph.ed.ac.uk westwood.ph.ed.ac.uk evans.ph.ed.ac.uk whiley.ph.ed.ac.uk blackburn.ph.ed.ac.uk peterson.ph.ed.ac.uk tong.ph.ed.ac.uk anniemac.ph.ed.ac.uk fearne.ph.ed.ac.uk mills.ph.ed.ac.uk goldfinger.ph.ed.ac.uk reggie.ph.ed.ac.uk bethan.ph.ed.ac.uk lamacq.ph.ed.ac.uk nelson.ph.ed.ac.uk pearce.ph.ed.ac.uk kwame.ph.ed.ac.uk marsden.ph.ed.ac.uk davies.ph.ed.ac.uk green.ph.ed.ac.uk nightingale.ph.ed.ac.uk kay.ph.ed.ac.uk freeman.ph.ed.ac.uk powell.ph.ed.ac.uk smashie.ph.ed.ac.uk wright.ph.ed.ac.uk saville.ph.ed.ac.uk chegwin.ph.ed.ac.uk brookes.ph.ed.ac.uk denning.ph.ed.ac.uk bowman.ph.ed.ac.uk bolt.ph.ed.ac.uk punjabihitsquad.ph.ed.ac.uk ugh1.ph.ed.ac.uk ugh2.ph.ed.ac.uk ugh3.ph.ed.ac.uk ugh4.ph.ed.ac.uk ugh5.ph.ed.ac.uk florence.ph.ed.ac.uk djsemtex.ph.ed.ac.uk nihal.ph.ed.ac.uk benjib.ph.ed.ac.uk sarahlove.ph.ed.ac.uk harkin.ph.ed.ac.uk titchmarsh.ph.ed.ac.uk knight.ph.ed.ac.uk dedicoat.ph.ed.ac.uk sloth.ph.ed.ac.uk grimshaw.ph.ed.ac.uk jones.ph.ed.ac.uk gingerexplosion.ph.ed.ac.uk oldjock.ph.ed.ac.uk' ,
       memory = 400 ,
       virtual_memory = 400 ,
       arch = 'INTEL' 
       ) 
    ) 
 ) 

