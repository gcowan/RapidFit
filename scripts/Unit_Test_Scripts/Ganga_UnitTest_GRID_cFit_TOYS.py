def TOY_Splitter( XML='XML.xml', total_toys=1000, toys_per_core=10, pull_name='pullPlots.root', SEED='12345'):
        args = []
	real_SEED = SEED
        for i in range(0,total_toys,toys_per_core):
                temp = []
                temp.append( str( XML ) )
                temp.append( str( toys_per_core ) )
                temp.append( "pullPlots.root" )
		temp.append( str( real_SEED ) )
		real_SEED = real_SEED + 1
                args.append( temp )
	print args
	return args

j = Job( application = Root( version = '5.26.00b' ) )
j.name = 'CFIT-Unit_TOYS'
j.application.script = File( name='./RapidFit_TOYS.C' )
j.inputsandbox = ['../unittest/Xdata/goldstandard_unbiased_timecut3.root','../unittest/Xdata/goldstandard_biased_timecut3PELC.root','./RapidFit.C','../unittest/cFit/DataFile_Tagged_cFit_MOMP_TOYS.xml','../lib/libRapidRun.so']
j.outputsandbox = [ 'pullPlots.root' ]
j.merger = RootMerger( files=[ 'pullPlots.root' ] )
j.backend = Dirac()
j.splitter=ArgSplitter( args = TOY_Splitter('DataFile_Tagged_cFit_MOMP_TOYS.xml', 1000, 10, 'pullPlots.root', 12345 ) )
j.submit()

