from math import *

def grid_points( par1min=0., par1max=0., par1res=1, par2min=0., par2max=0., par2res=1 ):
	grid=[]
	step1 = abs(par1min-par1max)/float(par1res)
	step2 = abs(par2min-par2max)/float(par2res)
	for i in range(0,par1res+1,1):
		row=[]
		par1val = par1min + float(i)*step1
		for j in range(0,par2res+1,1):
			par2val = par2min + float(j)*step2
			point=[]
			point.append( str(par1val) )
			point.append( str(par2val) )
			row.append(point)
		grid.append(row)
		#print row
	return grid

def cell_maker( grid_points, cellsize=2 ):
	len_x = len(grid_points)-1
	len_y = len(grid_points[0])-1
	grid_step = float(grid_points[0][1][1])-float(grid_points[0][0][0])
	split_grid=[]
	for i in range(0,len(grid_points),cellsize):
		for j in range(0,len(grid_points[0]),cellsize):
			min_i = i
			min_j = j
			max_i = i+cellsize-1
			max_j = j+cellsize-1
			points_x = cellsize
			points_y = cellsize
			if max_i > len_x :
				points_x-= int((max_i-len_x)/grid_step)
				max_i = len_x
			if max_j > len_y :
				points_y-= int((max_j-len_y)/grid_step)
				max_j = len_y
			grid_limits = []
			grid_limits.append( [ grid_points[min_i][min_j][0], grid_points[max_i][max_j][0], str(points_x) ] )
			grid_limits.append( [ grid_points[min_i][min_j][1], grid_points[max_i][max_j][1], str(points_y) ] )
			split_grid.append( grid_limits )
	return split_grid


def LL_Splitter(XML='XML.xml',par1='x',par1min=0.,par1max=0.,par1res=1,par2='y',par2min=0.,par2max=0.,par2res=1,cellsize=2):
	args = []
	par1res -= 1
	par2res -= 1
	grid_array = grid_points( par1min, par1max, par1res, par2min, par2max, par2res )
	grid_cells = cell_maker( grid_array, cellsize )
	for cell in grid_cells:
		par1str = par1+","+str(cell[0][0])+","+str(cell[0][1])+","+str(cell[0][2])+' '
		par2str = par2+","+str(cell[1][0])+","+str(cell[1][1])+","+str(cell[1][2])+' '
		temp_arg = [ XML, par1str, par2str ]
		args.append( temp_arg )
	print args
	return args

j = Job( application = Root( version = '5.26.00b' ) )
j.inputsandbox = ['../unittest/Xdata/goldstandard_unbiased_timecut3.root','../unittest/Xdata/goldstandard_biased_timecut3PELC.root','./RapidFit.C','../unittest/sFit/DataFile_Tagged_sFit_MOMP.xml','../lib/libRapidRun.so']
j.name = 'CFIT-Unit'
j.application.script = File( name='./RapidFit.C' )
j.outputsandbox = ['LLcontourScanData.root']
j.merger = RootMerger( files=['LLcontourScanData.root'] )
j.backend = Dirac()
#j.splitter=ArgSplitter( args = LL_Splitter('DataFile_Tagged_sFit_MOMP.xml','Phi_s',-pi,pi,40,'deltaGamma',-0.7,0.7,40,5) )
j.splitter=ArgSplitter( args = LL_Splitter('DataFile_Tagged_sFit_MOMP.xml','Phi_s',-3.2,3.2,40,'deltaGamma',-0.7,0.7,40,8) )
j.submit()
