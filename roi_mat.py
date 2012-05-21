#!/usr/bin/python

import csv, os, sys, argparse
import numpy as np, matplotlib.pyplot as plt, nibabel as nib

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Generate adjacency matrix from launch_prob run')
#Positional arguments
arg_parse.add_argument('dir',help='Tractography results directory',nargs=1)
arg_parse.add_argument('seeds',help='CSV seeds file used as input to launch_prob',nargs=1)
arg_parse.add_argument('out',help='Root for outupts',nargs=1)
#Optional arguments
arg_parse.add_argument('-nv',action='store_const',const=1,help='Normalize by number of voxels/verticies')
arg_parse.add_argument('-nw',action='store_const',const=1,help='Normalize by waytotals')
args = arg_parse.parse_args()

#Get ROI list
try:
	roi_file = open(args.seeds[0],'r')
except (IOError):
	print 'Cannot load seeds file at %s. Exiting...'%(args.seeds[0])
	sys.exit()
reader = csv.reader(roi_file)

#Create empty adjacency matrix
adj_mat = np.zeros((82,82))

#Loop through each ROI and get the connectivity values
i = 0
for roi in reader:
	
	#Set proper roi results directory
	dir = os.path.join(args.dir[0],roi[1])
	if os.path.exists(dir) is False: 
		print 'Cannot find roi dir for %s. Exiting...'%(roi[1])
		sys.exit()
	
	#Load matrix_seeds_to_all_targets for roi
	mat = os.path.join(dir,'matrix_seeds_to_all_targets')
	try:
		roi_mat = np.loadtxt(mat)
	except (IOError):
		print 'Cannot load matrix_seeds_to_all_targets for %s. Exiting...'%(roi[1])
		sys.exit()
	
	#Sum all the columns and reshape to get a 1,81 row vector.
	row = np.sum(roi_mat,axis=0).reshape(1,81)
	
	#Insert zero-diagonal element
	row = np.insert(row,i,0,axis=1)
	
	#Add summed row to adjacency matrix
	adj_mat[i,:]  = row

	i += 1

#Save adjacency matrix as nifti file
img = nib.Nifti1Image(adj_mat,np.eye(4))
img.to_filename('%s.nii.gz'%(args.out[0]))

#Plot the adjacency matrix
im = plt.imshow(adj_mat); plt.colorbar(im); plt.savefig('%s.png'%(args.out[0]))
