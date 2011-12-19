#!/usr/bin/python

"""

fc_scrub.py: Python script to scrub high movment frames from functional connectivty data.

Method described in Power et. al 2011

"""

#Import system modules
import sys, argparse, csv

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Motion Scrub fcMRI data')
#Positional arguments
arg_parse.add_argument('fc',help='Path to nifti fcMRI image',nargs=1)
arg_parse.add_argument('mcf',help='MCFLIRT estimated rotations and translations file',nargs=1)
arg_parse.add_argument('outroot',help='Root for outputed files',nargs=1)
#Optional arguments
arg_parse.add_argument('-fcut',help='Frame displacement cutoff. Default is 0.3.',type=float,default=[0.5],nargs=1)
arg_parse.add_argument('-dcut',help='DVAR cutoff. Default is 500. (Assumes 10,000 normalized).',type=float,default=[500.0],nargs=1)
arg_parse.add_argument('-mask',help='BOLD defined voxel mask.',default=[1.0],nargs=1)
args = arg_parse.parse_args()

#Import external modules
import numpy as np, nibabel as nib

#Load BOLD image
try:
	fc = nib.load(args.fc[0])
except (IOError,nib.spatialimages.ImageFileError):
	print 'Cannot find fcMRI image at %s'%(args.perf[0])
	sys.exit()	
fc_data = np.float32(fc.get_data())

#Load Mask
if args.mask[0] != 1.0:	
	try:
		fc_mask_3d = nib.load(args.mask[0])
	except (IOError,nib.spatialimages.ImageFileError):
		print 'Cannot find mask at %s'%(args.mask[0])
		sys.exit()
	fc_mask_3d_data = fc_mask_3d.get_data()
	fc_mask_3d_data = ( fc_mask_3d_data - 1 ) * -1
else:
	fc_mask_3d_data = np.less_equal(np.min(fc_data,axis=3),0.0)
	
print 'Determining root mean square outliers'

#Generate 4D mask
fc_mask_4d_data = np.empty_like(fc_data)
fc_mask_4d_data[:,:,:,:] = np.expand_dims(fc_mask_3d_data,axis=3)
fc_masked_data = np.ma.array(fc_data,mask=fc_mask_4d_data)

#Get root mean square error
fc_rms_data = np.ma.sqrt(np.ma.average(np.ma.average(np.ma.average(np.ma.diff(fc_masked_data)**2.0,axis=0),axis=0),axis=0))

#Add 0 back to bring back frame count
fc_rms_data = np.insert(fc_rms_data,0,0.0)

#Determine root mean square outliers
rms_outliers = fc_rms_data > 11.0

#Write out root mean square mask
frame_rms_mask_data = np.empty_like(fc_data)
frame_rms_mask_data[:,:,:,:] = rms_outliers
frame_rms_mask = nib.Nifti1Image(frame_rms_mask_data,fc.get_affine())
frame_rms_mask.to_filename(args.outroot[0] + '_rms_mask.nii.gz')

#Write out filtered image
frame_rms_masked_data = np.compress(np.invert(rms_outliers),fc_data,axis=3)
frame_rms_masked = nib.Nifti1Image(frame_rms_masked_data,fc.get_affine())
frame_rms_masked.to_filename(args.outroot[0]  + '_rms_masked.nii.gz')

print 'Determining frame displacement outliers'

#Read in motion paramaters
motion_params = np.loadtxt(args.mcf[0])

#Calculate framewise differences
motion_diff = np.diff(motion_params,axis=0)

#Split into rotation and arrays
rot_diff,trans_diff = np.hsplit(motion_diff,2)

#Calculate rotational displacements
rot_diff_dis = np.absolute(rot_diff) * 50.0

#Calculate total frame displacements
frame_dis = np.sum(np.absolute(trans_diff),axis=1) + np.sum(rot_diff_dis,axis=1)

#Add a zero to top to bring back frame count
frame_dis = np.insert(frame_dis,0,0.0)

#Write out frame displacements
np.savetxt(args.outroot[0] + '_frame_displacements.txt',frame_dis)

#Calculate temporal mask
outliers = frame_dis > .3
#Mask out one back
outlier_back = np.roll(outliers,-1); outlier_back[0] = False
outlier_one = np.roll(

#Write out displacement mask
frame_dis_mask_data = np.empty_like(fc_data)
frame_dis_mask_data[:,:,:,:] = outliers
frame_dis_mask = nib.Nifti1Image(frame_dis_mask_data,fc.get_affine())
frame_dis_mask.to_filename(args.outroot[0] + '_dis_mask.nii.gz')

#Write out filtered image
frame_dis_masked_data = np.compress(np.invert(outliers),fc_data,axis=3)
frame_dis_masked = nib.Nifti1Image(frame_dis_masked_data,fc.get_affine())
frame_dis_masked.to_filename(args.outroot[0]  + '_dis_masked.nii.gz')

print 'Determining outlier intersection'

outlier_intersection = np.logical_and(rms_outliers,outliers)

#Write out root mean square mask
frame_inter_mask_data = np.empty_like(fc_data)
frame_inter_mask_data[:,:,:,:] = outlier_intersection
frame_inter_mask = nib.Nifti1Image(frame_inter_mask_data,fc.get_affine())
frame_inter_mask.to_filename(args.outroot[0] + '_inter_mask.nii.gz')

#Write out filtered image
frame_inter_masked_data = np.compress(np.invert(outlier_intersection),fc_data,axis=3)
frame_inter_masked = nib.Nifti1Image(frame_inter_masked_data,fc.get_affine())
frame_inter_masked.to_filename(args.outroot[0]  + '_inter_masked.nii.gz')






	


