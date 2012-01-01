#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

"""

fc_scrub.py: Python script to scrub high movment frames from functional connectivty data.

Method described in Power et. al 2011 and following comment by (ADD)

"""

#Import system modules
import sys, argparse, csv, subprocess

#Parse arguments
arg_parse = argparse.ArgumentParser(description='Motion Scrub fcMRI data')
#Positional arguments
arg_parse.add_argument('fc',help='Path to nifti fcMRI image',nargs=1)
arg_parse.add_argument('mcf',help='MCFLIRT estimated rotations and translations file',nargs=1)
arg_parse.add_argument('outroot',help='Root for outputed files',nargs=1)
#Optional arguments
arg_parse.add_argument('-fcut',help='Frame displacement cutoff. Default is 0.3.',type=float,default=[0.1],nargs=1)
arg_parse.add_argument('-dcut',help='DVAR cutoff. Default is 500.',type=float,default=[300.0],nargs=1)
arg_parse.add_argument('-nonorm',help='Do not do mean 10000 intensity normilization.',action='store_const',default=[0.0],const=1)
arg_parse.add_argument('-mask',help='BOLD defined voxel mask.',default=[1.0],nargs=1)
arg_parse.add_argument('-inter',help='Write out intermediate masks and filtered images.',action='store_const',default=[0.0],const=1)
args = arg_parse.parse_args()

#Import external modules
import numpy as np, nibabel as nib
from scipy.interpolate import interp1d

print 'Loading image data'

#Load BOLD image. Use float 32 to avoid precision problems when calculating RMS error.
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
	
print 'Determining RMS outliers'

#Generate 4D mask
fc_mask_4d_data = np.empty_like(fc_data)
fc_mask_4d_data[:,:,:,:] = np.expand_dims(fc_mask_3d_data,axis=3)
fc_masked_data = np.ma.array(fc_data,mask=fc_mask_4d_data)

#Unless user overrides, do mean 10,000 intensity normilization
if args.nonorm[0] != 1.0:
	scale = 10000 / np.ma.average(fc_masked_data) 
	fc_masked_data = np.multiply(scale,fc_masked_data)

#Get RMS error
fc_rms_data = np.ma.sqrt(np.apply_over_axes(np.average,np.ma.diff(fc_masked_data)**2.0,[0,1,2]))

#Add 0 back to bring back frame count
fc_rms_data = np.insert(fc_rms_data,0,0.0)
print fc_rms_data

#Determine root mean square outliers
rms_outliers = fc_rms_data > 11.0

if args.inter[0] == 1.0:

	#Write out root mean square mask
	frame_rms_mask_data = np.empty_like(fc_data)
	frame_rms_mask_data[:,:,:,:] = rms_outliers
	frame_rms_mask = nib.Nifti1Image(frame_rms_mask_data,fc.get_affine())
	frame_rms_mask.to_filename(args.outroot[0] + '_rms_mask.nii.gz')

	#Write out masked image	
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
disp_outliers = frame_dis > .3

if args.inter[0] == 1.0:
	
	#Write out displacement mask
	frame_dis_mask_data = np.empty_like(fc_data)
	frame_dis_mask_data[:,:,:,:] = disp_outliers
	frame_dis_mask = nib.Nifti1Image(frame_dis_mask_data,fc.get_affine())
	frame_dis_mask.to_filename(args.outroot[0] + '_dis_mask.nii.gz')

	#Write out masked image
	frame_dis_masked_data = np.compress(np.invert(disp_outliers),fc_data,axis=3)
	frame_dis_masked = nib.Nifti1Image(frame_dis_masked_data,fc.get_affine())
	frame_dis_masked.to_filename(args.outroot[0]  + '_dis_masked.nii.gz')

print 'Determining outlier intersection'

outlier_intersection = np.logical_and(rms_outliers,disp_outliers)
print outlier_intersection

#Write Intersection Mask
frame_inter_mask_data = np.empty_like(fc_data)
frame_inter_mask_data[:,:,:,:] = outlier_intersection
frame_inter_mask = nib.Nifti1Image(frame_inter_mask_data,fc.get_affine())
frame_inter_mask.to_filename(args.outroot[0] + '_inter_mask.nii.gz')

print 'Interpolating bad frames'

#Get index for later use in intepolation
index = np.arange(0,outlier_intersection.shape[0],1)

#Delete the masked elements from both the index and the data
filtered_index = np.compress(np.invert(outlier_intersection),index)
filtered_data = np.compress(np.invert(outlier_intersection),fc_data,axis=3)

#Interpolate the removed values with nearest neighbor
interpolate = interp1d(filtered_index,filtered_data,kind='nearest')

#Use the original index to get a full series back. Note need to check if this changes the 
#non-removed data. If it does, will have to put interpolated values back into the data.
interpolated_data = interpolate(index)

#Write out filtered image
frame_inter_masked = nib.Nifti1Image(interpolated_data,fc.get_affine())
frame_inter_masked.to_filename(args.outroot[0] + '_interp.nii.gz')

print 'Using FSL to bandpass filter the data'

#Call FSL maths to temporal bandpass filter
fsl_bptf = subprocess.call('fslmaths ' + args.outroot[0] + '_interp' + ' -bptf 2.0 1.0 ' + args.outroot[0] + '_interp_bptf',shell=True)

#Import filtered data
try:
	bptf = nib.load(args.outroot[0] + '_interp_bptf.nii.gz')
except (IOError,nib.spatialimages.ImageFileError):
	print 'Error loading FSL filtered image'
	sys.exit()	
bptf_data = np.float32(bptf.get_data())

print 'Removing outliers after band pass filtering'

#Removed bad volumes using intersected mask
bptf_filtered_data = np.compress(np.invert(outlier_intersection),bptf_data,axis=3)

#Write out final cleaned image.
bptf_filtered = nib.Nifti1Image(bptf_filtered_data,fc.get_affine())
bptf_filtered.to_filename(args.outroot[0] + '_bptf_filtered.nii.gz')



	


