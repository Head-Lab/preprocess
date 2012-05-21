#!/bin/tcsh -f
#Launch bedpostx using results of a dti_preproc run

if ( $#argv != 1 ) then
	echo "Usage: ${0:t} <dti_preproc directory>"
	exit 1
endif
set pdir = $1; pushd $pdir; set subj = `awk '/Outroot/{print $2}' dti_preproc.log | head -1`

#Create links to make bedpostx happy
set ddir = `pwd`/dti
ln -s ${ddir}/${subj}_dti_bvals.txt ${ddir}/bvals; ln -s ${ddir}/${subj}_dti_bvecs_rot.txt ${ddir}/bvecs
ln -s ${ddir}/${subj}_dti_b0_brain_mask.nii.gz ${ddir}/nodif_brain_mask.nii.gz
ln -s ${ddir}/${subj}_dti_eddy.nii.gz ${ddir}/data.nii.gz

#Run bedpostx
bedpostx dti

