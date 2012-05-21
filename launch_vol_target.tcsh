#!/bin/tcsh -f
#Run probtrackx2 on each ROI with everyother ROI as a target

if ( $#argv != 2 ) then
	echo "Usage: $0:t <subject_id> <subjects_dir>"
	exit 1
endif
set subj = $1; set outdir = vol_struct; setenv SUBJECTS_DIR $2; mkdir $outdir; pushd $outdir

#Use FreeSurfer to get a stop mask
mri_binarize --i $SUBJECTS_DIR/$subj/mri/aparc+aseg.mgz --o ${subj}_stop.nii.gz --ventricles --match 0 15 24 30 62

#Use FreeSurfer to get a bunch of targets
mkdir vol_targets; set targets = (); echo -n "" > vol_targets.txt
foreach line ( `cat ~/preprocess/roi_lookups.txt` )
	set match = `echo $line | cut -d, -f1`; set target = `echo $line | cut -d, -f2`
	mri_binarize --i $SUBJECTS_DIR/$subj/mri/aparc+aseg.mgz --o vol_targets/${target}.nii.gz --match $match
	set targets = ( $targets vol_targets/$target ); echo "vol_targets/${target}.nii.gz" >> vol_targets.txt
end

#Run probtrackx2 on each target. 
foreach target ( $targets )
	echo "probtrackx2 --samples=../dti.bedpostX/merged --mask=../dti.bedpostX/nodif_brain_mask --seed=$target --dir=$target:t:r --forcedir --opd --os2t --s2tastext --targetmasks=vol_targets.txt --stop=${subj}_stop.nii.gz --xfm=../atlas/${subj}_free_to_dti.mat --seedref=../atlas/${subj}_T1.nii.gz --loopcheck --cthr=0.2 --nsteps=2000 --nsamples=5000" > command.txt 
end





