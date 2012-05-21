#!/bin/tcsh -f
#Run FSL's first to get subcortical surfaces and then convert to FreeSurfer space

if ( $#argv != 2 ) then
	echo "Usage: $0:t <subject_id> <subjects_dir>"
	exit 1
endif
set subj = $1; set outdir = first; setenv SUBJECTS_DIR $2; mkdir $outdir; pushd $outdir

#Convert .mgz files to .nii.gz
foreach img ( $SUBJECTS_DIR/$subj/mri/rawavg.mgz $SUBJECTS_DIR/$subj/mri/T1.mgz )
	mri_convert $img ${subj}_${img:t:r}.nii.gz
end

#Reorient rawavg to standard orieintation so that first_flirt doesn't fail
$FSLDIR/bin/fslreorient2std ${subj}_rawavg ${subj}_rawavg_std

#Run first to get all subcortical surfaces. Using modified script that does not use fsl_sub.
~/scripts/my_run_first_all -i ${subj}_rawavg_std -o ${subj}_rawavg_std -v

#Get a transformation from reoriented rawavg to conformed FreeSurfer T1
$FREESURFER_HOME/bin/tkregister2_cmdl --targ ${subj}_T1.nii.gz --mov ${subj}_rawavg_std.nii.gz --regheader --reg ${subj}_raw_std_to_free.reg --fslregout ${subj}_raw_std_to_free.mat --noedit

#Transform each subcortical surface to FreeSurfer space
foreach vtk ( *vtk )
	~/bin/surf2surf --surfin=$vtk --surfout=${vtk:r}.asc --convin=first --convout=freesurfer --volin=${subj}_rawavg_std.nii.gz --volout=${subj}_T1.nii.gz --xfm=${subj}_raw_std_to_free.mat
end

