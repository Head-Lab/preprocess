#!/bin/tcsh -f
#Run probtrackx2 on each ROI with everyother ROI as a target

if ( $#argv != 2 ) then
	echo "Usage: $0:t <subject_id> <subjects_dir>"
	exit 1
endif
set subj = $1; set outdir = struct; setenv SUBJECTS_DIR $2; mkdir $outdir; pushd $outdir

#Use FreeSurfer to get a stop mask
mri_binarize --i $SUBJECTS_DIR/$subj/mri/aparc+aseg.mgz --o ${subj}_stop.nii.gz --ventricles --match 0 15 24 30 62

#Create links to all the nodes
mkdir nodes; set nodes = ()
echo -n "" > nodes/nodes.txt
foreach snode ( `cat ~/preprocess/sub_nodes.txt` )
	if ( ! -e ../first/${subj}_rawavg_std-${snode} ) then
		echo "Error: Cannot find subcortical node $snode"; exit 1
	endif
	set nodes = ( $nodes nodes/$snode )
	ln -s ../../first/${subj}_rawavg_std-${snode} nodes/$snode; echo nodes/$snode >> nodes.txt 
end
foreach cnode ( `cat ~/preprocess/cort_nodes.txt` )
	if ( ! -e ../label/$cnode ) then
		echo "Error: Cannot find cortical node $cnode"; exit 1
	endif
	set nodes = ( $nodes nodes/$cnode )
	ln -s ../../label/$cnode nodes/$cnode; echo nodes/$cnode >> nodes.txt 
end

#Run probtrackx2 on each node. 
foreach node ( $nodes )
	echo "probtrackx2 --samples=../dti.bedpostX/merged --mask=../dti.bedpostX/nodif_brain_mask --seed=$node --dir=$node:t:r --forcedir --opd --os2t --s2tastext --targetmasks=nodes.txt --stop=${subj}_stop.nii.gz --xfm=../atlas/${subj}_free_to_dti.mat --seedref=../atlas/${subj}_T1.nii.gz --meshspace=freesurfer --loopcheck --cthr=0.2 --nsteps=2000 --nsamples=5000" > command.txt 
end





