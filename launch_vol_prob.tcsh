#!/bin/tcsh -f
#Run probtrackx2 on each ROI with everyother ROI as a target

if ( $#argv != 2 ) then
	echo "Usage: $0:t <subject_id> <subjects_dir>"
	exit 1
endif
set subj = $1; set outdir = `pwd`/vol_prob; setenv SUBJECTS_DIR $2; mkdir $outdir; pushd $outdir

#Use FreeSurfer to get a stop mask
mri_binarize --i $SUBJECTS_DIR/$subj/mri/aparc+aseg.mgz --o ${subj}_stop.nii.gz --ventricles --match 0 15 24 30 62

#Use FreeSurfer to get a bunch of seeds
mkdir vol_seeds; set seeds = ()
foreach line ( `cat ~/preprocess/roi_lookups.txt` )
	set match = `echo $line | cut -d, -f1`; set seed = `echo $line | cut -d, -f2`
	mri_binarize --i $SUBJECTS_DIR/$subj/mri/aparc+aseg.mgz --o vol_seeds/${seed}.nii.gz --match $match
	set seeds = ( $seeds vol_seeds/$seed )
end

#Run probtrackx2 on each seed. 
foreach seed ( $seeds )
	echo -n "" > ${seed}_targets.txt
	foreach target ( $seeds )
		if ( $seed != $target ) echo $target >> ${seed}_targets.txt
	end
	echo '#\!/bin/tcsh' > ${seed}_command.txt
	echo "#PBS -N ${subj}_${seed:t:r} -l nodes=1:ppn=1,walltime=24:00:00,vmem=6gb -q dque" >> ${seed}_command.txt
	echo 'source /home/tblazey/.tcshrc' >> ${seed}_command.txt
	echo "cd $outdir" >> ${seed}_command.txt
	echo "probtrackx --samples=../dti.bedpostX/merged --mask=../dti.bedpostX/nodif_brain_mask --seed=$seed --dir=$seed:t:r --forcedir --opd --os2t --s2tastext --targetmasks=${seed}_targets.txt --stop=${subj}_stop.nii.gz --xfm=../atlas/${subj}_free_to_dti.mat --seedref=../atlas/${subj}_T1.nii.gz --loopcheck --cthr=0.2 --nsteps=2000 --nsamples=5000" >> ${seed}_command.txt 
	#sleep 2
	#qsub ${seed}_command.txt
end





