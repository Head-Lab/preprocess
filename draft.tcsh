#!/usr/bin/tcsh

cd ~/Desktop
set overwrite = 0
set in_file = "~/Desktop/f"
setenv SUBJECTS_DIR /Applications/freesurfer/subjects
set subject = "110214_cc32903"


#Use FSL to extract reference frame
set ref_frame = "${in_file}_ref"
set vols = `$FSLDIR/bin/fslval $in_file dim4`; set ref_vol = `echo "$vols / 2 - 1" | bc`
if ( $overwrite == 1 || ! -e ${in_file}_ref.nii.gz ) then
	$FSLDIR/bin/fslroi $in_file $ref_frame $ref_vol 1
endif

#Use FSL for slicetime correction
set out_file = "${in_file}_stc"
if ( $overwrite == 1 || ! -e ${in_file}_stc.nii.gz ) then
	set stc_file = "/Users/Tyler/Desktop/siemens_even_interleaved.txt"
	$FSLDIR/bin/slicetimer -i $in_file -o $out_file --ocustom=$stc_file
endif
set in_file = $out_file

#Use FSL to get motion outliers
set motion_outliers = "~/Desktop/motion_outliers.txt"
if ( $overwrite == 1 || ! -e ~/Desktop/motion_outliers.txt ) then
	$FSLDIR/bin/fsl_motion_outliers $in_file 0 $motion_outliers
endif 

#Use FSL for motion correction
set out_file = "${in_file}_mc"
if ( $overwrite == 1 || ! -e ${in_file}_mc.nii.gz ) then
	$FSLDIR/bin/mcflirt -in $in_file -refvol $ref_vol -o $out_file -mats -plots -rmsrel -rmsabs
endif
set in_file = $out_file

#Use FSL to create motion plots
if ( $overwrite == 1 || ! -e rot.png || ! -e trans.png || ! -e disp.png ) then
	${FSLDIR}/bin/fsl_tsplot -i ${out_file}.par -t 'Estimated rotations (radians)' -u 1 --start=1 \
						 	 --finish=3 -a x,y,z -w 640 -h 144 -o rot.png
	${FSLDIR}/bin/fsl_tsplot -i ${out_file}.par -t 'Estimated translations (mm)' -u 1 --start=4 \
                         	 --finish=6 -a x,y,z -w 640 -h 144 -o trans.png
	${FSLDIR}/bin/fsl_tsplot -i ${out_file}_abs.rms,${out_file}_rel.rms \
                        	 -t 'Estimated mean displacement (mm)' -u 1 -w 640 -h 144 \
                        	 -a absolute,relative -o disp.png
endif

set regressors = ()
foreach tissue ( "ventricles" "wm" )
	
	#setup regressor file name
	if ( $tissue == "wm" ) then
		set regressor = "wm_regressor.txt"
	else
		set regressor = "ventrical_regressor.txt"
	endif
	set regressors = ( $regressors $regressor ${regressor:r}_dt.txt )
	
	#Use FreeSurfer to get white matter and csf masks
	if ( $overwrite == 1 || ! -e ${tissue}_mask.nii.gz ) then
		$FREESURFER_HOME/bin/mri_binarize --i $SUBJECTS_DIR/$subject/mri/aseg.mgz --$tissue \
										  --erode 1 --o ${tissue}_mask.nii.gz 
	endif
	
	#Use FreeSurfer to transform masks to bold space
	if ( $overwrite == 1 || ! -e ${tissue}_mask_to_fmri.nii.gz ) then
		$FREESURFER_HOME/bin/mri_vol2vol --mov ${ref_frame}.nii.gz --targ ${tissue}_mask.nii.gz \
										 --o ${tissue}_mask_to_fmri.nii.gz --reg ${fmri_to_t1}.dat \
										 --inv
	endif
	
	#Use FSL to rebinarize interpolated masks
	if ( $overwrite == 1 || ! -e ${tissue}_mask_to_fmri.nii.gz ) then
		$FSLDIR/bin/fslmaths ${tissue}_mask_to_fmri.nii.gz -thr .5 \
		                     -bin ${tissue}_mask_to_fmri.nii.gz
	endif		
	
	#Use FSL to generate regressor timecourses
	if (  $overwrite == 1 || ! -e $regressor ) then
		$FSLDIR/bin/fslmeants -i $in_file --eig -m ${tissue}_mask_to_fmri.nii.gz \
							  -o $regressor
	endif
	
	#Generate temporal derivative of regressor timecourse
	if ( $overwrite == 1 || ! -e ${regressor:r}_dt.txt ) then
		echo "0.0" > ${regressor:r}_dt.txt
		awk 'p{print $0-p}{p=$0}' $regressor >> ${regressor:r}_dt.txt
	endif	
	
end

#Use FSL to get brain mask
set fmri_mask = "${ref_frame}_brain_mask"
if ( $overwrite == 1 || ! -e ${fmri_mask}.nii.gz ) then
	$FSLDIR/bin/bet $ref_frame ${ref_frame}_brain -F
endif

#Use FSL to erode brain mask for use as a whole-brain regressor
set reg_brain_mask = "${fmri_mask}_eroded.nii.gz"
if ( $overwrite == 1 || ! -e ${fmri_mask}_eroded.nii.gz ) then
	$FSLDIR/bin/fslmaths $fmri_mask -ero -ero ${fmri_mask}_eroded.nii.gz
endif

#Use FSL to generate whole-brain regressor timecourse
set whole_brain_regressor = "whole_brain_regressor.txt"
set regressors = ( $regressors $whole_brain_regressor ${whole_brain_regressor:r}_dt.txt )
if ( $overwrite == 1 || ! -e $whole_brain_regressor ) then
	$FSLDIR/bin/fslmeants -i $in_file --eig -m $reg_brain_mask -o $whole_brain_regressor
endif

#Generate temporal derivative of whole brain regressor timecourse
if ( $overwrite == 1 || ! -e ${whole_brain_regressor:r}_dt.txt ) then
		echo "0.0" > ${whole_brain_regressor:r}_dt.txt
		awk 'p{print $0-p}{p=$0}' $whole_brain_regressor >> ${whole_brain_regressor:r}_dt.txt
endif	

#Create nuisance design matrix
set nuisance_matrix = "nuisance_matrix.txt"
if ( $overwrite == 1 || ! -e $nuisance_matrix ) then
	paste $regressors $motion_outliers > $nuisance_matrix
endif

#Use FSL to generate residuals
set res = "${in_file}_res"
if ( $overwrite == 1 || ! -e $res ) then
	$FSLDIR/bin/fsl_glm --in=$in_file --design=$nuisance_matrix --out_res=$res
endif
	
#Use FreeSurfer to register reference frame to T1
set fmri_to_t1 = "${ref_frame}_brain_to_t1"
if ( $overwrite == 1 || ! -e ${ref_frame}_brain_to_t1.nii.gz ) then
	$FREESURFER_HOME/bin/bbregister --s $subject --mov ${ref_frame}.nii.gz --reg ${fmri_to_t1}.dat \
									--init-fsl --bold --o ${fmri_to_t1}.nii.gz \
									--rms ${fmri_to_t1}.rms --fslmat ${fmri_to_t1}.mat
endif


