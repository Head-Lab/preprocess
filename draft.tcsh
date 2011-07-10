#!/usr/bin/tcsh -f

cd ~/Desktop/practice
set overwrite = 0
set echo
set in_file = "~/Desktop/practice/f"
setenv SUBJECTS_DIR /Applications/freesurfer/subjects
set subject = "110214_cc32903"
set trep = 2.2
set hpass = .009
set lpass = .08

#Use FSL to extract reference frame
set ref_frame = "${in_file}_ref"
set vols = `$FSLDIR/bin/fslval $in_file dim4`; set ref_vol = `echo "$vols / 2 - 1" | bc`
if ( $overwrite == 1 || ! -e ${in_file}_ref.nii.gz ) then
	$FSLDIR/bin/fslroi $in_file $ref_frame $ref_vol 1
endif

#Use FSL to get brain mask
set fmri_mask = "${ref_frame}_brain_mask"
if ( $overwrite == 1 || ! -e ${fmri_mask}.nii.gz ) then
	$FSLDIR/bin/bet $ref_frame ${ref_frame}_brain -F
endif

#Use FreeSurfer to register reference frame to T1
set fmri_to_t1 = "${ref_frame}_brain_to_t1"
if ( $overwrite == 1 || ! -e ${ref_frame}_brain_to_t1.nii.gz ) then
	$FREESURFER_HOME/bin/bbregister --s $subject --mov ${ref_frame}.nii.gz --reg ${fmri_to_t1}.dat \
									--init-fsl --bold --o ${fmri_to_t1}.nii.gz \
									--rms ${fmri_to_t1}.rms --fslmat ${fmri_to_t1}.mat
endif

#Use FSL for slicetime correction
set out_file = "${in_file}_stc"
if ( $overwrite == 1 || ! -e ${out_file}.nii.gz ) then
	set stc_file = "/Users/Tyler/Desktop/practice/siemens_even_interleaved.txt"
	$FSLDIR/bin/slicetimer -i $in_file -o $out_file --ocustom=$stc_file
endif
set in_file = $out_file

#Use FSL to get motion outliers
set motion_outliers = "motion_outliers.txt"
if ( $overwrite == 1 || ! -e $motion_outliers ) then
	$FSLDIR/bin/fsl_motion_outliers $in_file 0 $motion_outliers
endif 

#Use FSL for motion correction
set out_file = "${in_file}_mc"
if ( $overwrite == 1 || ! -e ${out_file}.nii.gz ) then
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

#Create a motion regressor file with temporal derivatives
set motion_regressor = "motion_regressor.txt"
if ( $overwrite == 1 || ! -e $motion_regressor ) then
	@ count = 1; set tmp_regressors = ()
	while ( $count <= 6 )
		sed 's/  / /g' ${out_file}.par | cut -d " " -f $count > tmp_${count}_regressor.txt
		awk 'p{print $0-p}{p=$0}' tmp_${count}_regressor.txt > tmp_${count}_regressor_dt.txt
		set tmp_regressors = ( $tmp_regressors tmp_${count}_regressor.txt \
							   tmp_${count}_regressor_dt.txt )
		@ count ++
	end
	paste $tmp_regressors > $motion_regressor
	rm $tmp_regressors
endif

#Use FSL to do 4D intensity normalization
set out_file = "${in_file}_inorm"
if ( $overwrite == 1 || ! -e ${out_file}.nii.gz ) then
	set fifty_intensity = `$FSLDIR/bin/fslstats $in_file -k $fmri_mask -p 50`
	set scale_factor = `echo "scale=5; 1000 / $fifty_intensity" | bc -l`
	$FSLDIR/bin/fslmaths $in_file -mul $scale_factor $out_file
endif
set in_file = $out_file

#Use FSL to do temporal bandpass filtering
set out_file = "${in_file}_bptf"
if ( $overwrite == 1 || ! -e ${out_file}.nii.gz ) then
	set lp_vol = `echo "scale=5; 1 / $lpass / $trep" | bc -l`
	set hp_vol = `echo "scale=5; 1 / $hpass / $trep" | bc -l`
	$FSLDIR/bin/fslmaths $in_file -bptf $hp_vol $lp_vol $out_file
endif
set in_file = $out_file

set regressors = ()
foreach tissue ( "ventricles" "wm" )
	
	set regressor = "${tissue}_regressor.txt"
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
	if ( $overwrite == 1 || ! -e ${tissue}_mask_to_fmri_bin.nii.gz ) then
		$FSLDIR/bin/fslmaths ${tissue}_mask_to_fmri.nii.gz -thr .5 \
		                     -bin ${tissue}_mask_to_fmri_bin.nii.gz
	endif		
	
	#Use FSL to generate regressor timecourses
	if (  $overwrite == 1 || ! -e $regressor ) then
		$FSLDIR/bin/fslmeants -i $in_file --eig -m ${tissue}_mask_to_fmri_bin.nii.gz \
							  -o $regressor
	endif
	
	#Generate temporal derivative of regressor timecourse
	if ( $overwrite == 1 || ! -e ${regressor:r}_dt.txt ) then
		echo "0.0" > ${regressor:r}_dt.txt
		awk 'p{print $0-p}{p=$0}' $regressor >> ${regressor:r}_dt.txt
	endif	
	
end

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
	paste $regressors $motion_outliers $motion_regressor > $nuisance_matrix
endif

#Use FSL to generate residuals
set out_file = "${in_file}_res"
if ( $overwrite == 1 || ! -e ${out_file}.nii.gz ) then
	$FSLDIR/bin/fsl_glm --in=$in_file --design=$nuisance_matrix --out_res=$out_file
endif
set in_file = $out_file

set out_files = ()
foreach hemi ( lh rh )
	set out_file = "${in_file}_res_${hemi}_fsaverage_fwhm5"
	set out_files = ( $out_files $out_file )
	#Use FreeSurfer to transform residuals to FreeSurfer surface space. Smooth by 5mm on surface
	if ( $overwrite == 1 || ! -e ${out_file}.nii.gz ) then
		$FREESURFER_HOME/bin/mri_vol2surf --mov ${in_file}.nii.gz --reg ${fmri_to_t1}.dat \
										  --trgsubject fsaverage --hemi $hemi \
										  --o ${out_file}.nii.gz --surf-fwhm 5
	endif
end

	
	
	
	
	
	
	
	
	
