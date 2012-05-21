#!/bin/tcsh -f
#Get a set of labeled surfaces from FreeSurfer segmentation

if ( $#argv != 2 ) then
	echo "Usage: $0:t <subject_id> <subjects_dir>"
	exit 1
endif
set subj = $1; set outdir = label; setenv SUBJECTS_DIR $2; mkdir $outdir; pushd $outdir

#First convert white matter surfaces to ascii
$FREESURFER_HOME/bin/mris_convert -a $SUBJECTS_DIR/$subj/surf/lh.white ./${subj}_lh.white.asc
$FREESURFER_HOME/bin/mris_convert -a $SUBJECTS_DIR/$subj/surf/rh.white ./${subj}_rh.white.asc

#Get label files for each hemis annotation
$FREESURFER_HOME/bin/mri_annotation2label --subject $subj --hemi lh --outdir ./
$FREESURFER_HOME/bin/mri_annotation2label --subject $subj --hemi rh --outdir ./

#Get a label surface for each parcellation. Keeping hemisphere's seperate in order to keep the file size down
foreach hemi ( lh rh )
	foreach label ( ${hemi}*label )
		echo "$label" > to_convert.txt
		~/bin/label2surf --surf=${subj}_${hemi}.white.asc --out=${label:r}.asc --labels=to_convert.txt
	end
end
rm to_convert.txt
