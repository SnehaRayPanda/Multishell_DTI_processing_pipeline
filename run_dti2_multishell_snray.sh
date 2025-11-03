#! /bin/bash
# -----------------------------------------------------------------------------
# Multishell DTI preprocessing with AP/PA fieldmaps
# Batch processing script using MRtrix3 + FSL + ANTs
#
# To run: 
#   chmod +x run_dti.sh
#   ./run_dti.sh
#
# Required per subject folder:
#   DWI data:      data.nii.gz, bvecs, bvals
#   Fieldmaps:     AP.nii.gz, AP.bvec, AP.bval
#                  PA.nii.gz, PA.bvec, PA.bval
#   T1-weighted:   T1.nii.gz
#   Atlas:         aal.nii (or dk68.nii) in MNI space

#   Author: Sneha Ray
#   Modified from MRTRIX pipeline
#   Guided and supervised by Dr. Neeraj Upadhyay
# -----------------------------------------------------------------------------

for j in ANAL ;do   # change 'ANAL' to your subject folders (e.g., sub-01 sub-02 ...)
echo "Processing subject: $j"
cd $j

d=data.nii.gz
bvc=bvecs
bvl=bvals
t=T1.nii.gz

ap=AP.nii.gz
apbvec=AP.bvec
apbval=AP.bval

pa=PA.nii.gz
pabvec=PA.bvec
pabval=PA.bval

#### Step 1 - Preprocessing of DTI data ######

### Step 1.1 - Conversion of DTI file formats into mif (MRtrix Image Format)
mrconvert $d dwi.mif -fslgrad $bvc $bvl -force
mrconvert $ap AP.mif -fslgrad $apbvec $apbval -force
mrconvert $pa PA.mif -fslgrad $pabvec $pabval -force

### Step 1.2 - Denoising and Gibbs ringing removal
dwidenoise dwi.mif dwi_den.mif -noise noise.mif -force
mrdegibbs dwi_den.mif dwi_den_degibbs.mif -force

### Step 1.3 - Create AP/PA b0 pair for topup/eddy
dwiextract AP.mif - -bzero | mrmath - mean b0_AP.mif -axis 3 -force
dwiextract PA.mif - -bzero | mrmath - mean b0_PA.mif -axis 3 -force
mrcat b0_AP.mif b0_PA.mif -axis 3 b0_pair.mif -force 

### Step 1.4 - Motion and distortion correction (TOPUP + EDDY)
dwifslpreproc dwi_den_degibbs.mif dwi_preproc.mif \
  -pe_dir ap \
  -rpe_pair -se_epi b0_pair.mif \
  -readout_time 0.038 \
  -eddy_options " --slm=linear --repol --data_is_shelled " \
  -force

### Step 1.5 - Bias Field correction
dwibiascorrect ants dwi_preproc.mif dwi_unbiased.mif -force

### Step 1.6 - Brain mask estimation 
dwi2mask dwi_unbiased.mif mask.mif -force

#### Step 2 - T1w registration in DTI space ####
# Brain extraction
mri_synthstrip -i $t -o T1_brain.nii.gz

# 5TT segmentation
5ttgen fsl $t 5tt.mif -force

# Create mean b0 for registration
dwiextract dwi_unbiased.mif - -bzero | mrmath - mean b0_mean.mif -axis 3 -force
mrconvert b0_mean.mif b0_mean.nii.gz -force

# Register b0 to T1
flirt -ref T1_brain.nii.gz -in b0_mean.nii.gz -omat b02t1.mat -dof 6 -interp nearestneighbour
transformconvert b02t1.mat b0_mean.nii.gz T1_brain.nii.gz flirt_import diff2struct_mrtrix.txt
convert_xfm -omat t12b0.mat -inverse b02t1.mat
mrtransform 5tt.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif
5tt2gmwmi 5tt_coreg.mif gmwmseed_coreg.mif

#### Step 3 - Fiber Orientation Distribution (FOD) ####
### Step 3.1 - Response function estimation (Dhollander for multishell)
dwi2response dhollander dwi_unbiased.mif wm.txt gm.txt csf.txt -mask mask.mif -voxels voxels.mif -force

### Step 3.2 - Estimation of Fiber Orientation Distribution (MSMT-CSD)
dwi2fod msmt_csd dwi_unbiased.mif \
  -mask mask.mif \
  wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif -force

### Step 3.3 - Intensity Normalization (optional but recommended)
mtnormalise wmfod.mif wmfod_norm.mif \
            gmfod.mif gmfod_norm.mif \
            csffod.mif csffod_norm.mif \
            -mask mask.mif -force

#### Step 4 - Tractography ####
### Step 4.1 - Whole-brain tractography with ACT
tckgen -act 5tt_coreg.mif -backtrack -seed_gmwmi gmwmseed_coreg.mif \
       -maxlength 250 -cutoff 0.06 \
       -select 10000000 wmfod_norm.mif tracks_10M.tck -force

### Step 4.2 - Reduce biases using SIFT2
tcksift2 -act 5tt_coreg.mif tracks_10M.tck wmfod_norm.mif sift2_weights.txt -force

#### Step 5 - Atlas registration ## Transform the atlas (e.g., aal, Shen, schaefer) parcellation to the T1w space
#applywarp --ref=T1w_preproc.nii.gz --in= Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz --warp=mni2t1_transf.nii.gz  --premat=mni2t1.mat --out=Schaefer400_on_T1.nii.gz --interp=nn
#flirt -in aal.nii -ref T1_brain.nii.gz -omat aal2t1.mat -out aal2t1_affine.nii.gz -dof 12 -interp nearestneighbour
#flirt -in aal2t1_affine.nii.gz -ref b0_mean.nii.gz -applyxfm -init t12b0.mat -interp nearestneighbour -out aal2dwi.nii.gz

# first convert the Schaefer atlas from MNI space to subject space
flirt -in Schaefer2018_100_1mm.nii.gz -ref T1_brain.nii.gz -omat Schaefer2t1.mat -out Schaefer2t1_affine.nii.gz -dof 12 -interp nearestneighbour
flirt -in Schaefer2t1_affine.nii.gz -ref b0_brain.nii.gz -applyxfm -init t12b0.mat -interp nearestneighbour -out Schaefer2t1_b0.nii.gz


#### Step 6 - Generate connectome ####
tck2connectome -symmetric -zero_diagonal -scale_invnodevol \
               -tck_weights_in sift2_weights.txt \
               tracks_10M.tck aal2dwi.nii.gz \
               SC_AAL116.csv -out_assignment assignments_AAL116.csv -force

# Generating connectome (ROI-to-ROI connectivity matrix) for whole brain
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift_1M.txt tracks_10M.tck Schaefer2t1_b0.nii.gz SC_Schaefer_116.csv -out_assignment assignments.csv -force



### Step 7- Create specific region (e.g., Thalamus) tractography
fslsplit T1_first_all_fast_origsegs.nii.gz
fslmaths vol0006.nii.gz -add vol0013.nii.gz -bin thalamus_mask # specifay spesified ROIs .nii file here
#fsleyes thalamus_mask.nii.gz 
flirt -in thalamus_mask.nii.gz -ref b0_brain.nii.gz -applyxfm -init t12b0.mat -interp nearestneighbour -out thal2b0.nii.gz 
tckedit -include thal2b0.nii.gz tracks_10M.tck thalamus.tck

##To limit the number of streamlines through a seeded region; N=number of streamlines eg.500
tckedit -include ROFC_rois_combined.nii -number 2000 tracks_top20per_density.tck output_subset.tck -force
#tckedit -include /Users/snray/Documents/data_lab/DTI_analysis/DP01/ANAL/combined_LACC_seeds.nii -include /Users/snray/Documents/data_lab/DTI_analysis/DP01/ANAL/combined_RACC_seeds.nii  -number 2000 tracks_top20per_density.tck output_subset_ACC.tck -force

### View the results in MRview or TrackVis or DTIstudio or Dive
cd ..
done

