# Multishell DTI Processing Pipeline

This repository contains a fully automated **multishell diffusion MRI preprocessing and connectome-generation pipeline** built using **MRtrix3**, **FSL**, and **ANTs**.

The pipeline performs:
- DWI denoising, Gibbs ringing correction, and motion/distortion correction using **TOPUP + EDDY**  
- Bias-field correction (ANTs) and brain extraction (SynthStrip)  
- Multi-tissue **FOD estimation (MSMT-CSD)** and normalization  
- Whole-brain and ROI-based **tractography (ACT + SIFT2)**  
- Atlas registration (AAL / Schaefer) to diffusion space  
- Generation of **structural connectivity matrices (SC.csv)**

### Usage
1. Organize your subject folders as described in the script.
2. Make the script executable:
   ```bash
   chmod +x run_dti.sh
   ./run_dti.sh


**Dependencies**
MRtrix3 ≥ 3.0.4
FSL ≥ 6.0
ANTs ≥ 2.3
FreeSurfer SynthStrip ≥ 1.0

**Citation**
1. Tournier et al., NeuroImage, 2019 (MRtrix3)
2. Andersson & Sotiropoulos, NeuroImage, 2016 (EDDY)
3. Dhollander et al., ISMRM, 2016 (MSMT-CSD)
4. Avants et al., NeuroImage, 2011 (ANTs)
