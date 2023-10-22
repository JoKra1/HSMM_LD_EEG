# HSMM_LD_EEG

## Description:

Code to perform the HsMM-MVPA model selection and trial-level GAMM analysis described in *Word type and frequency effects on lexical decisions are process-dependent and start early* (Krause, van Rij, & Borst; submitted). The code
makes use of the HsMM functions provided by Berberyan et al. (2021) that can be downloaded from [here](https://osf.io/z49me/files/).

## Performing the analysis:

 - Clone the repository
 - Copy the HsMM functions folder into the repository
 - Open Matlab and make sure the repository is the working/current directory
 - Install eeglab, the parallel computing, and wavelet toolboxes if they are not yet installed
 - Run ``HSMMMVPA_forward_selection.mat``
 - Open R and make sure the repository is the working/current directory, (optionally) install the dependencies via renv
 - Run ``trial_level_analyis.Rmd``