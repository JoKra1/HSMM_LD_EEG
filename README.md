# HSMM_LD_EEG

## Description:

Code to perform the HsMM-MVPA model selection and trial-level GAMM analysis described in *Word type and frequency effects on lexical decisions are process-dependent and start early* (Krause, van Rij, & Borst; submitted). The code
makes use of the HsMM functions provided by Berberyan et al. (2021) that can be downloaded from [here](https://osf.io/z49me/files/). Additionally, the code requires a copy of the item-level DLP data (Keuleers et al., 2010), which can for example be downloaded [here](https://lib.ugent.be/catalog/pug01:1076200). The ``vwr`` R package is required for parts of the trial-level analysis but appears to no longer be on CRAN - the most recent version can still be downloaded from [here](https://CRAN.R-project.org/package=vwr) and then installed locally.

## Performing the analysis:

 - Clone the repository
 - Copy the HsMM functions folder into the repository
 - Copy the ``dlp-items.Rdata`` and ``dlp-stimuli.Rdata`` files into the data folder
 - Open Matlab and make sure the repository is the working/current directory
 - Install eeglab, the parallel computing, and wavelet toolboxes if they are not yet installed
 - Run ``HSMMMVPA_forward_selection.mat``
 - Open R and make sure the repository is the working/current directory, (optionally) install the dependencies via renv
 - Run ``trial_level_analyis.Rmd``