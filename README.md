# EEGTMSfMRI
Code repository for simultaneous EEG-TMS-fMRI project

Created Jordan Muraskin @ Columbia University

Scripts and functions to be used on EEG data collected with custom-built MRI compatible amplifier (EP-link). 
Requires Matlab and EEGLAB, data can be made available upon request (jsm2112@columbia.edu)

ScriptToLoadandFilter.m 
  -- Main script to preprocess EEG data collected outside of MRI scanner with amplifier. 
  
import_eplink_NoGa.m
  -- function that imports raw EEG data into EEGLAB format **does not run** gradient artifact removal. 
  
read_amp_binaryv3.m
  -- function called by import_eplink_NoGa that reads in the raw binary EEG data files in the custom format from EP-link system.
  
shortestpath.m
  -- rereferencing function for custom EEG-fMRI cap
  
efmri34.ced
  -- electrode locations for custom EEG-fMRI cap. 
