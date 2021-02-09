# preclinical-reconstruction
Preclinical reconstruction software for the CCHMC CPIR.

## Installation
Make sure all files are on the Matlab path.

## Usage
Call reconstruction via the recon.m function, typically using the Matlab command line. If you provide a path to the folder containing data to be reconstructed, the reconstruction will do the rest, otherwise it will prompt you to select a folder in a ui. Data folder must contain method and fid files. To run through many reconstructions hands off, make sure to generate configuration files in each data folder as well.

## Reconstruction Pathway
The steps taken by the reconstruction code are as follows:
* Make sure supplied path contains method and fid files
* Check whether a configuration file is present - if not there, prompt user to generate file. Then, read in file
* Read in method file which also calculates trajectories
* Based on whether imaging a mouse or phantom (info in configuration file), set the triggering parameter so that phantom reconstructions will work properly
* Based on information from the method file, pass path, trajectories, method parameters, and configuration parameters to either the radial_recon function or the spiral_recon function
* In recon functions, first apply trajectory delay correction - ultimately want to have these calibrated to proper values for both xenon and proton, theoretical and measured
* Read in and reshape data based on information from the method file
* Delete Zero-filled points and acquisition shift points
* If config file says it's a mouse image and trigger is off, run retrospective gating - would love for this to be totally hands off and generate both inspiratory and expiratory images automatically
* If config file says it's a phantom or it's a mouse with triggered acquisition, run through a normal reconstruction
* Loop through all slices, Echo Times, b values, and repetitions and reconstruct each image, storing the images and k-space data for everything
* Display reconstructed image
* Write out Reconstructed and Raw Data to file in structures
* Back in original recon function, write out when and by whom the reconstruction was done along with important imaging and reconstruction parameters.
* Finish by closing all open files so as not to spam up matlab with a bunch of open method, fid, and text files in the background.
