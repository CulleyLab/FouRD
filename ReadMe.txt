***FouRD - Fourier Ring Depolarisation***
MATLAB implimentation v1.0 01/11/23

*Program to measure, calibrate and correct fluorescence polarisation effect in microscopy images of fiberous structures. Submitted with the paper "Characterisation and correctio of polarisation effects in fluorescently labelled fibres" to the Journal of Microscopy.*

***FouRD implemented in the file FouRD_PolCor.m***

***All other files are a modified version of "Alignment by Fourier Transform" for the paper "A Workflow for Rapid Unbiased Quantification of Fibrillar Feature Alignment in Biological Images" Front. Comput. Sci 14 October 2021 Stefania Marcotti et.al.
DOI: https://doi.org/10.3389/fcomp.2021.745831
Original version avalible at https://github.com/OakesLab/AFT-Alignment_by_Fourier_Transform ***

***FouRD operation***
To lanch FouRD run the scrit FouRD_PolCor in MATLAB
This will show the Main user interface with several buttons that launch the component functions as follows

1) "Anisotropy from image pair" :
This will generate an anisotropy spectrum from two polarisation resolved images of the same structure.

2) "Batch AFT"
This launches the AFT analysis on all the image files in the selected folder. In addition to the standard output two 3D MATLAB matricies are saved in the directory - Anglemat.mat & Eccentricitymat.mat. These containe the calculated angles and excentricities where the 3 dmention coresponds to each file in the directory

3) "Recursive Batch AFT"
This is the same as Batch AFT but will recurse sub directories in the chosen directory. Each subdirectory is treated as a seperate measurment and statistics on the results of each subdirectory are produced.

4) "Make Calibration"
This produces a calibration spectrum from an unaligned sample. The resultant calibration will be averaged over all the files in the chosen directory. This will normally corespond to a Z-stack or a single image. Calibration data is automaticaly saved in this directory.

5) "Apply Corection"
This will correct polarisation bias in all the images in a selected folder (usually a z stack) using previously save calibration data from (4) above. All parameters are automaticaly saved and corected images stored in a seperate folder, the name of wwhich also contains the parameters.
 
6) "Spectrum from image"
This will produce an anisotropy specterum from a single non-polarisation image. This is usefull to see if an unaligned sample has polarisation bias.

7) "Spectrum from Calibration"
This will produce an anisotropy spectrum from previous calibration data showing the degree of polarisation bias should this be of interest. If a tilled calibration is performed it will show a spectrum for each image tile.