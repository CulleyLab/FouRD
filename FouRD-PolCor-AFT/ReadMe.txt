***FouRD - Fourier Ring Depolarisation***
MATLAB implementation v1.0 01/11/23

*Program to measure, calibrate and correct fluorescence polarisation effect in microscopy images of fibrous structures. Submitted with the paper "Characterisation and correction of polarisation effects in fluorescently labelled fibres" to the Journal of Microscopy 3/11/23.*

***FouRD implemented in the file FouRD_PolCor.m***	Richard J. Marsh

***Recursive Batch AFT implemented in AFT_superbatch.m***	Richard J. Marsh

***All other files are a modified version of "Alignment by Fourier Transform (AFT)" ***
The originals are taken from the paper "A Workflow for Rapid Unbiased Quantification of Fibrillar Feature Alignment in Biological Images" Front. Comput. Sci 14 October 2021, Stefania Marcotti et.al.

DOI: https://doi.org/10.3389/fcomp.2021.745831
Original software available at https://github.com/OakesLab/AFT-Alignment_by_Fourier_Transform
See AFT_README.md for instructions for use and how to run as an independent program.

For the AFT components all rights and permissions belong to: Patrick Oakes poakes@gmail.com

***FouRD installation***
You will Need MATLAB to run this application
It was developed and tested on MATLAB version 2019a
Download the FouRD-PolCor-AFT folder to a local directory
Optionally download the DemoData folder for some data to try.
Open MATLAB, then open the script FouRD_PolCor.m
if prompted, click 'add to path'

***FouRD operation***
To launch FouRD open and run the script FouRD_PolCor.m in MATLAB
This will show the Main user interface with several buttons that launch the component functions as follows:

1) "Anisotropy from image pair" :
This will generate an anisotropy spectrum from two polarisation resolved images of the same structure.
Parameters:
G-factor : number to multiply second polarisation image by to account for any differences in detection efficiency.

2) "Batch AFT"
This launches the AFT analysis on all the image files in the selected folder. In addition to the standard output, two 3D MATLAB matrices are saved in the directory - Anglemat.mat & Eccentricitymat.mat. These contain the calculated angles and eccentricities where the 3rd dimension corresponds to each file in the directory.
Parameters:
Window size: size of the boxed sub region of the image that is used for each AFT calculation.
Window overlap: how much each box overlaps each other.
Neighbourhood radius: Neighbourhood size on which the order parameter is calculated.
For a more detailed description of AFT parameters and options see the AFT_README.md file

3) "Recursive Batch AFT"
This is the same as Batch AFT but will recurse sub directories in the chosen directory. Each subdirectory is treated as a separate measurement and statistics on the results of each subdirectory are produced.
Parameters:
Same as Batch AFT above

4) "Make Calibration"
This produces a calibration spectrum from an unaligned sample. The resultant calibration will be averaged over all the files in the chosen directory. This will normally correspond to a Z-stack or a single image. Calibration data is automatically saved in this directory.
Parameters:
Max pixels X&Y: The size of the smallest image in the folder. Normally all images should be the same size but if there are small differences (eg. due to stitching) all images will be cropped to this size.
Number of tiles X&Y: If performing a tiled correction, this is the number of tiles the calibration image is divided into in each dimension. (1 for an untiled correction)
StartRing: the radius of the first Fourier ring (spatial frequency) to use in the correction.

5) "Apply Correction"
This will correct polarisation bias in all the images in a selected folder (usually a z stack) using some previously saved calibration data from (4) above. All parameters are automatically saved and corrected images stored in a separate folder, the name of which also contains the parameters.
Parameters:
Max pixels X&Y: Same as for Make Calibration above. Use the same numbers here.
StartRing: the radius of the first Fourier ring (spatial frequency) to use in the correction.
Rings to smooth: Applies a sliding average smoothing of this many points to the calibration data anisotropy spectra before performing the correction. Useful if very noisy data produce large spikes in the spectra.
 
6) "Spectrum from image"
This will produce an anisotropy spectrum from a single non-polarisation image. This is useful to see if an unaligned sample has polarisation bias.
Parameters: None

7) "Spectrum from Calibration"
This will produce an anisotropy spectrum from previous calibration data showing the degree of polarisation bias should this be of interest. If a tilled calibration is performed it can show a spectrum for each image tile.
Parameters:
Rings to smooth: As in Apply Correction above but this time only for display purposes.
Tile? : Y or N . whether to display a separate spectrum for each tile in a tiled calibration data
