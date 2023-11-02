***Demo data to test operation of FouRD image polarisation software***

1) Demo Correction Sims
This folder contains the simulated images used for the demonstration corection in the assosiated publication "Characterisation and correction of polarisation effects in fluorescently labelled fibres".

The sub folder 'Calibration' contains the unaligned calibration sim. just select this folder when perfoming a calibration. This will leave a calibration data .mat file in this folder

For this data we use No tiles X and No tiles Y = 1 ie there is no need to do a tiled corection. All other parameters cn be default

The sub folder 'test' contains the alined & polarised image for correction. When performing a correction select this folder for the data and the calibration.mat file from above for the calibration.

Also in the main folder are the corected and Ground Truth images from the paper for comparison.


2) Example SoRa Frames
This folder contains a couple of example fromes form the SoRa microscope. Space prevents the upload of the full stacks.

The calibration folder again contains an example frame fro one stack to use as a calibration, select this folder

Parametes
the image size is 1408x1408 (for max image size)
For this data we recomend using No tiles X & No tile Y = 7.
For the start ring the default value is 2 which will leave some edge artefct for this data. We sugest a value of 5 but also to experiment with altering this value to gain an understanding of its efect

The Test folder contains to frames. One is an image from the same z stack as the calibration image but a different z slice. The other is a frame from a seperate z stack (independent measurment). Selecting this folder when applying the corection will automaticaly correct both files using the same calibration.


AFT casn then be used on each folder to analyse the degree of alignment and asses the performance of the correction.

3) Demo Polarisation Data
This folder contains two polarisation resolved images of the same structure for generating the anisotropy spectra from an image pair. They have been corrected for G factor so just use 1 for this value.