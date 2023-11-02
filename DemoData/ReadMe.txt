***Demo data to test operation of FouRD image polarisation software***

1) Demo Correction Sims
This folder contains the simulated images used for the demonstration corection in the assosiated publication "Characterisation and correction of polarisation effects in fluorescently labelled fibres".

The sub folder 'Calibration' contains the unaligned calibration sim. just select this folder when perfoming a calibration. This will leave a calibration data .mat file in this folder

The sub folder 'test' contains the alined & polarised image for correction. When performing a correction select this folder for the data and the calibration.mat file from above for the calibration.

Also in the main folder are the corected and Ground Truth images from the paper for comparison.


2) Example SoRa Frames
This folder contains a couple of example fromes form the SoRa microscope. Space prevents the upload of the full stacks.

The calibration folder again contains an example frame fro one stack to use as a calibration, select this folder

The Test folder contains to frames. One is an image from the same z stack as the calibration image but a different z slice. The other is a frame from a seperate z stack (independent measurment). Selecting this folder when applying the corection will automaticaly correct both files using the same calibration.


AFT casn then be used on each folder to analyse the degree of alignment and asses the performance of the correction.

