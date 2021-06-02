# Analysis

## Summary
This folder contains scripts for analysing bidirectionality, spatial correlations, manifold structure, and dimensionality. The folder `Utilities/` includes additional helper functions.

## Data requirements
Three data files per experiment are required to run the code in this folder. All files can be found at (Figshare link).

The first file (`ExperimentName.mat`) contains metadata and behavioural data. The relevant variables used in the analysis are:
* `aquisition_rate` is a scalar describing the acquisition rate in Hz.
* `Pixel_size` is a scalar describing the size of each pixel in microns.
* `Numb_patches` is a scalar describing the number of patches in the experiment.
* `Patch_coordinates` is a structure containing coordinate information about each patch. `Patch_coordinates.data` is a matrix in which each row represents a patch, and columns 5, 6, and 7 represent the X, Y, and Z positions (respectively) of that patch.
* `SpeedDataMatrix` and `SpeedTimeMatrix` are vectors containing the wheel speed time series and times from the wheel encoder.
* `dlc_whisk_angle` and `dlc_whisk_time` are vectors containing the whisking angle time series and times as determined via DeepLabCut.
* `wheel_MI` is a matrix whose second column contains the wheel motion index time series as determined from the wheel cameras and whose second column contains the corresponding times.
Note that this file may also contain variables extracted by now obsolete methods which were not included by the analysis in the paper (e.g., `Whiskers_angle_0` for old whisker position detection,`Axon_dFF` for old grouping procedure). You can ignore these.

The second file (`ExperimentName_GroupedData.mat`) contains the functional data after grouping.
* `Cn` is a cell of length `Numb_patches`. Each element of the cell is a d1 x d2 matrix of the correlation image of a different patch in the experiment, where d1 and d2 are the dimensions of the patch. The correlation image is calculated as the correlation of each pixel with its neighbours.
* `Ain_axons` is a cell of length `Numb_patches`. Each element of the cell is a (d1 * d2) x M matrix corresponding to the spatial filters (masks) of each putative axon after grouping, where M is the number of identified axons in that patch. To visualize, use `reshape` to convert a column of the matrix `Ain_axons{patch_no}` into a d1 x d2 matrix corresponding to the patch dimensions, e.g., `imagesc(reshape(Ain_axons{patch_number}(:,axon_number),d1,d2))`.
* `dFF_axons` is a cell of length `Numb_patches`. Each element of the cell is an M x T matrix corresponding to the dFF of each putative axon after grouping, where T is the number of timepoints in the experiment. The mth row of `dFF_axons{patch_no}` and the mth column of `Ain_axons{patch_no}` corresponds to the same putative axon.
* `time_axons` is a cell corresponding to the times for each axon, same structure as `dFF_axons`.
* `dFF_rois`, `Ain_rois`, and `time_rois` as above, but for individual granule cell varicosities before grouping. 
* `ix_axons_to_rois` is a cell of length `Numb_patches`. Each element of the cell is another cell of length M corresponding the the ROI id numbers for each varicosity associated with that putative axon. The mth element of this cell is a vector containing the ROI id numbers of all varicosities associated with the putative axon corresponding to axon id number m.
* `axon_ids` is a vector whose length is the number of varicosities, which performs the inverse mapping (from ROI id number to axon id number). The ith element of this vector is the axon id number for the varicosity corresponding to ROI id number i.

The final file (`ExperimentName_fibre_direction.mat` in Figshare, but in the code it is just called `fibre_direction.mat` within the folder for that experiment) contains information about the overall fibre orientation as determined from tracing small secctions from structural data.
* `vector_mean` is a normalized vector containing the average axon direction over all traced segments in the patch.
* `angle` is a vector listing the deflection of each individual traced axon segment from `vector_mean` (in radians).
* `angle_std` is the standard deviation of angle (in degrees).

## Overview of scripts

`define_dirs.m` defines the path and folder names for each experiment.

`bidirectional.m` identifies positively and negatively modulated parallel fibres (PFs), correlates to behaviour, and plots example traces.  This script generates Figures 1e, 2a,b and Extended Data Figures 3, 4b, 5, 6.

`bidirectional_soma.m` identifies positively and negatively modulated PFs based on somatic data for Supplementary Figure 2c.

`spatial_structure.m` calculates correlations between PFs, and plots as as a function of distance to look for spatial structure. This script generates Figures 2c,d and Extended Data Figure 4c,d, 7.

`manifold.m` contains scripts to calculate manifold distances and angles, generating Figure 3b,c,e,f.

`orthogonal_prctle.m` calculates angle between manifolds while removing increasing percentiles of positively and negatively modulated PFs, to generate Figure 3g.

`plot_different_behavioural_spaces_byhand.m` shows example manifolds with hand-labelled periods of isolated whisking, to generate Supplementary Figure 3.

`regression.m`contains scripts to regress PF activity against behavioural variables. This code generates Figures 4, 5 and Extended Data Figre 9.

`dimensionality.m` contains scripts to estimate the dimensionality of PF activity. This is used to generate Figure 6 and Extended Data Figure 10.

`behaviour_pca_video.m` contains the script to generate video S2.

`rotate_manifold_video.m` contains the script to generate video S3.
