# Analysis

## Summary
This folder contains scripts for analysing bidirectionality, spatial correlations, manifold structure, and dimensionality. The folder `Utilities/` includes additional helper functions.

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
