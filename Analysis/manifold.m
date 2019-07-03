% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;

%%

dataset_ix = 8;
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[coeff, score] = pca(dFF');

figure, plot(zscore(whisk_angle))
hold on, plot(zscore(score(:,1))+3)
hold on, plot(zscore(score(:,2))+6)
hold on, plot(zscore(score(:,3))+9)
