% This script plots spatial information that wasn't used for regrouping, as
% a sanity check to compare with existing literature on parallel fibre
% bouton information (Pichitpornchai et al., 1994; Shepherd et al., 2002)

%% First plot distribution of varicosity sizes 
%(before rejecting on basis of SNR etc)

clear all; clc
addpath('From CNMF_E/')
addpath('Utilities/')

basedir = '~/Documents/ParallelFibres/Data/';
datasets = {'FL87_180501_11_03_09',...  1
            'FL87_180501_10_47_25',...  2
            'FL87_180501_10_36_14',...  3 
            'FL87_180220_10_38_55',...  4
            'FL77_180213_10_46_41',...  5
            'FL_S_170906_11_26_25',...  6
            'FL_S_170905_10_40_52',...  7            
            'FL45_170125_14_47_04'}; %  8

% Choose dataset
dataset_ix = 1;

% Choose patch number
patch_no = 1;

fname = datasets{dataset_ix};
disp(fname)

while exist([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'],'file') == 2
    
    load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])

    Y = double(Y);

    [Ain,Cn] = detect_ROIs(Y, [d1,d2], 2, 5, 0); close; close;

end
    
%%
d = get_interbouton_dist(Ain,[d1,d2],ix_axons_to_rois,Pixel_size);
[mean(d),std(d)]