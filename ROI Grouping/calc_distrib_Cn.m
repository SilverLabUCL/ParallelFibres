% This script calculates a histogram of the distribution of Cn values after
% smoothing (over all datasets) to find a threshold for ROI rejection

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
            'FL45_170125_14_47_04',...  8
             ...%%
             'FL95_180425_10_53_40',...  9
             'FL87_180413_11_00_55',...  10
             'FL104_180725_10_42_37',... 11
             'FL87_180117_11_23_20',...  12
             'FL106_180807_10_52_25',... 13
             'FL_S_171109_14_54_34',...  14
             'FL_S_171109_15_19_52',...  15
             ...%
             'FL75_170912_10_33_29',... 16
             'FL76_170913_10_57_06',... 17
             'FL77_180113_10_58_50'};  %18
       		% 'FL92_180228_11_10_48'
		    % 'FL92_180228_11_18_24'

        
% Iterate over all datasets and all patches
        
Cn_all = [];

for dataset_ix = 1:8
    fname = datasets{dataset_ix}; disp(fname)
    load([basedir,fname,'/',fname,'.mat'],'Numb_patches')
    for patch_no = 1:Numb_patches
        [dataset_ix, patch_no]

        % Load data
        load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
        fname_fibreangles = [basedir,fname,'/processed/fibre_direction.mat'];
        Y = double(Y);

        % Calculate smoothed correlation image and save pixel values
        [~,Cn] = detect_ROIs(Y, [d1,d2], 2, 5);
        Cn_all = [Cn_all; Cn(:)];
    end
end

%% Plot distribution
h = histogram(Cn_all);
figure, plot(h.BinEdges(1:end-1),h.Values,'k','LineWidth',2)
set(gca,'FontSize',18)
ylabel('Number of pixels')
xlabel('Smoothed correlation image value')
savefig([basedir,'/figs/histogram_Cn.fig'])
