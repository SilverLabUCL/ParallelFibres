% This script calculates a histogram of the distribution of Cn values 
% (values of the correlation image) over all datasets to find a threshold 
% for identifying candidate varicosities

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
             'FL95_180425_10_53_40',... 8
             'FL87_180413_11_00_55',... 9
             'FL87_180117_11_23_20',... 10
             'FL_S_171109_14_54_34',... 11
             'FL_S_171109_15_19_52',... 12
             'FL76_170913_10_57_06'}; % 13
        
% Iterate over all datasets and all patches
        
Cn_all = [];

for dataset_ix = 1:13
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
