% This script compares the distribution of SNRs of ROIs with the SNRs of
% spurious ROIs in the neuropil to determine a threshold for ROI detection
% (set to the 95th percentile of the distribution of SNRs of spurious ROIs)

% This script generates the following variables for each experiment:
%    SNR_thresh     threshold SNR value for the SNR for that specific experimen
%    SNR            vector of SNRs of all detected ROIs
%    SNR_neuropil   vector of SNRs of spurious ROIs in the neuropil

clear all; clc; close all

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

SNR_all = cell(size(datasets,2),1);
SNR_np_all = cell(size(datasets,2),1);

%%
for dataset_ix = 13
    fname = datasets{dataset_ix};
    disp(fname)

    load([basedir,fname,'/',fname,'.mat'],'Numb_patches')
    SNR = cell(Numb_patches,1);
    SNR_np = cell(Numb_patches,1);

    for patch_no = 1:Numb_patches
        patch_no
        
        % Load data
        load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
        Y = double(Y);

        % Use CNMF initialization to estimate initial spatial filters
        [Ain,Cn] = detect_ROIs(Y, [d1,d2]); 
        
        % Get SNR for all ROIs
        [~,SNR{patch_no},~,~] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,0); 
        
        % Calculate fake ROIs in neuropil
        A_neuropil = get_neuropil_filters(Ain,Cn,[d1,d2]);
        
        % Get SNR for fake ROIs
        [~,SNR_np{patch_no},~,~] = remove_bad_cells(A_neuropil,Y,[d1,d2],acquisition_rate,0); 
        
    end
    SNR_all{dataset_ix} = SNR;
    SNR_np_all{dataset_ix} = SNR_np;
end

%%  Plot distributions and calculate thresholds, and save data

dbin = .5;
bins = 0:dbin:30;
bin_c = bins(1:end-1) + dbin/2;

for dataset_ix = 13

    fname = datasets{dataset_ix};
    
    SNR = SNR_all{dataset_ix};
    SNR_neuropil = SNR_np_all{dataset_ix};
    
    figure, hold on,
    histogram(vertcat(SNR{:}),bins,'Normalization','probability');
    histogram(vertcat(SNR_neuropil{:}),bins,'Normalization','probability');
    set(gca,'FontSize',18), xlabel('SNR'), ylabel('probability')
        
    title(fname,'Interpreter','None')
    
    % Threshold for accepting candidate varicosities
    SNR_thresh = prctile(vertcat(SNR_neuropil{:}),95);
    
    save([basedir,fname,'/processed/SNR.mat'],'SNR_thresh','SNR','SNR_neuropil');
end
