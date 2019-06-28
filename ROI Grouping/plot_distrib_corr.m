% This function compares the distribution of SNRs of ROIs with the SNRs of
% spurious ROIs in the neuropil to determine a threshold for ROI detection

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
            'FL45_170125_14_47_04'}; %  8

C_intra_patch = cell(8,1);
C_inter_patch = cell(8,1);

%%
for dataset_ix = 1:8
    fname = datasets{dataset_ix};
    fname_SNR = [basedir,fname,'/processed/SNR.mat'];

    disp(fname)

    load([basedir,fname,'/',fname,'.mat'],'Numb_patches')
    load(fname_SNR,'SNR_thresh');
    
    % smoothing window
    if ismember(dataset_ix,[4,5,6,7])
        smooth_win_s = 0.35;
    else
        smooth_win_s = [];
    end
    
    dFF = cell(Numb_patches,1);
    
    for patch_no = 1:Numb_patches
        patch_no

        % Load data 
        load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
        Y = double(Y);

        % Use CNMF initialization to estimate initial spatial filters
        [Ain,Cn] = detect_ROIs(Y, [d1,d2]); 
        
        % Remove low SNR ROIs and trim
        [Ain,~,~,~] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,SNR_thresh); 
        Ain = trim_ROIs(Ain,[d1,d2]);
        [Ain,~,~,~] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,SNR_thresh); 

        % Calculate dFF
        dFF{patch_no} = get_dFF(Ain,Y,acquisition_rate,smooth_win_s);
        
        % Get correlation within patch
        C = corrcoef(dFF{patch_no}');
        C = C(triu(true(size(C)),1));
         
        % Save to correlations within each patch
        C_intra_patch{dataset_ix} = [C_intra_patch{dataset_ix}; C];
    end
    
    % Compute cross correlations between patches
    for patch_no1 = 1:Numb_patches
        for patch_no2 = (patch_no1+1):Numb_patches
            
            % Get correlation between patches
            C = corr(dFF{patch_no1}',dFF{patch_no2}');    
            
            % Save to correlations between patches
            C_inter_patch{dataset_ix} = [C_inter_patch{dataset_ix}; C(:)];
            
            [patch_no1,patch_no2,max(C(:))]
        end
    end
end

%%  Plot distributions and calculate thresholds, and save data

dbin = .02;
bins = -1:dbin:1;
bin_c = bins(1:end-1) + dbin/2;

SNR_thresh = zeros(8,1);

for dataset_ix = 1:8

    fname = datasets{dataset_ix};
    
    C_intra = C_intra_patch{dataset_ix};
    C_inter = C_inter_patch{dataset_ix};
    
    [dataset_ix,prctile(C_inter,95)]
    
    figure, hold on,
    histogram(C_intra,bins,'Normalization','probability');
    histogram(C_inter,bins,'Normalization','probability');
    set(gca,'FontSize',18), xlabel('Correlation'), ylabel('probability')
        
    title(fname,'Interpreter','None')
    
    save([basedir,fname,'/processed/corr_histograms.mat'],'C_intra','C_inter');
end
