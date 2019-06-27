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

varratio_intra_patch = cell(8,1);
varratio_inter_patch = cell(8,1);

%%
for dataset_ix = 1:8
    fname = datasets{dataset_ix};
    fname_SNR = [basedir,fname,'/processed/SNR.mat'];
    fname_corr = [basedir,fname,'/processed/corr_histograms.mat'];

    disp(fname)

    load([basedir,fname,'/',fname,'.mat'],'Numb_patches')
    load(fname_SNR,'SNR_thresh');
    
    load(fname_corr,'C_inter');
    rho_min = prctile(C_inter,95);
    
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
        dFF{patch_no} = get_dFF(Ain,Y,acquisition_rate);
        
        % Get correlation within patch
        C = get_corrs(dFF{patch_no},rho_min);
        [a,b] = find(C>0);
        ratio_temp = zeros(length(a),1);
        for k = 1:length(a)
            ratio_temp(k) = get_var_ratio(dFF{patch_no}(a(k),:),dFF{patch_no}(b(k),:));
        end
         
        % Save to correlations within each patch
        varratio_intra_patch{dataset_ix} = [varratio_intra_patch{dataset_ix}; ratio_temp];
    end
    
    % Compute cross correlations between patches
    for patch_no1 = 1:Numb_patches
        for patch_no2 = (patch_no1+1):Numb_patches
            
            % Get correlation between patches
            C = corr(dFF{patch_no1}',dFF{patch_no2}');    
            C = max(C,rho_min);  % remove any corerlations less than rho_thresh
            C(C==rho_min) = 0; % set min corrs to 0 

            [a,b] = find(C>0);
            ratio_temp = zeros(length(a),1);
            for k = 1:length(a)
                ratio_temp(k) = get_var_ratio(dFF{patch_no1}(a(k),:),dFF{patch_no2}(b(k),:));
            end

            % Save to correlations between patches
            varratio_inter_patch{dataset_ix} = [varratio_inter_patch{dataset_ix}; ratio_temp];
            
            [patch_no1,patch_no2]
        end
    end
end

%%  Plot distributions and calculate thresholds, and save data

dbin = .1;
bins = 0:dbin:10;
bin_c = bins(1:end-1) + dbin/2;

for dataset_ix = 1:8

    fname = datasets{dataset_ix};
    
    varratio_intra = varratio_intra_patch{dataset_ix};
    varratio_inter = varratio_inter_patch{dataset_ix};
    
    figure, hold on,
    histogram(varratio_intra,bins,'Normalization','probability');
    histogram(varratio_inter,bins,'Normalization','probability');
    set(gca,'FontSize',18), xlabel('Var ratio'), ylabel('probability')
        
    title(fname,'Interpreter','None')
    
   save([basedir,fname,'/processed/varratio_histograms.mat'],'varratio_intra','varratio_inter');
end
