% This function plots a histogram of the SNR of all ROIs identified in the
% correlation image of each plot.

clear all; clc; close all
addpath('CNMFE_code/')

basedir = '~/Documents/FredPF/raw/offMC/';
datasets = {'FL87_180501_11_03_09',...  1
            'FL87_180501_10_47_25',...  2
            'FL87_180501_10_36_14',...  3 
            'FL87_180220_10_38_55',...  4
            'FL77_180213_10_46_41',...  5
            'FL_S_170906_11_26_25',...  6
            'FL_S_170905_10_40_52',...  7            
            'FL45_170125_14_47_04'}; %  8

% Choose which dataset
dataset_ix = 8;

fname = datasets{dataset_ix};
disp(fname)

load([basedir,fname,'/',fname,'.mat'],'Numb_patches')

SNR = cell(Numb_patches,1);

for patch_no = 1:Numb_patches
    patch_no

    load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])

    Y = double(Y);

    % Use CNMF initialization to estimate initial spatial filters
    [Ain,Cn] = detect_ROIs(Y, [d1,d2], 1, 5); 
    close(1) % Closes thresholded ROI figure

    % Remove low SNR ROIs
    [~,SNR_good,~,SNR_bad] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,6); 
    close % Closes thresholded high SNR ROI figure
    
    SNR{patch_no} = [SNR_good; SNR_bad];
    
    pause
    close % Closes correlation image
end

figure, histogram(vertcat(SNR{:}))
set(gca,'FontSize',15)
xlabel('SNR'), ylabel('Number')
title(fname,'Interpreter','None')

saveas(gcf,[basedir,fname,'/figs/SNR.png'])

fileID = fopen([basedir,fname,'/processed/SNR.txt'],'w');
fprintf(fileID,'%d ',vertcat(SNR{:}));
