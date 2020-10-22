%%%% This script removes neuropil from all varicosities
% And resaves 

clear all; clc; 
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
             'FL87_180117_11_23_20',...  12
             'FL_S_171109_14_54_34',...  14
             'FL_S_171109_15_19_52',...  15
             ...%
             'FL76_170913_10_57_06',... 17
             'FL77_180113_10_58_50'};  %18

dataset_ix = 6; 

fname = datasets{dataset_ix};
disp(fname)

patch_no = 1;

load([basedir,fname,'/',fname,'.mat'],'Numb_patches');
load([basedir,fname,'/processed/',fname,'_GroupedData.mat']);

% Load data
disp([num2str(patch_no),' / ',num2str(Numb_patches)])
load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
Y = double(Y); T = size(Y,2);

% Plot mask for a particular neuropil
figure, imagesc(Cn{patch_no})
colormap(gray), axis equal, axis([0,d2,0,d1])
set(gca,'XTick',[],'YTick',[])
hold on, plot([20,20+5 / Pixel_size],[d1-5,d1-5],'w','LineWidth',2)

Ain_npl = get_neuropil_masks(Ain_rois{patch_no},Cn{patch_no},[d1,d2],Pixel_size);

[~,F_raw,F_neuropil] = compare_npl(Ain_rois{patch_no},Ain_npl,Y,acquisition_rate);

% Plot mask for a particular neuropil
roi = 2;
figure, imagesc(1-reshape(Ain_npl(:,roi),d1,d2))
colormap(gray), axis equal, axis([0,d2,0,d1])
set(gca,'XTick',[],'YTick',[])

%% Plot transients
N = size(F_raw,1);
figure, hold on
for k = 1:5
    plot((1:T)/acquisition_rate,F_raw(k,:)+600*k,'k')
    plot((1:T)/acquisition_rate,F_neuropil(k,:)+600*k-100,'b')
end
set(gca,'FontSize',15)
xlabel('Time (s)'), xlim([-20,430])
ylabel('Fluorescence')

%%
figure, plot(F_raw(:),F_raw(:)-F_neuropil(:),'.k')
xlabel('Raw fluorescence')
ylabel('Raw fluorescence - Neuropil')
set(gca,'FontSize',15)

corr(F_raw(:),F_raw(:)-F_neuropil(:))

%%

C = [];
V_raw = [];
V_npl = [];

C_npl = [];

ix = 1;

for dataset_ix = 1:15 

    fname = datasets{dataset_ix};
    disp(fname)

    load([basedir,fname,'/',fname,'.mat'],'Numb_patches');
    load([basedir,fname,'/processed/',fname,'_GroupedData.mat']);

    for patch_no = 1:Numb_patches

        % Load data
        disp([num2str(patch_no),' / ',num2str(Numb_patches)])
        load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
        Y = double(Y); T = size(Y,2);

        Ain_npl = get_neuropil_masks(Ain_rois{patch_no},Cn{patch_no},[d1,d2],Pixel_size);

        [~,F_raw,F_neuropil] = compare_npl(Ain_rois{patch_no},Ain_npl,Y,acquisition_rate);

        for roi = 1:size(F_raw,1)
            C(ix) =  corr(F_raw(roi,:)',F_raw(roi,:)'-F_neuropil(roi,:)');
            V_raw(ix) = var(F_raw(roi,:));
            V_npl(ix) = var(F_neuropil(roi,:));
            ix = ix+1;
        end

    end

end

figure, histogram(C)
set(gca,'FontSize',15)
xlabel('Correlation r')
ylabel('Number')
set(gca,'Box','off')

figure,  histogram(V_npl./V_raw)
set(gca,'FontSize',15)
xlabel('Neuropil variance / Raw fluorescence variance')
ylabel('Number')
set(gca,'Box','off')