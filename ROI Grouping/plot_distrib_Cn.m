%% Plot histogram of correlation image value for each pixel,
% to determine threshold 'min_Cn' for ROI detection

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


Cn = cell(1,1);

for dataset_ix = 1:8
    fname = datasets{dataset_ix};
    disp(fname)

    Cn_temp_all = [];
    patch_no = 1;

    while exist([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'],'file') == 2

        patch_no

        load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])

        Y = double(Y);

        [~,Cn_temp] = get_spat_filtered_data(Y,[d1,d2],5,2);

        Cn_temp_all= [Cn_temp_all; Cn_temp(:)];

        patch_no = patch_no + 1;

    end
    figure, histogram(Cn_temp_all), 
    title(fname,'Interpreter','None'), drawnow
    Cn{dataset_ix} = Cn_temp_all;
end 

figure, histogram(vertcat(Cn{:}))
set(gca,'FontSize',15)
xlabel('Cn'), ylabel('Number')
title('All data')


