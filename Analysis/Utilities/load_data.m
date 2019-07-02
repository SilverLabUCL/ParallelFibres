%
% This function loads data from all patches of a given experiment
% and concatenates it 
%
% Input:
%    dataset_ix       Dataset number
% 
% Output:
%    Ain          
%    dFF          

function [dFF_all,acquisition_rate] = load_data(dataset_ix,grouped)

    basedir = '~/Documents/ParallelFibres/Data/';

    datasets = {'FL87_180501_11_03_09',...  1
                'FL87_180501_10_47_25',...  2
                'FL87_180501_10_36_14',...  3
                'FL87_180220_10_38_55',...  4
                'FL77_180213_10_46_41',...  5
                'FL_S_170906_11_26_25',...  6
                'FL_S_170905_10_40_52',...  7
                'FL45_170125_14_47_04'}; %  8

    fname = datasets{dataset_ix};
    disp(fname)
    
    load([basedir,fname,'/processed/',fname,'_GroupedData.mat'])
    
    Numb_patches = size(Ain_axons,1);
    T = size(dFF_axons{1},2);

    if grouped
        dFF_all = vertcat(dFF_axons{:});
    elseif ~grouped
        dFF_all = vertcat(dFF_rois{:});
    end
    
%     
%     % Get total number of axons after grouping
%     Num_axons = 0;
%     for patch_no = 1:Numb_patches
%         Num_axons = Num_axons + size(Ain_axons{patch_no},2);
%     end
%     
%     % Add all axons to list of dFF
%     dFF_all = zeros(Num_axons,T);
%     count = 1;
%     for patch_no = 1:Numb_patches
%         dFF_all(count:count + size(dFF_axons{patch_no},1) - 1,:) ...
%             = dFF_axons{patch_no};
%         count = count + size(dFF_axons{patch_no},1);
%     end
    

