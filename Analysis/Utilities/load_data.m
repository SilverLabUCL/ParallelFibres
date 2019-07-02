%
% This function loads data from all patches of a given experiment
% and concatenates it 
%
% Input:
%    dataset_ix       Dataset number
% 
% Output:       
%    dFF          

function [dFF_all,time,acquisition_rate] = load_data(dataset_ix,grouped)

    if nargin < 2 || isempty(grouped)
        grouped = 1;
    end

    define_dirs;

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
    
    load([basedir,fname,'/',fname,'.mat'],'TimeAxon');
    time = mean(TimeAxon,2) / 1000;
    
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
    

