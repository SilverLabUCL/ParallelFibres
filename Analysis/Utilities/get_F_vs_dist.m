%
% For each dataset, this script finds each ROI and calculates average F
% over time & over space centered around ROI center, to get spatial decay
%
% Input:
%    dataset_ix       Dataset number
% 
% Output:       
%    rho_all
%    distances_all

function [F_all,distances_all] = get_F_vs_dist(dataset_ix)

    define_dirs;

    fname = datasets{dataset_ix};
    disp(fname)
    
    load([basedir,fname,'/processed/',fname,'_GroupedData.mat'])
    load([basedir,fname,'/',fname,'.mat'],'Pixel_size','Numb_patches');
    
    distances_all = [];
    F_all = [];
    
    % Check a 10 um x 10 um square surrounding the centroid
    num_pix_to_check = ceil(5 / Pixel_size);

    for p = 1:Numb_patches
        p
        
        % Load patch data
        load([basedir,fname,'/raw/Patch',sprintf('%03d',p),'.mat'])
        Y = double(Y);
        
        % Get average F for each patch
        Y = mean(Y,2);
        Y = reshape(Y,d1,d2);
        
        % Number of ROIs in this patch
        N = size(Ain_rois{p},2);
        
        % Go through each roi in the patch
        for roi = 1:N
            
            % Get centroid of the ROI
            Ain = Ain_rois{p}(:,roi);
            c = regionprops(reshape(Ain,d1,d2),'centroid'); 
            c = round(c.Centroid); 
            c2 = c(1); c1 = c(2);
            
            % Normalize by avg. F at centroid
            Y0 = Y(c1,c2);
            
            min_pix1 = max(1,c1-num_pix_to_check);
            max_pix1 = min(d1,c1+num_pix_to_check);
            min_pix2 = max(1,c2-num_pix_to_check);
            max_pix2 = min(d2,c2+num_pix_to_check);
            
            for pix1 = min_pix1 : max_pix1
                for pix2 =  min_pix2 : max_pix2
                    dist = sqrt((pix1-c1)^2 + (pix2-c2)^2) * Pixel_size;
                    if dist <= 5
                        F_norm = Y(pix1,pix2)/Y0;
                        F_all = [F_all; F_norm];
                        distances_all = [distances_all; dist];
                    end
                end
            end

        end
        
        
    end
    
    
   