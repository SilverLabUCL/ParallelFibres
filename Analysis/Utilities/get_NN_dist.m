%
% This function loads data from all patches of a given experiment
% and concatenates it 
%
% Input:
%    dataset_ix       Dataset number
% 
% Output:       
%    dFF_all
%    time
%    acquisition_rate
%    distances

function [NN_dist,NN_dist_ON,NN_dist_OFF,NN_dist_ON_shuff,NN_dist_OFF_shuff] = get_NN_dist(dataset_ix,grouped,euclid_dist)

    % If euclid_dist = 0, uses distance after projecting onto plane
    % orthogonal to average fiber direction
    if nargin < 3 || isempty(euclid_dist)
        euclid_dist = 0;
    end
    
    if nargin < 2 || isempty(grouped)
        grouped = 1;
    end
    
    define_dirs;

    fname = datasets{dataset_ix};
    disp(fname)
    
    % Load functional data - keeping specific to
    load([basedir,fname,'/processed/',fname,'_GroupedData.mat'])
    
    % Restrict only to the same patch
    if grouped
        dFF = dFF_axons;
        Ain = Ain_axons;
        time = mean(vertcat(time_axons{:}))'/1000;
    elseif ~grouped
        dFF = dFF_rois;
        Ain = Ain_rois;
        time = mean(vertcat(time_rois{:}))'/1000;
    end

    
    % Load behavioural data and get A/QW periods
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,0);
    
    
    % Load fibre direction, Z in pix, pix size, number of patches
    load([basedir,fname,'/processed/fibre_direction.mat'])
    load([basedir,fname,'/',fname,'.mat'],'Pixel_size','Numb_patches');

    % Patch size
    [d1,d2] = size(Cn{1});

    NN_dist = [];
    NN_dist_ON = [];
    NN_dist_OFF = [];
    NN_dist_ON_shuff = [];
    NN_dist_OFF_shuff = [];
    
    for p = 1:Numb_patches
        N = size(dFF{p},1);
        
        % Find ON and OFF GCs in this patch
        [change_dFF,p_val] = change_dFF_sig(dFF{p},A,QW,acquisition_rate);
        ix_ON = find(change_dFF > 0 & p_val < .05);
        ix_OFF = find(change_dFF < 0 & p_val < .05);
        N_ON = size(ix_ON,1);
        N_OFF = size(ix_OFF,1);

        % If projecting distances, project centroids now
        if ~ euclid_dist
            % Get vector of projected XY positions
            XYproj_um = zeros(1,N);
            X_um = zeros(1,N);
            Y_um = zeros(1,N);
            for n1 = 1:N
                % Centroid of axon 1
                Ain1 = reshape(Ain{p}(:,n1),d1,d2);
                c1 = regionprops(Ain1,'centroid'); c1 = c1.Centroid; 

                % vector orthogonal to vector_mean
                vector_orth = [-vector_mean(2), vector_mean(1)];
                if round(norm(vector_orth),10) ~=1 || dot(vector_mean,vector_orth)~=0
                    error('Problem with vector_orth.')
                end

                % Project axon 1 onto vector_orth and convert to um
                XYproj_um(n1) = dot(vector_orth,c1) * Pixel_size;

                X_um(n1) = c1(1) * Pixel_size;
                Y_um(n1) = c1(2) * Pixel_size;
            end
        end

        distances = zeros(N,N);
        for n1 = 1:N
            if euclid_dist
                Ain1 = reshape(Ain{p}(:,n1),d1,d2);
                c1 = regionprops(Ain1,'centroid'); c1 = c1.Centroid * Pixel_size; 
            end
            for n2 = 1:N
                % If Euclidean distances, calculate distances of
                % centroids
                if euclid_dist
                    Ain2 = reshape(Ain{p}(:,n2),d1,d2);
                    c2 = regionprops(Ain2,'centroid'); c2 = c2.Centroid * Pixel_size; 
                    distances(n1,n2) = sqrt(sum((c1-c2).^2));
                else
                    distances(n1,n2) = sqrt((XYproj_um(n1)-XYproj_um(n2))^2);
                end
            end
        end
        
        % Nearest neighbours
        NN_patch = nan(N,1);
        for n = 1:N
            NN_patch(n) = min(distances(n,[1:(n-1),(n+1):N]));
        end
        
        % Nearest neighbours ON
        NN_patch_ON = nan(N_ON,1);
        distances_ON = distances(ix_ON,ix_ON);
        if N_ON > 1
            for n = 1:N_ON
                NN_patch_ON(n) = min(distances_ON(n,[1:(n-1),(n+1):N_ON]));
            end
        end
        
        % Nearest neighbours OFF
        NN_patch_OFF = nan(N_OFF,1);
        distances_OFF = distances(ix_OFF,ix_OFF);
        if N_OFF > 1
            for n = 1:N_OFF
                NN_patch_OFF(n) = min(distances_OFF(n,[1:(n-1),(n+1):N_OFF]));
            end
        end
        
        % Randomly relable
        ix_ON = randsample(1:N,N_ON);
        ix_OFF = randsample(1:N,N_OFF);
        
        % Nearest neighbours ON
        NN_patch_ON_sh = nan(N_ON,1);
        distances_ON = distances(ix_ON,ix_ON);
        if N_ON >1
            for n = 1:N_ON
                NN_patch_ON_sh(n) = min(distances_ON(n,[1:(n-1),(n+1):N_ON]));
            end
        end
        
        % Nearest neighbours OFF
        NN_patch_OFF_sh = nan(N_OFF,1);
        distances_OFF = distances(ix_OFF,ix_OFF);
        if N_OFF >1
            for n = 1:N_OFF
                NN_patch_OFF_sh(n) = min(distances_OFF(n,[1:(n-1),(n+1):N_OFF]));
            end
        end

        % Save
        NN_dist = [NN_dist; NN_patch];
        NN_dist_ON = [NN_dist_ON; NN_patch_ON];
        NN_dist_OFF = [NN_dist_OFF; NN_patch_OFF];
        NN_dist_ON_shuff = [NN_dist_ON_shuff; NN_patch_ON_sh];
        NN_dist_OFF_shuff = [NN_dist_OFF_shuff; NN_patch_OFF_sh];

    end
        
    
    