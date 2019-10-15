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

function [rho_all,distances_all] = get_corr_vs_dist(dataset_ix,grouped,all_patches,euclid_dist)

    if nargin < 4 || isempty(euclid_dist)
        euclid_dist = 0;
    end
    
    if nargin < 3 || isempty(all_patches)
        all_patches = 0;
    end

    if nargin < 2 || isempty(grouped)
        grouped = 1;
    end

    define_dirs;

    fname = datasets{dataset_ix};
    disp(fname)
    
    load([basedir,fname,'/processed/',fname,'_GroupedData.mat'])
    
    % Calculate distances
    distances = [];
    if all_patches == 1
        % Calculate distances & correlations including in PFs in different
        % patches

        if grouped
            dFF_all = vertcat(dFF_axons{:});
        elseif ~grouped
            dFF_all = vertcat(dFF_rois{:});
        end

        N = size(dFF_all,1);
        
        % Load fibre direction, Z in pix, pix size, number of patches
        load([basedir,fname,'/processed/fibre_direction.mat'])
        load([basedir,fname,'/',fname,'.mat'],'Patch_coordinates','Pixel_size','Numb_patches');

        % Patch size
        [d1,d2] = size(Cn{1});
        
        % Concatenate Ain and get Z (um) for all rois
        if grouped
            Ain_all = horzcat(Ain_axons{:});
            Z_um = cell(Numb_patches,1);
            for p = 1:Numb_patches
                Z_um{p} = ones(1,size(Ain_axons{p},2)) * ...
                    Patch_coordinates.data(p,7) * Pixel_size;
            end
            Z_um = horzcat(Z_um{:});
        elseif ~grouped
            Ain_all = horzcat(Ain_rois{:});
            Z_um = cell(Numb_patches,1);
            for p = 1:Numb_patches
                Z_um{p} = ones(1,size(Ain_rois{p},2)) * ...
                    Patch_coordinates.data(p,7) * Pixel_size;
            end
            Z_um = horzcat(Z_um{:});
        end
        
        if euclid_dist
            error('Euclidean distances not yet supported for option to do all patches at once')
        end
        
        % Get vector of projected XY positions
        XYproj_um = zeros(1,N);
        X_um = zeros(1,N);
        Y_um = zeros(1,N);
        for n1 = 1:N
            % Centroid of axon 1
            Ain1 = reshape(Ain_all(:,n1),d1,d2);
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
        
        distances = zeros(N,N);
        for n1 = 1:N
            for n2 = 1:N                
                distances(n1,n2) = sqrt( (XYproj_um(n1)-XYproj_um(n2))^2 + (Z_um(n1)-Z_um(n2))^2 );
            end
        end
        
        rho = corrcoef(dFF_all');

        % Remove doublecounting
        rho_all = rho(triu(ones(size(rho)),1)==1);
        distances_all = distances(triu(ones(size(distances)),1)==1);
    else
        % Restrict only to the same patch
        
        if grouped
            dFF = dFF_axons;
            Ain = Ain_axons;
        elseif ~grouped
            dFF = dFF_rois;
            Ain = Ain_rois;
        end

        % Load fibre direction, Z in pix, pix size, number of patches
        load([basedir,fname,'/processed/fibre_direction.mat'])
        load([basedir,fname,'/',fname,'.mat'],'Pixel_size','Numb_patches');

        % Patch size
        [d1,d2] = size(Cn{1});

        rho_all = [];
        distances_all = [];

        for p = 1:Numb_patches
            N = size(dFF{p},1);

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

            rho = corrcoef(dFF{p}');

            % Remove doublecounting
            rho = rho(triu(ones(size(rho)),1)==1);
            distances = distances(triu(ones(size(distances)),1)==1);

            rho_all = [rho_all; rho];
            distances_all = [distances_all; distances];
        end
    end
        
    
    