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

function [rho,distances,rho_ON,distances_ON,rho_OFF,distances_OFF] = get_corr_vs_dist_3D(dataset_ix,grouped,euclid_dist)

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
    load([basedir,fname,'/',fname,'.mat'],'Pixel_size','Numb_patches','Patch_coordinates');
    
    if norm(diff(Patch_coordinates.data(:,[7,10])')) ~= 0
        error('Error with loading patch coordinates')
    end
    
    % vector orthogonal to vector_mean
    vector_orth = [-vector_mean(2), vector_mean(1)];
    if round(norm(vector_orth),10) ~=1 || dot(vector_mean,vector_orth)~=0
        error('Problem with vector_orth.')
    end

    % Patch size
    [d1,d2] = size(Cn{1});

    rho_all = [];
    distances_all = [];
    rho_ON_all = [];
    distances_ON_all = [];
    rho_OFF_all = [];
    distances_OFF_all = [];
    
    XYproj_um_all = cell(Numb_patches,1);
    Z_um_all = cell(Numb_patches,1);

    change_dFF = cell(Numb_patches,1);
    p_val = cell(Numb_patches,1);
    
    for p = 1:Numb_patches
        N = size(dFF{p},1);
        
        patch_X = Patch_coordinates.data(p,5);
        patch_Y = Patch_coordinates.data(p,6);
        
        % Find ON and OFF GCs in this patch
        [change_dFF{p},p_val{p}] = change_dFF_sig(dFF{p},A,QW,acquisition_rate);

        % Get vector of projected XY positions
        XYproj_um = zeros(1,N);
        for n1 = 1:N
            % Centroid of axon 1
            Ain1 = reshape(Ain{p}(:,n1),d1,d2);
            c1 = regionprops(Ain1,'centroid'); c1 = c1.Centroid; 
            
            XY_um = (c1 + [patch_X,patch_Y]) * Pixel_size;

            % Project axon 1 onto vector_orth and convert to um
            XYproj_um(n1) = dot(vector_orth,XY_um) ;
        end
        
        XYproj_um_all{p} = XYproj_um';
        Z_um_all{p} = ones(N,1) * Patch_coordinates.data(p,7);
    end

    XYproj_um_all = vertcat(XYproj_um_all{:});
    Z_um_all = vertcat(Z_um_all{:});
    
    dFF_all = vertcat(dFF{:});
    rho = corrcoef(dFF_all');
    
    N = size(dFF_all,1);
    
    distances = zeros(N,N);
    for n1 = 1:N
        for n2 = 1:N
            distances(n1,n2) = sqrt((XYproj_um_all(n1)-XYproj_um_all(n2))^2+(Z_um_all(n1)-Z_um_all(n2))^2);
        end
    end

    % Remove doublecounting
    change_dFF = vertcat(change_dFF{:});
    p_val = vertcat(p_val{:});

    ix_ON = find(change_dFF > 0 & p_val < .05);
    J_ON = zeros(N,N); J_ON(ix_ON,ix_ON) = 1;
    ix_OFF = find(change_dFF < 0 & p_val < .05);
    J_OFF = zeros(N,N); J_OFF(ix_OFF,ix_OFF) = 1;
        
    rho_ON = rho(triu(J_ON,1)==1);
    rho_OFF = rho(triu(J_OFF,1)==1);
    rho = rho(triu(ones(size(rho)),1)==1);

    distances_ON = distances(triu(J_ON,1)==1);
    distances_OFF = distances(triu(J_OFF,1)==1);
    distances = distances(triu(ones(size(distances)),1)==1);

    