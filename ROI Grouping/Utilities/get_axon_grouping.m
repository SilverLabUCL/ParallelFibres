%
% This function performs the main axon regrouping based on:
%   (1) low-frequency correlations, and 
%   (2) fibre direction
%
%
% Input:
%    Ain              Spatial filters matrix (num pixels x num ROIs)
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    dims             Vector of pixel dims of patch: [d1, d2]
%    acq_rate         Acquisition rate (Hz)
%    um_per_pix       Microns per pixel
%    fibre_vec        Vector of overall fibre direction (est. from z stack)
%    angle_std        Standard deviation of angle from overall fibre direction (default 10 degrees)
%    rho_min          Minimum correlation between ROIs
%    rho_std_max      Max standard deviation of correlation between ROIs
%    pix_defl_max     Maximum pixel deflection from existing fibre  
%                      (for grouping 3rd+ ROI onto already regrouped fibre) 
%    manual_check     Set to 1 if want to check all cases, 0 if completely automated, 
%                      e if want to check border cases within e% of rho_min and ang_max 
% 
% Output:
%    Ain_new            Spatial filters matrix after ROI grouping 
%                        (num pixels x num putative axons)
%    ix_axons_to_rois   Cell (1 x num putative axons) in which each
%                        element is the ROIs for that axon
%    axon_ids           Vector (1 x num ROIs) in which each element is the
%                        axon id for that ROI - note that axon_ids and
%                        ix_axons_to_rois are opposite transformations
% 

function [Ain_new,ix_axons_to_rois,axon_ids_new] = get_axon_grouping(Ain,Y,dims,acq_rate,um_per_pix,fibre_vec,angle_std,rho_min,rho_std_max,pix_defl_max,manual_check)
    %% Default values and load parameters
    
    if nargin < 11 || isempty(manual_check)
        manual_check = 0;
    end
    
    if nargin < 10 || isempty(pix_defl_max)
        pix_defl_max = 2;
    end
    
    if nargin < 9 || isempty(rho_std_max)
        rho_std_max = 0.25;
    end
    
    if nargin < 8 || isempty(rho_min)
        rho_min = 0.5;
    end
    
    % Default angle taken from distribution of firbres from z stacks
    if nargin < 7 || isempty(angle_std)
        angle_std = 13.45; 
    end

    % for estimated parameters
    T = size(Y,2);
    N_ROIs = size(Ain,2);
    
    d1 = dims(1);
    d2 = dims(2);

    %% First, calculate average dFF per ROI and correlations 
    
    
    dFF = get_dFF(Ain,Y,acq_rate);
    Cn = correlation_image(Y, [1,2], d1,d2);

    % Calculate correlations (initialisation)
    C = get_corrs(dFF,rho_min);
    
    % Group identities: loners vs notloners
    % Loners are ROIs that are not significantly correlated with anyone
    ix_loners = find(sum(C,1)==0 & sum(C,2)'==0 ); 
    
    % Notloners are ROIs that have at least min_corr correlation with
    % another ROI
    ix_notloners = setdiff(1:N_ROIs,ix_loners);
    
    % Index for axons, starting with loners
    axon_ids = zeros(1,N_ROIs);
    axon_ids(ix_loners) = 1:numel(ix_loners);
    axon_ids(ix_notloners) = (1+numel(ix_loners)):N_ROIs;

    %% Reorder rois so that loner axons are at the top
    temp = [ix_loners,ix_notloners];
    dFF_axons = dFF(temp,:);
    dFF_smooth_axons = dFF(temp,:);
    Ain_axons = Ain(:,temp);

    % indices for grouping
    ix_axons_to_rois = cell(1,N_ROIs);
    for axon = 1:N_ROIs
        % create a new axon for each loner roi, containing only that roi
        ix_axons_to_rois{axon} = temp(axon);
    end

    % Calculate correlations
    C = get_corrs(dFF_smooth_axons,rho_min);
    
    % Reference matrix to determine if have already checked if those rois are
    % on same axon. 1 indicates already checked, 0 otherwise. Value of 1
    % generally means that the correlation is high but the rois are too far
    % apart to be on the same axon.
    rois_have_checked = eye(N_ROIs,N_ROIs);
    rois_have_checked(ix_loners,:) = 1;
    rois_have_checked(:,ix_loners) = 1;

    %% Go through correlated ROIs for grouping

    if manual_check
        figure, 
    end
    
    % While there are correlated axons (above min correlation), try merging 
    while max(C(:)) > 0

        % Find pair with max correlation
        C_max = max(C(:));
        [axon_1,axon_2]=find(C==C_max);
        A_merged = Ain_axons(:,axon_1)+Ain_axons(:,axon_2);
        A_merged(A_merged>1) = 1; % in case there is slight overlap in spatial filters
        
        % Initialize variable
        merge_rois = 0;
        
        % Calculate spatial information and use this to determine whether
        % to merge ROIs - this relies on assumption 
        if numel(ix_axons_to_rois{axon_1}) == 1 && numel(ix_axons_to_rois{axon_2}) == 1
            % Case 1: 2 ROIs
            % First infer putative axon if we merge the ROIs
            [vector_axon,~,pts] = get_axon_vector(A_merged,[d1,d2]);
            % Calculate angle between mean vector and vector for this axon
            ang = rad2deg(subspace(vector_axon',fibre_vec'));
            % Group if angle is within 95% of angle distribution
            if ang < 2*angle_std
                merge_rois = 1;
            end
            if manual_check && merge_rois
                disp(' ')
                disp(['Angle between putative axon and overall axon: ',num2str(ang),' deg.'])
            end
        elseif numel(ix_axons_to_rois{axon_1}) == 1 && numel(ix_axons_to_rois{axon_2}) > 1
            % Case 2: 1 axon, 1 ROI (axon_2 has multiple ROIs)
            [vector_axon,x_int,pts] = get_axon_vector(Ain_axons(:,axon_2),[d1,d2]);
            vorth = [-vector_axon(2),vector_axon(1)];
            c = regionprops(reshape(Ain_axons(:,axon_1),d1,d2),'centroid'); 
            c = c.Centroid;
            dist = abs((c-[x_int,0])*vorth') * um_per_pix;
            % Consider grouping if roi is sufficiently close to putative fibre
            if dist < pix_defl_max
                merge_rois = 1;
            end
            if manual_check && merge_rois
                disp(' ')
                disp(['Distance between axon and ROI: ',num2str(dist),' um.'])
            end
        elseif numel(ix_axons_to_rois{axon_1}) > 1 && numel(ix_axons_to_rois{axon_2}) == 1
            % Same as previous case (now axon_1 has multiple ROIs)
            [vector_axon,x_int,pts] = get_axon_vector(Ain_axons(:,axon_1),[d1,d2]);
            vorth = [-vector_axon(2),vector_axon(1)];
            c = regionprops(reshape(Ain_axons(:,axon_2),d1,d2),'centroid'); 
            c = c.Centroid;
            dist = abs((c-[x_int,0])*vorth') * um_per_pix;
            % Consider grouping if roi is sufficiently close to putative fibre
            if dist < pix_defl_max
                merge_rois = 1;
            end
            if manual_check && merge_rois
                disp(' ')
                disp(['Distance between axon and ROI: ',num2str(dist),' um.'])
            end
        elseif numel(ix_axons_to_rois{axon_1}) > 1 && numel(ix_axons_to_rois{axon_2}) > 1
            % Case 3: 2 axons
            [vector_axon1,x_int1,pts] = get_axon_vector(Ain_axons(:,axon_1),[d1,d2]);
            [~,x_int2,~] = get_axon_vector(Ain_axons(:,axon_2),[d1,d2]);
            vorth1 = [-vector_axon1(2),vector_axon1(1)];
            dist = abs(([x_int2-x_int1,0])*vorth1') * um_per_pix;
            % Consider grouping if roi is sufficiently close to putative fibre
            if dist < pix_defl_max
                merge_rois = 1;
            end
            if manual_check && merge_rois
                disp(' ')
                disp(['Distance between axon and ROI: ',num2str(dist),' um.'])
            end
        else
            error('Cannot have empty axons')
        end

        % If merge_rois = 1 based on correlations and spatial location,
        % check that error for large events is lower than baseline error 
        if 1 %merge_rois
            
            % Get distribution of correlations
            C_hist = get_corr_dist(dFF_axons(axon_1,:),dFF_axons(axon_2,:),acq_rate);
            
            % If standard deviation of correlations is large, don't group
            % them anyway
            if nanstd(C_hist) >= rho_std_max
                merge_rois = 0;
            end
        
            % Chance to overide previous decision if using manual_check
            if manual_check
                
                % plot it if option chosen to manually check merging
                % should be changed in the future so that red (blue) trace has 
                % red (blue) spatial field contour
                subplot(2,5,1:4), plot_contours(A_merged,Cn,0.9,false,[],[],2); colormap(gray), hold on
                % plot axon 1 ROI centres in red
                for roi = ix_axons_to_rois{axon_1}
                    c = regionprops(reshape(Ain(:,roi),d1,d2),'centroid'); c = c.Centroid; 
                    plot(c(1),c(2),'*r')
                end
                % plot axon 2 ROI centres in red
                for roi = ix_axons_to_rois{axon_2}
                    c = regionprops(reshape(Ain(:,roi),d1,d2),'centroid'); c = c.Centroid; 
                    plot(c(1),c(2),'*b')
                end
                % plot vector connecting
                plot(pts(:,1),pts(:,2),':w')
                % plot traces
                subplot(2,5,6:9), hold off, plot((1:T)/acq_rate, dFF_axons(axon_1,:),'r'), 
                hold on, plot((1:T)/acq_rate, 1+dFF_axons(axon_2,:),'b'), xlim([0,T/acq_rate])

                rho = corrcoef(dFF_axons(axon_1,:),dFF_axons(axon_2,:)); 
                rho = rho(1,2);
                
                title(num2str(rho)), set(gca, 'FontSize',15), ylim([-.5,3])
                subplot(2,5,5), plot(dFF_axons(axon_1,:),dFF_axons(axon_2,:),'.k')
                hold on, plot([-.5,2],[-.5,2],'k'), hold off, 
                set(gca,'FontSize',15)
                axis tight, axis equal, axis([-.5,2,-.5,2])
                
                subplot(2,5,10), histogram(C_hist), xlim([-.4,1])
                title(['M: ', num2str(mean(C_hist)),' S: ', num2str(nanstd(C_hist))])
                set(gca, 'FontSize',15)

                % Ask for input
                if merge_rois == 1
                    x = input('ROIs will be merged. Enter R to override: ','s');
                    if  strcmp(x,'r') || strcmp(x,'R')
                        merge_rois = 0;
                    end
                elseif merge_rois == 0
                    x = input('ROIs will NOT be merged. Enter M to override: ','s');
                    if  strcmp(x,'r') || strcmp(x,'R')
                        merge_rois = 1;
                    end
                end
            end
        end

        % If still decide to merge ROIs
        if merge_rois

            % Merge spatial fields for the two axons
            Ain_axons_new = Ain_axons;
            Ain_axons_new(:,axon_1) = A_merged;
            Ain_axons_new(:,axon_2) = []; 

            % Update axon identities
            ix_axons_to_rois_new = ix_axons_to_rois;
            ix_axons_to_rois_new{axon_1} = [ix_axons_to_rois{axon_1},ix_axons_to_rois{axon_2}];
            ix_axons_to_rois_new(axon_2) = []; % delete that cell element
            num_axons = length(ix_axons_to_rois_new);

            % update axon ids for the rois
            axon_ids_new = zeros(1,N_ROIs);
            for axon = 1:num_axons
                rois = ix_axons_to_rois_new{axon};
                axon_ids_new(rois) = axon;
            end

            % calculate dFF for new axons            
            dFF_axons_new = dFF_axons;
            dFF_axons_new(axon_1,:) = get_dFF(A_merged,Y,acq_rate);
            dFF_axons_new(axon_2,:) = [];
            
            dFF_smooth_axons_new = dFF_axons_new;

            % recompute correlations
            C = get_corrs(dFF_smooth_axons_new,rho_min);

            % update everything
            Ain_axons = Ain_axons_new;
            ix_axons_to_rois = ix_axons_to_rois_new;
            axon_ids = axon_ids_new;
            dFF_axons = dFF_axons_new;
            dFF_smooth_axons= dFF_smooth_axons_new;
        else
            % Merge has been rejected despite high correlations. This means the
            % rois or axons are too far apart. Flag that this has arleady been
            % checked, or else the program gets stuck
            for roi_1 = ix_axons_to_rois{axon_1}
                for roi_2 = ix_axons_to_rois{axon_2}
                    rois_have_checked(roi_1,roi_2) = 1;
                    rois_have_checked(roi_2,roi_1) = 1;
                end
            end
        end

        % set correlations to -inf if the rois have previously been determined
        % to be too far to be on the same axon
        [axon_r,axon_c] = find(C>0);
        for a = 1:length(axon_r)
            for roi_1 = ix_axons_to_rois{axon_r(a)}
                for roi_2 = ix_axons_to_rois{axon_c(a)}
                    if rois_have_checked(roi_1,roi_2) > 0
                        C(axon_r(a),axon_c(a)) = -inf;
                    end
                end
            end
        end
       
    end

    %% Now reorder axons so that all the loners are 1st

    num_axons = length(ix_axons_to_rois);
    ix_axons_to_rois_new = cell(1,num_axons);
    ix_loners = cell(1); i = 1; j = 1;
    for axon = 1:num_axons
        if numel(ix_axons_to_rois{axon}) == 1
            ix_loners{i} = ix_axons_to_rois{axon};
            i = i+1;
        else
            ix_axons_to_rois_new{j}  = ix_axons_to_rois{axon};
            j = j+1;
        end
    end
    ix_axons_to_rois_new(j:end) = ix_loners;

    axon_ids_new = zeros(1,N_ROIs);
    for axon = 1:num_axons
        rois = ix_axons_to_rois_new{axon};
        axon_ids_new(rois) = axon;
    end

    Ain_new = zeros(size(Ain_axons));
    for axon = 1:num_axons
        for roi = ix_axons_to_rois_new{axon}    
            Ain_new(:,axon) = Ain_new(:,axon) + Ain(:,roi);
        end
    end

ix_axons_to_rois = ix_axons_to_rois_new;
axon_ids = axon_ids_new;


