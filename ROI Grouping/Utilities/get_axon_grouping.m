%
% This function performs the main axon grouping procedure
%
% Input:
%    Ain              Spatial filters matrix (num pixels x num ROIs)
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    dims             Vector of pixel dims of patch: [d1, d2]
%    acq_rate         Acquisition rate (Hz)
%    smooth_win_s     Temporal window for smoothing data
%    um_per_pix       Microns per pixel
%    fibre_vec        Vector of overall fibre direction (est. from z stack)
%    angle_std        Standard deviation of angle from overall fibre direction (default 10 degrees)
%    rho_min          Minimum correlation between ROIs
%    var_ratio_max      Max standard deviation of correlation between ROIs
%    pix_defl_max     Maximum pixel deflection from existing fibre  
%                      (for grouping 3rd+ ROI onto already regrouped fibre) 
%    manual_check     Set to 1 if want to check all cases, 0 if completely automated, 
%                      e if want to check border cases within e% of rho_min and ang_max 
% 
% Output:
%    Ain_axons           Spatial filters matrix after ROI grouping 
%                        (num pixels x num putative axons)
%    ix_axons_to_rois   Cell (1 x num putative axons) in which each
%                        element is the ROIs for that axon
%    axon_ids           Vector (1 x num ROIs) in which each element is the
%                        axon id for that ROI - note that axon_ids and
%                        ix_axons_to_rois are opposite transformations
% 

function [Ain_axons,ix_axons_to_rois,axon_ids] = get_axon_grouping(Ain,Y,dims,acq_rate,smooth_win_s,um_per_pix,fibre_vec,angle_std,rho_min,var_ratio_max,pix_defl_max,manual_check)
    %% Default values and load parameters
    
    close all
    
    if nargin < 12 || isempty(manual_check)
        manual_check = 0;
    end
    
    if nargin < 11 || isempty(pix_defl_max)
        pix_defl_max = 1;
    end
    
    if nargin < 10 || isempty(var_ratio_max)
        var_ratio_max = 1.5;
    end
    
    if nargin < 9 || isempty(rho_min)
        rho_min = 0.5;
    end
    
    % Default angle taken from distribution of firbres from z stacks
    if nargin < 8 || isempty(angle_std)
        angle_std = 13.45; 
    end

    % for estimated parameters
    T = size(Y,2);
    N_ROIs = size(Ain,2);
    
    d1 = dims(1);
    d2 = dims(2);

    %% First, calculate average dFF per ROI and correlations 
    
    dFF = get_dFF(Ain,Y,acq_rate,smooth_win_s);
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
    %axon_ids = zeros(1,N_ROIs);
    %axon_ids(ix_loners) = 1:numel(ix_loners);
    %axon_ids(ix_notloners) = (1+numel(ix_loners)):N_ROIs;

    %% Reorder rois so that loner axons are at the top
    temp = [ix_loners,ix_notloners];
    dFF_axons = dFF(temp,:);
    Ain_axons = Ain(:,temp);

    % indices for grouping
    ix_axons_to_rois = cell(1,N_ROIs);
    for axon = 1:N_ROIs
        % create a new axon for each loner roi, containing only that roi
        ix_axons_to_rois{axon} = temp(axon);
    end

    % Calculate correlations
    C = get_corrs(dFF_axons,rho_min);
    
    % Reference matrix to determine if have already checked if those rois are
    % on same axon. 1 indicates already checked, 0 otherwise. Value of 1
    % generally means that the correlation is high but the rois are too far
    % apart to be on the same axon.
    rois_have_checked = eye(N_ROIs,N_ROIs);
    rois_have_checked(ix_loners,:) = 1;
    rois_have_checked(:,ix_loners) = 1;

    %% Go through correlated ROIs for grouping

    if manual_check
        figure(1), 
    end
    
    % Count how many merges (count_merge), and how many automatic merges 
    % are manually prevented (count_FA), and how many automatic non-merges
    % are manually forced (count_MISS)
    count_merge = 0;
    count_FA = 0;
    count_MISS = 0;
    
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
            
        % First infer putative axon if we merge the ROIs
        [vector_axon,x_int,pts] = get_axon_vector(A_merged,[d1,d2]);
        
        % Calculate angle between mean vector and vector for this axon
        ang = rad2deg(subspace(vector_axon',fibre_vec'));
        
        % Vector to project onto to get distances of boutons
        vorth = [-vector_axon(2),vector_axon(1)];
        
        % Get distances of boutons to axon
        dist = zeros(numel(ix_axons_to_rois{axon_1}) + numel(ix_axons_to_rois{axon_2}),1);
        ix = 1;
        for roi = [ix_axons_to_rois{axon_1}, ix_axons_to_rois{axon_2}]
            c = regionprops(reshape(Ain(:,roi),d1,d2),'centroid'); c = c.Centroid; 
            dist(ix) = abs((c-[x_int,0])*vorth') * um_per_pix;
            ix = ix+1;
        end
        
        % Group if angle is within 95% of angle distribution
        if ang < 2*angle_std && sum(dist < pix_defl_max) == numel(dist)
            merge_rois = 1;
        end
        
        if manual_check
            disp(' ')
            disp(['Angle between putative axon and overall axon: ',num2str(ang),' deg.'])
            for ix = 1:numel(dist)
                disp(['Distance of ROI from axon: ',num2str(dist(ix)),' um.'])
            end
        end

        % If merge_rois = 1 based on correlations and spatial location,
        % check that error for large events is lower than baseline error 
        if merge_rois || manual_check
            
            var_ratio = get_var_ratio(dFF_axons(axon_1,:),dFF_axons(axon_2,:),manual_check);
            
            % If var_ratio is too high, don't group them anyway
            if var_ratio >= var_ratio_max
                merge_rois = 0;
            end
        
            % Chance to overide previous decision if using manual_check
            if manual_check
                
                % plot it if option chosen to manually check merging
                % should be changed in the future so that red (blue) trace has 
                % red (blue) spatial field contour
                figure(1), subplot(2,1,1), hold off
                plot_contours(A_merged,Cn,0.9,false,[],[],2); 
                colormap(gray), hold on
                
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
                subplot(2,1,2), hold off, plot((1:T)/acq_rate, dFF_axons(axon_1,:),'r'), 
                hold on, plot((1:T)/acq_rate, mean(dFF_axons(axon_1,:)) * 3 + .5 + dFF_axons(axon_2,:),'b'), xlim([0,T/acq_rate])

                rho = corrcoef(dFF_axons(axon_1,:),dFF_axons(axon_2,:)); 
                rho = rho(1,2);
                
                title(num2str(rho)), set(gca, 'FontSize',15), temp = mean(dFF_axons(axon_1,:))*3 + mean(dFF_axons(axon_2,:)) * 3 +3; ylim([-.5,temp])%ylim([-.5,mean(dFF_axons(axon_1,:))*3 + mean(dFF_axons(axon_2,:)) * 3 +3 ])
                
                % Ask for input
                if merge_rois == 1
                    x = input('ROIs will be merged. Enter R to override: ','s');
                    if  strcmp(x,'r') || strcmp(x,'R')
                        merge_rois = 0;
                        count_FA = count_FA + 1;
                    end
                elseif merge_rois == 0
                    set(gca,'Color',[1,.8,.8])
                    title(num2str(rho),'Color','r')
                    x = input('ROIs will NOT be merged. Enter M to override: ','s');
                    if  strcmp(x,'m') || strcmp(x,'M')
                        merge_rois = 1;
                        count_MISS = count_MISS + 1;
                    end
                end
            end
        end

        % If still decide to merge ROIs
        if merge_rois
            
            % Update count
            count_merge = count_merge + 1;
            
            % Merge spatial fields for the two axons
            Ain_axons_new = Ain_axons;
            Ain_axons_new(:,axon_1) = A_merged;
            Ain_axons_new(:,axon_2) = []; 

            % Update axon identities
            ix_axons_to_rois_new = ix_axons_to_rois;
            ix_axons_to_rois_new{axon_1} = [ix_axons_to_rois{axon_1},ix_axons_to_rois{axon_2}];
            ix_axons_to_rois_new(axon_2) = []; % delete that cell element
            num_axons = length(ix_axons_to_rois_new);

            % calculate dFF for new axons            
            dFF_axons_new = dFF_axons;
            dFF_axons_new(axon_1,:) = get_dFF(A_merged,Y,acq_rate,smooth_win_s);
            dFF_axons_new(axon_2,:) = [];
            
            % recompute correlations
            C = get_corrs(dFF_axons_new,rho_min);

            % update everything
            Ain_axons = Ain_axons_new;
            ix_axons_to_rois = ix_axons_to_rois_new;
            %axon_ids = axon_ids_new;
            dFF_axons = dFF_axons_new;
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

    %% Now reorder axons so that all the loners are last
    % After reordering have to recompute axon_ids, Ain_Axons

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

    Ain_axons_new = zeros(size(Ain_axons));
    dFF_axons_new = zeros(size(dFF_axons));
    for axon = 1:num_axons
        for roi = ix_axons_to_rois_new{axon}    
            Ain_axons_new(:,axon) = Ain_axons_new(:,axon) + Ain(:,roi);
        end
        dFF_axons_new(axon,:) = mean(dFF(ix_axons_to_rois_new{axon},:),1);
    end

    ix_axons_to_rois = ix_axons_to_rois_new;
    axon_ids = axon_ids_new;
    Ain_axons = Ain_axons_new;
    dFF_axons = dFF_axons_new;

    disp(['Number of merges:',num2str(count_merge)])
    disp(['Number of false alarms:',num2str(count_FA)])
    disp(['Number of misses:',num2str(count_MISS)])

    %% Next ask for any possible manual merges
    
    if manual_check
        plot_grouped_rois(Ain_axons,Cn,dFF,ix_axons_to_rois,acq_rate);

        % Loop for adding manually
        in = [];
        while ~(strcmp(in,'A') || strcmp(in,'a') || strcmp(in,'M') || strcmp(in,'m'))
            in = input('Enter A to accept grouping or M to merge a pair by hand: ','s');
        end    
        stop = 0;
        while strcmp(in,'m') || strcmp(in,'M')
            stop = 0;
            while stop == 0
                axon_1 = str2double(input('Enter Axon 1 ID#: ','s'));
                axon_2 = str2double(input('Enter Axon 2 ID#: ','s'));
                if numel(axon_1) == 1 && numel(axon_2) ==1
                    if (axon_1 <= num_axons && axon_2 <= num_axons && axon_1 > 0 && axon_2 > 0 )
                        stop = 1;
                    end
                end
            end
            
            close all;
            
            A_merged =  Ain_axons(:,axon_1) + Ain_axons(:,axon_2);
            
            % First infer putative axon if we merge the ROIs
            [vector_axon,x_int,pts] = get_axon_vector(A_merged,[d1,d2]);
        
             % Calculate angle between mean vector and vector for this axon
            ang = rad2deg(subspace(vector_axon',fibre_vec'));

            % Vector to project onto to get distances of boutons
            vorth = [-vector_axon(2),vector_axon(1)];

            % Get distances of boutons to axon
            dist = zeros(numel(ix_axons_to_rois{axon_1}) + numel(ix_axons_to_rois{axon_2}),1);
            ix = 1;
            for roi = [ix_axons_to_rois{axon_1}, ix_axons_to_rois{axon_2}]
                c = regionprops(reshape(Ain(:,roi),d1,d2),'centroid'); c = c.Centroid; 
                dist(ix) = abs((c-[x_int,0])*vorth') * um_per_pix;
                ix = ix+1;
            end
            disp(' ')
            disp(['Angle between putative axon and overall axon: ',num2str(ang),' deg.'])
            for ix = 1:numel(dist)
                disp(['Distance of ROI from axon: ',num2str(dist(ix)),' um.'])
            end
            
            
            get_var_ratio(dFF_axons(axon_1,:),dFF_axons(axon_2,:),manual_check);
            
            
            % plot it if option chosen to manually check merging
                % should be changed in the future so that red (blue) trace has 
                % red (blue) spatial field contour
                figure(1), subplot(2,1,1), hold off
                plot_contours(A_merged,Cn,0.9,false,[],[],2); 
                colormap(gray), hold on
                
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
                subplot(2,1,2), hold off, plot((1:T)/acq_rate, dFF_axons(axon_1,:),'r'), 
                hold on, plot((1:T)/acq_rate, mean(dFF_axons(axon_1,:)) * 3 + .5 + dFF_axons(axon_2,:),'b'), xlim([0,T/acq_rate])

                rho = corrcoef(dFF_axons(axon_1,:),dFF_axons(axon_2,:)); 
                rho = rho(1,2);
                
                title(num2str(rho)), set(gca, 'FontSize',15), temp = mean(dFF_axons(axon_1,:))*3 + mean(dFF_axons(axon_2,:)) * 3 +3; ylim([-.5,temp])%ylim([-.5,mean(dFF_axons(axon_1,:))*3 + mean(dFF_axons(axon_2,:)) * 3 +3 ])
                
            % Ask for input
            x = input('ROIs will be merged. Enter R to override: ','s');
            if ~ ( strcmp(x,'r') || strcmp(x,'R') )
            
                % Merge spatial fields for the two axons
                Ain_axons_new = Ain_axons;
                Ain_axons_new(:,num_axons+1) = A_merged;
                Ain_axons_new(:,[axon_1,axon_2]) = [];
                Ain_axons = Ain_axons_new;
                
                % Update dFF
                dFF_axons_new = dFF_axons;
                dFF_axons_new(num_axons+1,:) = get_dFF(A_merged,Y,acq_rate,smooth_win_s);
                dFF_axons_new([axon_1,axon_2],:) = [];
                dFF_axons = dFF_axons_new;
                
                % Update axon identities
                ix_axons_to_rois_new = ix_axons_to_rois;
                ix_axons_to_rois_new{num_axons+1} = [ix_axons_to_rois{axon_1},ix_axons_to_rois{axon_2}];
                ix_axons_to_rois_new([axon_1,axon_2]) = []; % delete that cell element
                ix_axons_to_rois = ix_axons_to_rois_new;
                
                num_axons = length(ix_axons_to_rois);

                % Update axon_ids
                axon_ids_new = axon_ids;
                for axon = 1:num_axons
                    rois = ix_axons_to_rois{axon};
                    axon_ids_new(rois) = axon;
                end
                axon_ids = axon_ids_new;

                count_merge = count_merge+1;
                count_MISS = count_MISS+1;
            end

            % Plot new
            plot_grouped_rois(Ain_axons,Cn,dFF,ix_axons_to_rois,acq_rate);

            in = [];
            while ~(strcmp(in,'A') || strcmp(in,'a') || strcmp(in,'M') || strcmp(in,'m'))
                in = input('Enter A to accept grouping or M to merge a pair by hand: ','s');
            end
        end
        
        disp(['Number of merges:',num2str(count_merge)])
        disp(['Number of false alarms:',num2str(count_FA)])
        disp(['Number of misses:',num2str(count_MISS)])
    end
    
