
clear all
    define_dirs;

    dataset_ix = 15
   
    grouped = 0;
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

        %%
for p = 1:Numb_patches
        N = size(dFF{p},1);
        p
        
        % Find ON and OFF GCs in this patch
        [change_dFF,p_val] = change_dFF_sig(dFF{p},A,QW,acquisition_rate);
        ix_ON = find(change_dFF > 0 & p_val < .05);
        J_ON = zeros(N,N); J_ON(ix_ON,ix_ON) = 1;
        ix_OFF = find(change_dFF < 0 & p_val < .05);
        J_OFF = zeros(N,N); J_OFF(ix_OFF,ix_OFF) = 1;
        
        
        rho = corrcoef(dFF{p}');
        rho = rho - diag(diag(rho));
        [a,b] = find(rho > .8);
        if ~isempty(intersect(a,ix_OFF))
            ix_QW = [];
            for k = 1:length(QW)
                ix_QW = [ix_QW, QW(k,1):QW(k,2)];
            end
            temp=intersect(a,ix_OFF);
            % [corr(dFF{p}(j,:)',dFF{p}(k,:)'),corr(dFF{p}(j,ix_QW)',dFF{p}(k,ix_QW)')]
            % figure, plot(time_rois{p}(j,:)/1000,dFF{p}(j,:)), hold on, plot(time_rois{p}(k,:)/1000,1+dFF{p}(k,:))
            figure(1), subplot(1,2,1)
            imagesc(corrcoef(dFF{p}(temp,:)')), caxis([.5,1])
            subplot(1,2,2)
            imagesc(corrcoef(dFF{p}(temp,ix_QW)')), caxis([.5,1])
            disp('Press enter.')
            pause
            close(1)

        end
end
        disp('done')
    %%

        % If projecting distances, project centroids now
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

        distances = zeros(N,N);
        for n1 = 1:N
            for n2 = 1:N
                distances(n1,n2) = sqrt((XYproj_um(n1)-XYproj_um(n2))^2);
            end
        end

    