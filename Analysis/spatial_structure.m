% This script goes through all datasets and calculates the correlations
% between ROIs

clear all; clc

define_dirs;
%%
is_3D = 0;
grouped = 0;

rho_all = []; distances_all = []; 
rho_ON_all = []; distances_ON_all = []; 
rho_OFF_all = []; distances_OFF_all = [];
for dataset_ix = 1:13

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,~,loco,speed] = load_behav_data(dataset_ix,time);
    onset_indices = get_onsets(speed,acquisition_rate);
    
    % Choose between 2d (patch based) and 3d (all patches) versions
    if ~is_3D
        [rho,distances,rho_ON,distances_ON,rho_OFF,distances_OFF] = get_corr_vs_dist(dataset_ix,grouped);
    else
        [rho,distances,rho_ON,distances_ON,rho_OFF,distances_OFF] = get_corr_vs_dist_3D(dataset_ix);
        disp('Calculating in XYZ. Warning: always grouped.')
    end
    
    rho_all = [rho_all; rho];
    distances_all = [distances_all; distances];
    
    rho_ON_all = [rho_ON_all; rho_ON];
    distances_ON_all = [distances_ON_all; distances_ON];

    rho_OFF_all = [rho_OFF_all; rho_OFF];
    distances_OFF_all = [distances_OFF_all; distances_OFF];

end

clear rho distances

ranksum(rho_all(distances_all<2),rho_all(distances_all>=2))

if ~is_3D
    if grouped
        save([basedir,'processed/spatial_corr_grouped'],'rho_all','distances_all',...
        'rho_ON_all','distances_ON_all','rho_OFF_all','distances_OFF_all')
    else
        save([basedir,'processed/spatial_corr_ungrouped'],'rho_all','distances_all',...
        'rho_ON_all','distances_ON_all','rho_OFF_all','distances_OFF_all')
    end
else
    save([basedir,'processed/spatial_corr_3D_grouped'],'rho_all','distances_all',...
    'rho_ON_all','distances_ON_all','rho_OFF_all','distances_OFF_all')
end

%% Testing onset correlations

rho_all = []; distances_all = []; 
rho_ON_all = []; distances_ON_all = [];
rho_OFF_all = []; distances_OFF_all = [];
for dataset_ix = 1:13

    % Get A and QW states
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,~,loco,speed] = load_behav_data(dataset_ix,time);
    onset_indices = get_onsets(speed,acquisition_rate);
    
    time_ix = [];
    for k = 1:size(onset_indices,1)
       time_ix = [time_ix, onset_indices(k,1):onset_indices(k,2)];
    end
    
    if isempty(time_ix)
        [rho,distances,rho_ON,distances_ON,rho_OFF,distances_OFF] = deal([]);
    else
        [rho,distances,rho_ON,distances_ON,rho_OFF,distances_OFF] = get_corr_vs_dist(dataset_ix,1,time_ix);
    end
    
    rho_all = [rho_all; rho];
    distances_all = [distances_all; distances];
    
    rho_ON_all = [rho_ON_all; rho_ON];
    distances_ON_all = [distances_ON_all; distances_ON];

    rho_OFF_all = [rho_OFF_all; rho_OFF];
    distances_OFF_all = [distances_OFF_all; distances_OFF];
        
end

clear rho distances 

save([basedir,'processed/spatial_corr_onsets'],'rho_all','distances_all',...
    'rho_ON_all','distances_ON_all','rho_OFF_all','distances_OFF_all')

%% Shuffle case for onsets
num_its = 20;

dbin = 3;
dist_bins_e = 0:dbin:60;
dist_bins_c = dist_bins_e(2:end) - dbin/2;

rho_binned_sh = cell(size(dist_bins_c)); 
rho_ON_binned_sh = cell(size(dist_bins_c));
rho_OFF_binned_sh = cell(size(dist_bins_c));

for dataset_ix = 1:13

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,~,loco,speed] = load_behav_data(dataset_ix,time);
    onset_indices = get_onsets(speed,acquisition_rate);
    
    time_ix = [];
    for k = 1:size(onset_indices,1)
       time_ix = [time_ix, onset_indices(k,1):onset_indices(k,2)];
    end
    
    if ~isempty(time_ix)
        for it_ix = 1:num_its
            % Shuffle time_ixs
            T = size(dFF,2);
            for shuffle_ix = 1:10
                time_ix_sh = time_ix + randsample(T,1);
                time_ix_sh(time_ix_sh > T) = time_ix_sh(time_ix_sh > T) - T;
            end
            
            [rho_sh,distances_sh,rho_ON_sh,distances_ON_sh,rho_OFF_sh,distances_OFF_sh] = get_corr_vs_dist(dataset_ix,1,time_ix_sh);
            
            for k = 1:length(dist_bins_c)
                ix = distances_sh <= dist_bins_e(k+1) & distances_sh > dist_bins_e(k);
                rho_binned_sh{k} = [rho_binned_sh{k}; rho_sh(ix)];
                
                ix = distances_ON_sh <= dist_bins_e(k+1) & distances_ON_sh > dist_bins_e(k);
                rho_ON_binned_sh{k} = [rho_ON_binned_sh{k}; rho_ON_sh(ix)];
                
                ix = distances_OFF_sh <= dist_bins_e(k+1) & distances_OFF_sh > dist_bins_e(k);
                rho_OFF_binned_sh{k} = [rho_OFF_binned_sh{k}; rho_OFF_sh(ix)];
                
            end
            
        end
    end
        
end

clear rho distances rho_sh distances_sh

save([basedir,'processed/spatial_corr_shuffle'],'rho_binned_sh','rho_ON_binned_sh','rho_OFF_binned_sh','dist_bins_c')

%% Compare to random times in AS
num_its = 20;

dbin = 3;
dist_bins_e = 0:dbin:60;
dist_bins_c = dist_bins_e(2:end) - dbin/2;

rho_binned_AS = cell(size(dist_bins_c)); 
rho_ON_binned_AS = cell(size(dist_bins_c));
rho_OFF_binned_AS = cell(size(dist_bins_c));

for dataset_ix = 1:13

    % Get A and QW states
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,loco,speed] = load_behav_data(dataset_ix,time);
    onset_indices = get_onsets(speed,acquisition_rate);
    
    [A,~] = define_behav_periods(whisk_amp,loco,acquisition_rate);
    ix_A = [];
    for k = 1:length(A)
        ix_A = [ix_A,A(k,1):A(k,2)];
    end
    
    time_ix = [];
    for k = 1:size(onset_indices,1)
       time_ix = [time_ix, onset_indices(k,1):onset_indices(k,2)];
    end
    
    if ~isempty(time_ix)
        for it_ix = 1:num_its
            % Take ranomd
            time_ix_AS = randsample(ix_A,numel(time_ix));
            
            [rho_AS,distances_AS,rho_ON_AS,distances_ON_AS,rho_OFF_AS,distances_OFF_AS] = get_corr_vs_dist(dataset_ix,1,time_ix_AS);
            
            for k = 1:length(dist_bins_c)
                ix = distances_AS <= dist_bins_e(k+1) & distances_AS > dist_bins_e(k);
                rho_binned_AS{k} = [rho_binned_AS{k}; rho_AS(ix)];
                
                ix = distances_ON_AS <= dist_bins_e(k+1) & distances_ON_AS > dist_bins_e(k);
                rho_ON_binned_AS{k} = [rho_ON_binned_AS{k}; rho_ON_AS(ix)];
                
                ix = distances_OFF_AS <= dist_bins_e(k+1) & distances_OFF_AS > dist_bins_e(k);
                rho_OFF_binned_AS{k} = [rho_OFF_binned_AS{k}; rho_OFF_AS(ix)];
                
            end
            
        end
    end
        
end

clear rho distances rho_AS distances_AS

save([basedir,'processed/spatial_corr_AS_random'],'rho_binned_AS','rho_ON_binned_AS','rho_OFF_binned_AS','dist_bins_c')

%% Plot correlation vs distance
% Figure 1F, S6

% spatial_corr_3D_grouped
% spatial_corr_grouped
% spatial_corr_ungrouped

load([basedir,'processed/spatial_corr_onsets'])
figure, hold on
title('Axons XY Onsets')

dbin = 3;
dist_bins_e = 0:dbin:60;
dist_bins_c = dist_bins_e(2:end) - dbin/2;

for k = 1:length(dist_bins_c)
    ix = find(distances_all <= dist_bins_e(k+1) & distances_all > dist_bins_e(k));
    ste = std(rho_all(ix))/sqrt(numel(ix));
    plot(dist_bins_c(k),mean(rho_all(ix)),'sk');
    plot(dist_bins_c(k) * [1,1],mean(rho_all(ix)) + ste*[-1,1],'k');
    
    ix = find(distances_ON_all <= dist_bins_e(k+1) & distances_ON_all > dist_bins_e(k));
    ste = std(rho_ON_all(ix))/sqrt(numel(ix));
    plot(dist_bins_c(k),mean(rho_ON_all(ix)),'sm');
    plot(dist_bins_c(k) * [1,1],mean(rho_ON_all(ix)) + ste*[-1,1],'m');
    
    ix = find(distances_OFF_all <= dist_bins_e(k+1) & distances_OFF_all > dist_bins_e(k));
    ste = std(rho_OFF_all(ix))/sqrt(numel(ix));
    plot(dist_bins_c(k),mean(rho_OFF_all(ix)),'sc');
    plot(dist_bins_c(k) * [1,1],mean(rho_OFF_all(ix)) + ste*[-1,1],'c');
end

options = optimset('MaxFunEvals',500000);
expfun = @(params,x) params(1) + params(2) .* exp(-params(3) .* x)+ params(4) .* exp(-params(5) .* x); 
err = @(params) sum( ( expfun(params,distances_all) - rho_all ).^2 );
params_fit = fminsearch(err,[.1,.5,.5,.5,.1],options);
x = 1:.01:59;
hold on, plot(x,expfun(params_fit,x),'k','LineWidth',3)

err = @(params) sum( ( expfun(params,distances_ON_all) - rho_ON_all ).^2 );
params_fit = fminsearch(err,[.1,.5,.5,.5,.1],options);
hold on, plot(x,expfun(params_fit,x),'m','LineWidth',3)

err = @(params) sum( ( expfun(params,distances_OFF_all) - rho_OFF_all ).^2 );
params_fit = fminsearch(err,[.05,.5,.5,.5,.1],options);
hold on, plot(x,expfun(params_fit,x),'c','LineWidth',3)

set(gca,'FontSize',18)
xlabel('Distance between fibres (um)')
ylabel('Pairwise correlation')
set(gca,'Box','off')
ylim([.0,.4])

%% For shuffle

load([basedir,'processed/spatial_corr_shuffle'])

[y_mean,y_top,y_bot,y_ON_mean,y_ON_top,y_ON_bot,y_OFF_mean,y_OFF_top,y_OFF_bot] = deal(zeros(size(dist_bins_c)));

for k = 1:length(dist_bins_c)
    ste = std(rho_binned_sh{k})/sqrt(numel(rho_binned_sh{k}));
    y_mean(k) = mean(rho_binned_sh{k});
    y_top(k) = y_mean(k) + ste;
    y_bot(k) = y_mean(k) - ste;
    
    ste = std(rho_ON_binned_sh{k})/sqrt(numel(rho_ON_binned_sh{k}));
    y_ON_mean(k) = mean(rho_ON_binned_sh{k});
    y_ON_top(k) = y_ON_mean(k) + ste;
    y_ON_bot(k) = y_ON_mean(k) - ste;
    
    ste = std(rho_OFF_binned_sh{k})/sqrt(numel(rho_OFF_binned_sh{k}));
    y_OFF_mean(k) = nanmean(rho_OFF_binned_sh{k});
    y_OFF_top(k) = y_OFF_mean(k) + ste;
    y_OFF_bot(k) = y_OFF_mean(k) - ste;
end


hold on, plot(dist_bins_c,y_mean,'k');
hold on, plot(dist_bins_c,y_top,'k');
hold on, plot(dist_bins_c,y_bot,'k');

hold on, plot(dist_bins_c,y_ON_mean,'m');
hold on, plot(dist_bins_c,y_ON_top,'m');
hold on, plot(dist_bins_c,y_ON_bot,'m');

hold on, plot(dist_bins_c,y_OFF_mean,'c');
hold on, plot(dist_bins_c,y_OFF_top,'c');
hold on, plot(dist_bins_c,y_OFF_bot,'c');


%% Compare correlation during onsets to shuffle 

num_its = 20;

rho_all = cell(13,1);
rho_ON_all = cell(13,1);
rho_OFF_all = cell(13,1);

rho_all_sh = cell(13,1);
rho_ON_all_sh = cell(13,1);
rho_OFF_all_sh = cell(13,1);

rho_all_AS = cell(13,1);
rho_ON_all_AS = cell(13,1);
rho_OFF_all_AS = cell(13,1);

for dataset_ix = 1:13

    % Get A and QW states
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,loco,speed] = load_behav_data(dataset_ix,time);
    onset_indices = get_onsets(speed,acquisition_rate);
    
    N = size(dFF,1);
      
    % Find ON and OFF GCs in this patch
    [A,QW] = define_behav_periods(whisk_amp,loco,acquisition_rate,0);
    [change_dFF,p_val] = change_dFF_sig(dFF,A,QW,acquisition_rate);
    ix_ON = find(change_dFF > 0 & p_val < .05);
    J_ON = zeros(N,N); J_ON(ix_ON,ix_ON) = 1;
    ix_OFF = find(change_dFF < 0 & p_val < .05);
    J_OFF = zeros(N,N); J_OFF(ix_OFF,ix_OFF) = 1;
    
    ix_A = [];
    for k = 1:length(A)
        ix_A = [ix_A,A(k,1):A(k,2)];
    end

    % Time indices
    time_ix = [];
    for k = 1:size(onset_indices,1)
       time_ix = [time_ix, onset_indices(k,1):onset_indices(k,2)];
    end
    
    if ~isempty(time_ix)
        rho = corrcoef(dFF(:,time_ix)');

        % Remove doublecounting
        rho_ON = rho(triu(J_ON,1)==1);
        rho_OFF = rho(triu(J_OFF,1)==1);
        rho = rho(triu(ones(size(rho)),1)==1);

        rho_all{dataset_ix} = rho;
        rho_ON_all{dataset_ix} = rho_ON;
        rho_OFF_all{dataset_ix} = rho_OFF; 
        
        % Remove doublecounting
        [rho_ON_sh,rho_OFF_sh,rho_sh,rho_ON_AS,rho_OFF_AS,rho_AS] = deal([]);

        for it_ix = 1:num_its
            % Shuffle time_ixs
            T = size(dFF,2);
            time_ix_sh = time_ix + randsample(T,1);
            time_ix_sh(time_ix_sh > T) = time_ix_sh(time_ix_sh > T) - T;
            
            rho_sh_ = corrcoef(dFF(:,time_ix_sh)');
            
            % Remove doublecounting
            rho_ON_sh(:,it_ix) = rho_sh_(triu(J_ON,1)==1);
            rho_OFF_sh(:,it_ix) = rho_sh_(triu(J_OFF,1)==1);
            rho_sh(:,it_ix) = rho_sh_(triu(ones(size(rho_sh_)),1)==1);
            
            time_ix_AS = randsample(ix_A,numel(time_ix));
            rho_AS_ = corrcoef(dFF(:,time_ix_AS)');
            
            % Remove doublecounting
            rho_ON_AS(:,it_ix) = rho_AS_(triu(J_ON,1)==1);
            rho_OFF_AS(:,it_ix) = rho_AS_(triu(J_OFF,1)==1);
            rho_AS(:,it_ix) = rho_AS_(triu(ones(size(rho_AS_)),1)==1);
            
        end
        
        rho_all_sh{dataset_ix} = nanmean(rho_sh,2);
        rho_ON_all_sh{dataset_ix} = nanmean(rho_ON_sh,2);
        rho_OFF_all_sh{dataset_ix} = nanmean(rho_OFF_sh,2);
        
        rho_all_AS{dataset_ix} = nanmean(rho_AS,2);
        rho_ON_all_AS{dataset_ix} = nanmean(rho_ON_AS,2);
        rho_OFF_all_AS{dataset_ix} = nanmean(rho_OFF_AS,2); 
    else
        disp('This mouse had no running periods.');
    end
        
end

save([basedir,'processed/pw_corr_onsets'],'rho_all','rho_ON_all','rho_OFF_all',...
    'rho_all_sh','rho_ON_all_sh','rho_OFF_all_sh','rho_all_AS','rho_ON_all_AS','rho_OFF_all_AS');

%%
bins = -1:.01:1;

figure, histogram(vertcat(rho_ON_all_AS{:}),bins,'FaceColor','r','EdgeColor','r','Normalization','probability')
hold on, histogram(vertcat(rho_ON_all{:}),bins,'FaceColor','w','EdgeColor','r','Normalization','probability')
xlabel('Correlation'),ylabel('Probability'), title('PM axons')
set(gca,'FontSize',15)

hold on, plot(mean(vertcat(rho_ON_all{:})),.03,'vr')
hold on, plot(mean(vertcat(rho_ON_all_AS{:})),.03,'vr','MarkerFaceColor','r')

legend({'Random times (AS)','Onsets','',''})

figure, histogram(vertcat(rho_OFF_all_AS{:}),bins,'FaceColor','b','EdgeColor','b','Normalization','probability')
hold on, histogram(vertcat(rho_OFF_all{:}),bins,'FaceColor','w','EdgeColor','b','Normalization','probability')
xlabel('Correlation'),ylabel('Probability'), title('NM axons')
set(gca,'FontSize',15)

hold on, plot(mean(vertcat(rho_OFF_all{:})),.045,'vb')
hold on, plot(mean(vertcat(rho_OFF_all_AS{:})),.045,'vb','MarkerFaceColor','b')

legend({'Random times (AS)','Onsets','',''})

figure, histogram(vertcat(rho_all_AS{:}),bins,'FaceColor','k','EdgeColor','k','Normalization','probability')
hold on, histogram(vertcat(rho_all{:}),bins,'FaceColor','w','EdgeColor','k','Normalization','probability')
xlabel('Correlation'),ylabel('Probability'), title('All axons')
set(gca,'FontSize',15)

hold on, plot(mean(vertcat(rho_all{:})),.045,'vk')
hold on, plot(mean(vertcat(rho_all_AS{:})),.045,'vk','MarkerFaceColor','k')

legend({'Random times (AS)','Onsets','',''})
%%
c = [.5,.5,.5];

rho_all_mean = nan(13,1);
rho_all_AS_mean = nan(13,1);

rho_ON_all_mean = nan(13,1);
rho_ON_all_AS_mean = nan(13,1);

rho_OFF_all_mean = nan(13,1);
rho_OFF_all_AS_mean = nan(13,1);

figure,  hold on
for dataset_ix = 1:13
    rho_all_mean(dataset_ix) = mean(rho_all{dataset_ix});
    rho_all_AS_mean(dataset_ix) = mean(rho_all_AS{dataset_ix});
    plot([0,1],[rho_all_mean(dataset_ix),rho_all_AS_mean(dataset_ix)],'k')
    
    rho_ON_all_mean(dataset_ix) = mean(rho_ON_all{dataset_ix});
    rho_ON_all_AS_mean(dataset_ix) = mean(rho_ON_all_AS{dataset_ix});
    plot([2,3],[rho_ON_all_mean(dataset_ix),rho_ON_all_AS_mean(dataset_ix)],'k')
    
    rho_OFF_all_mean(dataset_ix) = mean(rho_OFF_all{dataset_ix});
    rho_OFF_all_AS_mean(dataset_ix) = mean(rho_OFF_all_AS{dataset_ix});
    plot([4,5],[rho_OFF_all_mean(dataset_ix),rho_OFF_all_AS_mean(dataset_ix)],'k')
end
plot(zeros,rho_all_mean,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones,rho_all_AS_mean,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(0+.2*[-1,1],nanmean(rho_all_mean)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],nanmean(rho_all_AS_mean)*[1,1],'k','LineWidth',3)

plot(2*ones,rho_ON_all_mean,'o','MarkerFaceColor','w','Color','m','MarkerSize',8)
plot(3*ones,rho_ON_all_AS_mean,'o','MarkerFaceColor','w','Color','m','MarkerSize',8)
plot(2+.2*[-1,1],nanmean(rho_ON_all_mean)*[1,1],'k','LineWidth',3)
plot(3+.2*[-1,1],nanmean(rho_ON_all_AS_mean)*[1,1],'k','LineWidth',3)

plot(4*ones,rho_OFF_all_mean,'o','MarkerFaceColor','w','Color','c','MarkerSize',8)
plot(5*ones,rho_OFF_all_AS_mean,'o','MarkerFaceColor','w','Color','c','MarkerSize',8)
plot(4+.2*[-1,1],nanmean(rho_OFF_all_mean)*[1,1],'k','LineWidth',3)
plot(5+.2*[-1,1],nanmean(rho_OFF_all_AS_mean)*[1,1],'k','LineWidth',3)

set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Angle (rad.)'))
xlim([-.5,5.5])

%% Plot example of speed onsets

dataset_ix = 1;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[~,~,~,~,speed] = load_behav_data(dataset_ix,time);
onset_indices = get_onsets(speed,acquisition_rate);

bins = onset_indices(1,2) - onset_indices(1,1) + 1;

figure, hold on
speed_onsets = zeros(bins,0);
for k = 1:size(onset_indices,1)
   speed_onsets(:,k) = speed(onset_indices(k,1):onset_indices(k,2));
   plot((1:bins)/acquisition_rate,speed_onsets(:,k),'Color',[.6,.6,.6],'LineWidth',1)
end
plot((1:bins)/acquisition_rate,mean(speed_onsets,2),'k','LineWidth',2)
set(gca,'FontSize',15)

xlabel('time (s)')
ylabel('speed')

%% Nearest neighbour distances
% between ROIs

clear all; clc

define_dirs;

NN_dist = cell(13,1);
NN_dist_ON = cell(13,1);
NN_dist_ON_shuff = cell(13,1);
NN_dist_OFF = cell(13,1);
NN_dist_OFF_shuff = cell(13,1);

for dataset_ix = 1:13
    
    [NN_dist{dataset_ix},NN_dist_ON{dataset_ix},NN_dist_OFF{dataset_ix},...
        NN_dist_ON_shuff{dataset_ix},NN_dist_OFF_shuff{dataset_ix}] = get_NN_dist(dataset_ix);

end

figure, hold on

temp = vertcat(NN_dist_ON{:}); temp = temp(~isnan(temp));
plot([1,1],mean(temp)+[-1,1]*std(temp)/sqrt(numel(temp)),'r')
plot(1,mean(temp),'sr','MarkerFaceColor','r')

temp_shuff = vertcat(NN_dist_ON_shuff{:}); temp_shuff = temp_shuff(~isnan(temp_shuff));
plot([2,2],mean(temp_shuff)+[-1,1]*std(temp_shuff)/sqrt(numel(temp_shuff)),'r')
plot(2,mean(temp_shuff),'sr')

signrank(temp,temp_shuff)

temp = vertcat(NN_dist_OFF{:}); temp = temp(~isnan(temp));
plot([3,3],mean(temp)+[-1,1]*std(temp)/sqrt(numel(temp)),'b')
plot(3,mean(temp),'sb','MarkerFaceColor','b')

temp_shuff = vertcat(NN_dist_OFF_shuff{:}); temp_shuff = temp_shuff(~isnan(temp_shuff));
plot([4,4],mean(temp_shuff)+[-1,1]*std(temp_shuff)/sqrt(numel(temp_shuff)),'b')
plot(4,mean(temp_shuff),'sb')

signrank(temp,temp_shuff)

set(gca,'FontSize',15)
ylabel('NN distances')
set(gca,'Box','off')
xlim([0,5])
ylim([0,6])