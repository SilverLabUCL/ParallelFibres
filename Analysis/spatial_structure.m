% This script goes through all datasets and calculates the correlations
% between ROIs

clear all; clc

define_dirs;
%%

rho_all = []; distances_all = []; %rho_shuff_all = [];
rho_ON_all = []; distances_ON_all = []; 
rho_OFF_all = []; distances_OFF_all = [];
for dataset_ix = 1:15

    % Get A and QW states
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Choose between 2d (patch based) and 3d (all patches) versions
    [rho,distances,rho_ON,distances_ON,rho_OFF,distances_OFF] = get_corr_vs_dist(dataset_ix,0);
    %[rho,distances,rho_ON,distances_ON,rho_OFF,distances_OFF] = get_corr_vs_dist_3D(dataset_ix);
    
    rho_all = [rho_all; rho];
    distances_all = [distances_all; distances];
    
    rho_ON_all = [rho_ON_all; rho_ON];
    distances_ON_all = [distances_ON_all; distances_ON];

    rho_OFF_all = [rho_OFF_all; rho_OFF];
    distances_OFF_all = [distances_OFF_all; distances_OFF];

    %N_pairs = size(rho,1);
    %rho_shuff = rho(randsample(1:N_pairs,N_pairs));
    %rho_shuff_all = [rho_shuff_all;rho_shuff];
        
end

clear rho distances

ranksum(rho_all(distances_all<2),rho_all(distances_all>=2))
%%
save([basedir,'processed/spatial_corr_grouped'],'rho_all','distances_all',...
    'rho_ON_all','distances_ON_all','rho_OFF_all','distances_OFF_all')

%% Plot correlation vs distance
% Figure 1F

load([basedir,'processed/spatial_corr_grouped'])

dbin = 3;
dist_bins_e = 0:dbin:60;
dist_bins_c = dist_bins_e(2:end) - dbin/2;

y_shuff_mean = zeros(size(dist_bins_c));
y_shuff_top = zeros(size(dist_bins_c));
y_shuff_bot = zeros(size(dist_bins_c));

[y_mean,y_top,y_bot,y_ON_mean,y_ON_top,y_ON_bot,y_OFF_mean,y_OFF_top,y_OFF_bot] = deal(zeros(size(dist_bins_c)));

figure, hold on

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
ylim([.1,.4])
%% Nearest neighbour distances
% between ROIs

clear all; clc

define_dirs;

NN_dist = cell(15,1);
NN_dist_ON = cell(15,1);
NN_dist_ON_shuff = cell(15,1);
NN_dist_OFF = cell(15,1);
NN_dist_OFF_shuff = cell(15,1);

for dataset_ix = 1:15
    
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