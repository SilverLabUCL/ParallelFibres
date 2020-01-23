% This script goes through all datasets and calculates the correlations
% between ROIs

clear all; clc

define_dirs;

rho_all = []; distances_all = []; %rho_shuff_all = [];
rho_ON_all = []; distances_ON_all = []; 
rho_OFF_all = []; distances_OFF_all = [];
for dataset_ix = 1:15

    % Get A and QW states
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    [rho,distances,rho_ON,distances_ON,rho_OFF,distances_OFF] = get_corr_vs_dist(dataset_ix);
    
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

%% Plot error snake of distances

dbin = 3;
dist_bins_e = 0:dbin:60;
dist_bins_c = dist_bins_e(2:end) - dbin/2;

y_shuff_mean = zeros(size(dist_bins_c));
y_shuff_top = zeros(size(dist_bins_c));
y_shuff_bot = zeros(size(dist_bins_c));

[y_mean,y_top,y_bot,y_ON_mean,y_ON_top,y_ON_bot,y_OFF_mean,y_OFF_top,y_OFF_bot] = deal(zeros(size(dist_bins_c)));

for k = 1:length(dist_bins_c)
    ix = find(distances_all <= dist_bins_e(k+1) & distances_all > dist_bins_e(k));
    y_mean(k) = median(rho_all(ix));
    ste = std(rho_all(ix))/sqrt(numel(ix));
    y_top(k) = y_mean(k) + ste;
    y_bot(k) = y_mean(k) - ste;
    
    ix = find(distances_ON_all <= dist_bins_e(k+1) & distances_ON_all > dist_bins_e(k));
    y_ON_mean(k) = median(rho_ON_all(ix));
    ste = std(rho_ON_all(ix))/sqrt(numel(ix));
    y_ON_top(k) = y_ON_mean(k) + ste;
    y_ON_bot(k) = y_ON_mean(k) - ste;
    
    ix = find(distances_OFF_all <= dist_bins_e(k+1) & distances_OFF_all > dist_bins_e(k));
    y_OFF_mean(k) = median(rho_OFF_all(ix));
    ste = std(rho_OFF_all(ix))/sqrt(numel(ix));
    y_OFF_top(k) = y_OFF_mean(k) + ste;
    y_OFF_bot(k) = y_OFF_mean(k) - ste;
end

figure
fill([dist_bins_c, fliplr(dist_bins_c)],[y_ON_top, fliplr(y_ON_bot)],'r','FaceColor','m','LineStyle','none','FaceAlpha',.3)
hold on, plot(dist_bins_c,y_ON_mean,'m')
fill([dist_bins_c, fliplr(dist_bins_c)],[y_OFF_top, fliplr(y_OFF_bot)],'r','FaceColor','c','LineStyle','none','FaceAlpha',.3)
hold on, plot(dist_bins_c,y_OFF_mean,'c')
fill([dist_bins_c, fliplr(dist_bins_c)],[y_top, fliplr(y_bot)],'k','FaceColor','k','LineStyle','none','FaceAlpha',.3)
plot(dist_bins_c,y_mean,'k')

set(gca,'FontSize',18)
xlabel('Distance between fibres (um)')
ylabel('Pairwise correlation')
set(gca,'Box','off')

%% Plot error snake of distances

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

expfun = @(params,x) params(1) + params(2) .* exp(-params(3) .* x); 
err = @(params) sum( ( expfun(params,distances_all) - rho_all ).^2 );
params_fit = fminsearch(err,[1,1,1]);
x = 1:.01:59;
hold on, plot(x,expfun(params_fit,x),'k','LineWidth',3)

err = @(params) sum( ( expfun(params,distances_ON_all) - rho_ON_all ).^2 );
params_fit = fminsearch(err,[1,1,1]);
hold on, plot(x,expfun(params_fit,x),'m','LineWidth',3)

err = @(params) sum( ( expfun(params,distances_OFF_all) - rho_OFF_all ).^2 );
params_fit = fminsearch(err,[1,1,1]);
hold on, plot(x,expfun(params_fit,x),'c','LineWidth',3)

set(gca,'FontSize',18)
xlabel('Distance between fibres (um)')
ylabel('Pairwise correlation')
set(gca,'Box','off')
ylim([.1,.4])