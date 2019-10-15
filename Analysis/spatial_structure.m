% This script goes through all datasets and calculates the correlations
% between ROIs

clear all; clc

define_dirs;

rho_all = []; distances_all = []; rho_shuff_all = [];
for dataset_ix = 1:8
    
    [rho,distances] = get_corr_vs_dist(dataset_ix,1);
    %[rho,distances] = get_corr_vs_dist(dataset_ix,0,[],1);
    
    
    rho_all = [rho_all; rho];
    distances_all = [distances_all; distances];
    
    N_pairs = size(rho,1);
    rho_shuff = rho(randsample(1:N_pairs,N_pairs));
    
    rho_shuff_all = [rho_shuff_all;rho_shuff];
        
end

clear rho distances

%% Plot density of all pw corr vs 
figure, dscatter(distances_all,rho_all)
colormap(gray)

set(gca,'FontSize',18)
xlabel('Distance between fibres (um)')
ylabel('Pairwise correlation')

%% Plot error snake of distances

dbin = .5;
dist_bins_e = 0:dbin:60;
dist_bins_c = dist_bins_e(2:end) - dbin/2;

y_shuff_mean = zeros(size(dist_bins_c));
y_shuff_top = zeros(size(dist_bins_c));
y_shuff_bot = zeros(size(dist_bins_c));

y_mean = zeros(size(dist_bins_c));
y_top = zeros(size(dist_bins_c));
y_bot = zeros(size(dist_bins_c));
for k = 1:length(dist_bins_c)
    ix = find(distances_all <= dist_bins_e(k+1) & distances_all > dist_bins_e(k));
    y_mean(k) = mean(rho_all(ix));
    ste = std(rho_all(ix))/sqrt(numel(ix));
    y_top(k) = y_mean(k) + ste;
    y_bot(k) = y_mean(k) - ste;
    
    y_shuff_mean(k) = mean(rho_shuff_all(ix));
    ste = std(rho_shuff_all(ix))/sqrt(numel(ix));
    y_shuff_top(k) = y_shuff_mean(k) + ste;
    y_shuff_bot(k) = y_shuff_mean(k) - ste;
end

figure
fill([dist_bins_c, fliplr(dist_bins_c)],[y_shuff_top, fliplr(y_shuff_bot)],'r','FaceColor','r','LineStyle','none','FaceAlpha',.3)
hold on, plot(dist_bins_c,y_shuff_mean,'r')
fill([dist_bins_c, fliplr(dist_bins_c)],[y_top, fliplr(y_bot)],'k','FaceColor','k','LineStyle','none','FaceAlpha',.3)
plot(dist_bins_c,y_mean,'k')

set(gca,'FontSize',18)
xlabel('Distance between fibres (um)')
ylabel('Pairwise correlation')
set(gca,'Box','off')

%%

expfun = @(params,x) params(1) + params(2) .* exp(-params(3) .* x); 
err = @(params) sum( ( expfun(params,distances_all) - rho_all ).^2 );
params_fit = fminsearch(err,[1,1,1]);
hold on, plot(0:.01:5,expfun(params_fit,0:.01:5),'r','LineWidth',3)

%%
% This script goes through all datasets and calculates the correlations
% between ROIs

clear all; clc

define_dirs;

rho_all = []; distances_all = []; rho_shuff_all = [];
for dataset_ix = 1:8
    
    [rho,distances] = get_F_vs_dist(dataset_ix);
    
    rho_all = [rho_all; rho];
    distances_all = [distances_all; distances];
    
        
end

clear rho distances

%% Plot density of all pw corr vs 
figure, dscatter(distances_all,rho_all)
colormap(gray)

set(gca,'FontSize',18)
xlabel('Distance between fibres (um)')
ylabel('Pairwise correlation')

%% Plot error snake of distances

dbin = .5;
dist_bins_e = 0:dbin:5;
dist_bins_c = dist_bins_e(2:end) - dbin/2;

y_shuff_mean = zeros(size(dist_bins_c));
y_shuff_top = zeros(size(dist_bins_c));
y_shuff_bot = zeros(size(dist_bins_c));

y_mean = zeros(size(dist_bins_c));
y_top = zeros(size(dist_bins_c));
y_bot = zeros(size(dist_bins_c));
for k = 1:length(dist_bins_c)
    ix = find(distances_all <= dist_bins_e(k+1) & distances_all > dist_bins_e(k));
    y_mean(k) = nanmean(rho_all(ix));
    ste =  std(rho_all(ix))/sqrt(numel(ix));
    y_top(k) = y_mean(k) + ste;
    y_bot(k) = y_mean(k) - ste;
end

figure, hold on
fill([dist_bins_c, fliplr(dist_bins_c)],[y_top, fliplr(y_bot)],'k','FaceColor','k','LineStyle','none','FaceAlpha',.3)
plot(dist_bins_c,y_mean,'k')

set(gca,'FontSize',18)
xlabel('Distance from center (um)')
ylabel('Avg norm F')
set(gca,'Box','off')

%%

expfun = @(params,x) params(1) + params(2) .* exp(-params(3) .* x); 
err = @(params) sum( ( expfun(params,distances_all) - rho_all ).^2 );
params_fit = fminsearch(err,[1,1,1]);
hold on, plot(0:.01:5,expfun(params_fit,0:.01:5),'r')
%%



expfun = @(params,x) 1 - params(1) + params(1) .* exp(-params(2) .* x); 
err = @(params) sum( ( expfun(params,distances_all) - rho_all ).^2 );
params_fit = fminsearch(err,[1,1]);
hold on, plot(0:.01:5,expfun(params_fit,0:.01:5),'r','LineWidth',3)