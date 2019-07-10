% This script calculates the dimensionality of the datasets


clear all; clc

% Common code to define base directory and datasets
define_dirs;

N_sub = 300;

varmax = nan(8,1);
dimmax = nan(8,1);
varexp = cell(8,1);

varmax_rois = nan(8,1);
dimmax_rois = nan(8,1);
varexp_rois = cell(8,1);

tic

for dataset_ix = 1:8
    
    dataset_ix, toc
    
    % Load data
    [dFF,~,acquisition_rate] = load_data(dataset_ix,1);

    if size(dFF,1) > N_sub

        % Calculate dimensionality for grouped axons
        [varexp{dataset_ix},dimmax(dataset_ix),varmax(dataset_ix)] = get_dim(dFF,N_sub,N_sub,acquisition_rate);

        %clear dFF

        % Calculate dimensionality for ungrouped ROIs
        [dFF,~,acquisition_rate] = load_data(dataset_ix,0);
        [varexp_rois{dataset_ix},dimmax_rois(dataset_ix),varmax_rois(dataset_ix)] = get_dim(dFF,N_sub,N_sub,acquisition_rate);
        
    end


end

toc

save([basedir,'processed/dimensionality_N',num2str(N_sub)],'varmax','dimmax','varexp','varmax_rois','dimmax_rois','varexp_rois')

%% Plot var explained vs. number components

load([basedir,'processed/dimensionality_N300'])

figure, hold on
for dataset_ix = 1:8
    dim_tested = find(~isnan(sum(varexp{dataset_ix},1)));
    plot_error_snake(dim_tested,varexp{dataset_ix}(:,dim_tested),[0,0,0])
end
plot(dimmax,varmax+.03,'vk','MarkerFaceColor','k')
xlabel('Number of components')
ylabel('Variance explained (cross-val)')
ylim([0,.7])

figure, hold on
for dataset_ix = 1:8
    dim_tested = find(~isnan(sum(varexp_rois{dataset_ix},1)));
    plot_error_snake(dim_tested,varexp_rois{dataset_ix}(:,dim_tested),[1,0,0])
end
plot(dimmax_rois,varmax_rois+.03,'vk','MarkerFaceColor','k')
xlabel('Number of components')
ylabel('Variance explained (cross-val)')
ylim([0,.7])

%% Extrapolate to maximum dimensionality

varmax_all = []; dimmax_all = [];
for dataset_ix = 1:8
    if ~isempty(varexp{dataset_ix})
    for k = 1:10
        [varmax_,dimmax_] = max(varexp{dataset_ix}(k,:));
        varmax_all = [varmax_all; varmax_];
        dimmax_all = [dimmax_all; dimmax_];
    end
    end
end

varmax_rois_all = []; dimmax_rois_all = [];
for dataset_ix = 1:8
    if ~isempty(varexp_rois{dataset_ix})
    for k = 1:10
        [varmax_,dimmax_] = max(varexp_rois{dataset_ix}(k,:));
        varmax_rois_all = [varmax_rois_all; varmax_];
        dimmax_rois_all = [dimmax_rois_all; dimmax_];
    end
    end
end

figure, 
ix = find(~isnan(varmax_all));
plot(varmax_all,dimmax_all,'vk','MarkerFaceColor',[.6,.6,.6],'MarkerEdgeColor',[.6,.6,.6])
hold on, plot([0,1],[0,1]*(varmax_all(ix)'*varmax_all(ix))\(varmax_all(ix)'*dimmax_all(ix)),'k')
ix = find(~isnan(varmax));
plot(varmax,dimmax,'vk','MarkerFaceColor','k')
xlabel('Max variance explained')
ylabel('Number of components')
set(gca,'FontSize',18)

figure,
ix = find(~isnan(varmax_rois_all));
plot(varmax_rois_all,dimmax_rois_all,'vr','MarkerFaceColor',[1,.6,.6],'MarkerEdgeColor',[1,.6,.6]);
hold on, plot([0,1],[0,1]*(varmax_rois_all(ix)'*varmax_rois_all(ix))\(varmax_rois_all(ix)'*dimmax_rois_all(ix)),'k')
ix = find(~isnan(varmax_rois));
plot(varmax_rois,dimmax_rois,'vr','MarkerFaceColor','r')
xlabel('Max variance explained')
ylabel('Number of components')
set(gca,'FontSize',18)

%% Extrapolate extrapolation

N_sub = 150:50:700;

slope = zeros(size(N_sub));
slope_rois = zeros(size(N_sub));
for k = 1:length(N_sub)
    
    load([basedir,'processed/dimensionality_N',num2str(N_sub(k))])
    
    ix = find(~isnan(varmax)); 
    slope(k) = (varmax(ix)'*varmax(ix))\(varmax(ix)'*dimmax(ix));

    
    ix = find(~isnan(varmax_rois)); 
    slope_rois(k) = (varmax_rois(ix)'*varmax_rois(ix))\(varmax_rois(ix)'*dimmax_rois(ix));
end

figure, bar(N_sub',slope./N_sub,'FaceColor','k','EdgeColor','k','LineWidth',2)
set(gca, 'FontSize',18)
xtickangle(45), xlim([100,750])
xlabel('Population size (N)')
ylabel('Dimensionaity / N')
%% Following fragments of code are for modelling what happens with code .. 

T = 5000;
N = 300;
D = 60;

varexp = cell(25,5);
dimmax = nan(25,5);
varmax = nan(25,5);

tic

for it = 1:25
    [it,toc]
    V_true = randn(D,T);
    V_true = orth(V_true')';

    S_vec = exprnd(10,D,1);
    S_true = diag(sort(S_vec,'descend'));

    SV_true = S_true*V_true;

    U_true = randn(N,D);
    U_true = orth(U_true);

    noise_ix = 1;
    for noise =  logspace(log10(.02),log10(.2),15)%logspace(log10(.005),log10(.15),5)
        
        F = U_true*SV_true + noise*randn(N,T);
        [varexp{it,noise_ix},dimmax(it,noise_ix),varmax(it,noise_ix)] = get_dim(F,N,N);
        noise_ix = noise_ix+1;

    end
end
%%
save([basedir,'processed/dimensionality_sim'],'varmax','dimmax','varexp','notes')

%% Plot

varmax = []; dimmax = [];
for it = 1:25
    for noise_ix = 1:15
        for k = 1:10
            [varmax_,dimmax_] = max(varexp{it,noise_ix}(k,:));
            varmax= [varmax;varmax_];
            dimmax = [dimmax;dimmax_];
        end
    end
end

ix = find(~isnan(varmax));

figure, plot(varmax,dimmax,'v','MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',1.2)
xlabel('Variance explained')
ylabel('Number of components')
set(gca,'FontSize',18)

