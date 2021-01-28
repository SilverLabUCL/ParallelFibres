% This script calculates the dimensionality of the datasets


clear all; clc

% Common code to define base directory and datasets
define_dirs;

%% Calculate dimensionality for chosen subpopulation size

N_sub = 650;

dimmax = nan(13,1);
varexp = cell(13,1);

tic

for dataset_ix = 1:13
    
    toc
    
    % Load data
    [dFF,~,acquisition_rate] = load_data(dataset_ix,1);

    if size(dFF,1) > N_sub

        % Calculate dimensionality for grouped axons
        [varexp{dataset_ix},dimmax(dataset_ix),~] = get_dim(dFF,N_sub,N_sub,acquisition_rate);


    end


end

toc

save([basedir,'Processed/dimensionality_N',num2str(N_sub)],'dimmax','varexp')

%% Plot var explained vs. number components
% Figure 6A

N = 300;

load([basedir,'Processed/dimensionality_N',num2str(N)])

c_87=[.1,.5,.8];
c_76=[1,.7,0];
c_77=[.7,.1,.7];
c_s=[1,.4,.4];
c_95 =[.6,.6,.6];

 dataset_colors = {c_87,c_87,c_87,c_87,c_77,c_s,c_s,c_95,c_87,c_87,c_s,c_s,c_76};

figure, hold on
varmax = nan(13,1);
dimmax = nan(13,1);
for dataset_ix = 1:13
    if ~isempty(varexp{dataset_ix})
        plot_error_snake(1:N,varexp{dataset_ix},dataset_colors{dataset_ix})
        [varmax(dataset_ix),dimmax(dataset_ix)] = max(nanmean(varexp{dataset_ix},1));
        plot(dimmax(dataset_ix),varmax(dataset_ix)+.03,'v','MarkerFaceColor',dataset_colors{dataset_ix},'MarkerEdgeColor',dataset_colors{dataset_ix})
    end
end
xlabel('Number of components')
ylabel('Variance explained (cross-val)')
ylim([0,.7])

%% Inset plot
varexp_1to5 = nan(13,5);
figure, hold on
for dataset_ix = 1:13
    if ~isempty(varexp{dataset_ix})
        varexp_1to5(dataset_ix,:) = mean(varexp{dataset_ix}(:,1:5),1);
        plot(1:5,varexp_1to5(dataset_ix,:),'-','Color',[.6,.6,.6],'LineWidth',2)
    end
end
bar(1:5,nanmean(varexp_1to5),'FaceAlpha',0,'LineWidth',2)
xlabel('Num. components')
ylabel('Var. exp.')
set(gca,'FontSize',15,'XTick',1:5)
ylim([0,.7])
%% Extrapolate to maximum dimensionality
% Figure 6B

varmax_all = []; dimmax_all = []; col_all = [];
for dataset_ix = 1:13
    if ~isempty(varexp{dataset_ix})
    for k = 1:10
        [varmax_,dimmax_] = max(varexp{dataset_ix}(k,:));
        varmax_all = [varmax_all; varmax_];
        dimmax_all = [dimmax_all; dimmax_];
    end
    end
end
figure, 
ix = find(~isnan(varmax_all));
plot(varmax_all,dimmax_all,'vk','MarkerFaceColor',[.6,.6,.6],'MarkerEdgeColor',[.6,.6,.6])
hold on, plot([0,1],[0,1]*(varmax_all(ix)'*varmax_all(ix))\(varmax_all(ix)'*dimmax_all(ix)),'k')

for dataset_ix = 1:13
    if ~isnan(varmax(dataset_ix))
        plot(varmax(dataset_ix),dimmax(dataset_ix),'vk','MarkerFaceColor',dataset_colors{dataset_ix},'MarkerEdgeColor',dataset_colors{dataset_ix})
    end
end
xlabel('Max variance explained')
ylabel('Number of components')
set(gca,'FontSize',18)


%% Extrapolate extrapolation for different subsample populations
% Figure 4B

N_sub = 100:50:650;

slope = zeros(size(N_sub));
slope_rois = zeros(size(N_sub));
for k = 1:length(N_sub)
    
    load([basedir,'Processed/dimensionality_N',num2str(N_sub(k))])

    % Recalculate max variance and dimensionality
    varmax = nan(13,1);
    dimmax = nan(13,1);
    for dataset_ix = 1:13
        if ~isempty(varexp{dataset_ix})
            [varmax(dataset_ix),dimmax(dataset_ix)] = max(nanmean(varexp{dataset_ix},1));
        end
    end
    
    ix = find(~isnan(varmax)); 
    slope(k) = (varmax(ix)'*varmax(ix))\(varmax(ix)'*dimmax(ix));
end

figure, bar(N_sub',N_sub./slope,'FaceColor',[.6,.6,.6],'EdgeColor','w','LineWidth',.8)
set(gca, 'FontSize',15, 'Box', 'off','XTick',100:100:700)
xtickangle(45), xlim([50,750])
xlabel('Number of neurons')
ylabel('Neurons per dimension')

%% Testing dimensionality for random matrices

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
    for noise =  logspace(log10(.02),log10(.2),15)
        
        F = U_true*SV_true + noise*randn(N,T);
        [varexp{it,noise_ix},dimmax(it,noise_ix),varmax(it,noise_ix)] = get_dim(F,N,N);
        noise_ix = noise_ix+1;

    end
end

save([basedir,'Processed/dimensionality_sim'],'varmax','dimmax','varexp','notes')

%% Plot extrapolation with noise
% Figure S6
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

figure, plot(varmax,dimmax,'.','MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',1.2)
xlabel('Variance explained')
ylabel('Number of components')
set(gca,'FontSize',18)

%% Example dimensionality
% Figure S6

it = 1
varexp{it,noise_ix}(k,:)


figure, hold on
varmax = nan(15,1);
dimmax = nan(15,1);
for noise_ix = 8:15
    plot_error_snake(1:N,varexp{it,noise_ix},'k')
end
xlabel('Number of components')
ylabel('Variance explained (cross-val)')
ylim([0,.7])


%% Dimensionality for PUFF data !

N_sub = 80;

dimmax = nan(4,1);
varexp = cell(4,1);

dimmax_shuff = nan(4,1);
varexp_shuff = cell(4,1);

tic

for dataset_ix = 1:4
    
    toc
    
    % Load data
    [dFF,~,acquisition_rate] = load_data(dataset_ix,1);
    [time_ix,time_ix_shuff] = get_puff_times(dataset_ix);

    if size(dFF,1) > N_sub

        % Calculate dimensionality for grouped axons
        [varexp{dataset_ix},dimmax(dataset_ix),~] = get_dim(dFF(:,time_ix),N_sub,N_sub,acquisition_rate,[],[],100);
        [varexp_shuff{dataset_ix},dimmax_shuff(dataset_ix),~] = get_dim(dFF(:,time_ix_shuff),N_sub,N_sub,acquisition_rate,[],[],100);

        %clear dFF

        % Calculate dimensionality for ungrouped ROIs
        %[dFF,~,acquisition_rate] = load_data(dataset_ix,0);
        %[varexp_rois{dataset_ix},dimmax_rois(dataset_ix),varmax_rois(dataset_ix)] = get_dim(dFF,N_sub,N_sub,acquisition_rate);
        
    end


end

toc
%%

c_87=[.1,.5,.8];
c_76=[1,.7,0];
c_77=[.7,.1,.7];
c_s=[1,.4,.4];

dataset_colors = {c_87,c_77,c_s,c_76};

figure, hold on
varmax = nan(4,1);
dimmax = nan(4,1);
for dataset_ix = 1:4
    if ~isempty(varexp{dataset_ix})
        plot_error_snake(1:80,varexp{dataset_ix},dataset_colors{dataset_ix})
        [varmax(dataset_ix),dimmax(dataset_ix)] = max(nanmean(varexp{dataset_ix},1));
        plot(dimmax(dataset_ix),varmax(dataset_ix)+.03,'v','MarkerFaceColor',dataset_colors{dataset_ix},'MarkerEdgeColor',dataset_colors{dataset_ix})
    end
end
% xlabel('Number of components')
ylabel('Variance explained (cross-val)')
ylim([0,.3])


figure, hold on
varmax = nan(4,1);
dimmax = nan(4,1);
for dataset_ix = 1:4
    if ~isempty(varexp_shuff{dataset_ix})
        plot_error_snake(1:80,varexp_shuff{dataset_ix},dataset_colors{dataset_ix})
        [varmax(dataset_ix),dimmax(dataset_ix)] = max(nanmean(varexp_shuff{dataset_ix},1));
        plot(dimmax(dataset_ix),varmax(dataset_ix)+.03,'v','MarkerFaceColor',dataset_colors{dataset_ix},'MarkerEdgeColor',dataset_colors{dataset_ix})
    end
end
xlabel('Number of components')
ylabel('Variance explained (cross-val)')
title('Shuffled times')
ylim([0,.3])

