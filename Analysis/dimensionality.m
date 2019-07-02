% This script calculates the dimensionality of the datasets


clear all; clc

% Common code to define base directory and datasets
define_dirs;

N_sub = 600;

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

        % Calculate dimensionality
        [varexp{dataset_ix},dimmax(dataset_ix),varmax(dataset_ix)] = get_dim(dFF,N_sub,N_sub,acquisition_rate);

        %clear dFF

        % Redo - ungrouped
        % Load data
        [dFF,~,acquisition_rate] = load_data(dataset_ix,0);

        % Calculate dimensionality
        [varexp_rois{dataset_ix},dimmax_rois(dataset_ix),varmax_rois(dataset_ix)] = get_dim(dFF,N_sub,N_sub,acquisition_rate);
        
    end


end

toc

save([basedir,'processed/dimensionality_N',num2str(N_sub)],'varmax','dimmax','varexp','varmax_rois','dimmax_rois','varexp_rois')

%% Plot var explained vs. number components

figure, hold on
for dataset_ix = 1:8
    dim_tested = find(~isnan(sum(varexp{dataset_ix},1)));
    plot_error_snake(dim_tested,varexp{dataset_ix}(:,dim_tested),[0,0,0])
end
xlabel('Number of components')
ylabel('Variance explained (cross-val)')

figure, hold on
for dataset_ix = 1:8
    dim_tested = find(~isnan(sum(varexp_rois{dataset_ix},1)));
    plot_error_snake(dim_tested,varexp_rois{dataset_ix}(:,dim_tested),[1,0,0])
end
xlabel('Number of components')
ylabel('Variance explained (cross-val)')

%% Extrapolate to maximum dimensionality

varmax = []; dimmax = [];
for dataset_ix = 1:8
    if ~isempty(varexp{dataset_ix})
    for k = 1:10
        [varmax_,dimmax_] = max(varexp{dataset_ix}(k,:));
        varmax= [varmax;varmax_];
        dimmax = [dimmax;dimmax_];
    end
    end
end

%%
ix = find(~isnan(varmax));

figure, plot(varmax,dimmax,'ok','MarkerFaceColor','w','LineWidth',2)
hold on, plot([0,1],[0,1]*(varmax(ix)'*varmax(ix))\(varmax(ix)'*dimmax(ix)),'k')
xlabel('Variance explained')
ylabel('Number of components')
set(gca,'FontSize',18)

figure, plot(varmax_rois,dimmax_rois,'or','MarkerFaceColor','w','LineWidth',2)
hold on, plot([0,1],[0,1]*(varmax_rois(ix)'*varmax_rois(ix))\(varmax_rois(ix)'*dimmax_rois(ix)),'k')
xlabel('Variance explained')
ylabel('Number of components')
set(gca,'FontSize',18)

%% Extrapolate extrapolation

N_sub = 150:50:700;

for k = 1:length(N_sub)
    
    load([basedir,'processed/dimensionality_N',num2str(N_sub(k))])
    
    ix = find(~isnan(varmax)); slope = (varmax(ix)'*varmax(ix))\(varmax(ix)'*dimmax(ix));
    disp('Axons:')
    [N_sub(k), slope, slope/N_sub(k)]
    
    ix = find(~isnan(varmax_rois)); slope = (varmax_rois(ix)'*varmax_rois(ix))\(varmax_rois(ix)'*dimmax_rois(ix));
    disp('ROIs:')
    [N_sub(k), slope, slope/N_sub(k)]
end

%% Calculate dimensionality and iterate over all 

%% Following fragments of code are for modelling what happens with code .. 

T = 10000;
N_sub = 300;
D = 6;

V_true = randn(D,T);
V_true = orth(V_true')';

S_vec = rand(D,1)*10;
S_true = diag(sort(S_vec,'descend'));

SV_true = S_true*V_true;

U_true = randn(N,D);
U_true = orth(U_true);
%%
dimmax = [];
varmax = [];
for noise = linspace(.005,.05,10)
    noise
    F = U_true*SV_true + noise*randn(N,T);
    [~,dimmax_,varmax_] = get_dim(F,N_sub,N_sub);
    dimmax=[dimmax,dimmax_];
    varmax=[varmax,varmax_];
end
%%
figure(1), hold on
plot(num_dim,mean(varexp,1),'k')
[maxval,argmax]=max(mean(varexp,1))
%hold on

z = [z; argmax, maxval];
%%

k = find(z(:,2)<.5);
z_ = z(k,:);

k_rnd = randsample(length(k),5);
z_ = z_(k_rnd,:);

