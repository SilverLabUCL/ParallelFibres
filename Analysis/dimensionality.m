% This script calculates the dimensionality of the datasets


clear all; clc

% Common code to define base directory and datasets
define_dirs;

N_sub = 700;

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
% 
% %% Following fragments of code are for modelling what happens with code .. 
% 
% T = 10000;
% N = 150;
% D = 40;
% 
% V_true = randn(D,T);
% 
% S_vec = rand(D,1)*10;%zeros(D,1);
% %for i = 1:D
% %    S_vec(i) = i;
% %end
% S_true = diag(sort(S_vec,'descend'));
% 
% SV_true = S_true*V_true;
% 
% U_true = randn(N,D);
% U_true = orth(U_true);
% %%
% noise = .01;
% 
% F = U_true*SV_true + noise*randn(150,10000);
% 
% 
% %%
% D_max = 50; k_max = 1;
% varexp = zeros(k_max,D_max);
% for k = 1:k_max
%    
%     ix_x = sort(randsample(N,round(.8*N)));
%     ix_y = setdiff(1:N,ix_x)';
% 
%     Neuron_x = F(ix_x,:);
% 
%     Neuron_y = F(ix_y,:);
% 
%     [ num_dim, Test_t, y_est, res, U_x, U_y ] = peer_predict_dim( Neuron_x, Neuron_y, 1:round(.8*T), 1:D_max);
%     varexp(k,:) = mean(res,1);
% end
% 
% figure(1), hold on
% plot(num_dim,mean(varexp,1),'k')
% [maxval,argmax]=max(mean(varexp,1))
% %hold on
% 
% z = [z; argmax, maxval];
% %%
% 
% k = find(z(:,2)<.5);
% z_ = z(k,:);
% 
% k_rnd = randsample(length(k),5);
% z_ = z_(k_rnd,:);

