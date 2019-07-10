% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;

%%

[C_wsp,C_amp,C_spd,C_pupil] = deal(nan(7,1));

for dataset_ix = 1:7
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get score
    [~, score] = pca(dFF');
    
    % Correlation with pc1
    pc1 = score(:,1);
    temp = corr(whisk_set_point,pc1);
    
    if temp < 0
        pc1 = -pc1;
        C_wsp(dataset_ix) = -temp;
    else
        C_wsp(dataset_ix) = temp;
    end
    
    figure, plot(pc1)
    
    
    
    C_amp(dataset_ix) = corr(whisk_amp,pc1);
    C_spd(dataset_ix) = corr(speed,pc1);
    if ~isempty(pupil)
        C_pupil(dataset_ix) = corr(pupil,pc1);
    end
    
   


end













%% ALL FOLLOWING CODE FOR variability
%%
[offset_indices, onset_indices] = get_onsets(whisk_amp,acquisition_rate);

figure, x = whisk_amp;
plot(x,'k'), hold on
for k = 1:size(onset_indices,2)
    t = offset_indices(k):onset_indices(k);
    plot(t,x(t),'c')
    plot(onset_indices(k),x(onset_indices(k)),'ok','MarkerFaceColor','c')
end
size(onset_indices,2)
%%
[N,T] = size(dFF);

% [U,S,V] = svd(dFF');
% S(1,1)=0;
% dFF = (U*S*V)';

[coeff, score] = pca(dFF');
dFF = score';

num_bins = round(.2 * acquisition_rate);

dFF_pre_onset = zeros(N,1);

for t = 1:size(onset_indices,2)
    t_ = ((onset_indices(t)-num_bins)-1) : (onset_indices(t)-1);
    dFF_pre_onset = [dFF_pre_onset,dFF(:,t_)];
end

E = eig(cov(dFF_pre_onset')); L = sqrt(sum(E.^2))


num_reps = 1000;

L_shuff = zeros(num_reps,1);

for i = 1:num_reps
    dFF_shuff = zeros(N,1);
    
    for t = 1:size(onset_indices,2)
        t_start = randsample(offset_indices(t): ((onset_indices(t)-8*num_bins)-1),1);
        t_ = t_start:(t_start+num_bins);
        dFF_shuff = [dFF_shuff,dFF(:,t_)];
    end
    
    E = eig(cov(dFF_shuff'));
    L_shuff(i) = sqrt(sum(E.^2));
end

figure, histogram(L_shuff,'Normalization','probability')
hold on, plot([L,L],ylim,':k','LineWidth',2)
set(gca,'FontSize',15)
xlabel('Total variance')
ylabel('Probability')

p = sum(L_shuff < L)/num_reps

%% Plot 3d manifold
[coeff, score] = pca(dFF');

zsp = zscore(whisk_set_point);

figure, hold on
for k = 1:length(score)-1
    c = (zsp(k)+1)/2;
    if c < 0
        c = 0;
    elseif c > 1
        c = 1;
    end
    c=c*[1,0,1]+(1-c)*[0,1,1];
   plot3(score(k:k+1,1),score(k:k+1,2),score(k:k+1,3),'-','Color',c)
end
view(-53,-52)

plot3(score(onset_indices,1),score(onset_indices,2),score(onset_indices,3),'or')
plot3(score(offset_indices,1),score(offset_indices,2),score(offset_indices,3),'ok')

%%

%%

figure, plot(zscore(whisk_angle))
hold on, plot(zscore(score(:,1))+3)
hold on, plot(zscore(score(:,2))+6)
hold on, plot(zscore(score(:,3))+9)