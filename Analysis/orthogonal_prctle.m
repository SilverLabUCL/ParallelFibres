%% Angle between subspaces
% Figure 2E,F

clear all; clc

define_dirs;

ix_ON = cell(13,1); C_ON = cell(13,1);
ix_OFF = cell(13,1); C_OFF = cell(13,1);
ix_fail = cell(13,1); C_fail = cell(13,1);


for dataset_ix = 1:13
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get A and QW states
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate); 
    [change_dFF,p_val] = change_dFF_sig(dFF,A,QW,acquisition_rate);
    
    ix_fail{dataset_ix} = find(p_val > 0.05);
    C_fail{dataset_ix} = change_dFF(ix_fail{dataset_ix});
    
    ix_ON{dataset_ix} = find(change_dFF > 0 & p_val < .05);
    C_ON{dataset_ix} = change_dFF(ix_ON{dataset_ix});
    
    ix_OFF{dataset_ix} = find(change_dFF < 0 & p_val < .05);
    C_OFF{dataset_ix} = change_dFF(ix_OFF{dataset_ix});
end

%%

angle_A_QW = nan(13,1);
angle_shuff = nan(13,1);
p_val = nan(13,1);

num_PCs = 2;
kmax = 1000;
kmax_rand = 50;

percentile = 100;

for dataset_ix = 1:13
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    
    if percentile > 0
        ix_delete_ON = ix_ON{dataset_ix}((C_ON{dataset_ix} > prctile(C_ON{dataset_ix},100 - percentile)));
        ix_delete_OFF = ix_OFF{dataset_ix}((C_OFF{dataset_ix} < prctile(C_OFF{dataset_ix}, percentile)));
        ix_delete = [ix_delete_ON; ix_delete_OFF];
        dFF(ix_delete,:) = [];
    end

    [N,T] = size(dFF);
    
    if N > 100

        [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);

        ix_A = [];
        for k = 1:length(A)
            ix_A = [ix_A, A(k,1):A(k,2)];
        end
        dFF_A = dFF(:,ix_A);
        coeff_A = pca(dFF_A');

        ix_QW = [];
        for k = 1:length(QW)
            ix_QW = [ix_QW, QW(k,1):QW(k,2)];
        end
        dFF_QW = dFF(:,ix_QW);
        coeff_QW = pca(dFF_QW');

        angle_A_QW(dataset_ix) = subspace(coeff_QW(:,1:num_PCs),coeff_A(:,1:num_PCs));

        angle_shuff_dist = zeros(kmax,1);
        for k = 1:kmax
            train_ixs = block_shuffle_time(T,acquisition_rate);
            test_ixs = train_ixs(1:round(T/2));
            train_ixs = setdiff(train_ixs,test_ixs); 

            coeff_1 = pca(dFF(:,test_ixs)');
            coeff_2 = pca(dFF(:,train_ixs)');
            angle_shuff_dist(k) = subspace(coeff_1(:,1:num_PCs),coeff_2(:,1:num_PCs));
        end

        angle_shuff(dataset_ix) = mean(angle_shuff_dist);
        p_val(dataset_ix) = sum(angle_shuff_dist > angle_A_QW(dataset_ix)) / kmax;

    end
end

c = [.5,.5,.5];

figure,  hold on
plot(zeros,angle_A_QW,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones,angle_shuff,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(0+.2*[-1,1],nanmean(angle_A_QW)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],nanmean(angle_shuff)*[1,1],'k','LineWidth',3)
set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Angle (rad.)'))
xlim([-.5,1.5])
ylim([0,1.6])
title(percentile)

signrank(angle_A_QW,angle_shuff)

save([basedir,'processed/orthogonal_',num2str(percentile),'thprctile'],'angle_A_QW','angle_shuff','p_val','ix_ON','C_ON','ix_OFF','C_OFF','ix_fail','C_fail');

%% Add random comparison - remove random neurons

percentile = 100;

filename = [basedir,'orthogonal_',num2str(percentile),'thprctile'];

load(filename)

angle_rand = nan(13,1);

num_PCs = 2;
kmax = 50;

for dataset_ix = 1:13
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    
    if percentile > 0
        ix_delete_ON = ix_ON{dataset_ix}((C_ON{dataset_ix} > prctile(C_ON{dataset_ix},100 - percentile)));
        ix_delete_OFF = ix_OFF{dataset_ix}((C_OFF{dataset_ix} < prctile(C_OFF{dataset_ix}, percentile)));
        num_to_delete = numel(ix_delete_ON) + numel(ix_delete_OFF);
    end

    [N,T] = size(dFF);
    if N > 100 + num_to_delete

        [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);

        ix_A = [];
        for k = 1:length(A)
            ix_A = [ix_A, A(k,1):A(k,2)];
        end
        dFF_A = dFF(:,ix_A);
        
        ix_QW = [];
        for k = 1:length(QW)
            ix_QW = [ix_QW, QW(k,1):QW(k,2)];
        end
        dFF_QW = dFF(:,ix_QW);
        
        angle_rand_dist = zeros(kmax,1);
        for k = 1:kmax
            ix_keep = randsample(N,N-num_to_delete);
            
            coeff_A = pca(dFF_A(ix_keep,:)');
            coeff_QW = pca(dFF_QW(ix_keep,:)');

            angle_rand_dist(k) = subspace(coeff_A(:,1:num_PCs),coeff_QW(:,1:num_PCs));
        end

        angle_rand(dataset_ix) = mean(angle_rand_dist);
    end
end


c = [.5,.5,.5];

figure,  hold on
plot(-ones,angle_A_QW,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(zeros,angle_rand,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones,angle_shuff,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(-1+.2*[-1,1],nanmean(angle_A_QW)*[1,1],'k','LineWidth',3)
plot(0+.2*[-1,1],nanmean(angle_rand)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],nanmean(angle_shuff)*[1,1],'k','LineWidth',3)
set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Angle (rad.)'))
xlim([-1.5,1.5])
ylim([0,1.6])
title(percentile)

save(filename,'angle_rand','-append')

%% Plot figure

clear all; clc

define_dirs;

c = [.5,.5,.5];

figure, hold on

angle_A_QW_mean = nan(1,8);
angle_A_QW_ste = nan(1,8);
angle_rand_mean = nan(1,8);
angle_rand_ste = nan(1,8);
angle_shuff_mean = nan(1,8);
angle_shuff_ste = nan(1,8);

for percentile = 0:10:100

    clear angle_rand angle_A_QW angle_shuff
    load([basedir,'orthogonal_',num2str(percentile),'thprctile'])
    
    k = 1+percentile/10;
    
    angle_A_QW_mean(k) = nanmean(angle_A_QW);
    angle_A_QW_ste(k) = nanstd(angle_A_QW) / sqrt(sum(~isnan(angle_A_QW)));
    
    if percentile > 0
        angle_rand_mean(k) = nanmean(angle_rand);
        angle_rand_ste(k) = nanstd(angle_rand) / sqrt(sum(~isnan(angle_rand)));
    elseif percentile == 0
        angle_rand_mean(k) = nanmean(angle_A_QW);
        angle_rand_ste(k) = nanstd(angle_A_QW) / sqrt(sum(~isnan(angle_A_QW)));
    end
    
    angle_shuff_mean(k) = nanmean(angle_shuff);
    angle_shuff_ste(k) = nanstd(angle_shuff) / sqrt(sum(~isnan(angle_shuff)));
    
    plot(percentile*[1,1],angle_A_QW_mean(k) + angle_A_QW_ste(k) * [-1,1],'k','LineWidth',1.5)
    plot(percentile*[1,1],angle_rand_mean(k) + angle_rand_ste(k) * [-1,1],'k','LineWidth',1.5)
    plot(percentile*[1,1],angle_shuff_mean(k) + angle_shuff_ste(k) * [-1,1],'Color',[.6,.6,.6],'LineWidth',1.5)

    find(isnan(angle_A_QW))
    [percentile,signrank(angle_A_QW,angle_shuff),signrank(angle_A_QW,angle_rand)]
end

plot(0:10:100, angle_A_QW_mean,'k','LineWidth',1.5)
plot(0:10:100, angle_rand_mean,':k','LineWidth',1.5)
plot(0:10:100, angle_shuff_mean,'Color',[.6,.6,.6],'LineWidth',1.5)

set(gca,'FontSize',15)
xlabel('Percentile')
ylabel('Angle (rad.)')
xlim([-10,110])
