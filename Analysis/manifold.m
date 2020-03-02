% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;


%% Plot 3d manifold
% Figure 2B

dataset_ix = 1;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,1);
cd = get_coding_dimension(dFF,A,QW);
A_or_QW = cd'*(dFF - nanmean(dFF,2));

[coeff, score] = pca(dFF');

A_or_QW = zscore(A_or_QW);
c = nan(size(A_or_QW,2),3);

figure, hold on
for k = 1:length(score)-1
    c_ = (A_or_QW(k)+1)/2;
    if c_ < 0
        c_ = 0;
    elseif c_ > 1
        c_ = 1;
    end
    c(k,:)=c_*[1,0,1]+(1-c_)*[0,1,1];
    
    plot3(score(k:k+1,1),score(k:k+1,2),score(k:k+1,3),'-','Color',c(k,:))
end
view(44,44)
set(gca,'FontSize',15,'XTick',[],'YTick',[],'ZTick',[])
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

%%  Distance between manifolds
% Figure 2C

angle_A_QW = nan(15,1);
angle_shuff = nan(15,1);

dist = @(x,y) sqrt(sum((x-y).^2));

dist_in_A_all = nan(15,1);
dist_in_QW_all = nan(15,1);
dist_A_QW_all = nan(15,1);

figure, hold on
for dataset_ix = 1:15
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get score
    
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
    
    dFF_A = [];
    for k = 1:length(A)
        ix = A(k,1):A(k,2);
        dFF_A = [dFF_A, dFF(:,ix)];
    end
    
    dFF_QW = [];
    for k = 1:length(QW)
        ix = QW(k,1):QW(k,2);
        dFF_QW = [dFF_QW, dFF(:,ix)];
    end
    
    T_A = size(dFF_A,2);
    dist_in_A = nan(T_A,1); ix = 1;
    for t1 = 1:T_A
        for t2 = t1:T_A
            dist_in_A(ix) = dist(dFF_A(:,t1),dFF_A(:,t2));
            ix = ix+1;
        end
    end
    
    T_QW = size(dFF_QW,2);
    dist_in_QW = nan(T_QW,1); ix = 1;
    for t1 = 1:T_QW
        for t2 = t1:T_QW
            dist_in_QW(ix) = dist(dFF_QW(:,t1),dFF_QW(:,t2));
            ix = ix+1;
        end
    end

    dist_A_QW = nan(T_QW,1); ix = 1;
    for t1 = 1:T_A
        for t2 = 1:T_QW
            dist_A_QW(ix) = dist(dFF_A(:,t1),dFF_QW(:,t2));
            ix = ix+1;
        end
    end
    
    dist_in_A_all(dataset_ix) = nanmean(dist_in_A);
    dist_in_QW_all(dataset_ix) = nanmean(dist_in_QW);
    dist_A_QW_all(dataset_ix) = nanmean(dist_A_QW);
    
    plot([0,1,2],[dist_in_QW_all(dataset_ix),dist_in_A_all(dataset_ix),dist_A_QW_all(dataset_ix)],'Color',[.8,.8,.8],'LineWidth',2)
    plot(0,dist_in_QW_all(dataset_ix),'oc','MarkerFaceColor','w','LineWidth',2,'MarkerSize',8)
    plot(1,dist_in_A_all(dataset_ix),'om','MarkerFaceColor','w','LineWidth',2,'MarkerSize',8)
    plot(2,dist_A_QW_all(dataset_ix),'ok','MarkerFaceColor','w','LineWidth',2,'MarkerSize',8)
    
end

set(gca,'FontSize',15,'XTick',[])
xlim([-.7,2.7])
ylabel('Avg. distance')

signrank(dist_A_QW_all,dist_in_A_all)
signrank(dist_A_QW_all,dist_in_QW_all)


%% Angle between subspaces
% Figure 2E,F

angle_A_QW = nan(15,1);
angle_shuff = nan(15,1);
p_val = nan(15,1);

num_PCs = 1;
kmax = 1000;

for dataset_ix = 1:15
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);

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
            train_ixs = block_shuffle_time(T,acquisition_rate);%,.5);
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

signrank(angle_A_QW,angle_shuff)

figure, hold on
histogram(angle_shuff_dist,'Normalization','probability','EdgeColor',c,'FaceColor',c)
plot(angle_A_QW(dataset_ix),.25,'vk','MarkerFaceColor','k')
xlabel('Angle (rad.)'), ylabel('Probability')
set(gca,'FontSize',15)
xlim([0,2])
%% Angle between subspaces & cartesian ON/OFF
% 

angle_A_ON = nan(15,1);
angle_A_OFF = nan(15,1);
angle_QW_ON = nan(15,1);
angle_QW_OFF = nan(15,1);

num_PCs = 2;
kmax = 1000;

for dataset_ix = 1:15
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);

    [N,T] = size(dFF);
    if N > 100
        [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
        [change_dFF,p_val] = change_dFF_sig(dFF,A,QW,acquisition_rate);
        ix_ON = find(p_val < 0.05 & change_dFF>0);
        ix_OFF = find(p_val < 0.05 & change_dFF<0);
        ix_fail = find(p_val > 0.05);

        dFF_A = [];
        for k = 1:length(A)
            ix = A(k,1):A(k,2);
            dFF_A = [dFF_A, dFF(:,ix)];
        end
        coeff_A = pca(dFF_A');

        dFF_QW = [];
        for k = 1:length(QW)
            ix = QW(k,1):QW(k,2);
            dFF_QW = [dFF_QW, dFF(:,ix)];
        end
        coeff_QW = pca(dFF_QW');

        cartesian_ON = zeros(N,N);
        cartesian_ON(ix_ON,ix_ON) = eye(numel(ix_ON));
        cartesian_OFF = zeros(N,N);
        cartesian_OFF(ix_OFF,ix_OFF) = eye(numel(ix_OFF));

        
        angle_A_ON(dataset_ix) = subspace(coeff_A(:,1:num_PCs),cartesian_ON);
        angle_A_OFF(dataset_ix) = subspace(coeff_A(:,1:num_PCs),cartesian_OFF);

        angle_QW_ON(dataset_ix) = subspace(coeff_QW(:,1:num_PCs),cartesian_ON);
        angle_QW_OFF(dataset_ix) = subspace(coeff_QW(:,1:num_PCs),cartesian_OFF);
    end
end

c = [.5, .5, .5];

figure,  hold on
plot(zeros,angle_QW_ON,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones,angle_QW_OFF,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(2*ones,angle_A_ON,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(3*ones,angle_A_OFF,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(0+.2*[-1,1],nanmean(angle_QW_ON)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],nanmean(angle_QW_OFF)*[1,1],'k','LineWidth',3)
plot(2+.2*[-1,1],nanmean(angle_A_ON)*[1,1],'k','LineWidth',3)
plot(3+.2*[-1,1],nanmean(angle_A_OFF)*[1,1],'k','LineWidth',3)
set(gca,'Xtick',0:3)
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Angle (rad.)'))
xlim([-.5,3.5])
ylim([0,2])

%%  Distance between manifolds
% Figure 2C

angle_A_QW = nan(15,1);
angle_shuff = nan(15,1);

dist = @(x,y) sqrt(sum((x-y).^2));

dist_in_A_all = nan(15,1);
dist_in_QW_all = nan(15,1);
dist_A_QW_all = nan(15,1);

dist_onset_all = nan(15,1);

figure, hold on
for dataset_ix = 1:15
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get score
    
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
    cd = get_coding_dimension(dFF,A,QW);
    A_or_QW = cd'*(dFF - mean(dFF,2));
    
    sig = std(diff(A_or_QW));
    ix = find(diff(A_or_QW) > 2*sig);
    
    var_A = diag(diag(cov(dFF_A')));
    var_QW = diag(diag(cov(dFF_QW')));
    % Need to get 
    figure, plot(A_or_QW)
    hold on, plot(ix,A_or_QW(ix),'.r')
    pause


    dFF_A = [];
    for k = 1:length(A)
        ix = A(k,1):A(k,2);
        dFF_A = [dFF_A, dFF(:,ix)];
    end
    
    dFF_QW = [];
    for k = 1:length(QW)
        ix = QW(k,1):QW(k,2);
        dFF_QW = [dFF_QW, dFF(:,ix)];
    end
    
    T_A = size(dFF_A,2);
    dist_in_A = nan(T_A,1); ix = 1;
    for t1 = 1:T_A
        for t2 = t1:T_A
            dist_in_A(ix) = dist(dFF_A(:,t1),dFF_A(:,t2));
            ix = ix+1;
        end
    end
    
    T_QW = size(dFF_QW,2);
    dist_in_QW = nan(T_QW,1); ix = 1;
    for t1 = 1:T_QW
        for t2 = t1:T_QW
            dist_in_QW(ix) = dist(dFF_QW(:,t1),dFF_QW(:,t2));
            ix = ix+1;
        end
    end

    dist_A_QW = nan(T_QW,1); ix = 1;
    for t1 = 1:T_A
        for t2 = 1:T_QW
            dist_A_QW(ix) = dist(dFF_A(:,t1),dFF_QW(:,t2));
            ix = ix+1;
        end
    end
    
    dist_in_A_all(dataset_ix) = nanmean(dist_in_A);
    dist_in_QW_all(dataset_ix) = nanmean(dist_in_QW);
    dist_A_QW_all(dataset_ix) = nanmean(dist_A_QW);
    
    plot([0,1,2],[dist_in_QW_all(dataset_ix),dist_in_A_all(dataset_ix),dist_A_QW_all(dataset_ix)],'Color',[.8,.8,.8],'LineWidth',2)
    plot(0,dist_in_QW_all(dataset_ix),'oc','MarkerFaceColor','w','LineWidth',2,'MarkerSize',8)
    plot(1,dist_in_A_all(dataset_ix),'om','MarkerFaceColor','w','LineWidth',2,'MarkerSize',8)
    plot(2,dist_A_QW_all(dataset_ix),'ok','MarkerFaceColor','w','LineWidth',2,'MarkerSize',8)
    
end

set(gca,'FontSize',15,'XTick',[])
xlim([-.7,2.7])
ylabel('Avg. distance')

signrank(dist_A_QW_all,dist_in_A_all)
signrank(dist_A_QW_all,dist_in_QW_all)




