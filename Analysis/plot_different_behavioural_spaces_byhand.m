% This script needs to be cleaned up
% This script reproduces Supplementary Figure X
% For this script, behavioural labels were chosen for 2 datasets by hand

clear all; clc
dataset_ix = 6; % set to either 5 or 6

%% Hand-chose periods of A and SW from behavioral data for mouse F77 (dataset 5) and FL_S (dataset 6)
if dataset_ix == 5
    A = [179,451;903,1319;1393,1755;7412,7780;9100,9119];
    SW = [1955,1997;3151,3229;3262,3453;3493,3573;3621,3731;4005,4068;4204,4300;4351,4397;4685,4737;4933,4989;5154,5198;5234,5282;5525,5568;5736,5812;5936,5975;6148,6233;6314,6377;6559,6589;6853,6917;7189,7232;7338,7397;7939,7994;8071,8121;8214,8245;8325,8349;8509,8539;8604,8647;8832,8891];
elseif dataset_ix == 6
    A = [529,772;3809,4200];
    SW = [60,81;169,203;370,400;481,510;818,841;869,886;905,935;977,1057;1093,1139;1210,1250;1528,1539;...
        1560,1592;1781,1807;1829,1863;2122,2154;2230,2248;2352,2376;2440,2467;2490,2501;...
        2686,2699;2691,2691;2831,2849;2868,2886;2972,2994;3051,3067;3323,3342;3383,3397;3475,3542;3603,3632;...
        3685,3721;3761,3782;3792,3801;4385,4405;4412,4439;4473,4491;4687,4714;4732,4760;4811,4850;4921,4942;5006,5031;5204,5235;5258,5299;5299,5329;5376,5404;...
        5459,5475;5530,5553;5568,5591;5715,5747;5758,5780;5819,5839;5861,5874;5923,5943;6059,6077;6159,6173;6193,6239;6252,6290;6369,6391;6409,6421;6440,6517;6568,6591;6641,6705;6752,6779];
end

%% Load data

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,loco,speed] = load_behav_data(dataset_ix,time);
pupil = load_pupil(dataset_ix,time);

% Get score
[coeff, score] = pca(dFF');

% Smooth by 250 ms for clarity
score(:,1) = smoothdata(score(:,1),'gaussian',round(.25 * acquisition_rate));
score(:,2) = smoothdata(score(:,2),'gaussian',round(.25 * acquisition_rate));
score(:,3) = smoothdata(score(:,3),'gaussian',round(.25 * acquisition_rate));

%% Plot behavioural data coloured by SW / A 

figure,
subplot(2,1,1), plot(time,whisk_set_point,'k','LineWidth',1.5), 
hold on, ylabel('WSP (deg.)'), set(gca,'FontSize',15)
subplot(2,1,2), plot(time,speed,'k','LineWidth',1.5), hold on
hold on, ylabel('Speed'), set(gca,'FontSize',15)
xlabel('Time (s)')

for k = 1:size(A,1)
    ix = A(k,1):A(k,2);
    subplot(2,1,1), plot(time(ix),whisk_set_point(ix),'Color',[.7,0,.8],'LineWidth',1.5)
    subplot(2,1,2), plot(time(ix),speed(ix),'Color',[.7,0,.8],'LineWidth',1.5)
end
for k = 1:size(SW,1)
    ix = SW(k,1):SW(k,2);
    subplot(2,1,1), plot(time(ix),whisk_set_point(ix),'Color',[1,.6,0],'LineWidth',1.5)
    subplot(2,1,2), plot(time(ix),speed(ix),'Color',[1,.6,0],'LineWidth',1.5)
end

%% Plot PC 1-3 coloured by SW / A

figure, plot3(score(:,1),score(:,2),score(:,3),'k','LineWidth',1.5)
hold on, xlabel('PC 1'), ylabel('PC 2'), zlabel('PC 3')
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'FontSize',15)

for k = 1:size(A,1)
    ix = A(k,1):A(k,2);
    plot3(score(ix,1),score(ix,2),score(ix,3),'Color',[.7,0,.8],'LineWidth',1.5)
end
for k = 1:size(SW,1)
    ix = SW(k,1):SW(k,2);
    plot3(score(ix,1),score(ix,2),score(ix,3),'Color',[1,.6,0],'LineWidth',1.5)
end

