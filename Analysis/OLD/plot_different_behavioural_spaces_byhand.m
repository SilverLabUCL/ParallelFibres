% This script needs to be cleaned up
% This script reproduces Supplementary Figure X
% For this script, behavioural labels were chosen for 2 datasets by hand

clear all; clc
dataset_ix = 5; % set to either 5 or 6

%% Hand-chose periods of A and SW from behavioral data for mouse F77 (dataset 5) and FL_S (dataset 6)
if dataset_ix == 5
    A = [179,451;903,1319;1393,1755;7412,7780;9100,9119];
    SW = [1955,1997;3151,3229;3262,3453;3493,3573;3621,3731;4005,4068;4204,4300;4351,4397;4685,4737;4933,4989;5154,5198;5234,5282;5525,5568;5736,5812;5936,5975;6148,6233;6314,6377;6559,6589;6853,6917;7189,7232;7338,7397;7939,7994;8071,8121;8214,8245;8325,8349;8509,8539;8604,8647;8832,8891];
elseif dataset_ix == 6;
    A = [529,772;3809,4200];
    SW = [60,81;169,203;370,400;481,510;818,841;869,886;905,935;977,1057;1093,1139;1210,1250;1528,1539;...
        1560,1592;1781,1807;1829,1863;2122,2154;2230,2248;2352,2376;2440,2467;2490,2501;...
        2686,2699;2691,2691;2831,2849;2868,2886;2972,2994;3051,3067;3323,3342;3383,3397;3475,3542;3603,3632;...
        3685,3721;3761,3782;3792,3801;4385,4405;4412,4439;4473,4491;4687,4714;4732,4760;4811,4850;4921,4942;5006,5031;5204,5235;5258,5299;5299,5329;5376,5404;...
        5459,5475;5530,5553;5568,5591;5715,5747;5758,5780;5819,5839;5861,5874;5923,5943;6059,6077;6159,6173;6193,6239;6252,6290;6369,6391;6409,6421;6440,6517;6568,6591;6641,6705;6752,6779];
end
%% Chosen from activity space 
if dataset_ix == 5
    A = [180,454;906,1306;1394,1763;7339,7782;9099,9119];
    SW = [1958,2031;3171,3336;3391,3461;3626,3704;4006,4078;4222,4290;4351,4399;4690,4762;4932,5010;5164,5227;5242,5295;5526,5581;5745,5824;5942,5999;6148,6237;6322,6376;6564,6612;6829,6917;7190,7249;7957,8030;8078,8124;8222,8264;8326,8369;8509,8557;8612,8659;8847,8901];
elseif dataset_ix == 6
    A = [528,777;3765,4195];
    SW = [49,102;170,214;367,398;481,507;6234,6251;...
    2830,2856;2789,2816;253,289;6633,6653;4871,4900;2554,2569;248,296;975,989;996,1011;...
        821,849;905,935;1097,1142;1238,1252;1324,1345;1395,1424;1492,1516;1530,1546;1562,1597;1777,1814;1829,1863;1893,1912;1995,2011;2134,2156;2233,2258;2357,2382;2438,2469;2489,2510;2686,2705;2876,2900;2979,2995;3053,3081;3181,3190;3226,3243;3326,3349;3383,3410;3526,3542;3576,3590;3605,3622;3682,3730;...
        4349,4371;4423,4442;4475,4505;4574,4597;4637,4653;4690,4717;4738,4755;4815,4858;4872,4869;4920,4940;4957,4973;5007,5044;5056,5091;5170,5232;5260,5330;5382,5416;5461,5486;5535,5598;5616,5670;5723,5752;5761,5780;5809,5834;5867,5885;5918,5951;6022,6040;6062,6087;6111,6127;6159,6177;6256,6296;6367,6396;6418,6476;6483,6512;6653,6710;6763,6779];
end

%% Load data

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
pupil = load_pupil(dataset_ix,time);

% Get score
[coeff, score] = pca(dFF');

% Smooth by 250 ms for clarity
score(:,1) = smoothdata(score(:,1),'gaussian',round(.25 * acquisition_rate));
score(:,2) = smoothdata(score(:,2),'gaussian',round(.25 * acquisition_rate));
score(:,3) = smoothdata(score(:,3),'gaussian',round(.25 * acquisition_rate));

smooth_win = round(acquisition_rate * .5);

x = whisk_amp(1:end-1);
x = smoothdata(x,'movmedian',smooth_win);
x = x - mode(round(x,4));
x = x / std(x);

% Running information based on wheel MI
y = abs(diff(speed));
y = smoothdata(y,'movmedian',smooth_win);
y = y - mode(round(y,4));
y = y / std(y);

%% Plot behavioural data coloured by SW / A 

figure,
subplot(2,1,1), plot((1:length(x))/acquisition_rate,x,'k','LineWidth',1.5), 
hold on, ylabel('Whisk'), set(gca,'FontSize',15)
subplot(2,1,2), plot((1:length(y))/acquisition_rate,y+5,'k','LineWidth',1.5), hold on
hold on, ylabel('Loc'), set(gca,'FontSize',15)
xlabel('Time (s)')

for k = 1:size(A,1)
    ix = A(k,1):A(k,2);
    subplot(2,1,1), plot(ix/acquisition_rate,x(ix),'Color',[1,0,1],'LineWidth',1.5)
    subplot(2,1,2), plot(ix/acquisition_rate,y(ix)+5,'Color',[1,0,1],'LineWidth',1.5)
end
for k = 1:size(SW,1)
    ix = SW(k,1):SW(k,2);
    subplot(2,1,1), plot(ix/acquisition_rate,x(ix),'Color',[1,.6,0],'LineWidth',1.5)
    subplot(2,1,2), plot(ix/acquisition_rate,y(ix)+5,'Color',[1,.6,0],'LineWidth',1.5)
end

%% Plot PC 1-3 coloured by SW / A

figure, plot3(score(:,1),score(:,2),score(:,3),'k','LineWidth',1.5)
hold on, xlabel('PC 1'), ylabel('PC 2'), zlabel('PC 3')
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'FontSize',15)

for k = 1:size(A,1)
    ix = A(k,1):A(k,2);
    plot3(score(ix,1),score(ix,2),score(ix,3),'Color',[1,0,1],'LineWidth',1.5)
end
for k = 1:size(SW,1)
    ix = SW(k,1):SW(k,2);
    plot3(score(ix,1),score(ix,2),score(ix,3),'Color',[1,.6,0],'LineWidth',1.5)
end

%%

% %%
% figure, plot(score(:,1),1:length(score(:,1)),'k','LineWidth',1.5), hold on
% for k = 1:size(A,1)
%     ix = A(k,1):A(k,2);
%     plot(score(ix,1),ix,'Color',[1,0,1],'LineWidth',1.5)
% end
% for k = 1:size(SW,1)
%     ix = SW(k,1):SW(k,2);
%     plot(score(ix,1),ix,'Color',[1,.6,0],'LineWidth',1.5)
% end