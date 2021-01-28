%% Calculate fibre direction from z stack of images
% Requirements:
%         Zstack Images/Axon_Results.csv
% i.e., csv files generated from Image J of estimates of axon fibres

% File > Import > Image Sequence
% Image > Adjust > Brightness/Contrast
% Process > Smooth
% Multi-point tool
% Analyze > Measure

% generated using multipoint tool. Generates the following variables:
%    vector_mean     2x1 vector of average fibre direction for that patch
%    angle_std       std of angles around vector_mean
% These values are saved in Processed/P_angle_given_S.mat

clear all; clc

basedir = '~/Documents/ParallelFibres/Data/';
datasets = {'FL87_180501_11_03_09',...  1
            'FL87_180501_10_47_25',...  2
            'FL87_180501_10_36_14',...  3 
            'FL87_180220_10_38_55',...  4
            'FL77_180213_10_46_41',...  5
            'FL_S_170906_11_26_25',...  6
            'FL_S_170905_10_40_52',...  7            
            'FL45_170125_14_47_04',...  8
             ...%%
             'FL95_180425_10_53_40',...  9
             'FL87_180413_11_00_55',...  10
             'FL87_180117_11_23_20',...  11
             'FL_S_171109_14_54_34',...  12
             'FL_S_171109_15_19_52',...  13
             'FL76_170913_10_57_06',... 14
             'FL77_180113_10_58_50'};  % 15
         
% Choose dataset
dataset_ix = 1;

% Choose patch number - just to plot Cn for reference
patch_no = 1;

fname = datasets{dataset_ix};
disp(fname)

load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'],'Y','d1','d2')
Cn = correlation_image(Y, [1,2], d1,d2); clear Y

%% Load coordinates (New formatting for saving coordinates)

coords = csvread([basedir,fname,'/Zstack Images/Axon_Results.csv'],1,1);
coords = coords(:,5:6);

init_coords = coords(1:2:end,:);
end_coords = coords(2:2:end,:);

%% Plot axon segments in xy plane for references

figure, hold on
vectors_norm = zeros(size(end_coords));
for k = 1 : size(end_coords,1)
    temp = [init_coords(k,:);end_coords(k,:)];
    plot(temp(:,1),temp(:,2),'k');
    n = diff(temp,[],1); n = n/norm(n);
    vectors_norm(k,:) = n;
end
set(gca, 'Ydir','reverse', 'FontSize',15) 
axis tight; axis equal;
set(gca,'XTick',[],'YTick',[]);
set(gca,'Box','on')

%% Correct sign of vector in case one is accidentally going in wrong direction

% Approximate average vector
vector_mean = median(vectors_norm);
vector_mean = vector_mean/norm(vector_mean);

for k = 1 : size(end_coords,1)
    if norm(vectors_norm(k,:)-vector_mean) > sqrt(2)
        vectors_norm(k,:) = -vectors_norm(k,:);
    end
end

%% Calculate angles of vectors

% Average vector after correction
vector_mean = mean(vectors_norm);
vector_mean = vector_mean/norm(vector_mean);

% Calculate angle from average vector
angle = zeros(size(vectors_norm,1),1);
for k = 1:size(vectors_norm,1)
    angle(k) = sign(vectors_norm(k,1) - vector_mean(1)) * subspace(vectors_norm(k,:)',vector_mean');
    % Below; less stable method
    % sign(vectors_norm(k,1) - vector_mean(1)) * acos(vectors_norm(k,:)*vector_mean');
end

%% Plot all vectors and average vector on top of correlation matrix

figure('papersize',[d2,d1]/40); 
imagesc(Cn); colormap(gray), hold on,
set(gca,'XTick',[],'YTick',[]);

o = [d2/2,d1];

for k = 1:size(vectors_norm,1)    
    temp = [o; o+30*vectors_norm(k,:)];
    plot(temp(:,1),temp(:,2),'-w','LineWidth',.5) 
end
temp = [o; o+30*vector_mean];
plot(temp(:,1),temp(:,2),'-m','LineWidth',3)

%% Plot histogram of angles with normal on top

angle_std = std(rad2deg(angle));

figure, histogram(rad2deg(angle),-180:10:180)
set(gca,'FontSize',15)
xlabel('Angle (deg)'), ylabel('Number')
xlim([-180,180])

x = -180:1:180;
hold on, plot(x,normpdf(x,0,angle_std)*numel(angle)*10,'b','LineWidth',2)
title(angle_std)
%% Save data
fmatname = strcat(basedir,fname,'/processed/fibre_direction.mat');

x = input('Enter Y if you want to save.   ','s');
if strcmp(x,'y') || strcmp(x,'Y')
    save(fmatname,'vector_mean','angle_std','angle')
end

%% Collate results - Do this only after all datasets are finished

angle_all = [];

for dataset_ix = 1:16
    
    fname = datasets{dataset_ix};
    disp(fname)

    load([basedir,fname,'/processed/fibre_direction.mat'],'angle','angle_std')
    disp(angle_std)
    
    angle_all = [angle_all; angle];

end

std(rad2deg(angle_all))
