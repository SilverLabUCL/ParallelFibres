patch_no = 2
%
disp([num2str(patch_no),' / ',num2str(Numb_patches)])

load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
Y = double(Y);

Y = reshape(Y,d1,d2,[]);
Y(1:10,:,:) = 0; Y(31:40,:,:) = 0;
Y(:,1:10,:) = 0; Y(:,290:299,:) = 0;
Y = reshape(Y,d1*d2,[]);

Cn = correlation_image(Y, [1,2], d1,d2);

% Use CNMF initialization to estimate initial spatial filters
[Ain,~] = detect_ROIs(Y, [d1,d2]);
close all


% Remove low SNR ROIs
[Ain,~,~,~] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,SNR_thresh); 

% Trim ROIs to remove 
Ain = trim_ROIs(Ain,[d1,d2]);

% Remove low SNR ROIs again
[Ain,SNR,A_bad,SNR_bad] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,SNR_thresh); 

% Calculate dFF
dFF = get_dFF(Ain,Y,acquisition_rate,smooth_win_s);

% Plot all dFF
figure, hold on
for k = 1:size(dFF,1)
    plot((1:size(dFF,2))/acquisition_rate,dFF(k,:)+k)
end

% Manually delete ROIs due to slow drift
%[Ain,dFF] = delete_ROIs_manually(Ain,dFF,acquisition_rate);

% Calculate time for rois
time_rois = zeros(size(dFF));
for roi = 1:size(dFF,1)
    Ain_temp = reshape(Ain(:,roi),d1,d2);
    c = regionprops(Ain_temp,'centroid'); c = round(c.Centroid);
    time_rois(roi,:) = get_times(MatrixTime_patch{patch_no}(c(2),c(1)),MatrixTime(end,end),Numb_cycle,Time_between_trial,Numb_trials);
end
%%
N = size(time_rois,1);
ix_axons_to_rois = cell(N,1);
for i  = 1:N
    ix_axons_to_rois{i} = i;
end

plot_grouped_rois(Ain,Cn,dFF,ix_axons_to_rois,acquisition_rate);