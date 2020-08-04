clear all, clc, close all


define_dirs
%%

dataset_ix = 2;
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[coeff, score] = pca(dFF');

[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
cd = get_coding_dimension(dFF,A,QW);
A_or_QW = cd'*(dFF - mean(dFF,2));
A_or_QW = zscore(A_or_QW);


%% Load eye cam frames

obj = VideoReader([basedir,'180501_10_49_21_VidRec/EyeCam-1.avi']);
all_frames = zeros(480,640,3,1250,'uint8');
for k = 1 : 1250  %fill in the appropriate number
  this_frame = readFrame(obj);
  all_frames(:,:,:,k) = this_frame(1:2:end,1:2:end,:);
end
%% Load eye cam times

TimeEyeCam = dlmread([basedir,'180501_10_49_21_VidRec/EyeCam-relative times.txt']);
TimeEyeCam = TimeEyeCam(:,2);
TimeEyeCam = TimeEyeCam/1000;

%% Resample everything to eye cam rate
pc1 = -score(:,1); % for dataset 5, anticorr

pc1_EyeCamFreq = interp1(time,pc1,TimeEyeCam);
A_or_QW_EyeCamFreq = interp1(time,A_or_QW,TimeEyeCam);
score1_EyeCamFreq = interp1(time,score(:,1),TimeEyeCam);
score2_EyeCamFreq = interp1(time,score(:,2),TimeEyeCam);
score3_EyeCamFreq = interp1(time,score(:,3),TimeEyeCam);

%%


h = plot3(score(:,1),score(:,2),score(:,3),'Color',[.7,.7,.7]);
pts = [h.XData; h.YData; h.ZData; ones(size(h.XData))];

mat = viewmtx(8,-38);% * makehgtform('scale',1./[20 1 20]);
pt2 = mat * pts;
pc_proj_x = pt2(1,:) ./ pt2(4,:);
pc_proj_y = pt2(2,:) ./ pt2(4,:);
figure, plot(pc_proj_x,pc_proj_y)

%% Resample everything to eye cam rate

pc_proj_x_EyeCamFreq = interp1(time,pc_proj_x,TimeEyeCam);
pc_proj_y_EyeCamFreq = interp1(time,pc_proj_y,TimeEyeCam);


A_or_QW_EyeCamFreq(A_or_QW_EyeCamFreq>1) = 1; 
A_or_QW_EyeCamFreq(A_or_QW_EyeCamFreq<0) = 0;
%% 2D projection version

k_start = 150;

v = VideoWriter([basedir,'behaviour_pca_corr_colour.avi']);
v.Quality = 20;
v.FrameRate = 24;
open(v);

f = figure('Position',[0,0,800,400]);
set(gcf,'color','w');

subtightplot(3,6,[4:6,10:12,16:18])

set(gca,'FontSize',15)
set(gca,'XTick',{})
set(gca,'YTick',{})
set(gca,'ZTick',{})
hold on,
s = scatter(pc_proj_x_EyeCamFreq(k_start),pc_proj_y_EyeCamFreq(k_start),50,A_or_QW_EyeCamFreq(k_start),'filled','MarkerEdgeColor','k'); caxis([0,1]); colormap('cool')
axis([-1.5,1.5,-1.5,1])
cb = colorbar('Location','northoutside','Position',[0.8 0.8349 0.1 0.0400],'Ticks',[0,1],'TickLabels',{'QW','AS'});

subtightplot(3,6,[1:3,7:9,13:15])
image(all_frames(:,:,:,1))
set(gca,'XTick',{})
set(gca,'YTick',{})


%%
for k = 250:size(TimeEyeCam,1)-1
    try   
        subtightplot(3,6,[1:3,7:9,13:15])
        image(all_frames(:,:,:,k)) 
        set(gca,'XTick',{})
        set(gca,'YTick',{})

        col = A_or_QW_EyeCamFreq(k)*[1,0,1]+(1-A_or_QW_EyeCamFreq(k))*[0,1,1];

        subtightplot(3,6,[4:6,10:12,16:18]), axis([-5,3,-3,4]), hold on
        plot(pc_proj_x_EyeCamFreq(k:k+1),pc_proj_y_EyeCamFreq(k:k+1),'Color',col,'LineWidth',1)
        s.XData = pc_proj_x_EyeCamFreq(k+1);
        s.YData = pc_proj_y_EyeCamFreq(k+1);
        s.CData = col;
        uistack(s,'top');

        if mod(k,3)==0
        drawnow;
        frame = getframe(f);
        writeVideo(v,frame);
        end
    end
end

close(v)
disp('Finished writing video')
