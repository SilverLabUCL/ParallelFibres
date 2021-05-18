% Makes video of rotating manifolds to show orthogonality
% Supplementary video S3

clear all; clc; close all
define_dirs;
dataset_ix=1;
%%
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
cd = get_coding_dimension(dFF,A,QW);
A_or_QW = cd'*(dFF - nanmean(dFF,2));

[coeff, score] = pca(dFF');

c_A_QW = [];
score_A_QW = [];
for k = 1:length(A)
    ix = A(k,1):A(k,2);
    score_A_QW = [score_A_QW; score(ix,1:3)];
    c_A_QW = [c_A_QW; repmat([1,0,1],[length(ix),1])];
end


for k = 1:length(QW)
    ix = QW(k,1):QW(k,2);
    score_A_QW = [score_A_QW; score(ix,1:3)];
    c_A_QW = [c_A_QW; repmat([0,1,1],[length(ix),1])];
end

%%
f = figure('Position',[0,0,400,400]);
axis tight;
   axis equal;
   axis vis3d
h = scatter3(score_A_QW(:,1),score_A_QW(:,2),score_A_QW(:,3),10,c_A_QW,'^','filled');
caxis([.05,.2]); colormap('cool')
  
  
   v = VideoWriter([basedir,'manifold_rotate.avi']);
   v.Quality = 20;
   v.FrameRate = 24;%n_cycles/duration;
   open(v);

   f = figure('Position',[0,0,400,400]);
   h = scatter3(score_A_QW(:,1),score_A_QW(:,2),score_A_QW(:,3),10,c_A_QW,'^','filled');
   axis tight;
   axis equal;
   axis vis3d
   set(gca,'nextplot','replacechildren');
   set(gca,'xtick',[]);
   set(gca,'xticklabel',[]);
   xlabel('PC 1');
   set(gca,'ytick',[]);
   set(gca,'yticklabel',[]);
   ylabel('PC 2');
   set(gca,'ztick',[]);
   set(gca,'zticklabel',[]);
   zlabel('PC 3');
   set(gcf,'closer',' ');
   set(gcf,'color','w');
   
   %timepoint = 50;
   %set(h,'CameraViewAngleMode','Manual')


   % Create video and save raw data
   caxis([.05,.2]); colormap('cool')
   for t = 60:270
       try
           view(t,15);axis([-3, 5, -3, 4, -2, 3]); 
           if mod(t,2) == 0
               drawnow;
               frame = getframe(f);
               writeVideo(v,frame);
           end
       %catch
       %    error_box('video not exported, but data backed up',1)
       end
   end

   % Make sure the file is always closed
  
       delete(f)
       close(v);


  