%
% Function to plot video of fluorescence (not actually dFF, should change
% the name at some point).
%
%
% Input:
%    Ain              Spatial filters matrix (num pixels x num ROIs)
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    dims             Vector of pixel dims of patch: [d1, d2]
%    acq_rate         Acquisition rate (Hz)
%    pathtosave       Path and filename if you want to save an avi
%                     (optional - otherwise does not save video)
%   

function make_video_dFF(Ain,Y,dims,acq_rate,pathtosave)

    d1 = dims(1);
    d2 = dims(2);
    T = size(Y,2);

    if nargin < 5 || isempty(pathtosave)
        save_video = 0;
    else
        save_video = 1;
    end

    %% Smooth in time (200ms) and in space (1 pixel)
    smooth_win = round(.2 * acq_rate);

    psf = fspecial('gaussian',5, 1);

    YY = Y;
    for d = 1:d1*d2
        YY(d,:) = smoothdata(Y(d,:),'Gaussian',smooth_win);
    end

    YY = imfilter(reshape(YY, d1,d2,[]), psf, 'replicate');
    YY = reshape(YY,d1*d2,T);

    %% Draw boundaries of spatial filters and prepare plot
    
    figure, hold on
    h=imagesc(reshape(YY(:,1),d1,d2));
    for kk = 1:size(Ain,2)
       B = bwboundaries(reshape(Ain(:,kk),d1,d2));
       if numel(B)>1
           error
       end
       boundary = B{1};
       plot(boundary(:,2), boundary(:,1),'w', 'LineWidth', 1.5)
    end
    axis([0,d2,0,d1]), axis equal, axis tight
    set(gca,'XTick',{}),set(gca,'YTick',{})
    set(gca,'FontSize',20)
    set(gca,'YDir','Reverse')

    min_ = 100;% prctile(YY(:),10);
    max_ = 500;%prctile(YY(:),99.99);
    caxis([min_,max_])
    
    %% Plot each frame
    for t=1:T
        h.CData = reshape(YY(:,t),d1,d2); caxis([min_,max_])
        title(['t=',num2str(t/acq_rate,'%.2f'),'s frame ',num2str(t)])
        F(t) = getframe(gcf);
        drawnow
    end

    %% Save video
    if save_video
        % create the video writer 
        writerObj = VideoWriter([pathtosave,'.avi']);
        writerObj.FrameRate = (acq_rate) * 5;
        writerObj.Quality=20;
        % set the seconds per image
        % open the video writer
        open(writerObj);
        % write the frames to the video
        for i=1:length(F)
            % convert the image to a frame
            frame = F(i) ;    
            writeVideo(writerObj, frame);
        end
        % close the writer object
        close(writerObj);
    end


