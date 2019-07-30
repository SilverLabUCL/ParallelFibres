


arousal_signal = pc1.*(pc1>0);

aroused = [];

isDone = 0;
aroused_end = 0;
while isDone == 0

    % Beginning of whisking period defined by rapid whisker movement
    aroused_start = find(time>=aroused_end & arousal_signal>0);

    if numel(aroused_start) == 0
        isDone = 1; % no whisking periods found      
    else
        aroused_start = time(aroused_start(1));
        
        aroused_end = find(time>aroused_start & arousal_signal==0);
        
        if numel(aroused_end) == 0
            aroused = [aroused;aroused_start,time(end)];
            isDone = 1;
        else
            aroused_end = time(aroused_end(1));
            aroused = [aroused;aroused_start,aroused_end];
        end
    end
end

quiescent = [];
if numel(aroused)==0
    quiescent = [time(1),time(end)];
else            
    if aroused(1,1) > time(1)
        quiescent = [quiescent; time(1), aroused(1,1)];
    end
    for k = 1:size(aroused,1)-1
        quiescent = [quiescent; aroused(k,2), aroused(k+1,1)];
    end 
    if aroused(end,2) < time(end)
        quiescent = [quiescent; aroused(end,2), time(end)];
    end
end


for k = 1:length(aroused)
    if aroused(k,1) <= time(1)
    else
        aroused(k,1) = aroused(k,1) + .300;
        if aroused(k,1) > aroused(k,2)
            aroused(k,:) = NaN;
        end
    end
end
aroused = aroused(~isnan(aroused(:,1)),:);
for k = 1:length(aroused)
    if aroused(k,2) >= time(end)
    else
        aroused(k,2) = aroused(k,2) - 1;
        if aroused(k,1) > aroused(k,2)
            aroused(k,:) = NaN;
        end
    end
end
aroused = aroused(~isnan(aroused(:,1)),:);

for k = 1:length(quiescent)
    if quiescent(k,1) <= quiescent(1)
    else
        quiescent(k,1) = quiescent(k,1) + 1;
        if quiescent(k,1) > quiescent(k,2)
            quiescent(k,:) = NaN;
        end
    end
end
quiescent = quiescent(~isnan(quiescent(:,1)),:);
for k = 1:length(quiescent)
    if quiescent(k,2) >= time(end)
    else
        quiescent(k,2) = quiescent(k,2) - .300;
        if quiescent(k,1) > quiescent(k,2)
            quiescent(k,:) = NaN;
        end
    end
end
quiescent = quiescent(~isnan(quiescent(:,1)),:);

figure, 

a_or_q = quiescent;
for k = 1:length(a_or_q)
    x = a_or_q(k,1);
    y = -1;
    width = a_or_q(k,2) - a_or_q(k,1);
    height = 4;
    hold on, rectangle('Position',[x y width height],'FaceColor',[.85 .85 .85],'EdgeColor','w')
end
plot(time,pc1)
