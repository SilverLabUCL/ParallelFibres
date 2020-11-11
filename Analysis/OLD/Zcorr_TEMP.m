
define_dirs;

fname = 'Z_FL87_180419_13_25_45';
disp(fname)
    
load([basedir,fname,'/processed/',fname,'_GroupedData.mat'])
load([basedir,fname,'/wheel_MI.mat'])
load([basedir,fname,'/',fname,'.mat'])
   
dFF = vertcat(dFF_axons{:});
time = mean(vertcat(time_axons{:}))'/1000;

smooth_win_s = 0.2;

% Correct double frames by taking average - wheel
ix_to_correct = find(diff(wheel_MI(:,2))==0);
for ix = ix_to_correct
    wheel_MI(ix,1) = (wheel_MI(ix,1)+wheel_MI(ix+1,1))/2;
    wheel_MI(ix+1,:) = [];
end

% Smooth speed
Speed_time = wheel_MI(:,2) / 1000;
dt_speed = mean(diff(Speed_time));
speed = smoothdata(wheel_MI(:,1),'gaussian',[round(smooth_win_s/dt_speed) 0] *2);

% Convert to units of standard deviations
speed = zscore(speed);

speed = interp1(Speed_time,speed,time);

%%
speed = speed(30:end);
dFF = dFF(:,30:end);
time = time(30:end);

%%

[~,~,~,~,C,p_val,~,~] ...
        = corr_sig(dFF,speed,speed,speed,[],acquisition_rate);
%%

bins = linspace(-1,1,50);
bins_c = bins(2:end)-mean(diff(bins))/2;

C_pass_shuff = C(p_val < 0.05);
C_fail_shuff = C(p_val >= 0.05);

C_up = C_pass_shuff(C_pass_shuff>0);
C_down = C_pass_shuff(C_pass_shuff<0);

h_up = histcounts(C_up,bins);
h_down = histcounts(C_down,bins);
h_fail = histcounts(C_fail_shuff,bins);

figure, bar(bins_c,h_fail,'FaceColor',[.7,.7,.7],'EdgeColor','k','BarWidth',1)
hold on, bar(bins_c,h_up,'FaceColor',[1,.3,1],'EdgeColor','k','BarWidth',1)
hold on, bar(bins_c,h_down,'FaceColor',[.2,1,1],'EdgeColor','k','BarWidth',1)
set(gca,'FontSize',18)
set(gca,'Box','off')
xlabel('Correlation')
ylabel('Number')

fracs = [sum(h_up), sum(h_down),sum(h_fail)];
fracs = fracs/sum(fracs);
figure, p = pie(fracs);
set(gca,'FontSize',18)
t = p(1); t.FaceColor = [1,.3,1]; t.EdgeColor = [1,.3,1];
t = p(2); t.FontSize=20;
t = p(3); t.FaceColor = [.2,1,1]; t.EdgeColor = [.2,1,1];
t = p(4); t.FontSize=20;
t = p(5); t.FaceColor = [.7,.7,.7]; t.EdgeColor = [.7,.7,.7];
t = p(6); t.FontSize=20;

%%
N = size(dFF,1);

ix_up = find(C > 0 & p_val < 0.05);
ix_down = find(C < 0 & p_val < 0.05);
ix_fail = find(p_val >= 0.05);
%%

figure, subplot(8,1,1), plot(time, speed,'k','LineWidth',1.5)
set(gca,'XTick',[],'Box','off')
%xlabel('Time (s)'), ylabel('Loc. (std)')
xlim([2,65])

subplot(8,1,2:8)
ix_down_show = [1, 21, 59, 95, 96, 141];
 hold on
for k = 1:numel(ix_down_show)
    plot(time, dFF(ix_down_show(k),:)+k*5,'LineWidth',1.5)
end
xlim([2,65])
xlabel('Time (s)')
%%

figure, subplot(8,1,1), plot(time, speed,'k','LineWidth',1.5)
set(gca,'XTick',[],'Box','off')
%xlabel('Time (s)'), ylabel('Loc. (std)')
xlim([2,65])

subplot(8,1,2:8)
ix_up_show = [4, 62, 75, 82, 125, 129];
hold on 
for k = 1:numel(ix_up_show)
    plot(time, dFF(ix_up_show(k),:)+k*5,'LineWidth',1.5)
end
xlim([2,65])
xlabel('Time (s)')


