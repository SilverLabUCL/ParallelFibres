
dataset_ix = 2;

[~,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);


dFF = nan(200,length(whisk_set_point));

bin_centers = linspace(-20,30,100);

for n = 1:100
    dFF(n,:) = normpdf(whisk_set_point,bin_centers(n),10);
end


 for n = 101:150
     dFF(n,:) = whisk_set_point.*(whisk_set_point>median(whisk_set_point));
 end
 
 for n = 151:200
     dFF(n,:) = whisk_set_point.*(whisk_set_point<median(whisk_set_point));
 end

dFF = zscore(dFF')';

dFF = dFF + .1 * randn(size(dFF));
%%

dataset_ix = 2;

[~,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);


dFF = nan(200,length(whisk_set_point));

bin_centers = linspace(-20,30,100);

for n = 1:100
    dFF(n,:) = normpdf(whisk_set_point,bin_centers(n),10);
end


 for n = 101:150
     dFF(n,:) = whisk_set_point.*(whisk_set_point>median(whisk_set_point));
 end
 
 for n = 151:200
     dFF(n,:) = whisk_set_point.*(whisk_set_point<median(whisk_set_point));
 end

dFF = zscore(dFF')';

dFF = dFF + .1 * randn(size(dFF));

%%
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[N,T] = size(dFF);

dwsp = 5;
wsp_bin = -20:dwsp:60;
tuningcurves = nan(N,size(wsp_bin,2));
for k = 1:size(wsp_bin,2)
    ix = find(whisk_set_point > wsp_bin(k) - dwsp/2 & whisk_set_point < wsp_bin(k) + dwsp/2);
    tuningcurves(:,k) = mean(dFF(:,ix),2);
end


dFF = nan(N,T);
for n = 1:N
    dFF(n,:) = interp1(wsp_bin,tuningcurves(n,:),whisk_set_point');
end

ix = find(isnan(sum(dFF,1)))
dFF(:,ix) = [];
whisk_amp(ix) = [];
speed(ix) = [];
whisk_set_point(ix) = [];
%%



[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
cd = get_coding_dimension(dFF,A,QW);
A_or_QW = cd'*(dFF - nanmean(dFF,2));

[coeff, score] = pca(dFF');

zsp = zscore(A_or_QW);

figure, hold on
for k = 1:length(score)-1
    c = (zsp(k)+1)/2;
    if c < 0
        c = 0;
    elseif c > 1
        c = 1;
    end
    c=c*[1,0,1]+(1-c)*[0,1,1];
    
    plot3(score(k:k+1,1),score(k:k+1,2),score(k:k+1,3),'-','Color',c)
end
view(44,44)
set(gca,'FontSize',15,'XTick',[],'YTick',[],'ZTick',[])
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

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

subspace(coeff_QW(:,1:num_PCs),coeff_A(:,1:num_PCs))