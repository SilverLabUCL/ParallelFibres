
clear all; clc

define_dirs;

for dataset_ix = 1:13
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,loco,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);
    
    onset_indices = get_onsets(speed,acquisition_rate);
    
    if ~isempty(onset_indices)

        % Get A and QW states
        [A,QW] = define_behav_periods(whisk_amp,loco,acquisition_rate); 
        [change_dFF,p_val] = change_dFF_sig(dFF,A,QW,acquisition_rate);

        ix_OFF = find(change_dFF < 0 & p_val < .05);

        dFF_NM = mean(dFF(ix_OFF,:),1);

        N_onsets = size(onset_indices,1);
        onset_bins = onset_indices(1,2)-onset_indices(1,1)+1;
        time_onsets = nan(N_onsets,onset_bins);
        dFF_NM_onsets = nan(N_onsets,onset_bins);

        for k = 1:N_onsets
            ix = onset_indices(k,1):onset_indices(k,2);
            time_onsets(k,:) = time(ix);
            dFF_NM_onsets(k,:) = dFF_NM(ix);
        end

        save([basedir,datasets{dataset_ix},'_avg_NM'],'time','dFF_NM','time_onsets','dFF_NM_onsets')
    end
end
