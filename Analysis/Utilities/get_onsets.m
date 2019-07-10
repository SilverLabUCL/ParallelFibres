


function [offset_indices, onset_indices] = get_onsets(whisk_amp,acquisition_rate)

    x = smoothdata(whisk_amp,'gaussian',[round(1 * acquisition_rate),0]);
    dx = abs(diff(x));

    % 
    ix = find(dx > .01);

    % Remove first 500 ms
    num_bins = round(acquisition_rate * 2);
    ix = ix(ix>num_bins);
    
    % Remove last one if there isnt 500 ms after onset
    ix = ix(ix <= numel(x) -num_bins);

    % Only keep indices that had dx < threshold for the 500 ms prior
    ix_ = [];
    for i = 1:length(ix)
        if dx((ix(i) - num_bins):(ix(i)-1)) < .01 & x((ix(i) - num_bins):(ix(i)-1)) < .03
            ix_ = [ix_,ix(i)];
        end
    end
    ix = ix_;

    % Remove indices with high value of whisker amplitude
    ix = ix(x(ix)<.03);

    % Refine timing using second derivative
    ddx = diff(dx);
    for i = 1:length(ix)
        t_ = (ix(i)-round(.2*acquisition_rate)):ix(i);
        t_ = t_(ddx(t_)>.001);
        ix(i)=t_(1);
    end
    
    % Only keep indices that have high whisker amplitude for 2 s
    ix_ = [];
    for i = 1:length(ix)
        if mean(x((ix(i)+1):(ix(i) + num_bins))) >= x((ix(i)))
            ix_ = [ix_,ix(i)];
        end
    end
    ix = ix_;
    onset_indices = ix-1;

    % For each onset, find the preceeding offset
    offset_indices = zeros(size(onset_indices));
    for i = 1:length(offset_indices)
        if i == 1
            ix = 1:(onset_indices(i));
        else
            ix = onset_indices(i-1):(onset_indices(i));
        end
        ix = ix(x(ix)<.03);
        ix_ = find(diff(ix)~=1);
        if isempty(ix_)
            offset_indices(i) = ix(1);
        else
            offset_indices(i) = ix(ix_(end)+1);
        end
    end
    
    i = 2;
    while i <= length(offset_indices)
        if offset_indices(i) == onset_indices(i-1)
            onset_indices(i-1) = onset_indices(i);
            onset_indices(i) = []; offset_indices(i) = [];
        end
        i = i+1;
    end
    