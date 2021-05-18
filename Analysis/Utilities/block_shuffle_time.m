%
% This function randomly shuffles timepoints in blocks
% Input:
%    T                    Number of time steps
%    acquisition_rate     Acquisition rate
%    block_s              Time of each block in s
% 
% Output:
%    T_shuffled           Vector of block-shuffled time points      

function T_shuffled = block_shuffle_time(T,acquisition_rate,block_s)

    if nargin < 3 || isempty(block_s)
        block_s = 1; % default 1 s blocks
    end

    % Length of blocks in bins
    block_length = round(block_s *acquisition_rate);
    num_blocks = ceil(T/block_length);
    
    % Random reordering of blocks
    block_order = randsample(num_blocks,num_blocks);
    
    % Vector of shuffled time indices
    T_shuffled = zeros(1,T);
    count = 1;
    for bl = 1:num_blocks
        if block_order(bl) == num_blocks
            t_shuff = (num_blocks-1)*block_length+1 : T;
        else
            t_shuff = (block_order(bl)-1)*block_length+1 : block_order(bl)*block_length ;
        end
        t = count : count + length(t_shuff)-1;
        T_shuffled(t) = t_shuff;
        count = count + length(t_shuff);
    end