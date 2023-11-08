function [input_smooth] = recursiveAvg(input, alpha, dim)

    % Smooth input vector by recursive averaging
    %
    % Input:
    %   input: input vector
    %   alpha: smoothing parameter
    %   dim:   along wich dimension to smooth
    %
    % Output:
    %   input_smooth: smoothed input vector
    %
    
    % Prepare Output
    input_smooth = zeros(size(input));
    
    % Initialize with first value
    if dim == 1
        input_smooth(1,:) = input(1,:);
    elseif dim == 2
        input_smooth(:,1) = input(:,1);
    else
        error('Currently only 2 dimensions are supported.');
    end
    
    % Apply recursive Avg
    end_idx = size(input, dim);
    for xx = 2:end_idx
        if dim == 1
            before_val = input_smooth(xx-1,:);
            input_smooth(xx,:) = alpha*before_val + (1-alpha)*input(xx,:);
        elseif dim == 2
            before_val = input_smooth(:,xx-1);
            input_smooth(:,xx) = alpha*before_val + (1-alpha)*input(:,xx);
        end
    end
    
end

    
    