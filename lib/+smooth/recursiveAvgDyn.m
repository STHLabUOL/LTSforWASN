function [input_smooth] = recursiveAvgDyn(input, alpha, dim)

    % Smooth input vector by recursive averaging with dynamic (time
    % variant) smoothing constant
    %
    % Input:
    %   input: (N x M) input vector
    %   alpha: (N x 1) or (M x 1) smoothing parameter vector
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
            input_smooth(xx,:) = alpha(xx)*before_val + (1-alpha(xx))*input(xx,:);
        elseif dim == 2            
            before_val = input_smooth(:,xx-1);
            input_smooth(:,xx) = alpha(xx)*before_val + (1-alpha(xx))*input(:,xx);
        end
    end
    
end

    
    