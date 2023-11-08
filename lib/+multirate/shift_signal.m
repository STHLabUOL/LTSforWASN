function y_shifted = shift_signal(y, smp)
    % shift 1D array, insert zeros
    
    if size(y, 1) == 1
        col_shift = 1;
    elseif size(y, 2) == 1
        col_shift = 0;
    else
        error('Cannot shift multidim array');
    end
     
    if col_shift
        y = y';
    end

    if smp > 0            
        y_shifted = [zeros(smp, 1); y(1:end-smp)];
    else
        smp = -smp;
        y_shifted = [y(smp+1:end); zeros(smp, 1)];
    end
    
    if col_shift
        y_shifted = y_shifted';
    end

end