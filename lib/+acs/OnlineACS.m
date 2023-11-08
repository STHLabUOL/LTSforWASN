% This implements the ACS algorithm as a class to enable
% elegant online ACS estimation in parallel to the main DXCPP algorithm.
%
% Outputs only the BINARY IACS

classdef OnlineACS < handle
    properties
        % Parameters
        alpha = 0.95; % smooth. const. for derivate
        dThresh = 9e-5; % max. increase rate 9e-5 is 70th percentile of all dchi
        mapping = "EXP";
        outputType = "bin";
        overwriteStart = 0; % num frames that initially are overwritten with "1" regardless of GCSD2 (to ensure CL-DXCPP convergence)
        % State:
        ell = 0; % frame count
        preACS = 0;
        preACS_unrestr = 0; % unrestricted
        dACS = 0;
        ACS = 0;
    end
    methods
        function obj = OnlineACS(mapping, outputType, overwriteStart)
            % set params
            obj.mapping = mapping;
            obj.outputType = outputType;
            obj.overwriteStart = overwriteStart;
        end
        function obj = process(obj, GCSD2_avg)
            % check for valid (non-zero) input
            if ~all(GCSD2_avg, 'all')
                return;
            end
            % Input Time averaged GCSD2, as computed within DXCPP
            preACS_unrestr_before = obj.preACS_unrestr;
            % calc preliminary iacs
            obj.preACS_unrestr = mean(abs(GCSD2_avg), 1);
            % control via derivative
            obj.dACS = obj.alpha*obj.dACS + (1-obj.alpha)*(obj.preACS_unrestr - preACS_unrestr_before);
            if obj.dACS <= obj.dThresh
                obj.preACS = obj.preACS_unrestr;
            end
            % calc ACS
            obj.ACS = acs.acs__mapping(obj.preACS, obj.mapping, obj.outputType);   
            % catch NaN during init
            if isnan(obj.ACS)
                obj.ACS = 0;
                warning('Encountered NaN during ACS calculation. Overwriting with zero.');
            end
            % overwrite?
            if obj.ell < obj.overwriteStart
                obj.ACS = 1;
            end
            % count
            obj.ell = obj.ell+1;
        end
    end
end