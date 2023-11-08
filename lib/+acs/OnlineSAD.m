% Class is instantiated with a vector containing SAD
% information for frames of the relevant signal(s)
% 
% Within CL-DXPP function, process() is called to return the SAD flag
% for the current frame. This method is convinient due to the similarity
% in implementation compared to OnlineIACS.

classdef OnlineSAD < handle
    properties
        % Parameters
        SAD_all;
        SAD = 0;
    end
    methods
        function obj = OnlineSAD(SAD_all)
            % Boolearn SAD vector for all frames
            obj.SAD_all = SAD_all;
        end
        function obj = process(obj, ell)
           obj.SAD = obj.SAD_all(ell);
        end
    end
end