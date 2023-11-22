classdef Config < handle
    properties
        pathTIMIT = '../audio/LibriSpeech/dev-clean/';
        pathNoise = '../audio/';
        pathDB_root = '../databases/';   
        pathDB_train = ''; % needs to be set in constructor
        pathDB_test = ''; % needs to be set in constructor
    end
    methods
        function obj = Config()
            obj.pathDB_train = [obj.pathDB_root '2023_11_15_13_38_48\'];
            obj.pathDB_test = [obj.pathDB_root '2023_11_16_12_19_41\'];
        end
    end
end