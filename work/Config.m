classdef Config < handle
    properties
        pathTIMIT = '../audio/';
        pathDB_root = '../databases/';
        pathDB_train = ''; % needs to be set in constructor
        pathDB_test = ''; % needs to be set in constructor
    end
    methods
        function obj = Config()
            obj.pathDB_train = [obj.pathDB_root '2023_11_2_14_13_35\'];
            obj.pathDB_test = [obj.pathDB_root '2023_11_2_14_13_35\'];
        end
    end
end