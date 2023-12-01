% Do not move or rename this file!
%
% This file is executed automatically, when you launch matlab from this directory.
% It will make sure that all required folders belong to your matlab path.

% Add required libraries to path
path(path, [pwd filesep 'lib']);

% Set Graphics Default Parameters
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultLegendFontSize', 12);
set(groot, 'DefaultLegendFontSizeMode', 'manual')