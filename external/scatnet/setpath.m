function setpath()
%   MAKE adds paths of the DeComp to Matlab.

%  Copyright (c) 2017 Haizhao Yang, National University of Singapore
%  This file is distributed under the terms of the MIT License.

global SetPath
global CSPT
global MSPT

type = computer;

if strcmp(type,'MAC2'),
    CSPT = ':';
    SetPath = [pwd, CSPT];
    MSPT = ';';
elseif isunix,
    % Mac OS X returns isunix=1
    CSPT = '/';
    SetPath = [pwd, CSPT];
    MSPT = ':';
elseif strcmp(type(1:2),'PC');
    CSPT = '\';
    SetPath = [pwd, CSPT];
    MSPT = ';';
end

disp('Begin to set MATLAB path...')

file_path = mfilename('fullpath');
tmp = strfind(file_path,'setpath');
file_path = file_path(1:(tmp(end)-1));

% Foulder for all soource files recursively
addpath(genpath([file_path 'scatnet']));
addpath(genpath([file_path 'compare']));

disp('Path set!');

clear tempPath front back
clear SetPath MATLABVERSION CSPT
clear type MSPT
end