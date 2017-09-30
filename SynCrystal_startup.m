function SynCrystal_startup()
%  Copyright (c) 2016 Haizhao Yang, Duke University 
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
tmp = strfind(file_path,'SynCrystal_startup');
file_path = file_path(1:(tmp(end)-1));

% Foulder for all soource files recursively
addpath(genpath([file_path 'external']));
addpath(genpath([file_path 'data']));
addpath(genpath([file_path 'results']));
addpath(genpath([file_path 'SSTmethod']));
addpath(genpath([file_path 'VarSSTmethod']));
addpath(genpath([file_path 'Sketching']));

disp('Begin to compile MEX files...');
rootDir = pwd;
cd(['external' CSPT 'SynLab' CSPT 'Source' CSPT 'SS_CT_2D' CSPT 'src' CSPT]);
mex SS_polar_v2.c;
mex SS_polar_v1.c;
mex SS_polar.c;
mex skeletonPolarMex.c;
cd(rootDir);
cd(['SSTmethod' CSPT 'src' CSPT]);
mex LocSmooth.c;
cd(rootDir);
cd(['VarSSTmethod' CSPT 'src' CSPT 'srcSST' CSPT]);
mex LocWeight.c;
mex LocWeight_v2.c;
cd(rootDir);

disp('Path set!');

clear tempPath front back
clear SetPath MATLABVERSION CSPT
clear type MSPT

