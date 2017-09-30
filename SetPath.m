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
tmp = strfind(file_path,'SetPath');
file_path = file_path(1:(tmp(end)-1));

% Foulder for all soource files recursively
addpath(genpath([file_path 'data']));
addpath(genpath([file_path 'results']));
addpath(genpath([file_path 'Sketching']));
addpath(genpath([file_path 'SSTmethod']));
addpath(genpath([file_path 'external']));
addpath(genpath([file_path 'VarSSTmethod']));

disp('Begin to compile MEX files...');
rootDir = pwd;
cd(['SynLab' CSPT 'Source' CSPT 'SS_CT_2D' CSPT 'src' CSPT]);
mex SS_polar.c;
mex SS_polar_v2.c;
mex SS_polar_v1.c;
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

