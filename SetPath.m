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

back = CSPT;
tempPath = path;
front = [MSPT SetPath];
tempPath = [tempPath front];

% Source folder
tempPath = [tempPath front 'data' back];

tempPath = [tempPath front 'results' back];

tempPath = [tempPath front 'SynLab' back];

tempPath = [tempPath front 'SynLab' back 'Source' back];

tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_CT_2D' back];
tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_CT_2D' back 'demo' back];
tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_CT_2D' back 'src' back];

tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_WP_1D' back];
tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_WP_1D' back 'demo' back];
tempPath = [tempPath front 'SynLab' back 'Source' back 'SS_WP_1D' back 'src' back];

% Applications folder
tempPath = [tempPath front 'SSTmethod' back];
tempPath = [tempPath front 'SSTmethod' back 'demo' back];
tempPath = [tempPath front 'SSTmethod' back 'src' back];
tempPath = [tempPath front 'VarSSTmethod' back];
tempPath = [tempPath front 'VarSSTmethod' back 'demo' back];
tempPath = [tempPath front 'VarSSTmethod' back 'src' back];
tempPath = [tempPath front 'VarSSTmethod' back 'src' back 'srcOptPart' back];
tempPath = [tempPath front 'VarSSTmethod' back 'src' back 'srcOther' back];
tempPath = [tempPath front 'VarSSTmethod' back 'src' back 'srcSST' back];

path(tempPath);

disp('Begin to compile MEX files...');
rootDir = pwd;
cd(['SynLab' CSPT 'Source' CSPT 'SS_CT_2D' CSPT 'src' CSPT]);
mex SS_polar.c;
mex SS_polar_v2.c;
mex SS_polar_v1.c;
cd(rootDir);
cd(['SSTmethod' CSPT 'src' CSPT]);
mex LocWeight.c;
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

