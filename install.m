% Get the full path of the current file
fullPath = mfilename('fullpath');

% Extract only the directory path
filePath = fileparts(fullPath);

addpath([filePath '\docs'])
addpath([filePath '\functions'])
addpath([filePath '\MATPOWER-FUBM\lib'])
addpath([filePath '\MATPOWER-FUBM\mips\lib'])
addpath([filePath '\MATPOWER-FUBM\most\lib'])
addpath([filePath '\MATPOWER-FUBM\mp-opt-model\lib'])
savepath