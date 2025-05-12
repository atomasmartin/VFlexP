% Get the full path of the current file
fullPath = mfilename('fullpath');

% Extract only the directory path
filePath = fileparts(fullPath);
addpath([filePath '\docs'])
addpath([filePath '\functions'])
savepath