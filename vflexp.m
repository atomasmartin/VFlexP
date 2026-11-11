function appOut = vflexp(casePath)
%VFLEXP Launch the VFlexP Power System Simulation & Analysis MATLAB App.
%
%   vflexp()
%   vflexp('bus_systems/4buses_hybrid.xlsx')
%   app = vflexp(...)
%
%   Opens the programmatic 6-tab VFlexP graphical interface backed by
%   VFlexPSession and the modular Simulink library architecture.

repoRoot = fileparts(mfilename('fullpath'));
if exist('VFlexPSession', 'class') ~= 8 || exist('VFlexPApp', 'class') ~= 8
    install();
end

if nargin < 1
    casePath = '';
end

% Reuse existing open VFlexPApp window if valid and no new case was requested,
% or load the requested case into it.
existingFig = findall(groot, 'Type', 'figure', 'Tag', 'VFlexPMainFigure');
if ~isempty(existingFig) && isvalid(existingFig(1))
    existingApp = getappdata(existingFig(1), 'VFlexPAppInstance');
    if ~isempty(existingApp) && isvalid(existingApp)
        figure(existingFig(1));
        if ~isempty(casePath)
            existingApp.loadCase(casePath);
        end
        if nargout > 0
            appOut = existingApp;
        end
        return;
    end
end

if isempty(casePath)
    defaultCase = fullfile(repoRoot, 'bus_systems', 'VSM_and_SG.xlsx');
    if isfile(defaultCase)
        casePath = defaultCase;
    end
end

app = VFlexPApp(casePath);
if nargout > 0
    appOut = app;
end
end
