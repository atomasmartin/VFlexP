function [oppoint,varargout] = fun_modified_findop(model,varargin)
%FINDOP Find operating points from specifications or simulation
%
%   [OP_POINT,OP_REPORT]=FINDOP('model',OP_SPEC) finds an operating point,
%   OP_POINT, of the model, 'model', from specifications given in OP_SPEC.
%
%   [OP_POINT,OP_REPORT]=FINDOP('model',OP_SPEC,OPTIONS) using several options
%   for the optimization are specified in the OPTIONS object, which you can
%   create with the function FINDOPOPTIONS.
%
%   [OP_POINT,OP_REPORT]=FINDOP('model',OP_SPEC, PARAMS, OPTIONS) specifies
%   parameters that you want to vary, PARAMS. PARAMS is a structure with
%   fields:
%      Name: String or MATLAB expression specifying the parameter name
%     Value: Double array containing the different values of the
%            parameter
%   For multiple parameters, specify PARAMS as a structure array.
%
%   The input to findop, OP_SPEC, is an operating point specification object.
%   Create this object with the function OPERSPEC. Specifications on the
%   operating points, such as minimum and maximum values, initial guesses,
%   and known values, are specified by editing OP_SPEC directly or by using
%   get and set. To find equilibrium, or steady-state, operating points, set
%   the SteadyState property of the states and inputs in OP_SPEC to 1. The
%   FINDOP function uses optimization to find operating points that closely
%   meet the specifications in OP_SPEC. By default, findop uses the optimizer
%   fmincon. To use a different optimizer, change the value of OptimizerType
%   in OPTIONS using the FINDOPOPTIONS function.
%
%   A report object, OP_REPORT, gives information on how closely FINDOP
%   meets the specifications. The function FINDOP displays the report
%   automatically, even if the output is suppressed with a semi-colon.
%   To turn off the display of the report, set DisplayReport to 'off' in
%   OPTIONS using the function FINDOPOPTIONS.
%
%   OP_POINT=FINDOP('model',SNAPSHOTTIMES) runs a simulation of the model,
%   'model', and extracts operating points from the simulation at the
%   snapshot times given in the vector, SNAPSHOTTIMES. An operating point
%   object, OP_POINT, is returned.
%
%   See also OPERSPEC, FINDOPOPTIONS

%  Author(s): John Glass
%  Revised:   Yu Jiang
%  Copyright 1986-2016 The MathWorks, Inc.

% If an options object is specified pass this to the operating condition
% spec object


% Parse inputs
narginchk(1,6)

% Error checking for model input and convert it to char array
model = linearize.linutil.validateModelName(model, 'findop');

% First input is the model name. create the model parameter manager and
% make sure models are open
ModelParameterMgr = linearize.ModelLinearizationParamMgr.getInstance(model);
try
	ModelParameterMgr.loadModels;
catch Me
    rethrow(Me);
end

% Second input should be either OperSpec or a time vector
if nargin >= 2
	op = varargin{1};
	if isa(op,'opcond.OperatingSpec')
		usesnapshot = false;
	elseif	isa(op,'double')
		usesnapshot = true;
	else
		ctrlMsgUtils.error('Slcontrol:findop:InvalidOperatingPoint');
	end
else
	% Return the operating point at the model initial condition
	usesnapshot = true;
	op = 0;
end

% By default, display to the workspace for trim.  Use an empty variable
% to optimize checking for this function during the optimization.
dispfcn = [];
stopfcn = [];
params = [];
options = [];

% Parse 3rd - 6th input arguments
if nargin > 2	
	if isa(varargin{2}, 'trim.TrimOptions') || ...
			isa(varargin{2}, 'linearize.OldLinoptions') % Legacy data
		% No params variation provided: findop(model, ops, opt)
		options = varargin{2};		
	else
		% With params variation provided: findop(model, ops, param, opt)
		params = varargin{2};	
		if 	nargin > 3
			options = varargin{3};
		end
	end	
	
	if nargin >= 5
		% Display in the output function using this function handle for trim
		% the stop function is a function to check to stop the optimization.
		% The display function should be in the form a vector cell array where
		% the first element is a function handle.  The updated string then
		% will be added as a last argument.
		dispfcn = varargin{end-1};
		stopfcn = varargin{end};
	end
end

% Set default options, if it is not provided.
if isempty(options)
	options = findopOptions;    
end

% Create a TrimObject
try
	this = trim.TrimObject(ModelParameterMgr, op, params, options, ...
		dispfcn, stopfcn, usesnapshot);
catch ExCreateTrimObject
	throwAsCaller(ExCreateTrimObject);
end

% Compute the grid and store original parameters
precompile(this);

if ~usesnapshot
	isarray = numel(op) > 1;
	% If op is an array, make sure model and time fields are the same
	if isarray
		t = get(op,'Time');
		opmodel = unique(get(op,'Model'));
		if numel(unique(cell2mat(t))) > 1
			ctrlMsgUtils.error('Slcontrol:findop:TimeDoesNotMatchInSpecificationArray');
		end
		if numel(opmodel) > 1
			ctrlMsgUtils.error('Slcontrol:findop:ModelDoesNotMatchSpecificationArray',model);
		end
	else
		opmodel = op.Model;
    end
    
    % Parallelization is not supported for optimization based trimming, for
    % details see g2062044
    if isfield(options.OptimizationOptions,'UseParallel') ...
            && ~isempty(options.OptimizationOptions.UseParallel)...
            && options.OptimizationOptions.UseParallel
        ctrlMsgUtils.error('Slcontrol:findop:InvalidUseParallelOption');
    end
	
	% Make sure that the model name matches the operating point
	% specification.
	if ~strcmp(model,opmodel)
		ctrlMsgUtils.error('Slcontrol:findop:ModelDoesNotMatchSpecification',model,opmodel);
	end
end

% Search for the operating point
try
	[oppoint,opreport] = fun_modified_trimModel(this);
catch TrimException
	throwAsCaller(TrimException)
end

% Export trimming report
if nargout == 2
	varargout{1} = opreport;
end
end


