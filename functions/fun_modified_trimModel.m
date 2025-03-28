function [oppoint,opreport] = fun_modified_trimModel(this)
%% TRIMMODEL: method to trim the model.

% Copyright 2016 The MathWorks, Inc.

% Prepare options
options = this.options;
%checkProjOptimizer(this);
% Trim the model
ModelParameterMgr = this.ParamMgr;
op = this.OperSpec;
% Make sure the model is loaded
model = ModelParameterMgr.Model;
iscompiled = isCompiled(ModelParameterMgr);
% get the compilation manager
cmgr = slcontrollib.internal.utils.getCompilationMgr(model);

if (this.OperSpecType == linearize.OpTypeEnum.Snapshot)

    % snapsot is not supported for fast restart mode
    if isFastRestartOn(cmgr) && iscompiled
        warning(message('Slcontrol:findop:FastRestartSnapshotCheck'));
    end

	if iscompiled
		ModelParameterMgr.term();
	end
	opreport = [];	% snapshot does not provide trimming report
	nump = prod(this.ParamGridSize);
	oppoint = [];
	% Loop over parameters
	for ct = 1:nump
		pushParameters(this,ct,false);
		try
			oppoint = [oppoint;
				runsnapshot(LinearizationObjects.OperPointSnapShotEvent(...
				ModelParameterMgr,op))]; %#ok<AGROW>

		catch SimSnapshotException
			ModelParameterMgr.closeModels;
			throwAsCaller(SimSnapshotException)
		end
	end
	restoreParameters(this);
	ModelParameterMgr.closeModels;
	% Get rid of the first dimension for snapshot times if there is only
	% one snapshot. We need to use numel(oppoint)/prod(this.ParamGridSize)
	% instead of numel(this.OperSpec), because user can introduce snapshot
	% block that takes extra snapshots.
	opsize = [numel(oppoint)/prod(this.ParamGridSize) this.ParamGridSize];
	if opsize(1) == 1 && numel(opsize) > 1
		opsize(1) = [];
	end
	oppoint = reshape(oppoint,opsize);

else % Branch for operspec objects

	% Get the desired Simulink set-param properties
	if numel(op(1).Inputs)
		useModelu = false;
	else
		useModelu = true;
	end

	[ConfigSetParameters,ModelParams] = linearize.linutil.createLinearizationParams(false,useModelu,op(1).Time,options,true);
	ModelParams.SimulationMode = 'normal';

	% Determine whether to use linearization points.  If all the output constraints
	% are outports then use the standard model output evaluation function.
	% Otherwise find all the source ports for the output constraints and mark
	% them for linear analysis.
	if isscalar(op)
		op(1:this.NModels) = op;
		this.OperSpec = op;
	end
	ios = cell(1,numel(op));
	% In scenario in the array case should be determined independently.
	portflag = false;

    allios = [];
    if strcmp(options.OptimizationOptions.Jacobian,'on')
        [ios,portflag] = slcontrollib.internal.utils.generateIOsFromOpspec(op);
    end
    % If any case required portflag, merge all I/Os and push them to the model
    if portflag
        for ctop = 1:numel(op)
            allios = [allios;ios{ctop}]; %#ok<AGROW>
        end
        models = getUniqueNormalModeModels(ModelParameterMgr);
        ModelParams.SCDPotentialLinearizationIOs = linearize.createSCDPotentialLinearizationIOsStructure(allios,models,[],false);
    end

    singlecompile = this.options.AreParamsTunable;

    % ----------- Begin of compilation -----------------------------------%

    % check to see if the model needs to be recompiled
    needsrecompile = fastRestartForLinearAnalysis(model,'check','AnalysisPoints',allios);
    if needsrecompile
        warning(message('Slcontrol:findop:FastRestartNeedsRecompilation',model));
    end

    if ~options.BypassCompilation || ~iscompiled
        if iscompiled && ~isFastRestartOn(cmgr)
            ModelParameterMgr.term();
        end
        ismodelcompiled = isCompiled(ModelParameterMgr);
        if ~isFastRestartOn(cmgr)
            setModelParameters(ModelParameterMgr,ModelParams,ConfigSetParameters);
            ModelParameterMgr.prepareModels;
            %setSimscapeSolverBlocks(this);
        end

        try
            if ~ismodelcompiled
                ModelParameterMgr.compile('lincompile');
            end
            ismodelcompiled = true;
            % Check to be sure that a single tasking solver is being used.
            linearize.linutil.checkSingleTaskingSolver(model);
            % Determine parameter usage and check if single compilation is possible
            if singlecompile
                singlecompile = LocalCacheParameterUsage(this);
            end
        catch CompileModelException
            if ~isFastRestartOn(cmgr)
                if ismodelcompiled
                    ModelParameterMgr.term();
                end
                ModelParameterMgr.restoreModels;
                restoreParameters(this);
                ModelParameterMgr.closeModels;
            end
            throwAsCaller(CompileModelException);
        end
    end
    % ----------- End of compilation -------------------------------------%

    % if fast restart is on, only tunable params are supported
    if ~singlecompile && isFastRestartOn(cmgr)
        error(message('Slcontrol:findop:FastRestartNonTunableParams'));
    end
    
    % create the opdata interface
    opdata = opcond.internal.createOPDataInterface(this.Model);
    
    % ignore input data type warnings when syncing from model
    ws = warning('off','SLControllib:opcond:ModelHasNonDoubleRootPortInputDataTypes');
    cln1 = onCleanup(@() warning(ws));

	if singlecompile
        % sync the XUY with the model once
        syncStates(opdata,true);
        syncInputs(opdata);
        syncOutputs(opdata);
        wb = waitbar(0, "Initialising");
		% Model must be compiled
		for ctop = 1:this.NModels
            waitbar((ctop-1)/this.NModels, wb, "Initialising")
			try
				pushParameters(this,ctop,true);
			catch E
				myE = MException(message('Slcontrol:sllinearizer:ParamVaryError'));
				myE = myE.addCause(E);
				throw(myE);
			end
			[oppoint(ctop),opreport(ctop)] = LocalTrimModel(this, ctop, ios, opdata); %#ok<AGROW>
        end
        close(wb)
		% Clean up
        if ~options.BypassTermination
            cleanup(this);
        end
	else % Multiple compilation


		% Model has been compiled, terminate it.
		term(ModelParameterMgr);

		for ctop = 1:this.NModels
			pushParameters(this,ctop,false);
            try
                ismodelcompiled = false;
                ModelParameterMgr.compile('lincompile');
                ismodelcompiled = true;
                % Check to be sure that a single tasking solver is being used.
                linearize.linutil.checkSingleTaskingSolver(model);
            catch CompileModelException
                if ismodelcompiled
                    ModelParameterMgr.term();
                end

                ModelParameterMgr.restoreModels;
                restoreParameters(this);
                ModelParameterMgr.closeModels;
                throwAsCaller(CompileModelException);
            end
            % sync the XUY each time the model is compiled
            syncStates(opdata,true);
            syncInputs(opdata);
            syncOutputs(opdata);
            
			[oppoint(ctop),opreport(ctop)] = LocalTrimModel(this, ctop, ios, opdata); %#ok<AGROW>
			term(ModelParameterMgr);
		end

        % Clean up
        restoreModels(ModelParameterMgr);
        restoreParameters(this);
        closeModels(ModelParameterMgr);
	end % End for single or multiple compilation

	% Reshape the output
	oppoint = reshape(oppoint, this.FullResultSize);
	opreport = reshape(opreport, this.FullResultSize);
	% Throw warning
	if ~isempty(this.OperSpecExpansion.WarningMessage)
		warning(this.OperSpecExpansion.WarningMessage)
	end
end % End for snapshot/optimization
end

function [oppoint,opreport] = LocalTrimModel(this, ctop, ios, opdata)
op = this.OperSpec;
options = this.options;
ModelParameterMgr = this.ParamMgr;
% Check to see that the operating condition is up to date
if ~any(strcmp('ByPassConsistencyCheck',properties(options))) || ...
		strcmp(options.ByPassConsistencyCheck,'off')
	try
        % update the spec with any model changes
        update(op(ctop),true,opdata);
        
		% If the operspec is expanded to match the params grid, the fixed
		% value of the Known states will be the same in all the expended
		% operspecs. They will not be reset by grabing the model initial
		% condition after applying new parameters. Therefore, a warning is
		% thrown in case one would expect the Known states to
		% automatically change.
		if this.OperSpecExpansion.isExpanded && isempty(this.OperSpecExpansion.WarningMessage)
			xstruct = feval(this.Model, 'get', 'state');
			this.OperSpecExpansion.WarningMessage = LocalDetectKnownStateChange(op(ctop), xstruct);
		end
		% update the fixed value of known states, if necessary

		% Update states
		% xstruct = feval(this.Model, 'get', 'state');
		% op(ctop).setxu(xstruct);
	catch UpdateOperatingPointException
		% Clean up
		% Release the compiled model
		cleanup(this);
		throwAsCaller(UpdateOperatingPointException);
	end
end

try
	optim = createOptimizationObject(this, ctop, ios);
catch CreateOptimizerObject
	% Release the compiled model
	cleanup(this);
	% Throw the error
	throwAsCaller(CreateOptimizerObject);
end

% tell the optimizer what iteration it is and if it is part of a batch
% process
optim.BatchIndex = ctop;
optim.IsBatchElement = numel(op) > 1;

% Set the display and stop functions
optim.dispfcn = this.dispfcn;
optim.stopfcn = this.stopfcn;

% Run the optimization
try
	[oppoint,opreport,exitflag,optimoutput] = optimize(optim);
catch SearchOperatingPointException
	% Release the compiled model
	cleanup(this);
	throwAsCaller(SearchOperatingPointException);
end

if isempty(exitflag)
    % This happens if no optimization is performed, e.g., validating an
    % existing op with given opspec.
    exitdata = [];
elseif exitflag > 0
	exitdata = ctrlMsgUtils.message('Slcontrol:findop:SuccessfulTermination');
elseif (exitflag == 0)
	exitdata = ctrlMsgUtils.message('Slcontrol:findop:MaximumFunctionEvalExceeded');
elseif (exitflag == -1)
	exitdata = ctrlMsgUtils.message('Slcontrol:findop:OptimizationTerminatedPrematurely');
elseif (exitflag == -2)
	exitdata = ctrlMsgUtils.message('Slcontrol:findop:CouldNotMeetConstraints');
else
	exitdata = ctrlMsgUtils.message('Slcontrol:findop:OptimizationDidNotConverge');
end

% Store the optimization data in the report
opreport.TerminationString = exitdata;
opreport.OptimizationOutput = optimoutput;

% Display the report of the optimization
if strcmp(options.DisplayReport,'on')
    if optim.IsBatchElement
        disp(ctrlMsgUtils.message('Slcontrol:findop:OptimizationOutputBatchDisplayTitle',optim.BatchIndex));
    else
        disp(ctrlMsgUtils.message('Slcontrol:findop:OptimizationOutputDisplayTitle'));
    end
	disp('---------------------------------');
	display(opreport);
elseif strcmp(options.DisplayReport,'iter') && isempty(this.dispfcn)
	fprintf(1,'\n%s\n\n',exitdata);
end

% Loop through all the states to be trimmed and check if it is inside a
% triggered/function-call subsystem.

% We need to temporarily load all models for this to be able to throw the
% same warning in accelerated model scenarios.
unloadedrefmdls = ModelParameterMgr.getUniqueUnLoadedRefModels;
for ct = 1:numel(unloadedrefmdls)
	load_system(unloadedrefmdls{ct});
end
% Loop states
for ct = 1:numel(op(ctop).States)
	if op(ctop).States(ct).SteadyState
		asyncsub = LocalGetAsyncSubsystem(op(ctop).States(ct).Block);
		if ~isempty(asyncsub)
			% Throw a warning about this state
			ctrlMsgUtils.warning('SLControllib:opcond:StateToBeTrimmedInAsyncSub',...
				op(ctop).states(ct).Block,asyncsub);
		end
	end
end
% Close them back
for ct = 1:numel(unloadedrefmdls)
	bdclose(unloadedrefmdls{ct});
end

end

%% -----------------------------------------------------------------------
% LocalCacheParameterUsage: Get parameter usage and check if a single
% compilation is possible.
function singlecompile = LocalCacheParameterUsage(this)
singlecompile = true;
% Early return for no parameter case
p = this.Parameters;
if isempty(p)
	return;
end

numparams = numel(p);
model = this.Model;

% Create objects
hwks = get_param(model,'ModelWorkspace');
modelvars = whos(hwks);
for ct = numparams:-1:1
	e = p(ct).Name;
	% Parse name of the variable
	v = slLinearizer.getVarNameFromExpression(e);

	% get variable object based on the workspace they live in
	if any(strcmp(v,{modelvars.name}))
		varobjs(ct) = Simulink.VariableUsage(v,model);
	else
		workspace = get_param(model, 'DataDictionary');
		if isempty(workspace)
			workspace = 'base workspace';
		end
		varobjs(ct) = Simulink.VariableUsage(v,workspace);
	end
end

% Find usage of everything together
% Assert that the model is compiled
assert(strcmp(get_param(model,'SimulationStatus'),'paused'));
varusage = Simulink.findVars(model,varobjs,'SearchMethod','cached');
this.ParameterUsage = varusage;

% Now check if we will be able to achieve single compilation
singlecompile = LocalCheckTunableParams(this);

end


%% -----------------------------------------------------------------------
% LocalCheckTunableParams: Check if parameters are tunable.
function tunableparams = LocalCheckTunableParams(this)
tunableparams = true;
params = this.Parameters;
% Check individual parameter tunability
for ct = 1:numel(this.ParameterUsage)
	ud = this.ParameterUsage(ct).DirectUsageDetails;
    for ctdd = 1:numel(ud)
        prop = ud(ctdd).Properties;
        expressions = ud(ctdd).Expressions;
        blk = ud(ctdd).Identifier;
        % Check if it is a model block
        if strcmp(get_param(blk,'BlockType'),'ModelReference')
            param = this.ParameterUsage(ct).Name;
            % Check if parameterized. If so, jump to the next param
            if linearize.linutil.checkParameterizedModelArguments(blk, param)
                continue
            end
            
            % Otherwise, this has to be a model workspace parameter
            srctype = this.ParameterUsage(ct).SourceType;
            assert(any(ismember({'base workspace','data dictionary'},srctype)));
            paramval = Simulink.data.evalinGlobal(this.Model,param);
            
            if ~isa(paramval,'Simulink.Parameter') || strcmp(paramval.StorageClass,'auto')
                warning(message('Slcontrol:findop:warnNoSingleCompilation'));
                tunableparams = false;
                return;
            end
        end
		dlg_param = get_param(blk,'DialogParameters');
		for ctp = 1:numel(prop)
			if isempty(prop{ctp})
				continue;
			end
			attribs = dlg_param.(prop{ctp}).Attributes;
			% Check if property is not tunable
			istunable = ~any(strcmp(attribs, 'read-only-if-compiled')) && ...
				~any(strcmp(attribs, 'read-only'));
			% If not tunable, check if any param.Name is contained in the expression
			if ~istunable
				for ctpp = 1:numel(params)
					% Strip out last indexing element
					if ~isempty(strfind(expressions{ctp},LocalStripIndexing(params(ctpp).Name)))
						% Warn and break
						warning(message('Slcontrol:findop:ParamsNotTunable', params(ctpp).Name));
						tunableparams = false;
						return;
					end
				end
			end
		end
	end
end
end

%% -----------------------------------------------------------------------
% LocalStripIndexing: process parameter names.
function strippedname = LocalStripIndexing(paramname)
ind = regexp(paramname,'[({]');
if isempty(ind)
	strippedname = paramname;
else
	strippedname = paramname(1:ind(end)-1);
end
end

%% -----------------------------------------------------------------------
function sub = LocalGetAsyncSubsystem(blockpath)
%% LocalGetAsyncSubsystem
% Return the conditional subsystem in the hierarchy of the block. Return
% empty if none.
sub = [];
% Handle blocks in model reference
cellpath = convertToCell(slcontrollib.internal.utils.getSimulinkBlockPath(blockpath));
for ct = numel(cellpath):-1:1
	thisparent = get_param(cellpath{ct},'Parent');
	% Check each level if there is a trigger or enable block at that level
	while ~isempty(thisparent)
		% Look for enable or trigger block at this level
		triggerblks = find_system(thisparent,'SearchDepth',1,'BlockType','TriggerPort');
		if ~isempty(triggerblks)
			sub = thisparent;
			return;
		end
		thisparent = get_param(thisparent,'Parent');
	end
end
end


%% -----------------------------------------------------------------------
function ios = LocalGetIOsForSpec(op)
% getIOsForSpec
% Get the I/O list for a given spec
ios = [];
% Get ready to create the I/O required for linearization
outs = op.Outputs;
bh = get_param(get(outs,{'Block'}),'Object');
bt = get([bh{:}],{'BlockType'});
h = linearize.IOPoint;
h.Type = 'output';
for ct = 1:length(bt)
	ios = [ios;h.copy];
	if strcmpi(bt(ct),'Outport')
		% Get the source block
		ios(ct).Block = [op.model,'/',get_param(bh{ct}.PortConnectivity.SrcBlock,'Name')];
		% Get the source port
		ios(ct).PortNumber = bh{ct}.PortConnectivity.SrcPort + 1;
	else
		% Set the Block and PortNumber properties
		ios(ct).Block = outs(ct).Block;
		ios(ct).PortNumber = outs(ct).PortNumber;
	end
end
% Set the Analysis I/O properties for the input ports
ins = op.Inputs;
h.Type = 'input';
for ct = 1:length(ins)
	ios = [ios;h.copy];
	ios(end).Block = ins(ct).Block;
	ios(end).PortNumber = 1;
end
end

%% LocalDetectKnownStateChange
% Detect if there is any state or input change,
% in case of scalar expansion of the operspec.
function msg = LocalDetectKnownStateChange(op,xstruct)
%% checkKnownStateChange
% Detect if there is any state or input change,
% in case of scalar expansion of the operspec.

% Copyright 2016 The MathWorks, Inc.
msg = [];
% Extract the states from the operating condition object
states = op.States;
% Get the state names to match up with in the struct
statenames = get(states,'Block');
% Loop over each element in the new state structure
xsignals = xstruct.signals;
for ct = 1:length(xsignals)
	% Find the index of States in op that matches the ct-th signal in
	% xstruct
	ind = find(strcmp(xsignals(ct).blockName,statenames), 1);

	if ~isempty(ind)
		% Detect changes of known elementes in the located state
		isKnown = states(ind).Known;
		% Nothing to detect if there is no known states
		if any(isKnown)
			oldx = states(ind).x;           % originally specified x
			newx = xsignals(ct).values(:);  % new x after pushing params
			% Get the indices of the state elements that are:
			% 1) known, and
			% 2) changed after applying new params
			ind_known_and_changed = isKnown(:) & (oldx ~= newx);

			% Once detected, prepare and throw the warning message
			if any(ind_known_and_changed)

				% 1. Prepare block/state name
				% repair the blockname and use it as statename, if
				% statename is empty.
				BlockName = regexprep(states(ind).Block,'\n',' ');
				if isempty(states(ind).StateName)
					StateName = BlockName;
				else
					StateName = regexprep(states(ind).StateName,'\n',' ');
				end
				% Add link to the block in the string
				str = sprintf('slcontrollib.internal.utils.dynamicHiliteSystem(''%s'');',BlockName);
				if usejava('Swing') && desktop('-inuse') && feature('hotlinks')
					strStateName = sprintf('<a href="matlab:%s">%s</a>',str,StateName);
				else
					strStateName = sprintf('%s',StateName);
				end

				% 2. Export warning info
				% Construct the warning message
				msg = message('Slcontrol:findop:OPExpansionWithKnownStates', ...
					strStateName);
				return
			end
		end

	end
end
end
