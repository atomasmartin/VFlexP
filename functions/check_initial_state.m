% Define the model name
modelName = bdroot;

% Load the model if it is not already loaded
if ~bdIsLoaded(modelName)
    load_system(modelName);
end

% Find all Integrator blocks in the model
integratorBlocks = find_system(modelName, 'LookUnderMasks', 'on', 'BlockType', 'Integrator');

% Loop through each Integrator block to enable signal logging on their input signals
for i = 1:length(integratorBlocks)
    blockPath = integratorBlocks{i}; % Full path of the Integrator block
    
    % Get the input port handle for the Integrator block
    portHandles = get_param(blockPath, 'PortHandles');
    inputPortHandle = portHandles.Inport;
    
    % Get the signal line connected to the input port
    signalLine = get_param(inputPortHandle, 'Line');
    % Ensure the signal line exists and is valid
    if iscell(signalLine)
        % Get a source port
            srcHandle = get_param(signalLine{1}, "SrcPortHandle");
            set_param(srcHandle(1), 'DataLogging', 'on')
    else
        if signalLine(1) ~= -1
            % Get a source port
            srcHandle = get_param(signalLine(1), "SrcPortHandle");
            set_param(srcHandle(1), 'DataLogging', 'on')
        end
    end
end

% Set up simulation options to log data
simIn = Simulink.SimulationInput(modelName);
simIn = simIn.setModelParameter('SaveOutput', 'on', 'SaveTime', 'on', 'StopTime', '0'); % Simulate 0 s

% Run the simulation
simOut = sim(simIn);

for i = 1:numElements(simOut.logsout)
    if any(abs(simOut.logsout{i}.Values.Data) > 1e-2)
        blockPath = simOut.logsout{i}.BlockPath.convertToCell{1};
        portHandles = get_param(blockPath, 'PortHandles');
        outputPortHandle = portHandles.Outport;
        signalLine = get_param(outputPortHandle, 'Line');
        dstHandle = get_param(signalLine, "DstPortHandle");
        dstBlock = string(get_param(dstHandle, "Parent"));
        dstBlock(~strcmp(get_param(dstBlock, 'BlockType'), 'Integrator')) = '';
        simOut.logsout{i}.Values.Data
        disp(strcat('<a href="matlab:hilite_system(''', dstBlock, ''')">', dstBlock, '</a>'));
    end
end

    % Loop through each Integrator block to disable signal logging on their input signals
for i = 1:length(integratorBlocks)
    blockPath = integratorBlocks{i}; % Full path of the Integrator block
    
    % Get the input port handle for the Integrator block
    portHandles = get_param(blockPath, 'PortHandles');
    inputPortHandle = portHandles.Inport;
    
    % Get the signal line connected to the input port
    signalLine = get_param(inputPortHandle, 'Line');
    % Ensure the signal line exists and is valid
    if iscell(signalLine)
        % Get a source port
            srcHandle = get_param(signalLine{1}, "SrcPortHandle");
            set_param(srcHandle(1), 'DataLogging', 'off')
    else
        if signalLine(1) ~= -1
            % Get a source port
            srcHandle = get_param(signalLine(1), "SrcPortHandle");
            set_param(srcHandle(1), 'DataLogging', 'off')
        end
    end
end
