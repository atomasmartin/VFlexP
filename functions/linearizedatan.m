function blocklin = linearizedatan(BlockData)
    blocklin = eye(length(BlockData.Inputs.Values)).*1/(1+BlockData.Inputs.Values.^2);
end
