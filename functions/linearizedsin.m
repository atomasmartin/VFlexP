function blocklin = linearizedsin(BlockData)
    blocklin = diag(cos(BlockData.Inputs.Values));
end
