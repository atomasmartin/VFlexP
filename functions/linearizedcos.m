function blocklin = linearizedcos(BlockData)
    blocklin = diag(-sin(BlockData.Inputs.Values));
end
