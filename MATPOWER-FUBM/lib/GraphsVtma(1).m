function GraphsVtma(i, fVtma)
%GRAPHSVTMA Plots the convergence of "ma" for Vt control.

%   GRAPHSVTMA(I, FVTMA)
%
%   Plots the convergence of "ma" for Vt control.
%   Inputs:
%       I : Iteration number 
x = 0:i;
figure
plot(x',fVtma')
title('Vt control ma convergence')
