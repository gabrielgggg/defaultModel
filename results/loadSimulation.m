clear; clc;
simData = dlmread('sim.tab'); %#ok<DLMRD> 
ix = 1;
ySimIx = simData(:, ix); ix = ix + 1;
bSimIx = simData(:, ix); ix = ix + 1;
bPrSimIx = simData(:, ix); ix = ix + 1;
dSimIx = simData(:, ix); ix = ix + 1;
spSim = simData(:, ix);
clear simData ix;

spSim = (1 + spSim).^4 - 1;

save simulation.mat;