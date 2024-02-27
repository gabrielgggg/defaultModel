clear; clc;
simData = dlmread('sim.tab'); %#ok<DLMRD> 
ix = 1;
ySimIx = simData(:, ix); ix = ix + 1;
bSimIx = simData(:, ix); ix = ix + 1;
bPrSimIx = simData(:, ix); ix = ix + 1;
dSimIx = simData(:, ix); ix = ix + 1;
spSim = simData(:, ix); ix = ix + 1;
cSim = simData(:, ix); ix = ix + 1;
gdpSim = simData(:, ix); ix = ix + 1;
tbSim = simData(:, ix);
clear simData ix;

spSim = (1 + spSim).^4 - 1;

sz = size(spSim, 1);
K = 20;
N = 20;

valid = false([sz, 1]);
for ix = K+N+1:sz
    valid(ix) = (sum(dSimIx(ix-N:ix)) == 0);
end

save simulation.mat;