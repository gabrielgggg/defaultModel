clear; clc;

params = dlmread('parameters.tab'); %#ok<DLMRD> 
ix = 1;
ySz = params(ix); ix = ix + 1;
bSz = params(ix); ix = ix + 1;
crra = params(ix); ix = ix + 1;
rf = params(ix); ix = ix + 1;
delta = params(ix); % ix = ix + 1;
clear params ix;

yGrid = loadBinary('yGrid.bin', 'float64', [ySz, 1]);
yPi = loadBinary('yPi.bin', 'float64', [ySz, ySz]);
bGrid = loadBinary('bGrid.bin', 'float64', [bSz, 1]);

V = loadBinary('V.bin', 'float64', [ySz, bSz]);
Vr = loadBinary('Vr.bin', 'float64', [ySz, bSz]);
Vd = loadBinary('Vd.bin', 'float64', [ySz, 1]);
q = loadBinary('q.bin', 'float64', [ySz, bSz]);
dPol = loadBinary('dPol.bin', 'float64', [ySz, bSz]);
bPol = loadBinary('bPol.bin', 'float64', [ySz, bSz, bSz]);

EbPr = loadBinary('EbPr.bin', 'float64', [ySz, bSz]);

save results.mat;
