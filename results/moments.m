clc; clear; 
load results.mat;
load simulation.mat;

bSim = bGrid(bSimIx);
bPrSim = bGrid(bPrSimIx);

loggdp = log(gdpSim(valid));
logc = log(cSim(valid));
tby = tbSim(valid) ./ gdpSim(valid);

fprintf("Mean Debt to GDP   %10.2f \n", 100.0 * mean(bSim(valid) ./ gdpSim(valid) ./ 4));
fprintf("Mean Spread        %10.2f \n", 100.0 * mean(spSim(valid)));
fprintf("Std Spread         %10.2f \n", 100.0 * std(spSim(valid)));

fprintf("Std log C          %10.2f \n", 100.0 * std(logc));
fprintf("Std log GDP        %10.2f \n", 100.0 * std(loggdp));

fprintf("Corr Sp, GDP       %10.2f \n", 100.0 * corr(spSim(valid), loggdp));
fprintf("Corr TB/GDP, GDP   %10.2f \n", 100.0 * corr(tby, loggdp));