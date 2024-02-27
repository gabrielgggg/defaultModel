clear; clc;
load results.mat;
load simulation.mat;

ryIx = round(ySz/2, 0);
fixYs = [ ryIx-10,  ryIx, ryIx+10 ];

figure;
histogram(spSim(valid), 50, Normalization="probability");
xline(mean(spSim(valid)), 'LineWidth', 1)
xlim([0 0.1]);
xlabel('Spread');
ylabel('Frequency');
formatFigure;
% exportgraphics(gca, "spreadDist.pdf", "ContentType", "vector");
% title('Spreads');

figure;
histogram(bGrid(bSimIx(valid)), 50, Normalization="probability");
xline(mean(bGrid(bSimIx(valid))), 'LineWidth', 1)
xlabel('Debt');
ylabel('Frequency');
formatFigure;
% exportgraphics(gca, "debtDist.pdf", "ContentType", "vector");
% title('Debt');

figure;
plot(bGrid, squeeze( q(fixYs, :) ));
title('q');

figure;
sp = (delta + rf) * (1 ./ squeeze( q(fixYs, :) ) - 1);
plot(bGrid, (1+sp).^4 - 1, 'LineWidth',2);
ylim([0.0 0.1]);
xlabel('Debt Next Period (B'')');
ylabel('Spread');
formatFigure;
% exportgraphics(gca, "spreadCurves.pdf", "ContentType", "vector");


figure;
plot(bGrid, squeeze(bPol(fixYs, 150, :)), 'o-', 'LineWidth',2);
xline(bGrid(150), 'LineWidth', 1);
xlim([0.17, 0.24]);
xlabel('Debt');
ylabel('Choice Probability');
formatFigure;


figure;
tmp = squeeze(EbPr(fixYs, :));
tmp(squeeze(dPol(fixYs, :)) > 0.75) = NaN;
plot(bGrid, tmp, 'LineWidth',2);
hold on;
plot(bGrid, bGrid, '--k', 'LineWidth',1);
xlabel('Current Debt');
ylabel('Next Period Debt');
xlim([0 0.5]);
formatFigure;