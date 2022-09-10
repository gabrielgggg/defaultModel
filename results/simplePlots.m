clear; clc;
load results.mat;
load simulation.mat;

figure;
histogram(spSim(dSimIx == 0));
xline(mean(spSim(dSimIx==0)), 'LineWidth', 2)
title('Spreads');

figure;
histogram(bGrid(bSimIx(dSimIx == 0)));
title('Debt');

figure;
plot(bGrid, squeeze( q(1:3:end, :) ));
title('q');

figure;
plot(bGrid, squeeze(EbPr(1:3:end, :)));
hold on;
plot(bGrid, bGrid, '--k', 'LineWidth',2);
title('b''');

figure;
sp = (delta + rf) * (1 ./ squeeze( q(1:3:end, :) ) - 1);
plot(bGrid, (1+sp).^4 - 1);
ylim([0.0 0.1]);
title('Spread');
