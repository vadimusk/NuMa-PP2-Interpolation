function log_achse(titel, ylab)
% LOG_ACHSE  Einheitliche Beschriftung für die Plots über n (y-Achse logarithmisch).

set(gca, 'YScale', 'log');
xlim([2 40]);
grid on;  box off;
set(gca, 'YMinorGrid', 'off');
xlabel('Polynomgrad n');
ylabel(ylab);
title(titel);
end
