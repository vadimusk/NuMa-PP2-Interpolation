function speichere_plot(dateiname)
% SPEICHERE_PLOT  Aktuelle Abbildung als PNG nach ../plots/matlab/ speichern.
% Immer im hellen Design (sonst übernimmt MATLAB ab R2025a den Dark Mode des Systems).

if ~exist('../plots/matlab', 'dir'), mkdir('../plots/matlab'); end
if ~isMATLABReleaseOlderThan("R2025a")
    theme(gcf, 'light');
end
exportgraphics(gcf, fullfile('../plots/matlab', dateiname), 'Resolution', 150);
end
