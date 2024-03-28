

bar(binPoints, histPoints*100, 'BarWidth', 1, 'FaceAlpha', 0.1); % Adjust BarWidth and FaceAlpha
hold on;

plot(xgrid, gauss1*100, 'LineWidth', 2);
hold on;

plot(xgrid, gauss2*100, 'LineWidth', 2);
hold on;

plot(xgrid, gauss3*100, 'LineWidth', 2);
hold on;

plot(xgrid, gauss_total*100, 'LineWidth', 2);
hold off;



%%-- Comment it if you dont want to show the range. 

% rectangle('Position', [range(1), 0, range(2)-range(1), max(histPoints)], ...
%             FaceColor=[0, 1, 1, 0.2], EdgeColor='k');
% rectangle('Position', [range(3), 0, range(4)-range(3), max(histPoints)], ...
%             FaceColor=[1, 1, 0, 0.2], EdgeColor='k');
% rectangle('Position', [range(5), 0, range(6)-range(5), max(histPoints)], ...
%             FaceColor=[0.4660, 0.6740, 0.1880, 0.2], EdgeColor='k');


xlabel("log_{10}(D) (\mu m^2 /s)")
ylabel("Frequency(%)")
axis square
set(gca, 'FontSize', 18);

% Add a textbox at the top right corner
annotation('textbox', [0.6, 0.80, 0.2, 0.1], 'String', 'mage1 Metaphase', 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');