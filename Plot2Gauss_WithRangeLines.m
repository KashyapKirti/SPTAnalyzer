% Define colors used for plotting
color1 = [0.4660 0.6740 0.1880]; % Color for additional plots or annotations
color2 = [0 0.4470 0.7410]; % Color for the total Gaussian fit

% Plot the histogram of data with a low transparency
bar(binPoints, histPoints, 'BarWidth', 1, 'FaceAlpha', 0.05, 'EdgeColor', 'k');
hold on;

% Plot the first Gaussian fit
plot(xgrid, gauss1, 'LineWidth', 2, 'Color', 'k');
hold on;

% Plot the second Gaussian fit
plot(xgrid, gauss2, 'LineWidth', 2, 'Color', 'k');
hold on;

% Plot the total Gaussian fit
plot(xgrid, gauss_total, 'LineWidth', 2, 'Color', color2);
hold off;

% Optional: Highlight specific ranges on the plot with rectangles
rectangle('Position', [range(1), 0, range(2)-range(1), max(histPoints)], ...
          'FaceColor', [0, 1, 1, 0.2], 'EdgeColor', 'k'); % Cyan range
rectangle('Position', [range(3), 0, range(4)-range(3), max(histPoints)], ...
          'FaceColor', [1, 1, 0, 0.2], 'EdgeColor', 'k'); % Yellow range

% Set axis labels and font size
xlabel("log_{10}(D) (\mu m^2 /s)");
ylabel("Frequency");
set(gca, 'FontSize', 20);

% Add annotation to the plot
annotation('textbox', [0.73, 0.80, 0.2, 0.1], 'String', 'Data: mage1 Metaphase', ...
           'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
