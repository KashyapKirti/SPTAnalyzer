
color1 = [0.4660 0.6740 0.1880];
color2 =  [0 0.4470 0.7410]	;   

bar(binPoints, histPoints , 'BarWidth', 1, 'FaceAlpha', 0.05, 'EdgeColor', 'k'); % Adjust BarWidth and FaceAlpha
hold on;

plot(xgrid, gauss1 , 'LineWidth', 2, 'Color', 'k');
hold on;

plot(xgrid, gauss2 , 'LineWidth', 2, 'Color', 'k');
hold on;

plot(xgrid, gauss_total , 'LineWidth', 2, 'Color',color2);
hold off;

%-- Comment it if you don't want to show the range. 

rectangle('Position', [range(1), 0, range(2)-range(1), max(histPoints)], ...
            'FaceColor', [0, 1, 1, 0.2], 'EdgeColor', 'k');
rectangle('Position', [range(3), 0, range(4)-range(3), max(histPoints)], ...
            'FaceColor', [1, 1, 0, 0.2], 'EdgeColor', 'k');

xlabel("log_{10}(D) (\mu m^2 /s)")
ylabel("Frequency")
%axis square
set(gca, 'FontSize', 20);


annotation('textbox', [0.73, 0.80, 0.2, 0.1], 'String', 'Data: mage1 Metaphase', 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none');
