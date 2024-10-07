% Define colors used for plotting different components
color1 = [0.4940 0.1840 0.5560]; % Anaphase main
color2 = [0.8500 0.3250 0.0980 0.5]; % Anaphase 1
color3 = [0.4660 0.6740 0.1880 0.5]; % Anaphase 2
color4 = [1 0 1]; % Unspecified
grey = [0.5 0.5 0.5]; % Grey color
color5 = [0.6350 0.0780 0.1840]; % Unspecified
bl = [0 0 1 0.5]; % Metaphase 1
pl = [1 0 1 0.5]; % Metaphase 2
rl = [1 0 0 0.6]; % Metaphase main
color6 = [21/255 20/255 216/255 ]; % Cystol 1
color7 = [140/255 57/255 69/255 0.5 ]; % Cystol 2
color8 = [108/255 155/255 3/255 0.5 ]; % Cystol 3
color9 = [0/255 163/255 96/255]; % SPB Main
color10 = [69/255 76/255 240/255 0.5]; % SPB 1
color11 = [100/255 18/255 91/255 0.5]; % SPB 2

% Plot histogram of data
bar(binPoints, histPoints * 100, 'BarWidth', 1, 'FaceAlpha', 0.0001, 'EdgeColor', 'k');
hold on;

% Plot the first Gaussian fit
plot(xgrid, gauss1 * 100, 'LineWidth', 5, 'Color', "#A2142F");
hold on;

% Plot the second Gaussian fit
plot(xgrid, gauss2 * 100, 'LineWidth', 5, 'Color', 'b');
hold on;

% Plot the total fit (sum of Gaussians)
plot(xgrid, gauss_total * 100, 'LineWidth', 5, 'Color', "#027148");
hold off;

% Set axis labels and appearance
xlabel("log_{10}(D) (\mu m^2 /s)");
ylabel("Frequency(%)");
set(gca, 'FontSize', 50,'LineWidth', 3);
ylim([0,40])
xlim([-7 0])
xticks([-7 -6 -5 -4 -3 -2 -1 0]);  % Add xticks

% Prepare strings for annotation
% str1 = strcat("F_{bound} = ", num2str(round(fBound(1)*100), 2), "%");
% str2 = strcat("F_{free} = ", num2str(round(fBound(2)*100), 2), "%");
% 
% % Annotate the plot with bound fraction
% annotation('textbox', [0.15, 0.80, 0.25, 0.1], 'String', str1, 'FontSize', 40, ...
%            'EdgeColor', 'none', 'Color', 'k');
% 
% % Annotate the plot with bound diffusion coefficient
% str1 = strcat("D_{bound} = ", num2str(round(10^(meanOut(1)), 4)));
% str2 = strcat("D_{free} = ", num2str(round(10^(meanOut(2)), 4)));
% annotation('textbox', [0.15, 0.70, 0.25, 0.1], 'String', str1, 'FontSize', 40, ...
%            'EdgeColor', 'none', 'Color', 'k');
% annotation('textbox', [0.15, 0.60, 0.25, 0.1], 'String', str2, 'FontSize', 40, ...
%            'EdgeColor', 'none', 'Color', 'k');
%        

       
    