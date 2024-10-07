boundLT = rangedLagtime(:,1);
boundMSD = rangedMSD(:,1);

boundLT = boundLT * 1e3; % Converting seconds to milliseconds

UnboundLT = rangedLagtime(:,2);
UnboundMSD = rangedMSD(:,2);

UnboundLT = UnboundLT * 1e3;  % Converting seconds to milliseconds

%-- Polynomial Fit

% For Bound
x = log(boundLT);
 %y = log(boundMSD);
y = log(boundMSD);

[p, S] = polyfit(x, y, 1);
[y_fit, ~] = polyval(p, x, S);

alpha_bound = p(1);

figure()

loglog(exp(x), exp(y_fit), '-', 'LineWidth', 10.0, 'DisplayName', 'Power Law Fit (t^{\alpha}) ',  "Color","#A2142F"	);
hold on; 
loglog(exp(x), exp(y), 'o', 'DisplayName', 'MSD (Bound)', 'MarkerSize', 24, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);

% text(30,1.00e-3, ...
%     ['\alpha = ' num2str(alpha_bound, '%.2f')], ...
%     'FontSize', 46, 'Color', '#A2142F', 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation',23);

text(30,0.13e-3, ...
    ['\alpha = ' num2str(alpha_bound, '%.2f')], ...
    'FontSize', 46, 'Color', '#A2142F', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation',8);


% For Unbound
x = log(UnboundLT);
% y = log(UnboundMSD);
y = log(UnboundMSD);
[p, S] = polyfit(x, y, 1);
[y_fit, ~] = polyval(p, x, S);

alpha_unbound = p(1);

loglog(exp(x), exp(y_fit), '-', 'LineWidth', 10.0, 'DisplayName', 'Power Law Fit (t^{\alpha})',"Color","b");
hold on;
loglog(exp(x), exp(y), 'o', 'DisplayName', 'MSD (Unbound)', 'MarkerSize', 24, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);


% text(30,2.28e-3, ...
%     ['\alpha = ' num2str(alpha_unbound, '%.2f')], ...
%     'FontSize', 46, 'Color', 'b', 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation',30);
text(28,1.2e-3, ...
    ['\alpha = ' num2str(alpha_unbound, '%.2f')], ...
    'FontSize', 46, 'Color', 'b', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation',15);


% xlim([12, 100])
  ylim([0 4.0e-3])

axis square;

 xlabel('Time Interval (ms)', 'FontSize', 32);
% ylabel('MSD (\mu m^2)', 'FontSize', 32);

legend('Location', 'northwest', 'FontSize', 23);

%annotation('textbox', [0.57, 0.05, 0.2, 0.1], 'String', {'Data : mage1 Anaphase'}, ...
%   'FontSize', 14, 'LineStyle', 'none');
% yticks([1e-3, 2e-3, 3e-3, 4e-3]);
% yticks([1, 2, 3, 4] * 1e-3);
% Change y-tick labels to 1, 2, 3, 4
% yticklabels({'1', '2', '3', '4'});
% text(0.5, max(y)*1.05, '\times 10^{-4}', 'Units', 'data', 'FontSize', 12);
% ytickformat('%.4e');
% yticks([1, 2, 3, 4] * 1e-3); % This changes the y-tick positions
% yticklabels({'1', '2', '3', '4'}); % This changes the labels to simple numbers

 yticks([1, 2, 3, 4] * 1e-3);  % Set y-tick positions
ytickformat('%.0e');  % Format y-ticks to scientific notation (e.g., 1e-3, 2e-3)
yticklabels({'1', '2', '3', '4'}); % This changes the labels to simple numbers
yl=ylabel({'MSD(t)\times 10^{-3}'});
yl.Position(1) = yl.Position(1) - 1.0;  % Adjust the offset (negative moves left)


set(gca, 'FontSize', 52, "LineWidth", 2); % Set font size and weight for axis ticks
