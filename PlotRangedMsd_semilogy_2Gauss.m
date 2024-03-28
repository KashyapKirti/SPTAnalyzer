boundLT = rangedLagtime(:,1);
boundMSD = rangedMSD(:,1);

boundLT = boundLT * 1e3; % Converting seconds to milliseconds

UnboundLT = rangedLagtime(:,2);
UnboundMSD = rangedMSD(:,2);

UnboundLT = UnboundLT * 1e3;  % Converting seconds to milliseconds

%-- Polynomial Fit

% For Bound
x = log(boundLT);
y = log(boundMSD);

[p, S] = polyfit(x, y, 1);
[y_fit, ~] = polyval(p, x, S);

figure()

semilogy(exp(x), exp(y_fit), '-', 'LineWidth', 5.0, 'DisplayName', 'Power Law Fit (Bound)');
hold on; 
semilogy(exp(x), exp(y), 'o', 'DisplayName', 'MSD (Bound)', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);

% For Unbound
x = log(UnboundLT);
y = log(UnboundMSD);
[p, S] = polyfit(x, y, 1);
[y_fit, ~] = polyval(p, x, S);


semilogy(exp(x), exp(y_fit), '-', 'LineWidth', 5.0, 'DisplayName', 'Power Law Fit (Unbound)');
hold on;
semilogy(exp(x), exp(y), 'o', 'DisplayName', 'MSD (Unbound)', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);


axis square;

xlabel('Time Interval (ms) semilogy scale', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('MSD (\mu m^2) (Semilogy scale)', 'FontSize', 14, 'FontWeight', 'bold');

legend('Location', 'northwest', 'FontSize', 12, 'FontWeight', 'bold');

annotation('textbox', [0.57, 0.05, 0.2, 0.1], 'String', {'Data : mage1 Anaphase'}, ...
   'FontSize', 14, 'LineStyle', 'none');

set(gca, 'FontSize', 22, 'FontWeight', 'bold'); % Set font size and weight for axis ticks
