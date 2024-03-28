boundLT = rangedLagtime(:,1);
boundMSD = rangedMSD(:,1);
boundLT = boundLT * 1e3; % Converting seconds to milliseconds

intMLT = rangedLagtime(:,2);
intMMSD = rangedMSD(:,2);
intMLT = intMLT * 1e3; % Converting seconds to milliseconds

UnboundLT = rangedLagtime(:,3);
UnboundMSD = rangedMSD(:,3);
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

% For Intermediate
x = log(intMLT);
y = log(intMMSD);

[p, S] = polyfit(x, y, 1);
[y_fit, ~] = polyval(p, x, S);

semilogy(exp(x), exp(y_fit), '-', 'LineWidth', 5.0, 'DisplayName', 'Power Law Fit (Intermediate)');
hold on; 
semilogy(exp(x), exp(y), 'o', 'DisplayName', 'MSD (Intermediate)', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);

% For Unbound
x = log(UnboundLT);
y = log(UnboundMSD);
[p, S] = polyfit(x, y, 1);
[y_fit, ~] = polyval(p, x, S);

semilogy(exp(x), exp(y_fit), '-', 'LineWidth', 5.0, 'DisplayName', 'Power Law Fit (Unbound)');
hold on;
semilogy(exp(x), exp(y), 'o', 'DisplayName', 'MSD (Unbound)', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);

axis square;

xlabel('Time Interval (ms)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('MSD (\mu m^2) semilogy scale', 'FontSize', 24, 'FontWeight', 'bold');

annotation('textbox', [0.57, 0.05, 0.2, 0.1], 'String', {'Data : WT TFA1 10ms REP3'}, ...
   'FontSize', 14, 'LineStyle', 'none');

legend('Location', 'northwest', 'FontSize', 10, 'FontWeight', 'bold');

set(gca, 'FontSize', 22, 'FontWeight', 'bold'); % Set font size and weight for axis ticks
