boundLT = rangedLagtime(:,1); % Extract the first column (Bound Lag Time) from the input data
boundMSD = rangedMSD(:,1); % Extract the first column (Bound MSD) from the input data
boundLT = boundLT * 1e3; % Convert bound lag times from seconds to milliseconds

intMLT = rangedLagtime(:,2); % Extract the second column (Intermediate Lag Time) from the input data
intMMSD = rangedMSD(:,2); % Extract the second column (Intermediate MSD) from the input data
intMLT = intMLT * 1e3; % Convert intermediate lag times from seconds to milliseconds

UnboundLT = rangedLagtime(:,3); % Extract the third column (Unbound Lag Time) from the input data
UnboundMSD = rangedMSD(:,3); % Extract the third column (Unbound MSD) from the input data
UnboundLT = UnboundLT * 1e3;  % Convert unbound lag times from seconds to milliseconds

%-- Polynomial Fit

% For Bound data
x = log(boundLT); % Take the natural logarithm of the bound lag times
y = log(boundMSD); % Take the natural logarithm of the bound MSD

[p, S] = polyfit(x, y, 1); % Perform a first-order polynomial fit (linear fit) on the log-log data
[y_fit, ~] = polyval(p, x, S); % Evaluate the fitted polynomial at the points in x

figure() % Create a new figure window

% Plot the power-law fit for the Bound data on a semilogarithmic scale
semilogy(exp(x), exp(y_fit), '-', 'LineWidth', 5.0, 'DisplayName', 'Power Law Fit (Bound)');
hold on; % Retain the current plot when adding new plots

% Plot the actual MSD (Bound) data on a semilogarithmic scale
semilogy(exp(x), exp(y), 'o', 'DisplayName', 'MSD (Bound)', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);

% For Intermediate data
x = log(intMLT); % Take the natural logarithm of the intermediate lag times
y = log(intMMSD); % Take the natural logarithm of the intermediate MSD

[p, S] = polyfit(x, y, 1); % Perform a first-order polynomial fit (linear fit) on the log-log data
[y_fit, ~] = polyval(p, x, S); % Evaluate the fitted polynomial at the points in x

% Plot the power-law fit for the Intermediate data on a semilogarithmic scale
semilogy(exp(x), exp(y_fit), '-', 'LineWidth', 5.0, 'DisplayName', 'Power Law Fit (Intermediate)');
hold on; % Retain the current plot when adding new plots

% Plot the actual MSD (Intermediate) data on a semilogarithmic scale
semilogy(exp(x), exp(y), 'o', 'DisplayName', 'MSD (Intermediate)', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);

% For Unbound data
x = log(UnboundLT); % Take the natural logarithm of the unbound lag times
y = log(UnboundMSD); % Take the natural logarithm of the unbound MSD

[p, S] = polyfit(x, y, 1); % Perform a first-order polynomial fit (linear fit) on the log-log data
[y_fit, ~] = polyval(p, x, S); % Evaluate the fitted polynomial at the points in x

% Plot the power-law fit for the Unbound data on a semilogarithmic scale
semilogy(exp(x), exp(y_fit), '-', 'LineWidth', 5.0, 'DisplayName', 'Power Law Fit (Unbound)');
hold on; % Retain the current plot when adding new plots

% Plot the actual MSD (Unbound) data on a semilogarithmic scale
semilogy(exp(x), exp(y), 'o', 'DisplayName', 'MSD (Unbound)', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);

axis square; % Set the aspect ratio so that the data units are the same in every direction

% Label the x-axis with a description and set font properties
xlabel('Time Interval (ms)', 'FontSize', 24, 'FontWeight', 'bold');

% Label the y-axis with a description and set font properties
ylabel('MSD (\mu m^2) semilogy scale', 'FontSize', 24, 'FontWeight', 'bold');

% Add an annotation (textbox) at the specified position with a custom message
annotation('textbox', [0.57, 0.05, 0.2, 0.1], 'String', {'Data : WT TFA1 10ms REP3'}, ...
   'FontSize', 14, 'LineStyle', 'none');

% Add a legend to the plot and set its properties
legend('Location', 'northwest', 'FontSize', 10, 'FontWeight', 'bold');

% Set the font size and weight for axis ticks, and adjust the line width of the plot box
set(gca, 'FontSize', 22, 'FontWeight', 'bold');
