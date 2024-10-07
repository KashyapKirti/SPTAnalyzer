% Convert lag times to milliseconds
boundLT = rangedLagtime(:,1) * 1e3;  % Bound Lag Time in milliseconds
boundMSD = rangedMSD(:,1);           % Bound MSD

UnboundLT = rangedLagtime(:,2) * 1e3;  % Unbound Lag Time in milliseconds
UnboundMSD = rangedMSD(:,2);           % Unbound MSD

%-- Polynomial Fit for Bound Data
x_bound = log(boundLT);
%y_bound = log(boundMSD);
y_bound = log(boundMSD-4*sigma^2);

[p_bound, S_bound] = polyfit(x_bound, y_bound, 1);  % Polynomial fit
[y_fit_bound, ~] = polyval(p_bound, x_bound, S_bound);

%-- Polynomial Fit for Unbound Data
x_unbound = log(UnboundLT);
y_unbound = log(UnboundMSD);

[p_unbound, S_unbound] = polyfit(x_unbound, y_unbound, 1);  % Polynomial fit
[y_fit_unbound, ~] = polyval(p_unbound, x_unbound, S_unbound);

%-- Plotting
figure();

% Plot for Bound Data
plot(exp(x_bound), (exp(y_fit_bound)), '-', 'LineWidth', 2, ...
    'DisplayName', 'Power Law Fit (Bound)');
hold on;
plot(exp(x_bound), (exp(y_bound)), 'o', 'DisplayName', 'MSD (Bound)', ...
    'LineWidth', 1.5, 'MarkerSize', 8);

% % Plot for Unbound Data
% plot(exp(x_unbound), log10(exp(y_fit_unbound)), '-', 'LineWidth', 2, ...
%     'DisplayName', 'Power Law Fit (Unbound)');
% plot(exp(x_unbound), log10(exp(y_unbound)), 'o', 'DisplayName', 'MSD (Unbound)', ...
%     'LineWidth', 1.5, 'MarkerSize', 8);

% Labels, Title, and Legend
xlabel('Lag Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('MSD (\mum^2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Mean Squared Displacement vs Lag Time', 'FontSize', 14, 'FontWeight', 'bold');
axis square;
legend('show', 'Location', 'Northwest', 'FontSize', 10, 'FontWeight', 'bold');

% Grid for better readability
grid on;
