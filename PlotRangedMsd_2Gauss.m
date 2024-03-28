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



plot(exp(x), log10(exp(y_fit)), '-', 'LineWidth', 2, 'DisplayName', ...
    'Power Law Fit (Bound)')
hold on;
plot(exp(x), log10(exp(y)), 'o', 'DisplayName', 'MSD (Bound)', 'LineWidth', 1.5,...
    'MarkerSize', 8);

% For Unbound
x = log(UnboundLT);
y = log(UnboundMSD);

[p, S] = polyfit(x, y, 1);
[y_fit, ~] = polyval(p, x, S);



plot(exp(x), log10(exp(y_fit)), '-', 'LineWidth', 2, 'DisplayName',...
    'Power Law Fit (Unbound)')
hold on;
plot(exp(x), log10(exp(y)), 'o', 'DisplayName', 'MSD (Unbound)', ...
    'LineWidth', 1.5, 'MarkerSize', 8);


xlabel('Lag Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('MSD (\mum^2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Mean Squared Displacement vs Lag Time', 'FontSize', 14, 'FontWeight', 'bold');
axis square;
legend('show', 'Location', 'Northwest', 'FontSize', 10, 'FontWeight', 'bold');


grid on;
