function [meanOut, stdOut, AmpOut, fBound, FinalFit] = Gaussian_Fit(Data, Ngauss, nbin, ...
                                                                    normalization, MaxIter, MaxFunEvals)
    % This function fits Gaussian distributions to the given data. Depending on the 
    % value of Ngauss, it fits either two or three Gaussian distributions.

    % Calculate histogram counts and scale factor
    [hcount, ~] = histcounts(Data, nbin, 'Normalization', normalization);
    scale = sum(hcount);  % Normalize the histogram to the pdf/probability

    if Ngauss == 2
        %% Two-Gaussian Fit
        % Define the Gaussian mixture model for fitting two Gaussian distributions
        pdf_normmixture = @(Data, p, mu1, mu2, sigma1, sigma2) ...
            p * normpdf(Data, mu1, sigma1) + (1 - p) * normpdf(Data, mu2, sigma2);

        % Initial guesses for the parameters
        pStart = 0.5;  % Initial mixing proportion
        muStart = quantile(Data, [.25 .75]);  % Initial means based on data quartiles
        sigmaStart = sqrt(var(Data) - 0.25 * diff(muStart).^2);  % Initial standard deviations
        start = [pStart muStart sigmaStart sigmaStart];  % Initial parameter vector

        % Define parameter bounds
        lb = [0 -Inf -Inf 0 0];  % Lower bounds
        ub = [1 Inf Inf Inf Inf];  % Upper bounds

        % Set optimization options
        options = statset('MaxIter', MaxIter, 'MaxFunEvals', MaxFunEvals);

        % Perform Maximum Likelihood Estimation (MLE) to estimate parameters
        paramEsts = mle(Data, normalization, pdf_normmixture, 'Start', start, ...
                        'LowerBound', lb, 'UpperBound', ub, 'Options', options);

        % Define x-axis range for plotting
        xgrid = linspace(-7, 0, 500);

        % Extract parameters for the first Gaussian
        mean1 = paramEsts(2);
        std1 = paramEsts(4);
        FitAmp(1) = (paramEsts(1) / (sqrt(2 * pi) * std1)) / scale;
        gauss1 = FitAmp(1) * exp(-(xgrid - mean1).^2 / (2 * std1^2));

        % Extract parameters for the second Gaussian
        mean2 = paramEsts(3);
        std2 = paramEsts(5);
        FitAmp(2) = ((1 - paramEsts(1)) / (sqrt(2 * pi) * std2)) / scale;
        gauss2 = FitAmp(2) * exp(-(xgrid - mean2).^2 / (2 * std2^2));

        % Store output parameters
        meanOut = [mean1, mean2];
        stdOut = [std1, std2];
        AmpOut = [FitAmp(1), FitAmp(2)];

        % Store the fractions of each Gaussian (bound and unbound)
        fBound = [paramEsts(1), (1 - paramEsts(1))];

        % Combine the Gaussian fits for plotting
        FinalFit = [xgrid', gauss1', gauss2', (gauss1 + gauss2)'];

    else
        %% Three-Gaussian Fit
        % Define the Gaussian mixture model for fitting three Gaussian distributions
        pdf_normmixture = @(Data, p1, p2, mu1, mu2, mu3, sigma1, sigma2, sigma3) ...
            p1 * normpdf(Data, mu1, sigma1) + p2 * normpdf(Data, mu2, sigma2) + ...
            (1 - (p1 + p2)) * normpdf(Data, mu3, sigma3);

        % Initial guesses for the parameters
        pStart = [1/3, 1/3];  % Initial mixing proportions
        muStart = quantile(Data, [1/6, 3/6, 5/6]);  % Initial means based on quantiles
        sigmaStart = sqrt(var(Data) - 0.25 * diff(muStart).^2);  % Initial standard deviations
        start = [pStart muStart sigmaStart sigmaStart(1)];  % Initial parameter vector

        % Define parameter bounds
        lb = [0 0 -Inf -Inf -Inf 0 0 0];  % Lower bounds
        ub = [1 1 Inf Inf Inf Inf Inf Inf];  % Upper bounds

        % Set optimization options
        options = statset('MaxIter', MaxIter, 'MaxFunEvals', MaxFunEvals);

        % Perform MLE to estimate parameters
        paramEsts = mle(Data, normalization, pdf_normmixture, 'Start', start, ...
                        'LowerBound', lb, 'UpperBound', ub, 'Options', options);

        % Define x-axis range for plotting
        xgrid = linspace(-7.0, 2.0, 200);

        % Extract parameters for the first Gaussian
        mean1 = paramEsts(3);
        std1 = paramEsts(6);
        FitAmp(1) = (paramEsts(1) / (sqrt(2 * pi) * std1)) / scale;
        gauss1 = FitAmp(1) * exp(-(xgrid - mean1).^2 / (2 * std1^2));

        % Extract parameters for the second Gaussian
        mean2 = paramEsts(4);
        std2 = paramEsts(7);
        FitAmp(2) = (paramEsts(2) / (sqrt(2 * pi) * std2)) / scale;
        gauss2 = FitAmp(2) * exp(-(xgrid - mean2).^2 / (2 * std2^2));

        % Extract parameters for the third Gaussian
        mean3 = paramEsts(5);
        std3 = paramEsts(8);
        FitAmp(3) = ((1 - (paramEsts(1) + paramEsts(2))) / (sqrt(2 * pi) * std3)) / scale;
        gauss3 = FitAmp(3) * exp(-(xgrid - mean3).^2 / (2 * std3^2));

        % Store output parameters
        meanOut = [mean1, mean2, mean3];
        stdOut = [std1, std2, std3];
        AmpOut = [FitAmp(1), FitAmp(2), FitAmp(3)];

        % Store the fractions of each Gaussian (bound, intermediate, and unbound)
        fBound = [paramEsts(1), paramEsts(2), 1 - (paramEsts(1) + paramEsts(2))];

        % Combine the Gaussian fits for plotting
        FinalFit = [xgrid', gauss1', gauss2', gauss3', (gauss1 + gauss2 + gauss3)'];
    end
end
