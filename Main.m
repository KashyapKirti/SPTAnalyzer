clc; close all; clear;

%% Load Data
% Load data from various CSV files. Uncomment the file you want to work with.

data = readtable("File_name.csv");

%% Data Filtration
% Define initial parameters for data filtering
n1 = 1;  % Starting data index
n2 = 5  % End data index

MinDataPoints =6;  % Minimum data points required for analysis
RsqTest = 0.8;  % Threshold for R-squared value in fitting
pixel_nm = 110;  % Conversion factor from pixels to nanometers
frameRate = 15; % Frame rate in milliseconds

% nbin = 10;  % Number of bins for histogram


%% Calculate Diffusion Coefficients
% Input parameters to calculate diffusion coefficients
inputParams = [n1,n2,MinDataPoints,RsqTest,pixel_nm, frameRate];
 
[D,Msd,LagTime,positiveD] = get_D_values(data, inputParams);
 
%% for getting sigma (localization error from fixed cells data)
% [MsdAll,TimeAll,MsdTime_AvgFinal] = FullMsd(data,inputParams);
% MSD_no_NaN = MsdAll(~any(isnan(MsdAll), 2), :);
% MSD_avg=mean(MSD_no_NaN,1);
% %
% figure(4)
% plot(TimeAll(1,:)*10^(3), MSD_avg, '-', 'LineWidth', 8.0, 'DisplayName', 'MSD fixed cells ',  "Color","#0072BD"	);
% hold on 
% plot( TimeAll(1,:)*10^(3), MSD_avg,'o', 'DisplayName', 'MSD (Bound)', 'MarkerSize', 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
% xlim([15 75])
% ylim([0.8e-4 5e-4])
% xlabel("Time interval(ms)");
% ylabel("MSD (\mu m^2)");
% set(gca, 'FontSize', 30);
%  set(gca, 'FontSize', 52, "LineWidth", 3); % Set font size and weight for axis ticks
%%
% sigma=sqrt(MSD_avg(1,1)/4)

sigma=0.0071; %% localization error


%% Gaussian Fitting
% Perform Gaussian fitting on the logarithm of diffusion coefficients
logD = log10(D);  % Take the logarithm of the diffusion coefficients
normalization = 'pdf';  % Set normalization method for histogram

binsize = 0.30%

nbin=int8(((max(logD)-min(logD))/binsize));  % Calculate bin size for histogram

% binsize =((max(logD)-min(logD))/nbin);

Ngauss = 2;  % Number of Gaussian distributions to fit

% Parameters for optimization
Data = logD;
MaxIter = 10000;  % Maximum iterations for fitting
MaxFunEvals = 2*MaxIter;  % Maximum function evaluations

% Get histogram data
[binPoints, histPoints] = get_Histogram(logD, nbin, normalization );

% Perform Gaussian fitting
[meanOut, stdOut, AmpOut, fBound,  FinalFit] = Gaussian_Fit(Data, Ngauss, nbin, ...
    normalization,  MaxIter,MaxFunEvals );


%% Calculating MSD for Bound and Unbound States
% Define range of indices for Mean Squared Displacement (MSD) calculation
nn1 = n1; nn2 =n2; 

% Calculate ranged MSD based on Gaussian fitting
[rangedLagtime, rangedMSD, range] = RangedMSD(Data, Msd,LagTime, Ngauss, meanOut, stdOut, nn1, nn2);

%% Plotting
% Plot results based on the number of Gaussian distributions fitted

if Ngauss==2
    xgrid = FinalFit(:,1);  % X-axis values for plotting
    gauss1 = FinalFit(:,2);  % First Gaussian component
    gauss2 = FinalFit(:,3);  % Second Gaussian component
    gauss_total = FinalFit(:,4);  % Total Gaussian fit
figure()
    Plot2Gauss  % Call function to plot 2-Gaussian fit
PlotRangedMsd_loglog_2Gauss


   
elseif Ngauss==3
    xgrid = FinalFit(:,1);
    gauss1 = FinalFit(:,2);
    gauss2 = FinalFit(:,3);
    gauss3 = FinalFit(:,4);
    gauss_total = FinalFit(:,5);

    Plot3Gauss  % Call function to plot 3-Gaussian fit
end


%% Functions
% Define helper functions used in the script

% Function to calculate ranged MSD based on Gaussian fitting
function [rangedLagtime, rangedMSD, range] = RangedMSD(Data, Msd, LagTime, Ngauss, meanOut, stdOut, nn1, nn2)
    if Ngauss==2
        % Parameters for 2-Gaussian fitting
        mean1 = meanOut(1);
        mean2 = meanOut(2);
        std1 = stdOut(1);
        std2 = stdOut(2);

        % Define range for bound state (first Gaussian)
        range1B = mean1 - std1/2;
        range2B = mean1 + std1/2;
        [rangedLagtime1, rangedMSD1] = FilteredMSD(Data, range1B, range2B, Msd, LagTime, nn1, nn2);

        % Define range for unbound state (second Gaussian)
        range1U = mean2 - std2/2;
        range2U = mean2 + std2/2;
        [rangedLagtime2, rangedMSD2] = FilteredMSD(Data, range1U, range2U, Msd, LagTime, nn1, nn2);

        % Combine results for bound and unbound states
        rangedLagtime = [rangedLagtime1', rangedLagtime2'];
        rangedMSD = [rangedMSD1', rangedMSD2'];
        range = [range1B, range2B, range1U, range2U ];

    elseif Ngauss==3
        % Parameters for 3-Gaussian fitting
        mean1 = meanOut(1);
        mean2 = meanOut(2);
        mean3 = meanOut(3);
        std1 = stdOut(1);
        std2 = stdOut(2);
        std3 = stdOut(2);

        % Define range for bound state (first Gaussian)
        range1B = mean1 - std1;
        range2B = mean1 + std1;
        [rangedLagtime1, rangedMSD1] = FilteredMSD(Data, range1B, range2B, Msd, LagTime, nn1, nn2);

        % Define range for intermediate state (second Gaussian)
        range1I = mean2 - std2;
        range2I = mean2 + std2;
        [rangedLagtime2, rangedMSD2] = FilteredMSD(Data, range1I, range2I, Msd, LagTime, nn1, nn2);

        % Define range for unbound state (third Gaussian)
        range1U = mean3 - std3;
        range2U = mean3 + std3;
        [rangedLagtime3, rangedMSD3] = FilteredMSD(Data, range1U, range2U, Msd, LagTime, nn1, nn2);

        % Combine results for all three states
        rangedLagtime = [rangedLagtime1', rangedLagtime2', rangedLagtime3'];
        rangedMSD = [rangedMSD1', rangedMSD2', rangedMSD3'];
        range = [range1B, range2B, range1I, range2I,  range1U, range2U ];
    end
end

% Function to filter MSD values based on defined range
function [rangedLagtime, rangedMSD] = FilteredMSD(Data, range1, range2, Msd, LagTime, nn1, nn2)
    % Find indices within the specified range
    indices = find(Data >= range1 & Data <= range2);
    filteredMSD = Msd(indices,:);
    filteredLagTime = LagTime(indices,:);
    
    % Calculate average MSD and lag time
    for i = 1:length(filteredLagTime(1,:))
        tempLT = filteredLagTime(:,i);
        tempMSD = filteredMSD(:,i);

        tempNZ = tempLT(tempLT~=0);
        tempMSDNZ = tempMSD(tempMSD~=0);

        avgLT(i) = mean(tempNZ);
        avgMSD(i) = mean(tempMSDNZ);
    end
    1;  % Starting data index
n2 = 3;  % End data index

MinDataPoints = 5;  % Minimum data points required for analysis
RsqTest = 0.8;  % 
    % Remove NaN values and store results in specified range
    avgLT = avgLT(~isnan(avgLT));
    avgMSD = avgMSD(~isnan(avgLT));
    rangedLagtime = avgLT(nn1:nn2);
    rangedMSD = avgMSD(nn1:nn2);
end

% Function to calculate histogram data
function [binPoints, histPoints] = get_Histogram(value,nbin, normalization)
    % Calculate histogram counts and bin edges
    [hcount, hedges] = histcounts(value, nbin,'Normalization',normalization );
    binPoints = (hedges(1:end-1) + hedges(2:end)) / 2.0;  % Calculate bin centers
    histPoints = hcount;  % Store histogram counts
    scale = sum(hcount);  % Normalize histogram to match pdf/probability
    histPoints = histPoints/scale;  % Scale histogram points
end

% Function to calculate diffusion coefficients
function [D,Msd,LagTime,positiveD] = get_D_values(data,inputParams)
    % Extract input parameters
    n1 = inputParams(1);
    n2 = inputParams(2);
    MinDataPoints = inputParams(3);
    RsqTest = inputParams(4);
    pixel_nm = inputParams(5);
    frameRate = inputParams(6);

    % Initialize variables
    positiveD = 0;
    negativeD = 0;
    numTrajectories = data.Trajectory(end);
    D = zeros(1, numTrajectories);

    i = 1;
    for n = 1:numTrajectories
        Data_points = sum(data.Trajectory == n);
        if Data_points >= MinDataPoints  % Filter trajectories based on minimum data points
            [MsdTime_AvgFinal, info] = msd(data,n,n1,n2,frameRate,pixel_nm);
            if info(1) > 0 && info(2) >= RsqTest  % Check if diffusion coefficient is positive and R-squared is sufficient
                D(i) = info(1);
                nn = length(MsdTime_AvgFinal(:,2));
                Msd(i,1:nn) =  MsdTime_AvgFinal(:,2)';
                LagTime(i,1:nn) = MsdTime_AvgFinal(:,1)';
                positiveD = positiveD + 1;
                i = i + 1;
            elseif info(1) < 0 && info(2) >= RsqTest  % Track negative diffusion coefficients
                negativeD = negativeD + 1;
            end
        end
    end
    D(i:end) = [];  % Remove unused entries
end

% Function to calculate Mean Squared Displacement (MSD) for a trajectory
function [MsdTime_AvgFinal, info] = msd(data,n,n1,n2,frameRate,pixel_nm )
    % Extract x, y coordinates, and frame numbers for a specific trajectory
    x = data.x(data.Trajectory == n);
    y = data.y(data.Trajectory == n);
    Frames = data.Frame(data.Trajectory == n);

    N = length(x);  % Number of data points in the trajectory

    deltaT_values = 1:N-1;  % Define range of time intervals (deltaT)

    msdT = zeros(length(deltaT_values), 1);  % Initialize MSD array

    lagtime_storage = zeros(length(deltaT_values),N);  % Initialize lag time storage
    msd_storage = zeros(length(deltaT_values),N);  % Initialize MSD storage

    for deltaT_index = 1:length(deltaT_values)
        deltaT = deltaT_values(deltaT_index);

        msd = zeros(N - deltaT, 1);  % Initialize MSD for each deltaT
        timeD = zeros(N - deltaT, 1);  % Initialize time differences for each deltaT
        for i = 1:(N - deltaT)
            % Find frames that are deltaT apart
            frame1 = Frames(i);
            frame2 = Frames(i + deltaT);

            % Calculate time difference between frames
            time_diff = abs(frame2 - frame1);

            % Calculate squared displacement using the time difference
            dx = x(i + deltaT) - x(i);
            dy = y(i + deltaT) - y(i);
            msd(i) = (dx^2 + dy^2);
            timeD(i) = time_diff;
        end

        % Store lag time and MSD for each deltaT
        lagtime_storage(deltaT,1:length(timeD)) = timeD;
        msd_storage(deltaT,1:length(msd)) = msd;
    end


    % Sort MSD and lag time based on lag time
    k_lim1 = max(lagtime_storage,[],'all');
    msdlagtime_sorted = zeros(k_lim1, N);
    check(:, 1) = (1:k_lim1)';
    for k1 = 1:k_lim1
        i1 = 1;
        for k = 1:size(lagtime_storage, 1)
            for i = 1:size(lagtime_storage, 2)
                if (lagtime_storage(k, i) == check(k1, 1))
                    msdlagtime_sorted(k1, i1) = msd_storage(k, i);
                    i1 = i1 + 1;
                end
            end
        end
    end

    % Calculate average MSD for each lag time
    for i=1:k_lim1
        MsdLagtime_nonZero = msdlagtime_sorted(i,msdlagtime_sorted(i,:)~=0);
        MsdTime_AvgFinal(i,2) = mean(MsdLagtime_nonZero)*(pixel_nm*10^-3)^2;
    end

    % Calculate lag time in seconds
    for i=1:k_lim1
        MsdTime_AvgFinal(i,1)=(i)*frameRate*10^(-3);
    end

    % Perform linear fitting on MSD vs. lag time
    [p, S] = polyfit(MsdTime_AvgFinal(n1:n2,1),MsdTime_AvgFinal(n1:n2,2),1);

  info(1) = (p(1)-4*0.0071^2)/4.0; % This is the fitted diffusion coefficient with localization error
 %info(1) = (p(1))/4.0; % This is the fitted diffusion coefficient
    info(2) = 1 - (S.normr/norm(MsdTime_AvgFinal(n1:n2,2) ...
        - mean(MsdTime_AvgFinal(n1:n2,2))))^2; % This is the R^2 value
end

function [MsdAll,TimeAll,MsdTime_AvgFinal] = FullMsd(data,inputParams)
n1 = inputParams(1);
    n2 = inputParams(2);
    MinDataPoints = inputParams(3);
    RsqTest = inputParams(4);
    pixel_nm = inputParams(5);
    frameRate = inputParams(6);
    numTrajectories = data.Trajectory(end);
    
    i = 1;
    for n = 1:numTrajectories
        Data_points = sum(data.Trajectory == n);
        if Data_points >= MinDataPoints  % Filter trajectories based on minimum data points
            [MsdTime_AvgFinal, info] = msd(data,n,n1,n2,frameRate,pixel_nm);
        end
        
        MsdAll(n,n1:n2)=MsdTime_AvgFinal(n1:n2,2)';
        TimeAll(n,n1:n2)=MsdTime_AvgFinal(n1:n2,1)';
    end
  
end