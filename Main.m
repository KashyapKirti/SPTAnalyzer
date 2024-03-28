clc; close all; clear;

%% Load Data

%data = readtable("H2B_10ms_REP1.csv");
%data = readtable("../WT_MED14_10ms_REP3.csv");
%data = readtable("../WT_SUA7_10ms_REP1.csv");
%data = readtable("../WT_TFA1_10ms_REP3.csv");


%data = readtable("TAF1AA_SUA7_10ms_DMSO.csv");
%data = readtable("RPB1AA_TBP_10ms_DMSO.csv");
%data = readtable("TAF1AA_TOA1_10ms_DMSO.csv");
%data = readtable("MED14AA_RPB1_10ms_RAPAMYCIN_REP2 (1).csv");

%data = readtable("../NewData/mage1_Anaphase_15ms.csv");
%data = readtable("../NewData/mage1_Metaphase_15ms.csv");
%data = readtable("../NewData/Atg11 at cytosol csv.csv");
%data = readtable("../NewData/Atg11 at SPB CSV.csv");
%data = readtable("../_ch00_23-09-28_15-55-51.csv");
data = readtable("../_ch00_23-10-31_18-13-55.csv");

%% Data Filteration

n1 = 2;  % starting data index
n2 = 5;  % end data index


MinDataPoints = 9;
RsqTest = 0.8;

pixel_nm = 110;
frameRate = 15; % in milisecond.

nbin = 10; 


%%

inputParams = [n1,n2,MinDataPoints,RsqTest,pixel_nm, frameRate];
[D,Msd,LagTime,positiveD] = get_D_values(data, inputParams);

%% Gaussian Fitting

logD = log10(D);
normalization = 'pdf';

binsize = ((max(logD)-min(logD))/nbin);


Ngauss = 2;

Data = logD;
MaxIter = 10000;
MaxFunEvals = 2*MaxIter;


[binPoints, histPoints] = get_Histogram(logD, nbin, normalization );
[meanOut, stdOut, AmpOut, fBound,  FinalFit] = Gaussian_Fit(Data, Ngauss, nbin, ...
    normalization,  MaxIter,MaxFunEvals );



%%-- Calculating MSD for bound and unbound
nn1 = 1; nn2 = 6; 

[rangedLagtime, rangedMSD, range] = RangedMSD(Data, Msd,LagTime, Ngauss, meanOut, stdOut, nn1, nn2);



%% Plotting



if Ngauss==2
    xgrid = FinalFit(:,1);
    gauss1 = FinalFit(:,2);
    gauss2 = FinalFit(:,3);
    gauss_total = FinalFit(:,4);

    Plot2Gauss
    

end
if Ngauss==3
    xgrid = FinalFit(:,1);
    gauss1 = FinalFit(:,2);
    gauss2 = FinalFit(:,3);
    gauss3 = FinalFit(:,4);
    gauss_total = FinalFit(:,5);

    Plot3Gauss
end










%% Functions


function [rangedLagtime, rangedMSD, range] = RangedMSD(Data, Msd, LagTime, Ngauss, meanOut, stdOut, nn1, nn2)
if Ngauss==2
    mean1 = meanOut(1);
    mean2 = meanOut(2);
    std1 = stdOut(1);
    std2 = stdOut(2);

    range1B = mean1 - std1/2;
    range2B = mean1 + std1/2;
   
    
    [rangedLagtime1, rangedMSD1] = FilteredMSD(Data, range1B, range2B, Msd, LagTime, nn1, nn2);


    range1U = mean2 - std2/2;
    range2U = mean2 + std2/2;


    [rangedLagtime2, rangedMSD2] = FilteredMSD(Data, range1U, range2U,...
                                                    Msd, LagTime, nn1, nn2);


    rangedLagtime = [rangedLagtime1', rangedLagtime2'];
    rangedMSD = [rangedMSD1', rangedMSD2'];
    range = [range1B, range2B, range1U, range2U ];


elseif Ngauss==3
    mean1 = meanOut(1);
    mean2 = meanOut(2);
    mean3 = meanOut(3);
    std1 = stdOut(1);
    std2 = stdOut(2);
    std3 = stdOut(2);

    range1B = mean1 - std1;
    range2B = mean1 + std1;

    
    [rangedLagtime1, rangedMSD1] = FilteredMSD(Data, range1B, range2B, Msd, LagTime, nn1, nn2);


    range1I = mean2 - std2;
    range2I = mean2 + std2;

    [rangedLagtime2, rangedMSD2] = FilteredMSD(Data, range1I, range2I,...
                                                    Msd, LagTime, nn1, nn2);

    range1U = mean3 - std3;
    range2U = mean3 + std3;

    [rangedLagtime3, rangedMSD3] = FilteredMSD(Data, range1U, range2U,...
                                                    Msd, LagTime, nn1, nn2);


    rangedLagtime = [rangedLagtime1', rangedLagtime2', rangedLagtime3'];
    rangedMSD = [rangedMSD1', rangedMSD2', rangedMSD3'];
    
    range = [range1B, range2B, range1I, range2I,  range1U, range2U ];
 


end
end
function [rangedLagtime, rangedMSD] = FilteredMSD(Data, range1, range2,...
                                                    Msd, LagTime, nn1, nn2)

%filteredData = Data(Data >= range1 & Data <= range2);

indices = find(Data >= range1 & Data <= range2);
filteredMSD = Msd(indices,:);
filteredLagTime = LagTime(indices,:);
for i = 1:length(filteredLagTime(1,:))

    tempLT = filteredLagTime(:,i);
    tempMSD = filteredMSD(:,i);

    tempNZ = tempLT(tempLT~=0);
    tempMSDNZ = tempMSD(tempMSD~=0);

    avgLT(i) = mean(tempNZ);
    avgMSD(i) = mean(tempMSDNZ);
end
avgLT = avgLT(~isnan(avgLT));
avgMSD = avgMSD(~isnan(avgLT));

rangedLagtime = avgLT(nn1:nn2);
rangedMSD = avgMSD(nn1:nn2);

end



function [binPoints, histPoints] = get_Histogram(value,nbin, normalization)

[hcount, hedges] = histcounts(value, nbin,'Normalization',normalization );
binPoints = (hedges(1:end-1) + hedges(2:end)) / 2.0;
histPoints = hcount;
scale = sum(hcount);
histPoints = histPoints/scale; % Scaling to match pdf and probability.
end


function [D,Msd,LagTime,positiveD] = get_D_values(data,inputParams)

n1 = inputParams(1);
n2 = inputParams(2);
MinDataPoints = inputParams(3);
RsqTest = inputParams(4);
pixel_nm = inputParams(5);
frameRate = inputParams(6);





positiveD = 0;
negativeD = 0;
numTrajectories = data.Trajectory(end);
D = zeros(1, numTrajectories);

i = 1;
for n = 1:numTrajectories
    Data_points = sum(data.Trajectory == n);
    if Data_points >= MinDataPoints  % For trajectories with at least 6 data points means 5 displacements
        [MsdTime_AvgFinal, info] = msd(data,n,n1,n2,frameRate,pixel_nm);
        if info(1) > 0 && info(2) >= RsqTest
            D(i) = info(1);
            nn = length(MsdTime_AvgFinal(:,2));
            Msd(i,1:nn) =  MsdTime_AvgFinal(:,2)';
            LagTime(i,1:nn) = MsdTime_AvgFinal(:,1)';
            positiveD = positiveD + 1;
            i = i + 1;
        elseif info(1) < 0 && info(2) >= RsqTest
            negativeD = negativeD + 1;
        end
    end
end

D(i:end) = [];


end



function [MsdTime_AvgFinal, info] = msd(data,n,n1,n2,frameRate,pixel_nm ) %Edit sahil
% Extract the x, y coordinates, and frames for a specific trajectory (e.g., Trajectory == 1)
x = data.x(data.Trajectory == n);
y = data.y(data.Trajectory == n);
Frames = data.Frame(data.Trajectory == n);

N = length(x);

deltaT_values = 1:N-1;  % Modify as needed


msdT = zeros(length(deltaT_values), 1);

lagtime_storage = zeros(length(deltaT_values),N);
msd_storage = zeros(length(deltaT_values),N);



for deltaT_index = 1:length(deltaT_values)
    deltaT = deltaT_values(deltaT_index);


    msd = zeros(N - deltaT, 1);
    timeD = zeros(N - deltaT, 1);
    for i = 1:(N - deltaT)
        % Find the indices of the frames that are deltaT apart
        frame1 = Frames(i);
        frame2 = Frames(i + deltaT);

        % Calculate the time difference between frames
        time_diff = abs(frame2 - frame1);

        % Calculate squared displacement using the time difference
        dx = x(i + deltaT) - x(i);
        dy = y(i + deltaT) - y(i);
        msd(i) = (dx^2 + dy^2) ;
        timeD(i) = time_diff;
    end

    lagtime_storage(deltaT,1:length(timeD)) = timeD;
    msd_storage(deltaT,1:length(msd)) = msd;
end

%%-- Modified by Kirti

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



for i=1:k_lim1
    MsdLagtime_nonZero = msdlagtime_sorted(i,msdlagtime_sorted(i,:)~=0);
    MsdTime_AvgFinal(i,2) = mean(MsdLagtime_nonZero)*(pixel_nm*10^-3)^2;
end


for i=1:k_lim1
    MsdTime_AvgFinal(i,1)=(i)*frameRate*10^(-3);
end


[p, S] = polyfit(MsdTime_AvgFinal(n1:n2,1),MsdTime_AvgFinal(n1:n2,2),1);

info(1) = p(1)/4.0; %This is the fitted D
info(2) = 1 - (S.normr/norm(MsdTime_AvgFinal(n1:n2,2) ...
    - mean(MsdTime_AvgFinal(n1:n2,2))))^2; % This is the R^2 value.

end



