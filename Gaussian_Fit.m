function [meanOut, stdOut, AmpOut, fBound, FinalFit] = Gaussian_Fit(Data, Ngauss,nbin,...
                                                normalization,  MaxIter, MaxFunEvals)

[hcount, ~] = histcounts(Data, nbin, 'Normalization', normalization); %%Normalizatin pdf
scale = sum(hcount);


if Ngauss==2
    
    pdf_normmixture = @(Data,p,mu1,mu2,sigma1,sigma2) ...
        p*normpdf(Data,mu1,sigma1) + (1-p)*normpdf(Data,mu2,sigma2);
    pStart = .5;
    muStart = quantile(Data,[.25 .75]);
    sigmaStart = sqrt(var(Data) - 0.25*diff(muStart).^2);
    start = [pStart muStart sigmaStart sigmaStart];
    lb = [0 -Inf -Inf 0 0];
    ub = [1 Inf Inf Inf Inf];
    
    options = statset('MaxIter',MaxIter,'MaxFunEvals',MaxFunEvals);
    paramEsts = mle(Data,normalization,pdf_normmixture,'Start',start, ...
        'LowerBound',lb,'UpperBound',ub,'Options',options);
    
    %xgrid = linspace(1.2*min(Data),1.5*max(Data),500);
    xgrid = linspace(-3.5,2,500);
    
    
    mean1 = paramEsts(2);
    std1 = paramEsts(4);
    FitAmp(1)=(paramEsts(1)/(sqrt(2*pi)*std1))/scale;
    gauss1 = FitAmp(1)* exp(-(xgrid - mean1).^2 / (2*std1^2));
    

    
    %2nd gaussian
    mean2 = paramEsts(3);
    std2 = paramEsts(5);
    FitAmp(2)=((1-(paramEsts(1)))/(sqrt(2*pi)*std2))/scale;
    gauss2 = FitAmp(2) * exp(-(xgrid - mean2).^2 / (2*std2^2));
    %

    meanOut = [mean1,mean2];
    stdOut = [std1, std2];
    AmpOut = [FitAmp(1), FitAmp(2)];
    
    fBound = [paramEsts(1), (1-paramEsts(1))];
    FinalFit = [xgrid', gauss1', gauss2', (gauss1 + gauss2)'];

else
    %% function handler for fitting gaussians
    pdf_normmixture = @(Data, p1,p2, mu1, mu2, mu3, sigma1, sigma2, sigma3) ...
        p1*normpdf(Data, mu1, sigma1) + p2*normpdf(Data, mu2, sigma2) + (1-(p1+p2))*normpdf(Data, mu3, sigma3);
    
    pStart = [1/3, 1/3];
    %pStart=[1/4 1/2]
    muStart = quantile(Data, [1/6, 3/6, 5/6]);
    sigmaStart = sqrt(var(Data) - .25 * diff(muStart).^2);
    start = [pStart muStart sigmaStart sigmaStart(1)];% sigmaStart];
    lb = [0 0 -Inf -Inf -Inf 0 0 0];
    ub = [1 1 Inf Inf Inf Inf Inf Inf];
    
    options = statset('MaxIter', MaxIter , 'MaxFunEvals',  MaxFunEvals);
    
    paramEsts = mle(Data, normalization , pdf_normmixture, 'Start', start, ...
        'LowerBound', lb, 'UpperBound', ub, 'Options', options);
    
    
    %histogram(Data, nbin, 'Normalization', 'probability');
    %hold on;


    %xgrid = linspace(1.1*min(Data), 1.1*max(Data), 200);
    xgrid = linspace(-3.5, 2.0, 200);
 
    % 1st gaussian
    mean1 = paramEsts(3);
    std1 = paramEsts(6);
    FitAmp(1)=(paramEsts(1)/(sqrt(2*pi)*std1))/scale;
    gauss1 = FitAmp(1)* exp(-(xgrid - mean1).^2 / (2*std1^2));


    
    %2nd gaussian
    mean2 = paramEsts(4);
    std2 = paramEsts(7);
    FitAmp(2)=(paramEsts(2)/(sqrt(2*pi)*std2))/scale;
    gauss2 = FitAmp(2) * exp(-(xgrid - mean2).^2 / (2*std2^2));
   
    %3rd gaussian
    
    mean3 = paramEsts(5);
    std3 = paramEsts(8);
    FitAmp(3)=((1-(paramEsts(1)+paramEsts(2)))/(sqrt(2*pi)*std3))/scale;
    gauss3 = FitAmp(3) * exp(-(xgrid - mean3).^2 / (2*std3^2));

    meanOut = [mean1,mean2, mean3];
    stdOut = [std1, std2, std3];
    AmpOut = [FitAmp(1), FitAmp(2), FitAmp(3)];
    
    fBound = [paramEsts(1), paramEsts(2), 1-(paramEsts(1)+paramEsts(2))];
    FinalFit = [xgrid', gauss1', gauss2',gauss3', (gauss1 + gauss2 + gauss3)'];

end
end