


aucT = diff(gauss_total).*diff(xgrid);

g1 = diff(gauss1).*diff(xgrid);
g2 = diff(gauss2).*diff(xgrid);


(sum(g1) + sum(g2));
ss = sum(aucT);


mean1 = meanOut(1);
mean2 = meanOut(2);

std1 = stdOut(1);
std2 = stdOut(2);



y1 = exp(-(xgrid - mean1).^2 / (2*std1^2));
y2 = exp(-(xgrid - mean2).^2 / (2*std2^2));


trapz(xgrid, gauss1)/trapz(xgrid, gauss_total)
trapz(xgrid, gauss2)/trapz(xgrid, gauss_total)

