


indices = find(xgrid >= range(1) & xgrid <= range(2));

xbound = xgrid(indices);
populationB = gauss1(indices);

pop1 = trapz(xbound,populationB);

indices = find(xgrid >= range(3) & xgrid <= range(4));

xUnbound = xgrid(indices);
populationUB = gauss2(indices);

pop2 = trapz(xUnbound, populationUB);


popt = pop1 + pop2;

pop1/popt
pop2/popt