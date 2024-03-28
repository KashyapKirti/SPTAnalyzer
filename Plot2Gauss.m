
color1 = [0.4940 0.1840 0.5560]; %ansphase main
color2 = [0.8500 0.3250 0.0980 0.5]	; %anaphase 1
color3 = [0.4660 0.6740 0.1880 0.5]; %anaphase 2
color4 = [1 0 1];
grey = [0.5 0.5 0.5];
color5 = [0.6350 0.0780 0.1840];
bl = [0 0 1 0.5]; %metaphase 1
pl = [1 0 1 0.5]; %metaphase 2
rl = [1 0 0 0.6]; %metaphase main
color6 = [21/255 20/255 216/255 ]; %cystol 1
color7 = [140/255 57/255 69/255 0.5 ]; %cystol 2
color8 = [108/255 155/255 3/255 0.5 ]; %cystol 3
color9 = [0/255 163/255 96/255]; %SPB Main
color10 = [69/255 76/255 240/255 0.5]; %SPB 1
color11 = [100/255 18/255 91/255 0.5]; % SPB 2



bar(binPoints, histPoints * 100, 'BarWidth', 1, 'FaceAlpha', 0.0001, 'EdgeColor', 'k'); % Adjust BarWidth and FaceAlpha
hold on;

plot(xgrid, gauss1 * 100, 'LineWidth', 2.5, 'Color', color2);
hold on;

plot(xgrid, gauss2 * 100, 'LineWidth', 2.5, 'Color',color3);
hold on;

plot(xgrid, gauss_total * 100, 'LineWidth', 3, 'Color',color6);
hold off;



xlabel("log_{10}(D) (\mu m^2 /s)")
ylabel("Frequency(%)")
%axis square
set(gca, 'FontSize', 28);
axis([-3.5 2 0 29])


str1 = strcat("Fbound = ", num2str(fBound(1)*100), "%");
str2 = strcat("Ffree = ", num2str(fBound(2)*100), "%");

annotation('textbox', [0.72, 0.80, 0.2, 0.1], 'String', 'Data: ',...
            'FontSize', 22, 'EdgeColor', 'none');
annotation('textbox', [0.72, 0.65, 0.2, 0.1], 'String',str1, 'FontSize', 22,...
            'EdgeColor', 'none','Color', color2);
annotation('textbox', [0.72, 0.6, 0.2, 0.1], 'String',str2, 'FontSize', 22,...
            'EdgeColor', 'none','Color', color3);

str1 = strcat("Dbound = ", num2str(10^(meanOut(1))));
str2 = strcat("Dfree = ", num2str(10^(meanOut(2))));
annotation('textbox', [0.72, 0.55, 0.2, 0.1], 'String',str1, 'FontSize', 22,...
            'EdgeColor', 'none','Color', color2);
annotation('textbox', [0.72, 0.5, 0.2, 0.1], 'String',str2, 'FontSize', 22,...
            'EdgeColor', 'none','Color', color3, "FontName",'Times');



