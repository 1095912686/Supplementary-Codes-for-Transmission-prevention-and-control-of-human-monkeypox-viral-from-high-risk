function [ax] = logScale3D(ax, x, y, M, titleName, labelNames, plotType)

[X, Y] = meshgrid(x, y);

if strcmp(plotType, 'surf')
    ax = surf(X, Y, M);
else %if strcmp(plotType, 'contour')
    x_refined = linspace(x(1), x(end), 1e3);
    y_refined = linspace(y(1), y(end), 1e3);
    [XX, YY] = meshgrid(x_refined, y_refined);
    MM = interp2(X,Y,M,XX,YY);
    [~, temp] = contour(XX,YY,MM); % default
    m1 = floor( log(min(MM, [], 'all')) / log(2) );
    if ~isfinite(m1)
        m1 = 0;
    end
    m2 = floor( log(max(MM, [], 'all')) / log(2) );
    if ~isfinite(m2)
        m2 = 0;
    end

    format short;
    contourf(XX, YY, MM,...
                'LevelList', 2.^(linspace(m1, m2, m2-m1+1)),...
                'TextList', 2.^(linspace(m1, m2-1, m2-m1)),...
                'ShowText','on');
    %contourf(XX, YY, MM, 'LevelList', temp.LevelList, 'ShowText','on');
    set(ax,'ColorScale','log')
    color = cbrewer2('RdBu', 256);
    colormap(ax, color(1:end-30, :));
    %temp = get(ax, 'colormap');
    %set(ax, 'colormap', flip(temp));
    axis equal
end

c = colorbar;

xlabel(labelNames{1}, 'FontSize', 14);
ylabel(labelNames{2}, 'FontSize', 14);
xticks(x);
yticks(y);
xticklabels(string(x*100) + "%");
yticklabels(string(y*100) + "%");
title(titleName, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
set(gca, 'FontName', 'Times New Roman');

end