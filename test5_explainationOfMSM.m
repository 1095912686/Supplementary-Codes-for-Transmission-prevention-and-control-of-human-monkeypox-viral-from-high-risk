clear; close all; clc;

n = 3;
modelSEIR = dynamicalModel_SEIRn;
modelSEIR.n = n;
modelSEIR.N = [3970238390 * 0.04;...
    3970238390 * 0.96;...
    3904727342]; % default global population size
modelSEIR.C = readmatrix('globalContactMatrix.xlsx');
modelSEIR.omega = 1/8.5;
modelSEIR.gamma = 1/11;
modelSEIR.VE = 0;



whos modelSEIR

idx = 1:n;
x0 = zeros(4*n, 1);
x0(idx + 0*n) = modelSEIR.N;
x0(idx + 2*n) = 1; % first case assigned proportional to population size
modelSEIR.xInit = x0;


%% setup (r1,r2) for simulation
r = [1e-2, 1e-3, 1e-4 1e-5];
Reff = [2, 3, 4];% [2, 3, 4];
record = zeros(numel(Reff) * numel(r)^2, 6);
temp = 1;

%% simulate for different Reff
for k = 1:numel(Reff)

    % initialize figure
    % curves
    fig1 = figure;
    fig1.WindowState = 'maximized';
    tile1 = tiledlayout(4,4, 'TileSpacing', 'tight', 'Padding', 'tight');

    % contours
    fig2 = figure;
    fig2.WindowState = 'maximized';
    tile2 = tiledlayout(2,3, 'TileSpacing', 'compact');

    % setup initial model
    modelSEIR.Reff = Reff(k);
    modelSEIR.q = Reff(k) * modelSEIR.gamma / max(eig(modelSEIR.C .* (1 - modelSEIR.VE)));


    %% Simulate with (r1,r2) and Plot
    M = zeros(numel(r));

    % ratio of contact decomposition
    for i = 1:numel(r)
        r1 = r(i);


        % ratio of relative risks
        for j = 1:numel(r)
            r2 = r(j);

            % setup model2 of given intervention
            model2 = modelSEIR;

            % r1: ratio of high/low risk contact
            % r2: ratio of relative risk
            model2.C = modelSEIR.C * (r1 +(1-r1) * r2);
            model2.C(1,1) = 0.22 * 4 + 0.08 * r2;
            model2.q = Reff(k) * model2.gamma / max(eig(model2.C .* (1-model2.VE)));

            % simulates
            tSpan = [0, 1000];
            [t, x] = predictModel(model2, tSpan, model2.xInit);

            dI.value = model2.omega * x(:,idx + 1*n);
            dI.lb = icdf('Poisson', 0.05, dI.value);
            dI.ub = icdf('Poisson', 0.95, dI.value);

            % visualize
            t1 = nexttile(tile1);
            plot(t1, t, dI.value); hold on;
            title(t1, "r1 = " + r1 + ", r2 = " + r2, 'FontWeight','normal');
            legend(t1, "group" + (1:3));



            % save record (incidence rate)
            tSpan2 = [0, 200];
            [t2, x2] = predictModel(model2, tSpan2, model2.xInit);


            record(temp, 1) = tSpan2(end);
            record(temp, 2) = r1;
            record(temp, 3) = r2;
            record(temp, 4) = Reff(k);
            record(temp, 5:7) = trapz(t2, model2.omega * x2(:,idx + 1*n)) ./ model2.N';
            temp = temp + 1;

        end
    end

    title(tile1, "\fontname{Times New Roman}\it R_{eff}\rm = " + Reff(k), 'FontWeight','normal');
    xlabel(tile1, "Time (in days)");
    ylabel(tile1, "Daily New Incidence Rate (person per day)");


    %% Contour plot
    plotTitle = {'Ratio of total number of cases in MSM and general population', 'Total number of cases in non-MSM males', 'Total number of cases in females',...
                 'Ratio of cumulative incidence rate between MSM and general population', 'Cumulative incidence rate of non-MSM males', 'Cumulative incidence rate of females'};
    for s = 1:6
        rr = 10.^(linspace(-2, -5, 20));
        for i = 1:numel(rr)
            rr1 = rr(i);

            for j = 1:numel(rr)
                rr2 = rr(j);

                % setup model2 of given intervention
                model2 = modelSEIR;

                % rr1: ratio of high/low risk contact
                % rr2: ratio of relative risk
                model2.C = modelSEIR.C * (rr1 +(1-rr1) * rr2);
                model2.C(1,1) = 0.22 * 4 + 0.08 * r2;
                model2.q = Reff(k) * model2.gamma / max(eig(model2.C .* (1-model2.VE)));

                % simulates
                tSpan = [0, 200];
                [t3, x3] = predictModel(model2, tSpan, model2.xInit);


                % for contour plot
                tSpan3 = [0, 1000];
                [t3, x3] = predictModel(model2, tSpan3, model2.xInit);
                cumMSM = trapz(t3, model2.omega * x3(:,n+1));
                cumMF = sum(trapz(t3, model2.omega.' .* x3(:,n+2 : n+3), 1), 2);
                cumM = trapz(t3, model2.omega.' .* x3(:,n+2), 1);
                cumF = trapz(t3, model2.omega.' .* x3(:,n+3), 1);
                switch s
                    case 1
                        M(i,j) = cumMSM/cumMF;
                    case 2
                        M(i,j) = cumM;
                    case 3
                        M(i,j) = cumF;
                    case 4
                        M(i,j) = (cumMSM / modelSEIR.N(1)) / (cumMF / sum(modelSEIR.N(2:3)));
                    case 5
                        M(i,j) = cumM / modelSEIR.N(2);
                    case 6
                        M(i,j) = cumF / modelSEIR.N(3);
                end

            end
        end

        ax = nexttile(tile2);
        x = linspace(1, numel(r), numel(rr));
        y = linspace(1, numel(r), numel(rr));
        labels = {'', ''};
        if s == 5
                labels{1} = "\fontname{Times New Roman}Contact degree ratios \it r_1 \rm of high-risk and low-risk behaviors";
                %labels{2} = "\fontname{Times New Roman}Relative risk ratios \it r_2 \rm of low-risk behaviors";
        end


        logScale3D(ax, x, y, M', plotTitle{s}, labels, 'contour');
    
        ax = gca;
        ax.XTick = 1:numel(r);
        ax.YTick = 1:numel(r);
        ax.XTickLabel = string(r);
        ax.YTickLabel = string(r);
        ax.Colormap = flip(ax.Colormap);
        ax.Children.LabelSpacing = 300;
        text(-0.1, 1.1, string(char(s+64)), 'Unit', 'Normalized', 'FontSize', 13, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
        if s == 1
            ytext = ylabel(ax, "\fontname{Times New Roman}Relative risk ratios \it r_2 \rm of low-risk behaviors");
            ytext.Position = [0.5, 1];
            ytext.Units = 'normalized';
            ytext.FontSize = 14;
        end
    end

    %temp = get(ax, 'colormap');
    %set(ax, 'colormap', flip(temp));


    title(tile2, "\fontname{Times New Roman}\it R_{eff}\rm = " + Reff(k));

    exportgraphics(fig1, "contactDecompositionCurvesWithReff" + Reff(k) + ".jpg", Resolution=300);
    exportgraphics(fig2, "contactDecompositionContourWithReff" + Reff(k) + ".jpg", 'Resolution', 300);
end



recordTable = array2table(record, 'VariableNames', {'timeSpan', 'r_1', 'r_2', 'Reff', 'MSM', 'Male', 'Female'});
writetable(recordTable, 'incidenceMaleFemaleMSM.xlsx');

