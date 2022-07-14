clear; close all; clc;

% n = 3;
% modelSEIR = dynamicalModel_SEIRn;
% modelSEIR.n = n;
% modelSEIR.N = [4030711645 * 0.04;...
%     4030711645 * 0.96;...
%     3957932033]; % default global population size
% modelSEIR.C = readmatrix('globalContactMatrix.xlsx');
% modelSEIR.omega = 1/12;
% modelSEIR.gamma = 1/11;
% modelSEIR.VE = 0;

% n = 1;
% modelSEIR = dynamicalModel_SEIRn;
% modelSEIR.n = n;
% modelSEIR.N = 4030711645 + 3957932033;
% modelSEIR.C = zeros(n); modelSEIR.C(1) = 13.45;
% modelSEIR.omega = 1/12;
% modelSEIR.gamma = 1/11;
% modelSEIR.VE = 0;

n = 2;
modelSEIR = dynamicalModel_SEIRn;
modelSEIR.n = n;
VC = 0.8 * 0.03; % current vaccination coverage
modelSEIR.N = 1.588e8 * [(1-VC); VC];
modelSEIR.C = repmat(13.45 * [(1-VC), VC], [2,1]);
modelSEIR.omega = 1/8.5;
modelSEIR.gamma = 1/11;
modelSEIR.VE = [0, 0.85] .* ones(2,1);


whos modelSEIR

idx = 1:n;
x0 = zeros(4*n, 1);
x0(idx + 0*n) = modelSEIR.N;
x0(idx + 2*n) = 1; %modelSEIR.N / sum(modelSEIR.N); % first case assigned proportional to population size
modelSEIR.xInit = x0;

Reff = [2, 3, 4];
Record = {};
for k = 1:numel(Reff)


    fig = figure;
    fig.WindowState = 'maximized';
    tile1 = tiledlayout('flow');
    tile1.TileSpacing = 'tight';
    tile1.Padding = 'tight';

    modelSEIR.Reff = Reff(k);
    modelSEIR.q = Reff(k) * modelSEIR.gamma / max(eig(modelSEIR.C .* (1 - modelSEIR.VE)));

    %% predict with (single & combined) intervention

    % interventions
    intensity = 0 : 0.1 : 0.9; % efficiency = 1 - intensity;
    newVC = 0 : 0.1 : 0.9; % numel(newVC) == numel(efficiency)
    m = numel(intensity);

    %interventionNames = {'reducing 1/\gamma', 'reducing Cq', 'vaccine coverage'};
    interventionNames1 = {'Reduction of infectious period', 'Reduction of Cq', 'Vaccine coverage'};
    interventionNames2 = {'\it \gamma = \gamma / \rm', '\it Cq = Cq\rm × ', '\it VC = \rm '};
    
    for i = 1:3  % 3 kinds of interventions
        
        %% single-intervention simulation
        t1 = nexttile;
        for j = 1:m % m = numel(newVC) == numel(intensity)

            % setup model2 of given intervention
            model2 = modelSEIR;
            switch i
                case 1
                    dgamma = 1 - intensity(j);
                    model2.gamma = model2.gamma / dgamma;
                case 2
                    dq = 1 - intensity(j);
                    model2.q = model2.q * dq;
                case 3
                    VC = newVC(j);
                    model2.N = 1.588e8 * [(1-VC); VC];
                    model2.C = repmat(13.45 * [(1-VC), VC], [2,1]);
                    model2.xInit((1:n) + 0*n) = model2.N; 
            end

            

            % simulates
            tSpan = [0, 10e2];
            [t, x] = predictModel(model2, tSpan, model2.xInit);

            dI.value = sum(model2.omega * x(:,(1:n) + 1*n), 2);
            dI.lb = icdf('Poisson', 0.05, dI.value);
            dI.ub = icdf('Poisson', 0.95, dI.value);

            % visualize
            %plot(t1, t, dI.value / sum(model2.N)); hold on;
            plot(t1, t, dI.value); hold on;

            % Record
            weights = ones(n,1);
            tSpan_longPeriod = [0,1000];
            tSpan_shortPeriod = [0,30];
            M_cumI_longPeriod = weightedAccumulativeCases(model2, model2.xInit, tSpan_longPeriod, weights);
            M_cumI_shortPeriod = weightedAccumulativeCases(model2, model2.xInit, tSpan_shortPeriod, weights);
            M_peak = peak(model2, model2.xInit, weights);
            M_duration = last(model2, model2.xInit, weights);
            Recordi = {"Model" + n, Reff(k), interventionNames1{i}, intensity(j), tSpan_longPeriod(end), M_cumI_longPeriod, M_peak, M_duration};
            Record = [Record; Recordi];
        end

        if i == 3
            % invervention of (3) vaccination
            legend(t1, interventionNames2(i) + string((intensity)*100) + "%");
        else
            % intervention of (1) and (2)
            legend(t1, interventionNames2(i) + string((1-intensity)*100) + "%");
        end
        ylabel('\fontname{Times New Roman}Daily new Cases', 'FontSize', 14)
        xlabel('\fontname{Times New Roman}Time (in days)', 'FontSize', 14);
        %t1.YLim = [0, 2.5e6];
        t1.YLim(1) = 0;;

        % set color
        color = cbrewer2('RdBu', 13, 'linear');
        color(6:8, :) = [];
        for j = 1:numel(t1.Children)
            t1.Children(j).Color = color(j, :);
        end

        text(-0.1, 1.05, string(char(i+64)), 'Unit', 'normalized', 'FontSize', 12, 'FontWeight', 'normal');

    end


    %% combined intervention simulation
    mm = 20;
    %intensityRefined = linspace(0, 0.9, mm); % efficiency = 1 - intensity;
    intensityRefined = intensity;
    newVCRefined = intensityRefined;
    for i = 1:3

        M_cumI = zeros(numel(intensityRefined));
        M_duration = zeros(numel(intensityRefined));
        M_peak = zeros(numel(intensityRefined));
        
        for ii = 1:numel(intensityRefined)
            for jj = 1:numel(intensityRefined)

                %idx = (ii-1)*m + jj;

                % initialize model
                model2 = modelSEIR;
                switch i
                    case 1
                    % intervention 23

                        % y-axis in contour plot
                        dq = 1 - intensityRefined(jj);
                        model2.q = model2.q * dq;

                        % x-axis in contour plot
                        VC = newVCRefined(ii);
                        model2.N = 1.588e8 * [(1-VC); VC];
                        model2.C = repmat(13.45 * [(1-VC), VC], [2,1]);
                        model2.xInit((1:n) + 0*n) = model2.N;
                        labels = interventionNames1([2,3]);

                    case 2
                    % intervention 13
                        % y-axis in contour plot
                        dgamma = 1 - intensityRefined(jj);
                        model2.gamma = model2.gamma / dgamma;
                        
                        % x-axis in contour plot
                        VC = newVCRefined(ii);
                        model2.N = 1.588e8 * [(1-VC); VC];
                        model2.C = repmat(13.45 * [(1-VC), VC], [2,1]);
                        model2.xInit((1:n) + 0*n) = model2.N;
                        labels = interventionNames1([1,3]);

                    case 3
                    % intervention 12

                        % base-line for comparison
                        VC = newVCRefined(1);
                        model2.N = 1.588e8 * [(1-VC); VC];
                        model2.C = repmat(13.45 * [(1-VC), VC], [2,1]);
                        model2.xInit((1:n) + 0*n) = model2.N;

                        % y-axis in contour plot
                        dgamma = 1 - intensityRefined(jj);
                        model2.gamma = model2.gamma / dgamma;
                        
                        % x-axis in contour plot
                        dq = 1 - intensityRefined(ii);
                        model2.q = model2.q * dq;
                        labels = interventionNames1([1,2]);

                end

                weights = ones(n,1);
                tSpan = [0,100];
                M_cumI(ii,jj) = weightedAccumulativeCases(model2, model2.xInit, tSpan, weights);
                M_peak(ii,jj) = peak(model2, model2.xInit, weights);
                %M_duration(ii,jj) = last(model2, model2.xInit, weights);

            end
        end

        x = intensityRefined;
        y = intensityRefined;


        nexttile;
        logScale3D(gca, x, y, M_cumI/sum(model2.N), 'Cumulative incidence rate Within 100-days', labels, 'contour');
        ax = gca;
        ax.XTick = intensity;
        ax.YTick = intensity;
        ax.XTickLabel = string(intensity * 100) + "%";
        ax.YTickLabel = string(intensity * 100) + "%";
        temp = get(ax, 'colormap');
        set(ax, 'colormap', flip(temp));

        text(-0.1, 1.05, string(char(i+67)), 'Unit', 'normalized', 'FontSize', 12, 'FontWeight', 'normal');


    end
    title(tile1, "\fontname{Times New Roman}\it R_{eff}\rm = " + Reff(k));

    exportgraphics(fig, "simulatInterventionReff" + Reff(k) + ".jpg", 'Resolution', 300);

end


    Recordi = {"Model" + n, Reff(k), interventionNames1{i}, intensity(j), tSpan_longPeriod(end), M_cumI_longPeriod, M_peak, M_duration};
    Record = cell2table(Record, 'VariableNames', {'Model', 'Reff', 'InterventionNames', 'InterventionIntensity', 'tSpan', 'Totalnumberofcases', 'Dailyincidencepeak', 'Epidemicduration'});
    writetable(Record, 'singleInterventionAnalysis.xlsx');




