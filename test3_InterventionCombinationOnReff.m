clear; close all; clc;
ImportDIR

%%% tabular for interventions
tableRecord = table();
for n = 2%[1,2,3]
    switch n
        case 1
            modelSEIR = dynamicalModel_SEIRn;
            modelSEIR.n = n;
            modelSEIR.N = 4030711645 + 3957932033;
            modelSEIR.C = zeros(n); modelSEIR.C(1) = 13.45;
            modelSEIR.omega = 1/12;
            modelSEIR.gamma = 1/11;
            modelSEIR.VE = 0;

        case 2
            modelSEIR = dynamicalModel_SEIRn;
            modelSEIR.n = n;
            VC = 0.8 * 0.03; % current vaccination coverage
            modelSEIR.N = 1.588e8 * [(1-VC); VC]; %7.9886e9 * [(1-VC); VC];
            modelSEIR.C = repmat(13.45 * [(1-VC), VC], [2,1])
            modelSEIR.omega = 1/8.5;
            modelSEIR.gamma = 1/11;
            modelSEIR.VE = [0, 0.85]';

        case 3
            modelSEIR = dynamicalModel_SEIRn;
            modelSEIR.n = n;
            modelSEIR.N = [4030711645 * 0.04;...
                4030711645 * 0.96;...
                3957932033]; % default global population size
            modelSEIR.C = readmatrix('globalContactMatrix.xlsx');
            modelSEIR.omega = 1/12;
            modelSEIR.gamma = 1/11;
            modelSEIR.VE = 0;
    end

    for Reff = 3%4:7

        modelSEIR.q = Reff * modelSEIR.gamma / max(eig(modelSEIR.C .* (1 - modelSEIR.VE)));

 



        %% predict with intervention

% interventions
intensity = (1 : -0.05 : 0.1)';
efficiency = 1 - intensity;
m = numel(intensity);
interventionNames = {'reducing 1/\gamma', 'reducing Cq'};

    % intervention intensity, initialized as 'none'
    M = zeros(m^2, 4); % {'Duration', 'dgamma', 'dCq', 'accumulativeI', 'peakLocation_E_plus_I', 'LocationWhenDailyIncidenceLesserThan1'}
    M_cumI_longPeriod = zeros(m);
    M_cumI_shortPeriod = zeros(m);
    M_peak = zeros(m);
    M_duration = zeros(m);

    for i = 1:m
        for j = 1:m
            idx = (i-1)*m + j;

            model2 = modelSEIR;
            model2.gamma = modelSEIR.gamma / intensity(i);
            model2.q = modelSEIR.q * intensity(j);

            weights = ones(n,1);
            xInit = zeros(4*n, 1);
            xInit((1:n) + 0*n) = modelSEIR.N; % S = N
            xInit(2*n+1) = 1; % initI_1 = 1
            tSpan_longPeriod = [0,100]; % for 100-days
            tSpan_shortPeriod = [0,30];
            M_cumI_longPeriod(i,j) = weightedAccumulativeCases(model2, xInit, tSpan_longPeriod, weights);
            M_cumI_shortPeriod(i,j) = weightedAccumulativeCases(model2, xInit, tSpan_shortPeriod, weights);
            M_peak(i,j) = peak(model2, xInit, weights);
            M_duration(i,j) = last(model2, xInit, weights);

            if idx == 1
                [t,x] = predictModel(model2, tSpan_longPeriod, xInit);
                figure;
                plot(t,sum(x(:,2*n+1:3*n), 2)); hold on;
            end
            M(2*idx-1, 1) = 30;
            M(2*idx-1, 2) = intensity(i);
            M(2*idx-1, 3) = intensity(j);
            M(2*idx-1, 4) = M_cumI_shortPeriod(i,j);
            M(2*idx-1, 5) = M_peak(i,j);
            M(2*idx-1, 6) = M_duration(i,j);

            M(2*idx, 1) = 100;
            M(2*idx, 2) = intensity(i);
            M(2*idx, 3) = intensity(j);
            M(2*idx, 4) = M_cumI_longPeriod(i,j);
            M(2*idx, 5) = M_peak(i,j);
            M(2*idx, 6) = M_duration(i,j);
        end
    end

    %% visualize for three matrices
    fig = figure;
    fig.WindowState = 'maximized';
    tile1 = tiledlayout('flow');
    
    
    x = efficiency;
    y = efficiency;
   

%     nexttile;
%     logScale3D(gca, x, y, M_cumI, 'Accumulative Cases Within 100-days', 'surf')
% 
%     nexttile;
%     logScale3D(gca, x, y, M_peak, 'Peak Instance of Epidemic', 'surf')
% 
%     nexttile;
%     logScale3D(gca, x, y, M_duration, 'Duration of Epidemic', 'surf')

    labelNames = {'Reduction of Infectious Period', 'Reduction of C*q'}

    nexttile;
    logScale3D(gca, x, y, M_cumI_longPeriod/sum(model2.N, 'all'), 'Accumulative Cases Within 100-days', labelNames, 'contour'); 
    
    nexttile;
    logScale3D(gca, x, y, M_cumI_shortPeriod/sum(model2.N, 'all'), 'Accumulative Cases Within 30-days', labelNames, 'contour'); 
    
%     nexttile;
%     logScale3D(gca, x, y, M_peak, 'Peak Instance of Epidemic', 'contour'); 
% 
%     nexttile;
%     logScale3D(gca, x,y, M_duration, 'Duartion of Epidemic', 'contour');
    sgtitle("Intervention Analysis (" + "Model" + n + ", " + "R_{eff} =" + Reff  + ")")
    exportgraphics(fig, "interventionAnalysis"+n+Reff+".pdf", 'Resolution', 600);
    1;

    end

    temp1 = array2table(repmat(Reff, [2*m^2, 1]), 'VariableNames', {'Reff'});
    temp2 = array2table(repmat({['SEIR', num2str(n)]}, [2*m^2, 1]), 'VariableNames', {'Model'});
    temp3 = array2table(M(:,1:6), 'VariableNames', {'Duration', 'dgamma', 'dCq', 'accumulativeI', 'peak', 'timeToEnd'});

    Table_modeln_regions = [temp1, temp2, temp3];
    tableRecord = [tableRecord; Table_modeln_regions];
end

writetable(tableRecord, 'interventionAnalysis.xlsx');




