clear; close all; clc;
ImportDIR

%%% tabular for interventions
tableRecord = table();
for dVC = 0 : 5e-2 : 1-0.8*0.03

%             n = 1;
%             modelSEIR = dynamicalModel_SEIRn;
%             modelSEIR.n = n;
%             modelSEIR.N = 4030711645 + 3957932033;
%             modelSEIR.C = zeros(n); modelSEIR.C(1) = 13.45;
%             modelSEIR.omega = 1/12;
%             modelSEIR.gamma = 1/11;
%             modelSEIR.VE = 0;
% 
            n = 2;
            modelSEIR = dynamicalModel_SEIRn;
            modelSEIR.n = n;
            VC = 0.8 * 0.03; % current vaccination coverage
            modelSEIR.N = 1.588e8 * [(1-VC); VC]; %7.9886e9 * [(1-VC); VC];
            modelSEIR.C = repmat(13.45 * [(1-VC), VC], [2,1])
            modelSEIR.omega = 1/12;
            modelSEIR.gamma = 1/11;
            modelSEIR.VE = [0, 0.85]';
% 
%             n = 3;
%             modelSEIR = dynamicalModel_SEIRn;
%             modelSEIR.n = n;
%             modelSEIR.N = [4030711645 * 0.04;...
%                 4030711645 * 0.96;...
%                 3957932033]; % default global population size
%             modelSEIR.C = readmatrix('globalContactMatrix.xlsx');
%             modelSEIR.omega = 1/12;
%             modelSEIR.gamma = 1/11;
%             modelSEIR.VE = 0;


    for Reff = 5%4:7

        modelSEIR.q = Reff * modelSEIR.gamma / max(eig(modelSEIR.C .* (1 - modelSEIR.VE)));

 



        %% predict with intervention

% interventions
intensity = (0.1 : 0.05 : 1)';
m = numel(intensity);
interventionNames = {'reducing 1/\gamma', 'reducing Cq'};

    % intervention intensity, initialized as 'none'
    M = zeros(m^2, 4); % {'Duration', 'dgamma', 'dCq', 'accumulativeI', 'peakLocation_E_plus_I', 'LocationWhenDailyIncidenceLesserThan1'}
    M_cumI100 = zeros(m);
    M_cumI14 = zeros(m);
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
            tSpan100 = [0,100]; % for 100-days
            tSpan14 = [0,28];
            M_cumI100(i,j) = weightedAccumulativeCases(model2, xInit, tSpan100, weights);
            M_cumI14(i,j) = weightedAccumulativeCases(model2, xInit, tSpan14, weights);
            M_peak(i,j) = peak(model2, xInit, weights);
            M_duration(i,j) = last(model2, xInit, weights);

            if idx == 1
                [t,x] = predictModel(model2, tSpan100, xInit);
                figure;
                plot(t,sum(x(:,2*n+1:3*n), 2)); hold on;
            end
            M(2*idx-1, 1) = 14;
            M(2*idx-1, 2) = intensity(i);
            M(2*idx-1, 3) = intensity(j);
            M(2*idx-1, 4) = M_cumI14(i,j);
            M(2*idx-1, 5) = M_peak(i,j);
            M(2*idx-1, 6) = M_duration(i,j);

            M(2*idx, 1) = 100;
            M(2*idx, 2) = intensity(i);
            M(2*idx, 3) = intensity(j);
            M(2*idx, 4) = M_cumI100(i,j);
            M(2*idx, 5) = M_peak(i,j);
            M(2*idx, 6) = M_duration(i,j);
        end
    end

    %% visualize for three matrices
    fig = figure;
    fig.WindowState = 'maximized';
    tile1 = tiledlayout('flow');
    
    
    x = intensity;
    y = intensity;
   

%     nexttile;
%     logScale3D(gca, x, y, M_cumI, 'Accumulative Cases Within 100-days', 'surf')
% 
%     nexttile;
%     logScale3D(gca, x, y, M_peak, 'Peak Instance of Epidemic', 'surf')
% 
%     nexttile;
%     logScale3D(gca, x, y, M_duration, 'Duration of Epidemic', 'surf')

    labels = {'intensity of reducing 1 / \gamma', 'intensity of reducing C*q'};
    nexttile;
    logScale3D(gca, x, y, M_cumI100/sum(model2.N, 'all'), 'Accumulative Cases Within 100-days', labels, 'contour'); 
    
    nexttile;
    logScale3D(gca, x, y, M_cumI14/sum(model2.N, 'all'), 'Accumulative Cases Within 14-days', labels, 'contour'); 
    
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




