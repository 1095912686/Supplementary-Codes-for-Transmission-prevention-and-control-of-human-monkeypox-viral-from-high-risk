clear; close all; clc;
ImportDIR

%%% tabular for interventions
tableRecord = table();
for n = [1,2,3]
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
            VC = 0.8 * 0.345; % current vaccination coverage
            modelSEIR.N = 7.9886e9 * [(1-VC); VC];
            modelSEIR.C = repmat(13.45 * [(1-VC), VC], [2,1])
            modelSEIR.omega = 1/12;
            modelSEIR.gamma = 1/11;
            modelSEIR.VE = [0, 0.85] .* ones(2,1);

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

            whos modelSEIR1
    end
%% read data
[~,sheetNames,~] = xlsfinfo('OtherCountries.xlsx');

for s = 2%1:6 % all regions
data = readtable('OtherCountries.xlsx', 'Sheet', sheetNames{s});

tData = days(data.Confirmation - data.Confirmation(1)) + 1; % time in days
dData = data.Confirmation;  % time in dates
IData = data.Cases;

%% fit
model1 = fitModel(modelSEIR, tData, IData)

%% predict with intervention

% interventions
intensity = (1 : -0.05 : 0.1)';
efficiency = 1 - intensity;
m = numel(intensity);
interventionNames = {'reducing 1/\gamma', 'reducing Cq'};

    % intervention intensity, initialized as 'none'
    M = zeros(m^2, 4); % {'Duration', 'dgamma', 'dCq', 'accumulativeI', 'peakLocation_E_plus_I', 'LocationWhenDailyIncidenceLesserThan1'}
    M_cumI100 = zeros(m);
    M_cumI30 = zeros(m);
    M_peak = zeros(m);
    M_duration = zeros(m);

    for i = 1:m
        for j = 1:m
            idx = (i-1)*m + j;

            model2 = model1;
            model2.gamma = model1.gamma / intensity(i);
            model2.q = model1.q * intensity(j);

            weights = ones(n,1);
            xInit = zeros(4*n, 1);
            xInit(1:n) = model2.xInit(1:n); % S = N
            xInit(2*n+1:3*n) = 1 / n; % initI_i = 1 / n
            tSpan100 = linspace(1,100,1e3); % for 100-days
            tSpan30 = linspace(1,30,1e3);
            M_cumI100(i,j) = weightedAccumulativeCases(model2, xInit, tSpan100, weights);
            M_cumI30(i,j) = weightedAccumulativeCases(model2, xInit, tSpan30, weights);
%             M_peak(i,j) = peak(model2, xInit, weights);
%             M_duration(i,j) = last(model2, xInit, weights);
            M(2*idx-1, 1) = 30;
            M(2*idx-1, 2) = intensity(i);
            M(2*idx-1, 3) = intensity(j);
            M(2*idx-1, 4) = M_cumI30(i,j);

            M(2*idx, 1) = 100;
            M(2*idx-1, 2) = intensity(i);
            M(2*idx-1, 3) = intensity(j);
            M(2*idx-1, 4) = M_cumI100(i,j);
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
    logScale3D(gca, x, y, M_cumI100, 'Accumulative Cases Within 100-days', labelNames, 'contour'); 
    
    nexttile;
    logScale3D(gca, x, y, M_cumI30, 'Accumulative Cases Within 30-days', labelNames, 'contour'); 
    
%     nexttile;
%     logScale3D(gca, x, y, M_peak, 'Peak Instance of Epidemic', 'contour'); 
% 
%     nexttile;
%     logScale3D(gca, x,y, M_duration, 'Duartion of Epidemic', 'contour');
    sgtitle("Intervention Analysis (" + "Model" + n + ", " + string(sheetNames{s} + ")"))
    exportgraphics(fig, "interventionAnalysis"+n+s+".pdf", 'Resolution', 600);
    1;

    end

    temp1 = array2table(repmat(sheetNames(s), [2*m^2, 1]), 'VariableNames', {'Region'});
    temp2 = array2table(repmat({['SEIR', num2str(n)]}, [2*m^2, 1]), 'VariableNames', {'Model'});
    temp3 = array2table(M(:,1:4), 'VariableNames', {'Duration', 'dgamma', 'dCq', 'accumulativeI'});

    Table_modeln_regions = [temp1, temp2, temp3];
    tableRecord = [tableRecord; Table_modeln_regions];
end

writetable(tableRecord, 'interventionAnalysis.xlsx');




