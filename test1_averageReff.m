      clear; close all; clc; format shortE
%%
% For the average Reff by now


% read data 
[~,sheetNames,~] = xlsfinfo('OtherCountries.xlsx');
opts = detectImportOptions('OtherCountries.xlsx');


%% for each region
allRegionRecord = []; 
recordR = []; % [q_medain, q_std, Reff_median, Reff_std, R0_median, R0_std]
shortRecord.region = sheetNames;

% initialize figures
fig1 = figure;
fig1.WindowState = 'maximized';
tile1 = tiledlayout('flow');
sgtitle('Fitted Accumulative Cases');

fig2 = figure;
fig2.WindowState = 'maximized';
tile2 = tiledlayout('flow');
sgtitle('Corresponding Daily New Incidence');

fig3 = figure;
fig3.WindowState = 'maximized';
tile3 = tiledlayout('flow');
sgtitle('Fitted Probability q of Infection via Single Contact');

fig4 = figure;
fig4.WindowState = 'maximized';
tile4 = tiledlayout('flow');
sgtitle('Distribution of Fitted q');

for k = 1:numel(sheetNames)
data = readtable('OtherCountries.xlsx', 'Sheet', sheetNames{k});


%%% preview daily incidence counts
tData = days(data.Confirmation - data.Confirmation(1)); % time in days
dData = data.Confirmation;  % time in dates
dIData = data.Cases;
if k == 6
    dIData(2:end) = dIData(2:end) * 0.78;
end




%% Fit Models
% %%% un-grouped Model
% n = 1;
% modelSEIR1 = dynamicalModel_SEIRn;
% modelSEIR1.n = n;
% modelSEIR1.N = 4030711645 + 3957932033;
% modelSEIR1.C = zeros(n); modelSEIR1.C(1) = 13.45;
% modelSEIR1.omega = 1/12;
% modelSEIR1.gamma = 1/11;
% modelSEIR1.VE = 0;

%%% vaccinated Model
n = 2;
modelSEIR1 = dynamicalModel_SEIRn;
modelSEIR1.n = n;
VC = 0.8 * 0.03; % current vaccination coverage
modelSEIR1.N = 1.588e8 * [(1-VC); VC];
modelSEIR1.C = repmat(13.45 * [(1-VC), VC], [2,1])
modelSEIR1.omega = 1/8.5;
modelSEIR1.gamma = 1/11;
modelSEIR1.VE = [0, 0.85] .* ones(2,1);


idx = 1:n;


    %t.start = numel(tData) - 14; % start by 
    t.start = 20; % start by 
    t.end = numel(tData); % end by 

    q = zeros(t.end-t.start+1, 1);
    E0 = zeros(t.end-t.start+1, n);
    Reff = zeros(t.end-t.start+1, 1);
    R0 = zeros(t.end-t.start+1, 1);
    Rsquare1 = zeros(t.end-t.start+1, 1);
    Rsquare2 = zeros(t.end-t.start+1, 1);
    MAP = zeros(t.end-t.start+1, 1);

    for i = t.start : t.end
        tDatai = tData(1:i);
        dDatai = dData(1:i);
        dIDatai = dIData(1:i);

        
        %fittedModel = fitModel(modelSEIR1, tDatai, dIDatai);
        breakPoints = [tDatai(1), tDatai(end)];
        %[fittedModel, ~, ~] = piecewiseFit(modelSEIR1, tDatai, dIDatai, breakPoints);
        [fittedModel, ~, ~] = fitModel(modelSEIR1, tDatai, dIDatai);
    
        q(i-t.start+1) = fittedModel.q;
        E0(i-t.start+1,:) = fittedModel.E0;
        Reff(i-t.start+1) = fittedModel.Reff;
        R0(i-t.start+1) = fittedModel.R0;

        tRefined = fittedModel.t;
        xPredicted = fittedModel.x;
        IPredicted = sum(fittedModel.gamma * xPredicted(:,idx + 2*n), 2);
        IData = dIDatai;
        cumPredict = cumtrapz(tRefined, IPredicted)  + dIDatai(1);
        cumPredict = interp1(tRefined, cumPredict, tDatai);
        cumData = cumsum(dIDatai);
        SSE1 = norm(cumData - cumPredict)^2;
        SST1 = norm(cumData - mean(cumData))^2;
        Rsquare1(i-t.start+1) = 1 - SSE1 / SST1;

        IPredicted = interp1(tRefined, IPredicted, tDatai);
        SSE2 = norm(IData - IPredicted) ^ 2;
        SST2 = norm(IData - mean(IData)) ^ 2;
        Rsquare2(i-t.start+1) = 1 - SSE2 / SST2;

        temp = abs((cumData - cumPredict) ./ cumData);
        %temp(isnan(temp)) = [];
        MAP(i-t.start+1) = sum(temp);

        if i == t.end

           
            t1 = nexttile(tile1);
            predictedAccumluativeCases.value = cumtrapz(tRefined, sum(fittedModel.gamma.' .* xPredicted(:,idx + 2*n), 2))  + dIDatai(1);
            predictedAccumluativeCases.lb = icdf('Poisson', 0.05, predictedAccumluativeCases.value);
            predictedAccumluativeCases.ub = icdf('Poisson', 0.95, predictedAccumluativeCases.value);
            hold on;
            plot(dDatai, cumsum(dIDatai),'o', 'MarkerSize', 3, 'MarkerEdgeColor', [79, 81, 140]/255, 'DateTimeTickFormat', 'yyyy-MM-dd');
            plot(tRefined, predictedAccumluativeCases.value, 'LineStyle', '-', 'Color', '#3E8914', 'LineWidth', 1);
            p2 = fill([tRefined; flip(tRefined)], [predictedAccumluativeCases.lb; flip(predictedAccumluativeCases.ub)],'red');
            p2.FaceColor = [1 0.8 0.8];
            p2.FaceAlpha = 0.6;
            p2.EdgeColor = 'none';
            uistack(p2, "bottom");
            legend('observation', 'model prediction', 'Location', 'northwest');
            title(t1, sheetNames{k});
            ylabel(tile1, 'accumlative case count');
            xlabel(tile1, 'Date')


            t2 = nexttile(tile2);
            predictedDailyReportedCases.value = sum(fittedModel.gamma.' .* xPredicted(:,idx + 2*n), 2);
            predictedDailyReportedCases.lb = icdf('Poisson', 0.05, predictedDailyReportedCases.value);
            predictedDailyReportedCases.ub = icdf('Poisson', 0.95, predictedDailyReportedCases.value);
            hold on;
            plot(dData, dIData, 'bo', 'DatetimeTickFormat', 'yyyy-MM-dd');
            plot(tRefined, predictedDailyReportedCases.value, 'r-', 'LineWidth', 1);
            p1 = fill([tRefined; flip(tRefined)], [predictedDailyReportedCases.lb; flip(predictedDailyReportedCases.ub)], 'red');
            p1.FaceColor = [1 0.8 0.8];
            p1.EdgeColor = 'none';
            uistack(p1, "bottom");
            legend('observation', 'model prediction', 'Location', 'northwest');
            title(t2, sheetNames{k});
            ylabel(tile2, 'Daily New Incidence Rate');
            xlabel(tile2, 'Date')

        end
        
    end


    %% Tabular and write
    % full record table
    dateTable = array2table(dData(t.start:t.end), 'VariableNames', {'Date'});
    countryNames = cell2table(repmat(sheetNames(k), [t.end-t.start+1, 1]), 'VariableNames', {'Region'});
    fittedParameterTable = array2table([q, Reff, R0, Rsquare1, Rsquare2 MAP], 'VariableNames', {'q', 'Reff', 'R0', 'Rsquare1', 'Rsquare2', 'MAPE'});
    countryRecord = [countryNames, dateTable, fittedParameterTable];
    allRegionRecord = [allRegionRecord; countryRecord];

    % short record for country-wise median R0/Reff
    recordR = [recordR; ...
        [median(fittedParameterTable.q), std(fittedParameterTable.q),...
         median(fittedParameterTable.Reff), std(fittedParameterTable.Reff),...
         median(fittedParameterTable.R0), std(fittedParameterTable.R0),...
         median(fittedParameterTable.Rsquare1), std(fittedParameterTable.Rsquare1)...
         median(fittedParameterTable.Rsquare2), std(fittedParameterTable.Rsquare2)...
         median(fittedParameterTable.MAPE), std(fittedParameterTable.MAPE)]];

    %% visualizing
    t3 = nexttile(tile3);
    ax1 = plot(dData(t.start:t.end), q,'r*-');
    ax1.DatetimeTickFormat = "yyyy-MM-dd";
    title(t3, sheetNames{k});
    xlabel(t3, "Interval Determined From " + string(dData(1)) + " to Date");
    ylabel(tile3, 'Fitted q');


end


recordR = array2table(recordR, 'VariableNames', {'q_median', 'q_std', 'Reff_median', 'Reff_std', 'R0_median', 'R0_std', 'Rsquare1_median', 'Rsquare1_std', 'Rsquare2_median', 'Rsquare2_std', 'MAPE_median', 'MAPE_std'});
recordR = [cell2table(sheetNames(:), 'VariableNames', {'Region'}), recordR]
writetable(allRegionRecord, 'allRegionRecord.xlsx');

% q and its 95% quantiles
Region = unique(allRegionRecord.Region);
q_median = groupsummary(allRegionRecord.q, allRegionRecord.Region, "median");
q_upper75 = groupsummary(allRegionRecord.q, allRegionRecord.Region, @(X)quantile(X,0.75));
q_lower25 = groupsummary(allRegionRecord.q, allRegionRecord.Region, @(X)quantile(X,0.25));
q = table(Region, q_median, q_lower25, q_upper75)


% Reff and its 95% quantiles
Region = unique(allRegionRecord.Region);
Reff_median = groupsummary(allRegionRecord.Reff, allRegionRecord.Region, "median");
Reff_upper75 = groupsummary(allRegionRecord.Reff, allRegionRecord.Region, @(X)quantile(X,0.75));
Reff_lower25 = groupsummary(allRegionRecord.Reff, allRegionRecord.Region, @(X)quantile(X,0.25));
Reff = table(Region, Reff_median, Reff_lower25, Reff_upper75)


% R0 and its 95% quantiles
% Region = unique(allRegionRecord.Region);
% R0_median = groupsummary(allRegionRecord.R0, allRegionRecord.Region, "median");
% R0_upper75 = groupsummary(allRegionRecord.R0, allRegionRecord.Region, @(X)quantile(X,0.75));
% R0_lower25 = groupsummary(allRegionRecord.R0, allRegionRecord.Region, @(X)quantile(X,0.25));
% R0 = table(Region, R0_median, R0_lower25, R0_upper75)





%%% Figure 4
t4_0 = nexttile(tile4);
namedRegion = categorical(allRegionRecord.Region, sheetNames, sheetNames);
b1 = boxchart(namedRegion, allRegionRecord.q, 'Notch', 'on');

%b1 = boxplot(allRegionRecord.q, namedRegion, 'Notch', 'on');
ylabel('q');


t4_1 = nexttile(tile4);
namedRegion = categorical(allRegionRecord.Region, sheetNames, sheetNames);
b2 = boxchart(namedRegion, allRegionRecord.Reff, 'Notch', 'on');
%b2 = boxplot(allRegionRecord.Reff, namedRegion, 'Notch', 'on');
ylabel('Reff');

% t4_2 = nexttile(tile4);
% b3 = boxchart(namedRegion, allRegionRecord.R0, 'Notch', 'on');
% %b3 = boxplot(allRegionRecord.R0, namedRegion, 'Notch', 'on');
% ylabel('R0');

exportgraphics(fig1,'FittedAccumulativeCases.jpg', 'Resolution', 600);
exportgraphics(fig2, 'CorrspondingDailyNewIncidence.jpg', 'Resolution', 600);
exportgraphics(fig3, 'FittedProbability_q.jpg', 'Resolution', 600);
% exportgraphics(fig4, 'DistributionOfReffR0.jpg', 'Resolution', 600);

exportgraphics(fig1,'FittedAccumulativeCases.pdf', 'Resolution', 600);
exportgraphics(fig2, 'CorrspondingDailyNewIncidence.pdf', 'Resolution', 600);
exportgraphics(fig3, 'FittedProbability_q.pdf', 'Resolution', 600);
% exportgraphics(fig4, 'DistributionOfReffR0.pdf', 'Resolution', 600);



