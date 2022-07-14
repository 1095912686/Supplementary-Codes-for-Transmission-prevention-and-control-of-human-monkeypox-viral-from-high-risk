clear; close all; clc;

% n = 3;
% modelSEIR1 = dynamicalModel_SEIRn;
% modelSEIR1.n = n;
% modelSEIR1.N = [4030711645 * 0.04;...
%     4030711645 * 0.96;...
%     3957932033]; % default global population size
% modelSEIR1.C = readmatrix('globalContactMatrix.xlsx');
% modelSEIR1.omega = 1/12;
% modelSEIR1.gamma = 1/11;
% modelSEIR1.VE = 0;

% n = 1;
% modelSEIR1 = dynamicalModel_SEIRn;
% modelSEIR1.n = n;
% modelSEIR1.N = 4030711645 + 3957932033;
% modelSEIR1.C = zeros(n); modelSEIR1.C(1) = 13.45;
% modelSEIR1.omega = 1/12;
% modelSEIR1.gamma = 1/11;
% modelSEIR1.VE = 0;

n = 2;
idx = 1:n;
modelSEIR1 = dynamicalModel_SEIRn;
modelSEIR1.n = n;
VC = 0.8 * 0.03; % current vaccination coverage
modelSEIR1.N = 1.588e8 * [(1-VC); VC];
modelSEIR1.C = repmat(13.45 * [(1-VC), VC], [2,1])
modelSEIR1.omega = 1/8.5;
modelSEIR1.gamma = 1/11;
modelSEIR1.VE = [0, 0.85] .* ones(2,1);


whos modelSEIR1

%% read data info
[~,sheetNames,~] = xlsfinfo('OtherCountries.xlsx');



%% predict with intervention

% interventions
intensity = (0 : 0.2 : 1)';
interventionNames = {'Reducing 1/\gamma', 'Reducing Cq'};


% initialize figure
fig = figure;
fig.WindowState = 'maximized';
tiledlayout1 = tiledlayout('flow');

% in each figure, fit and simulate intervention for all 6 regions
Record = [];
for k = 1:6

    % read data for region j
    data = readtable('OtherCountries.xlsx', 'Sheet', sheetNames{k});
    tData = days(data.Confirmation - data.Confirmation(1)) + 1; % time in days
    dData = data.Confirmation;  % time in dates
    dIData = data.Cases;

    % fit model
    model1 = fitModel(modelSEIR1, tData, dIData);


    %%% segment 1, fitted data and model validation
    tSpan1 = linspace(tData(1), tData(end), 1e3);
    [t1, x1] = predictModel(model1, tSpan1, model1.xInit);
    dI1.value = sum(model1.gamma * x1(:,idx + 2*n),2);
    dI1.lb = icdf('Poisson', 0.05, dI1.value);
    dI1.ub = icdf('Poisson', 0.95, dI1.value);


    % visualize for segment 1
    tile1 = nexttile;
    hold on;
    p1 = fill([tSpan1, flip(tSpan1)], [cumtrapz(tSpan1, dI1.lb); flip(cumtrapz(tSpan1, dI1.ub))], 'red');
    p1.FaceColor = [1 0.8 0.8];
    p1.FaceAlpha = 0.6;
    p1.EdgeColor = 'none';

    plot(tData, cumsum(dIData), '.', 'MarkerSize', 8, 'MarkerEdgeColor', [139, 95, 191]/255);  % training data
    plot(tSpan1, cumtrapz(tSpan1, dI1.value), 'LineStyle', '-', 'Color', '#BDD358', 'LineWidth', 1.25); % model validation

    for i = 1:2  % 2 kinds of interventions

        % intervention intensity, initialized as 'none'
        M = ones(numel(intensity),2);

        % columns (1,2) represent intervention of reducing (1/gamma, Cq)
        M(:, i) = 1 - intensity;



        %%% segment 2, predication of one intervention with series efficiencies
        for j = 1:size(M,1)

            % extract efficiency
            [dgamma, dCq] = deal(M(j,1), M(j,2));

            if dgamma == 0
                continue;
            end

            % setup model2 of given intervention
            model2 = model1;
            model2.gamma = model2.gamma / dgamma;
            model2.C = model2.C * dCq;

            % simulates
            predicationLength = 30;
            tSpan2 = model1.tFinal + linspace(0, predicationLength, 1e3);
            [t2, x2] = predictModel(model2, tSpan2, x1(end,:));
            dI2.value = sum(model2.gamma * x2(:,idx + 2*n),2);
            dI2.lb = icdf('Poisson', 0.05, dI2.value);
            dI2.ub = icdf('Poisson', 0.95, dI2.value);

            % visualize for segment 2
            %         p2 = fill([tSpan2, flip(tSpan2)], [dI2.lb; flip(dI2.ub)], 'red');
            %         p2.FaceColor = [1 0.8 0.8];
            %         p2.EdgeColor = 'none';

            cumPredicted = trapz(tSpan1, dI1.value) + cumtrapz(tSpan2, dI2.value);
            if j == 1
                plot(tSpan2, cumPredicted, 'LineStyle', '-.', 'Color', '#FF686B', 'LineWidth', 1.25); % model predication
            else
                if i == 1
                    plot(tSpan2, cumPredicted, 'LineStyle', '-.', 'Color', '#C3979F'); % model predication
                end

                if i == 2
                    plot(tSpan2, cumPredicted, 'LineStyle', '-.', 'Color', '#3FA7D6'); % model predication
                end

            end



            Recordi = {dData(end) + days(predicationLength),...
                sheetNames(k),...
                interventionNames(i),...
                1-M(j,i),...
                cumPredicted(end)};
            Record = [Record; Recordi];
        end

        YLim = get(gca, 'YLim');
        plot(tSpan1(end) * ones(1e3, 1), linspace(YLim(1),YLim(2), 1e3),'black:');

        %title([sheetNames{k} + ", reduction of (\gamma, q, c) set as (" + dgamma + ", " + dC + ", " + dq + ')']);
        %title([sheetNames{k}, ', ', interventionNames{i}]);
        title(sheetNames{k});

        timeTicks = tSpan1(1) : 7 : tSpan2(end);
        dayTicks = dData(1) + days(timeTicks - 1);
        xticks(timeTicks);
        xticklabels(string(dayTicks));

        if k == 1
            legend('95% CIs', 'Training Data', 'Model Validation',...
                'Predication without Intervention', '', '', '', 'Predication by Reducing 1/\gamma', '', ...
                '', '', '', '', '', 'Predication by Reducing Cq', ...
                'Start of Intervention',...
                'Location', 'northwest', 'FontSize', 10, 'Box', 'off');
        end

        ax = gca;
        set(ax, 'FontName', 'Times New Roman');
        set(ax, 'FontSize', 12);
        %ax.Title.FontWeight = 'bold';

        if i == 2
            text(-0.1, 1.1, string(char(k+64)), 'Unit', 'Normalized', 'FontSize', 12, 'FontWeight', 'normal');
        end
    end

end

RecordTable = cell2table(Record, 'VariableNames', {'until', 'Region', 'Intervention', 'Intensity', 'CumulativeReportedCases'})


ylabel(tiledlayout1, '\fontname{Times New Roman} Cumulative Reported Cases', 'FontSize', 14, 'FontWeight', 'bold');
xlabel(tiledlayout1, '\fontname{Times New Roman}  Date', 'FontSize', 14, 'FontWeight', 'bold');


exportgraphics(fig, "predicationWithIntervention.jpg", 'Resolution', 1200);

