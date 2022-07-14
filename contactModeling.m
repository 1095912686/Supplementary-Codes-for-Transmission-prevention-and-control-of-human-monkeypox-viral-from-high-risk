%% partition of High/Low risks
Data = readmatrix('12976_2009_207_MOESM1_ESM.csv');
Data = array2table(Data, 'VariableNames', {'Number', 'Duration', 'Conversation', 'PhysicalContact', 'Kissing'});
durationLegend = {'<5min', '5-15min', '15-60min', '1-2h', '2-4h', '>4h'};

sData = groupsummary(Data, "Duration");
sData.Properties.RowNames = durationLegend


%idx =  (Data.Duration >= 3 & Data.PhysicalContact == 1) | Data.Kissing == 1;
%idx = Data.Duration >= 4;
%idx = Data.PhysicalContact == 1 & Data.Conversation == 0;
% idx1 = Data.PhysicalContact == 1 | Data.Kissing == 1;   % by direct physical contact
% idx2 = (Data.Duration >=4 & Data.Conversation == 1) & ~idx1;      % only by air 
% idx = idx1 | idx2;
%idx = (Data.Duration >= 4 & Data.Conversation == 1 & Data.PhysicalContact == 1) | Data.Kissing == 1; 
%idx = Data.PhysicalContact == 1 & Data.Duration >= 4;

idx = Data.PhysicalContact == 1 | Data.Kissing == 1;
sampleSize = numel(idx)
specificCount = sum(idx)
percentage = specificCount / sampleSize





%% Full Contact Matrix
% groupLegend = {'MSM', 'Male not MSM', 'Female'};
c = [13.51; 13.51; 13.39];
r = 0.04;
rN =  [4030711645 * 0.04;...
            4030711645 * 0.96;...
            3957932033];

k = c .* rN / (c' * rN);
C1 = c * k'
writematrix(C1, 'globalContactMatrix.xlsx');


% % equation
% m = 1e3;
% x = [linspace(5,20,m)', ones(m,1)];
% c = [1, 13];
% q_fitted_x = params.N(1) * fittedParameterRecords(:,1)' ./ sum(x.*c, 2);
% 
% % t4 = nexttile;
% % plot(x(:,1), q_fitted, 'r', 'LineWidth', 1);
% % xlabel('relatively susceptibility x');
% % ylabel('q_{fitted}')

