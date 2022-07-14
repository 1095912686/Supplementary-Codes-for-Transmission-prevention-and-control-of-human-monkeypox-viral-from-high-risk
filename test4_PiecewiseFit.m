close all; clc; clear;

n = 2;
modelSEIR = dynamicalModel_SEIRn;
modelSEIR.n = n;
VC = 0.8 * 0.345; % current vaccination coverage
modelSEIR.N = 7.9886e9 * [(1-VC); VC];
modelSEIR.C = repmat(13.45 * [(1-VC), VC], [2,1])
modelSEIR.omega = 1/12;
modelSEIR.gamma = 1/11;
modelSEIR.VE = [0, 0.85] .* ones(2,1);

% n = 1;
% modelSEIR = dynamicalModel_SEIRn;
% modelSEIR.n = n;
% modelSEIR.N = 4030711645 + 3957932033;
% modelSEIR.C = zeros(n); modelSEIR.C(1) = 13.45;
% modelSEIR.omega = 1/12;
% modelSEIR.gamma = 1/11;
% modelSEIR.VE = 0;

%% read data
k = 1;  % number of sheet
[~,sheetNames,~] = xlsfinfo('OtherCountries.xlsx');
data = readtable('OtherCountries.xlsx', 'Sheet', sheetNames{k});

tData = days(data.Confirmation - data.Confirmation(1)) + 1; % time in days
dData = data.Confirmation;  % time in dates
dIData = data.Cases;

%% fit
breakPoints = [1; 40; numel(tData)];
[modelList, t, x] = piecewiseFit(modelSEIR, tData, dIData, breakPoints);

% show objects
for i = 1:numel(modelList)
    disp("Model of Segment " + i + " :")
    disp(modelList(i))
end

%% visualize
idx = 1:n;
 
fig1 = figure;
tile1 = tiledlayout('flow');

t1 = nexttile(tile1);
hold on;
plot(tData, dIData, 'bo'); 
plot(t, sum(modelSEIR.gamma' .* x(:, idx + 2*n), 2), 'r*');
ylabel(t1, "Daily New Incidence");

t2 = nexttile(tile1);
plot(tData, cumsum(dIData), 'bo'); hold on;
plot(t, cumtrapz(t, sum(modelSEIR.gamma' .* x(:, idx + 2*n), 2)), 'r*');
ylabel(t2, "Accumulative Cases");

xlabel(tile1, "Time (in days)");