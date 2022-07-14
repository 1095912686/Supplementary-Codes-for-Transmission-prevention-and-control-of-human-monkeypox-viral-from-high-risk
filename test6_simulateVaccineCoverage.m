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
modelSEIR.C = repmat(13.45 * [(1-VC), VC], [2,1])
modelSEIR.omega = 1/12;
modelSEIR.gamma = 1/11;
modelSEIR.VE = [0, 0.85] .* ones(2,1);


whos modelSEIR

idx = 1:n;
x0 = zeros(4*n, 1);
x0(idx + 0*n) = modelSEIR.N;
x0(idx + 2*n) = 1; %modelSEIR.N / sum(modelSEIR.N); % first case assigned proportional to population size
modelSEIR.xInit = x0;

Reff = 3;%[1,2,3,4,5,6,7,8];
for k = 1:numel(Reff)


    fig = figure;
    fig.WindowState = 'maximized';
    tile1 = tiledlayout('flow');


    modelSEIR.Reff = Reff(k);
    modelSEIR.q = Reff(k) * modelSEIR.gamma / max(eig(modelSEIR.C .* (1 - modelSEIR.VE)));

    %% predict with intervention

    % interventions
    efficiency = flip(0.1 : 0.1 : 1)';
    newVC = 0 : 0.05 : 0.9;
    interventionNames = {'reducing 1/\gamma by', 'reducing C*q by', 'set vaccine coverage as'};
    for i = 3%1:3  % 3 kinds of interventions

        t1 = nexttile;

        for j = 1:numel(newVC)

            % setup model2 of given intervention
            model2 = modelSEIR;
            switch i
                case 1
                    dgamma = efficiency(j);
                    model2.gamma = model2.gamma / dgamma;
                case 2
                    dq = efficiency(j);
                    model2.q = model2.q * dq;
                case 3
                    VC = newVC(j);
                    model2.N = 1.588e8 * [(1-VC); VC];
                    model2.C = repmat(13.45 * [(1-VC), VC], [2,1])
                    model2.xInit(idx + 0*n) = model2.N; 
            end

            

            % simulates
            tSpan = [0, 10e2];
            [t, x] = predictModel(model2, tSpan, model2.xInit);

            dI.value = sum(model2.omega * x(:,idx + 1*n), 2);
            dI.lb = icdf('Poisson', 0.05, dI.value);
            dI.ub = icdf('Poisson', 0.95, dI.value);


            %plot(t1, t, dI.value / sum(model2.N)); hold on;
%             plot(t1, t, dI.value); hold on;
            plot(t, x(:,idx + 1*n));
            legend('un-vaccinated', 'vaccinated')
            title(newVC(j))

            if j ~= numel(newVC)
                nexttile;
            end


        end


        if i == 3
            legend(t1, interventionNames(i) + " " + string(newVC*100) + "%")
        else
            legend(t1, interventionNames(i) + " " + string((1-efficiency)*100) + "%")
        end

        ylabel('Daily New Incidence')
        xlabel('Time (in days)');
    end
    title(tile1, "Reff = " + Reff(k));

end



