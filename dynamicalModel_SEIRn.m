classdef dynamicalModel_SEIRn
%dynamicalModel_SEIRn 
%
%   dynamicalModel_SEIRn methods:
%       dxdtSEIRn - compute derivatives of model ODEs
%       objfun - objective function of square loss or negative log-likelihood loss
%       fitModel - fit SEIRn model with data
%       predictModel - predict with initial values

%   dynamicalModel_SEIRn properties:
%       n - number of sub-groups
%       BetaN - normalized matrix of transmission coefficients
%       C - contact matrix
%       q - probability of infection via a single contact
%       VE - vaccine efficacy matrix/vector/scalar
%       N - population size of sub-groups
%       omega - incubation period of sub-groups
%       gamma - infectious period of sub-groups
%       f - mobility rate coefficient of sub-groups
%       dr - natual death rate coefficient of sub-groups
%       br - natual birth rate coefficient of sub-groups
%       Lambda - natual birth rate of sub-groups
%       xInit - initial values
%       xFinal - final values. [xInit, xFinal] are optional initial or boundary values
%       tInit - initial time instance
%       tFinal - final time instance. The model is defined on the time interval [tInit, tFinal] 
%



    properties
        n(1,1) double {mustBeNonnegative, mustBeInteger} = 1;
        N(:,1) double {mustBeNonnegative, mustBeFinite} = [];
        VE(:,:) double {mustBeNonnegative, mustBeFinite} = 0;
        C(:,:) double {mustBeNonnegative, mustBeFinite} = [];
        q(:,:) double {mustBeNonnegative, mustBeFinite} = [];
        omega(:,1) double {mustBeNonnegative, mustBeFinite} = [];
        gamma(:,1) double {mustBeNonnegative, mustBeFinite} = [];
        Reff(:,:) = [];
        R0(:,:) = [];
        f(:,1) double {mustBeNonnegative, mustBeFinite} = 0;
        dr(:,1) double {mustBeNonnegative, mustBeFinite} = 0;
        br(:,1) double {mustBeNonnegative, mustBeFinite} = 0;
        Lambda(:,1) double {mustBeNonnegative, mustBeFinite} = 0;
        tInit(:,1) double {mustBeNonnegative, mustBeFinite} = [];
        tFinal(:,1) double {mustBeNonnegative, mustBeFinite} = [];
        xInit(:,:) double {mustBeNonnegative, mustBeFinite} = [];
        xFinal(:,:) double {mustBeNonnegative, mustBeFinite} = [];
        E0(:,:) = [];
        t(:,:) = [];
        x(:,:) = [];

    end

    methods
        % ---------------------------------------------------------------------------------------------------------
        function dxdt = dxdt(model, t, x)
            %%% Vectorized derivatives for normalized ODEs of SEIRn model
            %
            % n*1 vector for state variable x
            % n*1 vector or scalar or n*n matrix for parameters
            % scalar for time variable t
            %
            % dxdt is derivative output of the same size as x
            %

            % number of sub-groups
            n = model.n;

            % extract variables
            idx = 1:n;
            S = x(idx + 0 * n);
            E = x(idx + 1 * n);
            I = x(idx + 2 * n);
            R = x(idx + 3 * n);

            % extract parameters
            N = model.N;
            BetaN = (model.C .* (1-model.VE)) * model.q;
            Lambda = model.Lambda(:);
            br = model.br(:);
            dr = model.dr(:);
            omega = model.omega(:);
            gamma = model.gamma(:);
            f = model.f(:);

            % derivatives
            newlyInfectionTerm = (BetaN.' * I ) .* (S ./ N); % BetaN(i,j) denotes beta_ij * N_i
            newlyInfectionTerm(isnan(newlyInfectionTerm)) = 0;
            dSdt = br .* N + Lambda - newlyInfectionTerm;
            dEdt = newlyInfectionTerm - omega.*E - dr.* E;
            dIdt = omega.*E - gamma.*I - dr.*I;
            dRdt = gamma.*I - dr.* R;

            dxdt = [dSdt; dEdt; dIdt; dRdt];
            dxdt = dxdt(:);

        end

        % ---------------------------------------------------------------------------------------------------------
        function [loss, t, x] = objfun(model, parameterValues, tData, dIData)

            % parameterValues = [q; E0], an (n+1)*1 vector

            n = model.n; % number of groups
            model.q = parameterValues(1); % initial guess of beta * N of unvaccinated


            % initial value of the ODE model
            idx = 1:n;
            if isempty(model.xInit)

                % treated as first segment (since the first case)
                % R0 = 0
                x0 = zeros(4*n, 1);

                % S0 = N
                x0(idx + 0*n) = model.N;

                % E0 to be fitted
                x0(idx + 1*n) = parameterValues(2:end);

                % I0 = ones
                x0(idx + 2*n) = max(1, dIData(1)) ./ model.gamma / n;
                %x0(idx + 2*n) = max(1, dIData(1)) / n;

            else

                % specified initial value
                x0 = model.xInit(:);

            end



            if isempty(model.tInit)
                model.tInit = min(tData);
            end
            if isempty(model.tFinal)
                model.tFinal = max(tData);
            end

            tSpan = [model.tInit, model.tFinal];
            odefun = @(t, x)dxdt(model, t, x);
            [t,x] = ode15s(odefun,tSpan,x0);

            accumulativeReprotedI = cumtrapz(t, sum(model.gamma' .* x(:,idx + 2*n),2));  % \sum gamma_i * I_i
            xPredicted = interp1(t, accumulativeReprotedI, tData) + dIData(1);

            cumsumData = cumsum(dIData);
            loss = sum((xPredicted - cumsumData) .^ 2); % total square error
            %loss = -sum(cumsumData .* log(xPredicted + eps) - xPredicted - log(cumsumData + eps)); % negative log-likelihood


        end

        % ---------------------------------------------------------------------------------------------------------
        function [fittedModel, t, x] = fitModel(model, tData, xData)

            % set default
            parametersToBeFitted = [0.1; 10 * ones(model.n, 1)]; % [q; E0];

            fun = @(parametersToBeFitted) objfun(model, parametersToBeFitted, tData, xData);


            A = [];
            b = [];
            Aeq = [];
            beq = [];
            lb = 0 * ones(size(parametersToBeFitted));
            ub = 1e3 * ones(size(parametersToBeFitted));

            fun(parametersToBeFitted)
            [fittedParameters] = fmincon(fun, parametersToBeFitted, A, b, Aeq, beq, lb, ub);
            [~, t, x] = objfun(model, fittedParameters, tData, xData);

            fittedModel = model;
            fittedModel.q = fittedParameters(1);
            fittedModel.E0 = fittedParameters(2:end);
            %fittedModel.(fieldNames{i}) = fittedParameters(i);


            fittedModel.Reff = max(eig(model.C .* (1-model.VE)) * fittedModel.q ./ fittedModel.gamma);
            fittedModel.R0 = max(eig(model.C) * fittedModel.q ./ fittedModel.gamma);
            fittedModel.t = t;
            fittedModel.x = x;
            fittedModel.tInit = t(1);
            fittedModel.tFinal = t(end);
            fittedModel.xInit = x(1, :);
            fittedModel.xFinal = x(end, :);

        end

        % ---------------------------------------------------------------------------------------------------------
        function [t,x] = predictModel(model, tSpan, xInit)
            odefun = @(t, x)dxdt(model, t, x);
            [t,x] = ode15s(odefun, tSpan, xInit);
        end

        % --------------------------------------------------------------------------------------------------------
        function [weightedAccumulativeCases] = weightedAccumulativeCases(model, xInit, tSpan, weights)
            [t,x] = predictModel(model, tSpan, xInit);
            accumulativeCases = trapz(t, model.gamma' .* x(:,2*model.n+1:3*model.n), 1);
            weightedAccumulativeCases = (accumulativeCases(:)).' * weights(:);

        end

        % --------------------------------------------------------------------------------------------------------

        function d = peak(model, xInit, weights)
            maxDuration = 1e5;
            tSpan = [0, maxDuration];
            [t,x] = predictModel(model, tSpan, xInit);

            idx = 1:model.n;
            omegaE = sum(model.omega .* x(:,idx + 1*model.n) .* weights(:)', 2);
            d = find(omegaE == max(omegaE));
            d = t(d(end));
        end

        % --------------------------------------------------------------------------------------------------------

        function d = last(model, xInit, weights)
            maxDuration = 1e4;
            tSpan = [0, maxDuration];
            [t,x] = predictModel(model, tSpan, xInit);

            idx = 1:model.n;
            omegaE = sum(model.omega.' .* x(:,idx + 1*model.n) .* weights(:)',2); % daily new incidence
            idx = find(omegaE > 1);
            if isempty(idx)
                d = maxDuration;
                return;
            end

            idx = idx(end);
            d = t(idx);
            if d == 0 
                d = maxDuration;
            end
        end

        % --------------------------------------------------------------------------------------------------------

        function [modelList, t, x] = piecewiseFit(model, tData, xData, breakPoints)

            tStart = min(tData);
            tEnd = max(tData);
            breakPoints(breakPoints <= tStart | breakPoints >= tEnd) = [];
            breakPoints = [tStart; breakPoints; tEnd];
            segmentCount = numel(breakPoints) - 1;

            modelList = repmat(model, [segmentCount, 1]);
            t = [];
            x = [];
            for i = 1:segmentCount
                % data for this segment
                idx = tData >= breakPoints(i) & tData <= breakPoints(i+1);
                tDatai = tData(idx);
                xDatai = xData(idx, :);

                % init values for this segment
                model0 = model;
                model0.tInit = breakPoints(i);
                model0.tFinal = breakPoints(i+1);
                if i > 1
                    % model0.xInit = model.xInit for i == 1
                    model0.xInit = modelList(i-1).xFinal;
                end
                [modeli, ti, xi] = fitModel(model0, tDatai, xDatai);

                modelList(i) = modeli;
                if i == 1
                    t = [t; ti];
                    x = [x; xi];
                else
                    t = [t; ti(2:end)];
                    x = [x; xi(2:end, :)];
                end
            end





        end

        % --------------------------------------------------------------------------------------------------------

        function dI = dailyNewIncidence(model, x)
            dI = model.omega' .* x(:, 1:model.n + model.n); % omega' .* E
        end


    end  % end methods





end % end classdef



