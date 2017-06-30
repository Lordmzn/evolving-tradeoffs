% AUTHOR: Emanuele Mason, Politecnico Di Milano

%% GLOBAL SETTINGS
set(0, 'DefaultFigureWindowStyle', 'docked') % remove when exporting

experimentId = 'SOP_inflowMitFeatures2_fxUtil_few';

param.system = 'Test';
param = feval([param.system '_parameters'], param);

% overwriting because i'm running other exp atm
param.policy.evaluate = @Test_SOP;
param.thetaLimits = [0 200; 0 200];

% choose a specific scenario combination
load SOP_inflowMitFeatures2_physicalScenario.mat
load SOP_inflowMitFeatures2_decisionScenario.mat
clear paramUsed

inflowScenario = inflows; clear inflows
param.inflows = inflowScenario;
param.horizon = length(param.inflows) - 1;
param.nBlocks = floor(size(inflowScenario,1) / param.blockSize);
nTheta = length(param.thetaLimits);
param.initialTheta = decisionScenario.policies(:, 1);

agentsNames = fieldnames(param.agents);
nAgents = length(agentsNames);
% calculate initial concession coefficients
% initial concession coefficients: those that in T steps and starting from
% the best possible state of the world (SOW) for the agent, allow to achieve
% the averageWeight SOW w.r.t. the full range of SOWs
T = param.averageT; % indicative number of negotiation steps
averageWeights = mean(decisionScenario.referenceWeights, 2); % target
for name = 1:nAgents
    param.agents.(agentsNames{name}).initialConcessionCoefficient = 1 / T * ...
        averageWeights(name) * range(decisionScenario.blockObjectives(name,:));
end
clear theBest T

% Create agents
agents = struct();
for name = 1:nAgents
    agent = CooperativeAgent(agentsNames{name}, param);
    agents.(agent.name) = agent;
end
disp(['flooding' 9 ' alpha = ' num2str(agents.flooding.concessionCoefficient)]);
disp(['irrigation' 9 ' alpha = ' num2str(agents.irrigation.concessionCoefficient)]);
if nAgents > 2
    disp(['reliability' 9 ' alpha = ' num2str(agents.reliability.concessionCoefficient)]);
end
clear agent

mediator = Mediator(param.system);

% create policy struct
policy.theta = nan(1, nTheta);
policy.evaluate = param.policy.evaluate;

% exhaustive search of alpha values
granularity = 3; %75;
values = [linspace(param.agents.flooding.concessionCoefficientLimits(1), ...
    param.agents.flooding.concessionCoefficientLimits(2), granularity); ...
    linspace(param.agents.irrigation.concessionCoefficientLimits(1), ...
    param.agents.irrigation.concessionCoefficientLimits(2), granularity)];
if nAgents > 2
    values = [values; ...
        linspace(param.agents.reliability.concessionCoefficientLimits(1), ...
        param.agents.reliability.concessionCoefficientLimits(2), granularity)];
end

concessionCoefficientsValues = [];
for x = 1 : granularity
    for y = 1 : granularity
        if nAgents > 2
            for z = 1 : granularity
                concessionCoefficientsValues(end+1, :) = [values(1,x) ...
                    values(2,y) values(3,z)];
            end
        else
            concessionCoefficientsValues(end+1, :) = [values(1,x) values(2,y)];
        end
    end
end
nAlpha = length(concessionCoefficientsValues);
clear values

%% ======================== BIG SIMULATION BLOCK ==========================
runTimes = nan(param.nBlocks, 1);

% resetting prior negotiation
concessionCoefficients = [];
for name = agentsNames'
    agents.(char(name)) = agents.(char(name)).resetState();
    concessionCoefficients.(char(name)) = [];
end
block = 1;

% ..................... FIRST SYSTEM SIMULATION ..........................
inflows = inflowScenario((block - 1) * param.blockSize + 1 : ...
    block * param.blockSize + 1,1);

% first system simulation 
policy.theta = param.initialTheta;
[initObj, initStates, initOther, initDecisions] = ...
    feval([param.system '_simulateSystem'], param.s0, inflows, policy, ...
    agents);

% make space
allAgreements = cell(nAlpha, 1);
allStates = cell(nAlpha, 1);
allDecisions = cell(nAlpha, 1);
allOtherVariables = cell(nAlpha, 1);
allT = nan(nAlpha, param.nBlocks);
allBest = nan(nAgents, param.nBlocks); % doesnt depend on alpha
MSE = nan(nAlpha, param.nBlocks);
r2 = nan(nAlpha, param.nBlocks);

for alphaValue = 1:nAlpha
    
    % make space for simulation
    allAgreements{alphaValue} = nan(param.nBlocks, nTheta + nAgents*2);
    allStates{alphaValue} = nan(size(inflowScenario, 1), 1);
    allDecisions{alphaValue} = nan(size(inflowScenario, 1), 1);
    allOtherVariables{alphaValue} = nan(size(inflowScenario, 1), 1 + nAgents);

    % copy first block
    % copy objective values as if they were agreed upon
    allAgreements{alphaValue}(1, 1:nTheta) = param.initialTheta;
    allAgreements{alphaValue}(1, nTheta+1:nTheta+nAgents) = initObj;
    allAgreements{alphaValue}(1, nTheta+nAgents+1:end) = initObj;
    allStates{alphaValue}(1:param.blockSize + 1,:) = initStates;
    allDecisions{alphaValue}(1:param.blockSize,:) = initDecisions(1:end-1);
    allOtherVariables{alphaValue}(1:param.blockSize + 1,:) = initOther;
    
end

clear initObj initStates initOther initDecisions

bestR2 = 1; % any value is ok
% start the actual cycle
while block < param.nBlocks
    
    blockTic = tic;

    block = block + 1;
    
    % re initialize the agents for negotiation
    % use bestR2' last block final state as next block initial state
    initialState = allStates{bestR2}((block-1) * param.blockSize+1, :);
    for name = 1:nAgents
        
        % check if utility function has to be updated
        if param.movingUtility
            agents.(agentsNames{name}) = agents.(agentsNames{name}). ...
                getStarted(inflows, initialState);
        elseif block == 2 
            % if it's the first negotiation, all inflows are passed to
            % construct the function
            agents.(agentsNames{name}) = agents.(agentsNames{name}). ...
                getStarted(inflowScenario, initialState);
        end
        
        % save the best
        allBest(name,block) = agents.(agentsNames{name}).objectives(1);
    end
    
    % select new block inflow values
    inflows = inflowScenario((block - 1) * param.blockSize + 1 : ...
        block * param.blockSize + 1,1);
    
    for alphaValue = 1:nAlpha
        
        alphaTic = tic;
        
        % log
        disp('__________________________________');
        disp([num2str(alphaValue) 9 'alpha'])
        
        % initial concession coefficients: pick one from list of exhaustive
        % search
        for name = 1:nAgents
            agents.(agentsNames{name}).concessionCoefficient = ...
                concessionCoefficientsValues(alphaValue, name);
        end
        
        % negotiate
        [allAgreements{alphaValue}(block, 1:nTheta + nAgents), allT(alphaValue, block)] = ...
            mediator.negotiate(agents);
        
        % use bestR2' last block final state as next block initial state for system sim
        policy.theta = allAgreements{alphaValue}(block, 1:nTheta);
        [allAgreements{alphaValue}(block, nTheta + nAgents + 1:end), states, otherVars, ...
            decisions] = feval([param.system '_simulateSystem'], initialState, ...
            inflows, policy, agents);
        
        % store values
        allStates{alphaValue}((block - 1) * param.blockSize + 2 : block * ...
            param.blockSize + 1, :) = states(2:end, :); % first state is already there
        allOtherVariables{alphaValue}((block - 1) * param.blockSize + 2 : ...
            block * param.blockSize + 1, :) = otherVars(2:end, :);
        allDecisions{alphaValue}((block - 1) * param.blockSize + 1 : block * ...
            param.blockSize, :) = decisions(1:end-1, :);
        
        % evaluate MSE
        decisionInTheScenario = decisionScenario.decisions((block - 1) * ...
            param.blockSize + 1 : block * param.blockSize)';
        MSE(alphaValue, block) = mean((decisions(1:end-1) - decisionInTheScenario).^2);
        r2(alphaValue, block) = 1 - MSE(alphaValue, block) / var(decisionInTheScenario);
        
        disp(['Run time of alpha #' num2str(alphaValue) ' [s]: ' ...
            num2str(toc(alphaTic)) ])
        
    end
    
    [~, bestR2] = min(r2(:,block)); % find best alpha set
    initialState = allStates{bestR2}(block * param.blockSize + 1);
    
    runTimes(block) = toc(blockTic);
    disp(['Run time for block =' num2str(block) ' [s]: ' ...
        num2str(runTimes(block)) ]);
    
end


%% saving

save(['A2 - set-based cooperative negotiation protocol/' ...
    param.system '/Data/' experimentId '_thetaRes' ...
    num2str(agents.flooding.thetaResolution) '_gran' ...
    num2str(granularity) '_negotiation.mat']);

