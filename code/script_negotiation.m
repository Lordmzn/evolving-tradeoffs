%% GLOBAL SETTINGS 
set(0, 'DefaultFigureWindowStyle', 'docked') % remove when exporting

experimentId = '';

param.system = 'Test';
param = feval([param.system '_parameters'], param);

% choose a specific scenario combination
load physicalScenario_long.mat
load decisionScenario_long.mat
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
averageWeights = nanmean(decisionScenario.referenceWeights, 2); % target
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
clear agent

mediator = Mediator(param.system);

% create policy struct
policy.theta = nan(1, nTheta);
policy.evaluate = param.policy.evaluate;

block = 1;
agreements = nan(param.nBlocks, nTheta + nAgents * 2);
agreements(1,1:nTheta) = param.initialTheta;

%% ======================== BIG SIMULATION BLOCK ==========================
timeToGoHome = false;
maxTrials    = 9;

MSE                       = nan(maxTrials, 1);
r2                        = nan(maxTrials, 1);
runTimes                  = nan(maxTrials, 1);

allAgreements             = cell(maxTrials, 1);
allStates                 = cell(maxTrials, 1);
allDecisions              = cell(maxTrials, 1);
allOtherVariables         = cell(maxTrials, 1);
allT                      = nan(maxTrials, param.nBlocks);
allConcessionCoefficients = cell(maxTrials, 1);

% for now, no search, just exhaustive search within the discretized domain
memory = linspace(1e-2, 1, maxTrials);

trial = 0;
while ~timeToGoHome
    
    memoryValueTic = tic;
    
    trial = trial + 1;
    
    % resetting prior negotiation
    concessionCoefficients = [];
    for name = agentsNames'
        agents.(char(name)) = agents.(char(name)).resetState();
        concessionCoefficients.(char(name)) = [];
    end
    block = 1;
    
    % set the parameter of alpha dynamic equation
    for name = agentsNames'
        agents.(char(name)).memory = memory(trial);
    end
    
    % ..................... FIRST SYSTEM SIMULATION ..........................
    inflows = inflowScenario((block - 1) * param.blockSize + 1 : ...
        block * param.blockSize + 1,1);
    
    % make space
    states         = nan(size(inflowScenario, 1), 1);
    otherVariables = nan(size(inflowScenario, 1), 1 + nAgents);
    decisions      = nan(size(inflowScenario, 1), 1);
    
    % first system simulation
    policy.theta = agreements(1, 1:nTheta);
    [agreements(1, nTheta + nAgents + 1:end), states(1:param.blockSize+1,:), ...
        otherVariables(1:param.blockSize+1,:), decisions(1:param.blockSize+1,:)] = ...
        feval([param.system '_simulateSystem'], param.s0, inflows, policy, ...
        agents);
    % copy objective values as if they were agreed upon
    agreements(1, nTheta + 1:nTheta + nAgents) = agreements(1, nTheta + nAgents + 1:end);
    
    while block < param.nBlocks
        
        block = block + 1;
        
        blockTic = tic;
        
        % log
        disp('__________________________________');
        disp([num2str(block) 9 'block'])
        
        % re initialize the agents for negotiation
        % use last block final state as next block initial state
        initialState = states((block-1) * param.blockSize+1, :);
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
        end
        
        % update alpha: f(last objective, last agreement) iff is not the
        % first block
        if block > 2
            for name = 1:nAgents
                agents.(agentsNames{name}) = agents.(agentsNames{name})...
                    .updateConcessionCoefficient(agreements(block - 1, nTheta + ...
                    nAgents + name), agreements(block - 1, nTheta + name));
            end
        end
        
        % save concessionCoefficients
        for name = agentsNames'
            concessionCoefficients.(char(name))(end+1) = ...
                agents.(char(name)).concessionCoefficient;
        end
        
        [agreements(block, 1:nTheta + nAgents), allT(trial, block)] = ...
            mediator.negotiate(agents);
        
        % select new block inflow values
        inflows = inflowScenario((block - 1) * param.blockSize + 1 : ...
            block * param.blockSize + 1,1);
        
        % use last block final state as next block initial state for system sim
        policy.theta = agreements(block, 1:nTheta);
        [thisPerformances, thisStates, thisOtherVars, thisDecisions] = ...
            feval([param.system '_simulateSystem'], initialState, inflows, ...
            policy, agents);
        
        % store values
        agreements(block, nTheta + nAgents + 1:end) = thisPerformances;
        states((block - 1) * param.blockSize + 2 : block * ...
            param.blockSize + 1, :) = thisStates(2:end, :); % first state is already there
        otherVariables((block - 1) * param.blockSize + 2 : block * ...
            param.blockSize + 1, :) = thisOtherVars(2:end, :);
        decisions((block - 1) * param.blockSize + 1 : block * ...
            param.blockSize, :) = thisDecisions(1:end-1, :);
        
        disp(['Run time of block #' num2str(block) ' [s]: ' ...
            num2str(toc(blockTic)) ])
        
    end
        
    % save last concessionCoefficients
    for name = agentsNames'
        concessionCoefficients.(char(name))(end+1) = ...
            agents.(char(name)).concessionCoefficient;
    end
    
    % evaluate MSE
    horizon = param.nBlocks * param.blockSize;
    MSE(trial) = mean((decisions(1:horizon) - decisionScenario.decisions(1:horizon)').^2); 
    r2(trial) = 1 - MSE(trial) / var(decisionScenario.decisions(1:horizon));
    
    % save stuff
    allAgreements{trial} = agreements;
    allStates{trial} = states;
    allOtherVariables{trial} = otherVariables;
    allDecisions{trial} = decisions;
    allConcessionCoefficients{trial} = concessionCoefficients;
    
    % check if it's beer time
    timeToGoHome = trial == maxTrials;
    runTimes(trial) = toc(memoryValueTic);
    
    disp(['Run time for memory=' num2str(trial) ' [s]: ' ...
        num2str(runTimes(trial)) ]);
    
end


%% saving

save(['data/negotiation_' experimentId '_long_mem' num2str(min(memory)) '_' ...
    num2str(max(memory)) '_thetaRes' ...
    num2str(agents.flooding.thetaResolution) '.mat']);

