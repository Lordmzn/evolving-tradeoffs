function [objectives, states, otherVars, decisions] = Test_simulateSystem(...
    initialState, disturbanceRealization, policy, agents)
%EVOLVESYSTEM simulates deterministically how the lake system evolves 
%__________________________________________________________________________
%   INPUTS: 
%   - initial_state:            starting lake state  
%   - disturbanceRealization:	inflows for the horizon
%   - policy:                   struct containing .evaluate(state)
%                               field which returns the action
%   - agents:                   struct array of objects representing ...
%                               operational objectives
%__________________________________________________________________________
%   OUTPUTS:
%   - objectives:   row vector of objective values, a column per each agent
%   - states:       array with the level of every day in the block
%   - otherVars:     matrix, each row contains <release, all stepcosts>
%   - decisions:      decisions taken
%__________________________________________________________________________
% =========================================================================
% DESCRIPTION:
%
% It simulates the system evolution, according to mass balance equation and
% given release policy.
%
% AUTHOR: Paolo Gazzotti, Emanuele Mason
% Politecnico Di Milano
%==========================================================================

verbose = false;

% some simulation parameter
env = Test_environment();

nSteps = length(disturbanceRealization); % first disturbance isn't used!
agentsNames = fieldnames(agents);

% row alignment indicates the timestep at which variable values is known
states = nan(nSteps, 1);
decisions = nan(nSteps, 1);
otherVars = nan(nSteps, 1 + length(agentsNames)); % releases and stepcosts

% initialization
states(1, :) = initialState;

for step = 1:nSteps - 1
    
    inflow = disturbanceRealization(step + 1);
    
    decision = policy.evaluate(policy, states(step));
    decisions(step) = decision; % save the action
    
    % Minimum and maximum release for current storage:
    min_release = max(states(step) - 100, 0);
    max_release = states(step);
    if min_release > decision || max_release < decision
        decision = max(min_release, min(max_release, decision));
    end
    otherVars(step + 1, 1) = decision; % save the release
    
    % Transition dynamic
    states(step + 1) = states(step) + env.stepLength / env.surface * ...
        (inflow - decision);
    
    % stepcosts
    for obj = 1:length(agentsNames)
        otherVars(step + 1, 1 + obj) = ...
            agents.(agentsNames{obj}).evaluateStepCost(states(step + 1), ...
            otherVars(step + 1, 1));
    end
    
    if verbose
        disp(['step: ' num2str(step) '; xt: ' num2str(states(step)) '; ut: ' ...
            num2str(decisions(step)) '; rt+1: ' num2str(otherVars(step+1,1)) ...
            '; xt+1: ' num2str(states(step+1)) ';gt+1: ' num2str(otherVars(step+1,2:end))]);
    end
end

% last call before the bar closes
objectives = nan(1, length(agentsNames));
for obj = 1:length(agentsNames)
    objectives(1, obj) = agents.(agentsNames{obj}).evaluateObjective(...
        otherVars(2:end, 1 + obj));
end

