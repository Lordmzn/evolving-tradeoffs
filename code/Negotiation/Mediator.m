classdef Mediator  
    properties
        nTheta = 1;
        verbose = true;
    end
    
    methods
        
        function obj = Mediator(domain)
            settings = feval([domain '_parameters']);
            obj.nTheta = size(settings.thetaLimits, 1);
            obj.verbose = settings.verboseMediator;
        end
        
        function [agreement, tau, allProposals] = negotiate(obj, agents)
            % reset agents for negotiation and retrieve first proposals
            agentsNames = fieldnames(agents);
            for i = 1:length(agentsNames)
                [agents.(agentsNames{i}), proposals.(agentsNames{i})] = ...
                    agents.(agentsNames{i}).getFirstProposal();
            end
            
            % store proposals in case they're asked
            if nargout > 2
                for i = 1:length(agentsNames)
                    allProposals.(agentsNames{i}){1} = proposals.(agentsNames{i});
                end
            end
            
            % proposal.(agentsName).theta; proposal.(agentsName).objective;
            agreement = obj.isAgreement(proposals);
            
            tau = 1;
            while isempty(agreement)
                tau = tau + 1;
                for i = 1:length(agentsNames)
                    [agents.(agentsNames{i}), proposals.(agentsNames{i})] = ...
                        agents.(agentsNames{i}).getNewProposal();
                end
                
                if nargout > 2
                    for i = 1:length(agentsNames)
                        allProposals.(agentsNames{i}){tau} = proposals.(agentsNames{i});
                    end
                end
                
                if obj.verbose
                    disp(['Negotiation #' num2str(tau) ': Flood=' ...
                        num2str(max(proposals.(agentsNames{1})(:, end))) ...
                        '; ' 9 'Irrigation=' num2str(max(proposals.(agentsNames{2})(:, end))) ]);
                end
                agreement = obj.isAgreement(proposals);
            end
            
            % select only one agreement
            agreement = obj.selectFinalAgreement(agreement);
            
            % logging
            disp(['Agreement at \tau = ' num2str(tau)]);
            
        end
    end
    
    methods (Access = private)
        
        function agreement = isAgreement(obj, proposals)
            % no guard against the case proposal has just one field,
            % corresponding to proposals from only one agent ...            
            agreement = [];
            agentsNames = fieldnames(proposals);
            if length(agentsNames) == 2
                % check for each theta of the first, there's a
                % correspondence in the second agent
                for row = 1:size(proposals.(agentsNames{1}), 1)
                    [isAnAgreement, atIndex] = ismember(...
                        proposals.(agentsNames{1})(row, 1:obj.nTheta), ...
                        proposals.(agentsNames{2})(:, 1:obj.nTheta), 'rows');
                    if isAnAgreement
                        % [parameters obj1 obj2]
                        agreement(end + 1, :) = [...
                            proposals.(agentsNames{1})(row, 1:obj.nTheta) ...
                            proposals.(agentsNames{1})(row, end) ...
                            proposals.(agentsNames{2})(atIndex, end)];
                    end
                end
            else
                % copy all agents but the first
                for agent=2:length(agentsNames)
                    subsetProposal.(agentsNames{agent}) = ...
                        proposals.(agentsNames{agent});
                end
                % recurse
                subsetAgreement = obj.isAgreement(subsetProposal);
                % check for each theta of the agreement, there's a
                % correspondence in the agent excluded
                for row = 1:size(subsetAgreement, 1)
                    [isAnAgreement, atIndex] = ismember(...
                        subsetAgreement(row, 1:obj.nTheta), ...
                        proposals.(agentsNames{1})(:, 1:obj.nTheta), 'rows');
                    if isAnAgreement
                        agreement(end + 1, :) = [subsetAgreement(row, :) ...
                            proposals.(agentsNames{1})(atIndex, end)];
                    end
                end
            end
        end
        
        function agreement = selectFinalAgreement(obj, agreement)
            % select pareto optimal solutions
            isParetoOptimal = paretofront(agreement(:, obj.nTheta+1:end));
            agreement = agreement(isParetoOptimal, :);
            productOfObjectives = agreement(:, obj.nTheta+1);
            for i = obj.nTheta+1 + 1:size(agreement,2) % skip 1st obj
                productOfObjectives = productOfObjectives .* agreement(:, i); 
            end
            [~, nashSolution] = max(productOfObjectives);
            agreement = agreement(nashSolution, :);
        end
        
    end
end