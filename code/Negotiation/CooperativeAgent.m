classdef CooperativeAgent
    
    properties
        name = 'generic';
        env; % = [];
        settings; % = [];
        
        thetaLimits;
        
        concessionCoefficient; % initial value % the alpha
        concessionCoefficientLimits;
        
        % parameters of the transition
        memory = 0.0;
        areObjectivesUpdated = false;
        
        % set at runtime
        stepcost;
        evaluateObjective;
        objectiveThreshold = 0.0; % the Gamma of gazzotti
        last_best;
        
        % the function J->theta
        theta;
        objectives;
        
        % parameters of the search for the function J -> theta
        thetaResolution =  36; % how many point for each theta
        thetaSteps;
        
        % plot 
        theAx;
        theAxIdx = 0;
    end
    
    methods
        
        function obj = CooperativeAgent(name, param)
            obj.name = name;
            obj.settings = param;
            obj.env = feval([param.system '_environment']);
            obj.thetaLimits = obj.settings.thetaLimits;
            obj.concessionCoefficient = obj.settings.agents.(name).initialConcessionCoefficient;
            obj.concessionCoefficientLimits = obj.settings.agents.(name).concessionCoefficientLimits;
            obj.stepcost = obj.settings.agents.(name).stepcost;
            obj.evaluateObjective = obj.settings.agents.(name).objective;
            
            ranges = abs(diff(obj.thetaLimits, 1, 2));
            obj.thetaSteps = ranges / obj.thetaResolution;
        end
            
        function obj = getStarted(obj, inflows, s0)
            % build policy struct
            policy.evaluate = obj.settings.policy.evaluate;
            
            % check if theta -> f(J) has to be updated
            if ~obj.settings.movingUtility
                if ~isempty(obj.theta) && ~isnan(obj.theta(1,1))% is not the first negotiation
                    disp(obj.theta)
                    if isempty(obj.theta) || isnan(obj.theta(1,1)); disp('true'); else; disp('false'); end
                    return % no need to recompute
                end
            end
            
            % first round: build theta sample
            for th = 1:size(obj.thetaLimits, 1)
                thetaValues{th} = obj.thetaLimits(th,1):obj.thetaSteps(th):obj.thetaLimits(th,2);
            end
            thetaGrid = cell(1, numel(thetaValues));
            [thetaGrid{:}] = ndgrid(thetaValues{:});
            thetaToTry = [];
            for th = 1:length(thetaGrid)
                thetaToTry = [thetaToTry thetaGrid{th}(:)];
            end
            objectivesToTry = nan(size(thetaToTry, 1), 1);
            
            % simulate all them
            agent.(obj.name) = obj;
            thetaTic = tic;
            for i = 1:size(thetaToTry, 1)
                policy.theta = thetaToTry(i, :);
                sub_inflows = inflows(1:min([10000 length(inflows)]));
                objectivesToTry(i) = Test_simulateSystem(s0, sub_inflows, ...
                    policy, agent);
                if mod(i, min(round(size(thetaToTry, 1)/10), 5000)) == 0
                    disp(['CooperativeAgent.getStarted: name = ' obj.name ...
                        '; thetaToTry n' num2str(i) ' = ' num2str(toc(thetaTic)) ' [s]'])
                end
            end
            
            % sort (ascending!!)
            [obj.objectives, indexes] = sort(objectivesToTry);
            obj.theta = thetaToTry(indexes, :);
            obj.areObjectivesUpdated = true;
            
            % plotting
            if obj.settings.plotJvsTheta == 2 || obj.settings.plotJvsTheta == 1
                plotJvsTheta(obj)
                if obj.settings.plotJvsTheta == 1
                   obj.settings.plotJvsTheta = 0; % 1 plot only 
                end
            end
        end
        
        function [obj, proposal] = getFirstProposal(obj)
            the_mins = find(obj.objectives == obj.objectives(1));
            proposal = [obj.theta(the_mins,:) obj.objectives(the_mins)];
            obj.objectiveThreshold = obj.objectives(1);
        end
        
        function [obj, proposal] = getNewProposal(obj)
            % update threshold = oldThreshold + alpha
            obj.objectiveThreshold = obj.objectiveThreshold + obj.concessionCoefficient;
            % look for the proposals
            toPropose = obj.objectives < obj.objectiveThreshold;
            proposal = [obj.theta(toPropose, :) obj.objectives(toPropose)];
        end
        
        function obj = updateConcessionCoefficient(obj, last, agreement)
            if ~obj.areObjectivesUpdated
                warning('objs not updated'); return;
            end % update alpha only if obj is updated
            T = obj.settings.averageT;
            if isempty(obj.last_best) % 1st negotiation
                obj.last_best = obj.objectives(1);
            end
            best = obj.last_best;
            obj.last_best = obj.objectives(1); % update
            % memory
            mem = obj.memory * obj.concessionCoefficient + (1 - obj.memory) ...
                * (last - best) / T;
            
            % update alpha
            obj.concessionCoefficient = mem;
            if obj.concessionCoefficient < obj.concessionCoefficientLimits(1)
                obj.concessionCoefficient = obj.concessionCoefficientLimits(1);
            end
            if isfield(obj.settings, 'movingUtility') && ...
                    obj.settings.movingUtility
                obj.areObjectivesUpdated = false;
            end
            
            % logging (ascii=9 => tab)
            disp([obj.name 9 ' alpha = ' num2str(obj.concessionCoefficient)]);
            
        end
        
        function obj = resetState(obj)
            obj.concessionCoefficient = obj.settings.agents.(obj.name).initialConcessionCoefficient;
            obj.theta = [];
            obj.objectives = [];
            obj.last_best = [];
            obj.areObjectivesUpdated = false;
            % parameters of the transition
            obj.memory = 0.0;
            obj.objectiveThreshold = 0.0; % the Gamma of gazzotti
            % plot
            obj.theAx = [];
            obj.theAxIdx = 0;
        end

        function cost = evaluateStepCost(obj, nextState, otherVars)
            cost = obj.stepcost(nextState, otherVars);
        end

    end
end
