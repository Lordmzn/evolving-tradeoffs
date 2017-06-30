%% build Pareto Front for inflows
set(0, 'DefaultFigureWindowStyle', 'docked') % remove when exporting

experimentId = '';
param.plotFigures = true;

param.system = 'Test';
param.dataFolder = 'data/';
if ~exist(param.dataFolder, 'dir')
    mkdir(param.dataFolder)
    INSTALL;
end
param = feval([param.system '_parameters'], param);

%% build physical scenario or load it
if exist(['physicalScenario' experimentId '.mat'], 'file')
    disp('Loading physical scenario...')
    load(['physicalScenario' experimentId '.mat'])
    param.physicalScenario = physicalScenario;
    param.inflows = inflows;
    clear physicalScenario inflows
else
    disp('Creating physical scenario...')
    param = Test_generatePhysicalScenario(param, true);
    physicalScenario = param.physicalScenario;
    inflows = param.inflows;
    save([param.dataFolder 'physicalScenario' experimentId '.mat'], ...
        'physicalScenario', 'inflows');
	clear physicalScenario inflows
end

% plot to check everything is working
if param.plotFigures
    figure;
    hold on; firstBlock = 19; lastBlock = 41;
    for b = firstBlock:lastBlock
        plot([b * param.blockSize b * param.blockSize], [0 250], '--', ...
            'Color', [.5 .5 .5]); 
        text((b - 1) * param.blockSize + 0.1 * param.blockSize, 245, ...
            param.physicalScenario.features{param.physicalScenario.featureScenario(b)}(1), ...
            'FontName', 'Economica', 'FontSize', 18);
    end
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
    thisAx.XLim = [(firstBlock-1) * param.blockSize lastBlock * param.blockSize];
    thisAx.XTick = (firstBlock-1) * param.blockSize:50:lastBlock * param.blockSize;
    thisAx.YLim = [0 250];
    plot(thisAx.XLim, [50 50], '--r');
    plot(thisAx.XLim, [110 110], '--r');
    plot((firstBlock-1) * param.blockSize + 1:lastBlock * param.blockSize, ...
        param.inflows((firstBlock-1) * param.blockSize + 1:lastBlock * param.blockSize));
    title('Snippet of inflows scenario')
    xlabel('Day in the scenario');
    ylabel('Inflow [m^3/s]');
    hold off;
end


%% Create agents for the stepcosts
agents = struct();
agentsNames = fieldnames(param.agents);
nAgents = length(agentsNames);
for name = 1:nAgents
    agent = CooperativeAgent(agentsNames{name}, param);
    agents.(agent.name) = agent;
end
param.agents = agents;
clear agent name

% create policy struct
nTheta = size(param.thetaLimits, 1); % NOTE: initial theta not used in this script
param.policy.theta = nan(1, nTheta);

% set up system initial state
param.initialState = param.s0;

%% optimize over the inflows
borgoptimset = {'frequency', 10000, 'rngstate', 1}; % seed
epsilons = [10^-1 10^-1 10^0]; % length = nAgents
nfe = 50000; % 400 = ~ 10' % 500 = 2h 30'?

files = dir(['data/' experimentId 'borgEMODPS*.mat']);
if isempty(files)
    error('Cant find any optimal policu... where are they?')
end

for file = 1:length(files)
    disp(files(file).name)
end
toload = input('which one? write just the _XXNFE part\n','s');
load([experimentId 'borgEMODPS' toload '.mat'])

[objectives, indexes] = sortrows(objectives);
theta = theta(indexes, :);

policy = param.policy;
for pol = 1:size(objectives, 1)
    decisionScenarioBuildingBlocks(pol).theta = theta(pol, :);
    policy.theta = theta(pol, :);
    [decisionScenarioBuildingBlocks(pol).objectives, ~, ...
        decisionScenarioBuildingBlocks(pol).otherVars] = Test_simulateSystem(...
        param.initialState, param.inflows, policy, param.agents);
    if decisionScenarioBuildingBlocks(pol).objectives ~= objectives(pol, :)
        disp(['Pol ' num2str(pol) ' error? Might be okay if pol were opt on a different scenario']);
    end
end

clear theta objectives indexes

if param.plotFigures
    % first is extremal
    policy.theta = decisionScenarioBuildingBlocks(1).theta;
    figure; hold on; grid on; title(['Policy [' num2str(policy.theta) ']']);
    x = 0:.1:200;
    plot(x, x, '--r'); plot(x, x-100, '--r'); % bounds
    plot(x, policy.evaluate(policy, x), '--');
    axis([0 max(x) 0 max(x)]); xlabel('level'); ylabel('decision');
    
    % mid one is a trade-off
    thePol = floor(mean(size(decisionScenarioBuildingBlocks))); % vector not matrix
    policy.theta = decisionScenarioBuildingBlocks(thePol).theta;
    figure; hold on; grid on; title(['Policy [' num2str(policy.theta) ']']);
    x = 0:.1:200;
    plot(x, x, '--r'); plot(x, x-100, '--r'); % bounds
    plot(x, policy.evaluate(policy, x), '--');
    axis([0 max(x) 0 max(x)]); xlabel('level'); ylabel('decision');
    
    % last is another extreme
    policy.theta = decisionScenarioBuildingBlocks(end).theta;
    figure; hold on; grid on; title(['Policy [' num2str(policy.theta) ']']);
    x = 0:.1:200;
    plot(x, x, '--r'); plot(x, x-100, '--r'); % bounds
    plot(x, policy.evaluate(policy, x), '--');
    axis([0 max(x) 0 max(x)]); xlabel('level'); ylabel('decision');
    
end

if param.plotFigures && nAgents == 2
    % plot the PF
    objToPlot = reshape([decisionScenarioBuildingBlocks.objectives], 2, ...
        length(decisionScenarioBuildingBlocks))';
    objToPlot = sortrows(objToPlot);
    figure; 
    plot(objToPlot(:,1), objToPlot(:,2), '*');
    grid on; title('Pareto front of optimal policies over all the inflows')
    xlabel('flooding [m]'); ylabel('irrigation [(m^3/s)^2]');
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
    
    % plot all the policies
    figure; hold on; grid on; 
    theTitle = title('Optimal policies over all inflows');
    theTitle.FontName = 'Economica'; theTitle.FontSize = 24;
    plot([30 30 70 70 30 12 104 70], [30 70 70 30 30 109 109 30], '--.', 'Color', [.6 .6 .6]);
    policy = param.policy;
    objToPlot = reshape([decisionScenarioBuildingBlocks.objectives], 2, ...
        length(decisionScenarioBuildingBlocks))';
    [sortedObjectives, indexes] = sortrows(objToPlot, 1);
    theta = reshape([decisionScenarioBuildingBlocks.theta], ...
        length(decisionScenarioBuildingBlocks(1).theta), ...
        length(decisionScenarioBuildingBlocks))';
    sortedTheta = theta(indexes, :);
    x = 0:.1:200; colors = colormap(parula(size(sortedTheta, 1)));
    plot(x, x, '--r'); plot(x, x-100, '--r'); % bounds
    for pol = 1:size(theta,1)
        policy.theta = sortedTheta(pol,:);
        plot(x, policy.evaluate(policy, x), 'Color', colors(pol,:), ...
            'LineWidth', 2);
    end
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
    plot(xlim, [50 50], '-.r'); plot([110 110], ylim, '-.r'); % thresholds
    axis([0 max(x) 0 max(x)]); xlabel('level'); ylabel('decision');
    theBar = colorbar;
    theBar.Ticks = linspace(0, 1, length(decisionScenarioBuildingBlocks));
    theBar.TickLabels = cellstr(num2str(objToPlot));
    theBar.FontName = 'Economica'; theBar.FontSize = 12;

    % zoom on the part releavant to irrigation
    axes(gcf, 'Position',[0.17,0.55,0.35,0.35]); hold on; grid on; box on;
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
    x = 30:.1:70; colors = colormap(parula(size(sortedTheta, 1)));
    plot(x, x, '--r'); plot(x, x-100, '--r'); % bounds
    for pol = 1:size(theta,1)
        policy.theta = sortedTheta(pol,:);
        plot(x, policy.evaluate(policy, x), 'Color', colors(pol,:), ...
            'LineWidth', 2);
    end
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
    plot(xlim, [50 50], '-.r'); plot([110 110], ylim, '-.r'); % thresholds
    xlabel('level'); ylabel('decision');
    xlim([30 70]); ylim([30 70])
    
end
clear objToPlot theta


%% compose the decision scenario

% a rule: fix the weights, then choose the policy that leads towards these
% weights over the last N blocks

% compute weights = direction in the objective space
% over the full scenario = average
nBlocks = floor(length(param.inflows) / param.blockSize);
nPol = length(decisionScenarioBuildingBlocks);
objs = reshape([decisionScenarioBuildingBlocks.objectives], nAgents, nPol)';
weights = 1 - (objs - repmat(min(objs), nPol, 1)) ./ (repmat(max(objs), ...
    nPol, 1) - repmat(min(objs), nPol, 1));
weights = weights ./ repmat(sum(weights, 2), 1, nAgents); % into the simplex

% calculate objective per block
for i = 1:nPol
    decisionScenarioBuildingBlocks(i).weights = weights(i, :)';
    for agent = 1:nAgents
        stepcostReshaped = reshape(...
            decisionScenarioBuildingBlocks(i).otherVars(2:end, 1 + agent), ...
            param.blockSize, []); % skip first stepcost (g_t+1) is a nan
        for block = 1:nBlocks
            decisionScenarioBuildingBlocks(i).blockObjectives(agent, block) = ...
                param.agents.(agentsNames{agent}). ...
                evaluateObjective(stepcostReshaped(:,block));
        end
    end
end

% then calculate weights. First find bounds, then normalize objs into
% weights
for block = 1:nBlocks
    thisObjs = [];
    for i = 1:nPol
        thisObjs(:, i) = decisionScenarioBuildingBlocks(i).blockObjectives(:, block);
    end
    maxBlockObjs(:, block) = max(thisObjs, [], 2);
    minBlockObjs(:, block) = min(thisObjs, [], 2);
end

clear thisObjs
for i = 1:nPol
    decisionScenarioBuildingBlocks(i).blockWeights = 1 - ...
        (decisionScenarioBuildingBlocks(i).blockObjectives - minBlockObjs) ...
        ./ (maxBlockObjs - minBlockObjs);
    decisionScenarioBuildingBlocks(i).blockWeights(isnan(decisionScenarioBuildingBlocks(i).blockWeights)) = 1;
    % into the simplex
    decisionScenarioBuildingBlocks(i).blockWeights = ...
        decisionScenarioBuildingBlocks(i).blockWeights ./ ...
        repmat(sum(decisionScenarioBuildingBlocks(i).blockWeights, 1), nAgents, 1);
end
clear maxBlockObjs minBlockObjs

if param.plotFigures    
    figure;
    nFigures = min(nPol, 14);
    for i = 1:nFigures
        subplot(nFigures, 1, i);
        hold on
        for agent = 1:nAgents
            plot(decisionScenarioBuildingBlocks(i * ...
                floor(nPol/nFigures)).blockWeights(agent,:));
        end
        ylimits = get(gca,'YLim'); xlimits = get(gca,'XLim');
        text(xlimits(1) - 0.05 * (diff(xlimits)), mean(ylimits), ...
            num2str(decisionScenarioBuildingBlocks(i * ...
            floor(nPol/nFigures)).weights), 'HorizontalAlignment', 'right', ...
            'FontName', 'Economica', 'FontSize', 14)
        thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 14;
    end
    subplot(nFigures, 1, 1); 
    title('Weights evolution for the building block policies')
    ylimits = get(gca,'YLim'); xlimits = get(gca,'XLim');
    text(xlimits(1) - 0.05 * (diff(xlimits)), ylimits(2) + 0.05 * (diff(ylimits)), ...
        'Policy weights:', 'HorizontalAlignment', 'right', 'FontName', ...
        'Economica', 'FontSize', 14)  
    theLeg = legend(agentsNames);
    theLeg.Position(1) = 0.92;
end

%% Here starts the real construction of the scenario
% weights = f( meanPerformance (MC) )
decisionScenario.movingWindow = 5;
decisionScenario.MonteCarloRuns = 5;
decisionScenario.referenceWeights(:, 1) = 0.5 * ones(nAgents,1);
% simulation vars
decisionScenario.states = [];
decisionScenario.otherVars = [];
decisionScenario.decisions = [];
state = param.initialState;
for block = 1:nBlocks
    
    disp([ 'Block #' num2str(block) ])
    
    % init stuff
    blockStart = (block - 1) * param.blockSize + 1;
    blockEnd = block * param.blockSize + 1;
    decisionScenario.blockBounds(:, block) = [blockStart; blockEnd];
    allWeights = []; allObjectives = [];
    if block ~= 1
        % update theReference weights (based on last 5 weighting over actual
        % performances)
        theWeights = decisionScenario.referenceWeights(:, block - ...
            min(decisionScenario(1).movingWindow, block - 1):block - 1);
        ARcoefficients = [.05 .07 .12 .2 .56]';
        if block < 6
            ARcoefficients = ARcoefficients(end-size(theWeights, 2)+1:end);
            ARcoefficients = ARcoefficients ./ sum(ARcoefficients);
        end
        theReference(:, block) = theWeights * ARcoefficients;
        % calculate the weights of the next policy to be chosen, given the 
        % last state and a biased montecarlo pool of inflows:
        selector.Dry = {'Dry', 'Dry', 'Normal', 'Wet'}; % + lastBlock
        selector.Normal = {'Dry', 'Normal', 'Normal', 'Wet'}; % + lastBlock
        selector.Wet = {'Dry', 'Normal', 'Wet', 'Wet'}; % + lastBlock
        for pol = 1:nPol
            policy.theta = decisionScenarioBuildingBlocks(pol).theta;
            % first simulate lastBlock
            [~, ~, temp] = Test_simulateSystem(decisionScenario.states(blockStart - 1), ...
                thisInflows, policy, param.agents);
            allStepcosts = temp(2:end, 2:end);
            % then simulate the others
            for item = 1:decisionScenario.MonteCarloRuns - 1
                thisBlockFeature = param.physicalScenario.features(...
                    param.physicalScenario.featureScenario(block));
                % stupid matlab cells
                thisBlockFeature = thisBlockFeature{1};
                mcInflowsFeature = selector.(thisBlockFeature)(item);
                mcInflowsFeature = mcInflowsFeature{1};
                % pick a random inflow block with that feature
                rBlock = randi(nBlocks, 1);
                rBlockFeature = param.physicalScenario.features(...
                    param.physicalScenario.featureScenario(rBlock));
                rBlockFeature = rBlockFeature{1};
                while ~strcmp(mcInflowsFeature, rBlockFeature)
                    rBlock = randi(nBlocks, 1);
                    rBlockFeature = param.physicalScenario.features(...
                        param.physicalScenario.featureScenario(rBlock));
                    rBlockFeature = rBlockFeature{1};
                end
                mcBlockStart = (rBlock - 1) * param.blockSize + 1;
                mcBlockEnd = rBlock * param.blockSize + 1;
                mcInflows = param.inflows(mcBlockStart:mcBlockEnd);
                [~, ~, temp] = Test_simulateSystem(decisionScenario.states(blockStart - 1), ...
                    mcInflows, policy, param.agents);
                allStepcosts = [allStepcosts; temp(2:end, 2:end)];
            end
            % then evaluate objectives over the stepcosts collected
            for agent = 1:nAgents
                allObjectives(agent, pol) = param.agents.(agentsNames{agent}). ...
                    evaluateObjective(allStepcosts(:, agent));
            end
        end
        % look for extremes
        maxBlockObjs = repmat(max(allObjectives, [], 2), 1, nPol);
        minBlockObjs = repmat(min(allObjectives, [], 2), 1, nPol);
        allWeights = 1 - (allObjectives - minBlockObjs) ./ (maxBlockObjs - minBlockObjs);
        % nan are asshoels
        allWeights(isnan(allWeights)) = 1;
        % into the simplex
        allWeights = allWeights ./ repmat(sum(allWeights, 1), nAgents, 1);
    else % block == 1
        theReference(:, block) = decisionScenario.referenceWeights(:, block);
        % the weights of each policy:
        for i = 1:nPol
            allWeights(:, i) = decisionScenarioBuildingBlocks(i).blockWeights(:, block);
            allObjectives(:, i) = decisionScenarioBuildingBlocks(i).blockObjectives(:, block);
        end
    end
    % calculate distances
    for i = 1:nPol
        distances(:,i) = abs(allWeights(:, i) - theReference(:, block));
    end
    % find the closest pol and save
    [~, thePol(block)] = min(sum(distances.^2, 1).^0.5);
    decisionScenario.policies(:, block) = decisionScenarioBuildingBlocks(thePol(block)).theta;
    theReferenceObjectives(:, block) = allObjectives(:, thePol(block));
    % simulate this policy
    thisInflows = param.inflows(blockStart:blockEnd);
    policy.theta = decisionScenario.policies(:, block);
    [decisionScenario.blockObjectives(:, block), state, others, decisions] = ...
        Test_simulateSystem(state(end), thisInflows, policy, param.agents);
    % store sim stuff
    decisionScenario.states = [decisionScenario.states state(1:end - 1, :)'];
    decisionScenario.otherVars = [decisionScenario.otherVars others(2:end, :)']; % first is nan
    decisionScenario.decisions = [decisionScenario.decisions decisions(1:end - 1, :)']; % last is nan
    % recalculate weights for future updates of the reference according to
    % the sim on the new inflows
    allWeights = [];
    if block ~= 1
        for pol = 1:nPol
            policy.theta = decisionScenarioBuildingBlocks(pol).theta;
            allWeights(:, pol) = Test_simulateSystem(...
                decisionScenario.states(blockStart - 1), thisInflows, policy, param.agents);
        end
    else
        for i = 1:nPol
            allWeights(:, i) = decisionScenarioBuildingBlocks(i).blockWeights(:, block);
        end
    end
    % allWeights now contains objs. look for extremes
    maxBlockObjs = repmat(max(allWeights, [], 2), 1, nPol);
    minBlockObjs = repmat(min(allWeights, [], 2), 1, nPol);
    allWeights = 1 - (allWeights - minBlockObjs) ./ (maxBlockObjs - minBlockObjs);
    % nan are annoying as hell
    allWeights(isnan(allWeights)) = 1;
    % into the simplex
    allWeights = allWeights ./ repmat(sum(allWeights, 1), nAgents, 1);
    decisionScenario.referenceWeights(:, block) = allWeights(:, thePol(block));
end
decisionScenario.theTargetReference = theReference;
decisionScenario.theTargetReferenceObjectives = theReferenceObjectives;

decisionScenario.objectives = mean(decisionScenario.blockObjectives, 2);

%% simulate'n'plot 
if param.plotFigures    
    figure; nPlots = 4;
    chunk = nBlocks / nPlots;
    for i = 1:nPlots
        subplot(nPlots,1,i); 
        plot((i-1)*chunk+1:i*chunk, decisionScenario.referenceWeights(:,(i-1)*chunk+1:i*chunk)');
        thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
        thisAx.XLim = [(i-1)*chunk+1 i*chunk];
        grid on
        if i == 1
            title('Weights evolution for the policy chosen in the scenario')
        end
        ylabel('Weights')
    end
    xlabel('# block')
    legend({'FLO', 'IRR', 'REL'}, 'Location', 'EastOutside')
    
    figure; nPlots = 4;
    chunk = nBlocks / nPlots;
    for i = 1:nPlots
        subplot(1,nPlots,i); 
        plot(flipud(decisionScenario.policies(:,(i-1)*chunk+1:i*chunk)'), ...
            flipud((i-1)*chunk+1:i*chunk));
        thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
        thisAx.YTickLabel = flipud(thisAx.YTickLabel);
        grid on
        if i == 1
            title('Parameter evolution over the blocks')
        end
        ylabel('# block')
    end
    xlabel('parameters')
    legend({'\theta_1', '\theta_2', '\theta_3', '\theta_4'}, ...
                'Location', 'EastOutside')

    figure; hold on; grid on
    title('Decisions VS state assumed in the scenario')
    x = 0:.1:200;
    plot(x, x, '--r'); plot(x, x-100, '--r'); % bounds
    plot(decisionScenario.states, decisionScenario.decisions, '.');
    text(10, max(x) - 10, {
        % ['Objective values with the reference policy ' num2str(decisionScenario.referenceObjectives)]
        ['Objective values after simulation ' num2str(decisionScenario.objectives')]
        });
    axis([0 max(x) 0 max(x)]); xlabel('level'); ylabel('decision');
    
    % policy sequence: find which policy are used and tag them
    [~, id, policyIndex] = unique(decisionScenario.policies', 'rows');
    % create the ordering
    indexes = unique(policyIndex);
    for dude = 1:length(indexes)
        objOfPol(:, dude) = median(decisionScenario.blockObjectives(:, ...
            policyIndex == dude), 2);
    end
    [~, reMap] = sortrows(objOfPol');
    figure;
    nPlots = 5;
    x = 1:length(policyIndex);
    for i = 1:length(x)
        y(i) = find(policyIndex(i) == reMap);
    end
    y = y(:)'; % in row
	nEl = length(policyIndex) / nPlots;
    colors = colormap(jet(length(indexes)));
    for plt = 1:nPlots
        subplot(nPlots, 1, plt);  grid on; hold on
        windowToPlot = false(size(x));
        windowToPlot((plt - 1) * nEl + 1:plt * nEl) = true;
        thePlot = plot(x(windowToPlot), y(windowToPlot));
        thePlot.LineWidth = 2; thePlot.Color = [.5 .5 .5];
        thePlot.Parent.YLim = [min(y) max(y)];
        % toSkip = false(length(indexes), 1); % if no block use this policy
        for i = 1:length(indexes)
            toPlot = y == indexes(i);
            if sum(windowToPlot & toPlot) == 0; toSkip(i) = true; continue; end
            h = plot(x(windowToPlot & toPlot), y(windowToPlot & toPlot), '.');
            h.MarkerSize = 18;
            h.Color = colors(i, :);
        end
        xlabel('# block')
        ylabel('Policy used [FLO IRR]')
        theColorbar = colorbar;
        theColorbar.Location = 'west';
        theColorbar.TickLabels = {};
        if plt == 1; title('Sequence of policy used in the scenario'); end
        h = gca;
        h.YTick = floor(linspace(min(y), max(y), range(y)/10));
        % thePols = decisionScenario.policies(:, id);
        % texts = num2str(thePols(:, ~toSkip)', '%3.0f');
        stuff = objOfPol(:, reMap);
        texts = num2str(stuff(:, h.YTick)', ',% 4.1f');
        texts = [repmat('[ ', length(h.YTick), 1) ...
            texts(:, 3:end) ...
            repmat(' ]', length(h.YTick), 1)];
        h.YTickLabel = mat2cell(texts, ones(size(texts, 1), 1));
        thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 14;
    end
    clear h x y indexes nPlots nEl plt id texts toSkip 
    
    % policy frequencies
    figure; hold on; grid on;
    thePlot = histogram(categorical(policyIndex, unique(policyIndex)));
    thePlot.Orientation = 'horizontal';
    ylabel('Median performance over scenario blocks per policy');
    title('Frequency of occurence of each policy in the scenario')
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 14;
    texts = num2str(objOfPol', ',% 4.1f');
    texts = [repmat('[ ', length(objOfPol), 1) ...
        texts(:, 3:end) ...
        repmat(' ]', length(objOfPol), 1)];
    thisAx.YTickLabel = mat2cell(texts, ones(size(texts, 1), 1));
    
    % plot all the policies
    figure; hold on; grid on; 
    theTitle = title('Policies in the scenario');
    plot([30 30 70 70 30 12 104 70], [30 70 70 30 30 109 109 30], '--.', 'Color', [.6 .6 .6]);
    policy = param.policy;
    [polToPlot, id, policyIndex] = unique(decisionScenario.policies', 'rows', 'stable');
    % create the ordering
    for dude = 1:length(id)
        objOfPol(:, dude) = median(decisionScenario.blockObjectives(:, ...
            policyIndex == dude), 2);
    end
    [~, reMap] = sortrows(objOfPol');
    objOfPol = objOfPol(:, reMap)';
    polToPlot = polToPlot(reMap, :);    
    x = 0:.1:200; colors = colormap(parula(length(id)));
    plot(x, x, '--r'); plot(x, x-100, '--r'); % bounds
    for pol = 1:length(id)
        policy.theta = polToPlot(pol,:);
        plot(x, policy.evaluate(policy, x), 'Color', colors(pol,:), ...
            'LineWidth', 2);
    end
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
    plot(xlim, [50 50], '-.r'); plot([110 110], ylim, '-.r'); % thresholds
    axis([0 max(x) 0 max(x)]); xlabel('level [m]'); ylabel('decision [m^3/s]');
    theBar = colorbar;
    theBar.Ticks = linspace(0, 1, length(id));
    theBar.TickLabels = cellstr(num2str(objOfPol));
    theBar.FontName = 'Economica'; theBar.FontSize = 12;
    ylabel(theBar, 'Median performance of policy across scenario blocks', ...
        'Rotation', 270, 'FontSize', 18, 'VerticalAlignment', 'bottom');

    % zoom on the part releavant to irrigation
    axes(gcf, 'Position',[0.17,0.55,0.35,0.35]); hold on; grid on; box on;
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
    x = 30:.1:70; colors = colormap(parula(length(id)));
    plot(x, x, '--r'); plot(x, x-100, '--r'); % bounds
    for pol = 1:length(id)
        policy.theta = polToPlot(pol,:);
        plot(x, policy.evaluate(policy, x), 'Color', colors(pol,:), ...
            'LineWidth', 2);
    end
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
    plot(xlim, [50 50], '-.r'); plot([110 110], ylim, '-.r'); % thresholds
    xlabel('level [m]'); ylabel('decision [m^3/s]');
    thisAx.XLim = [30 70]; thisAx.YLim = [30 70];
    
    clear polToPlot objOfPol
    
    % the ugly plot
    figure;
    param.horizon = length(param.inflows) - 1;
    subplot(3,1,1); 
    plot(param.inflows);
    ylabel('Inflow [m^3/s]')
    title('Inflow scenario')
    xlim([0 param.horizon]); ylim([0 400]);
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
    subplot(3,1,2); 
    plot(decisionScenario.decisions);
    title('Decision scenario')
    ylabel('Decision [m^3/s]')
    xlim([0 param.horizon]); ylim([0 400]);
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;    
    subplot(3,1,3); 
    plot(decisionScenario.states)
    title('Lake levels')
    xlabel('steps')
    xlim([0 param.horizon]); ylim([0 400]);
    ylabel('Lake level [m]')
    thisAx = gca; thisAx.FontName = 'Economica'; thisAx.FontSize = 18;
end


%% save scenario and details
paramUsed = param;

save([param.dataFolder experimentId 'decisionScenario3.mat'], 'decisionScenario', 'paramUsed');

figlist = findobj('type', 'figure');
savefig(figlist, [param.dataFolder experimentId 'allFiguresScenarioBuilding4.fig']);

