function param = Test_generatePhysicalScenario(param, climate_change)

% the inflows: feature based generation
param.inflows(1) = 50;
features = {'Wet', 'Normal', 'Dry'};
featureProbability = {0.1, 0.75, 0.15};
featureCDF = cumsum([featureProbability{:}]);
% b = log q_t+1 / log q_t 
% alpha = autocorrelation coefficient
% muq, sigmaq = average inflow value and (standard deviation in the log)
w.b = .5; w.alpha = .75; w.muq = 75; w.sigmaq = .75;
n.b = .7; n.alpha = .75; n.muq = 40; n.sigmaq = .65;
d.b = .9; d.alpha = .75; d.muq = 15; d.sigmaq = .3;
featureParams = {w, n, d};
clear w n d

if nargin > 1 && climate_change
    % the inflows: feature based generation
    param.inflows(1) = 50;
    features = {'Wet', 'Normal', 'Dry'};
    featureProbability = {0.2, 0.5, 0.3};
    featureCDF = cumsum([featureProbability{:}]);
    % b = log q_t+1 / log q_t
    % alpha = autocorrelation coefficient
    % muq, sigmaq = average inflow value and (standard deviation in the log)
    w.b = .75; w.alpha = .75; w.muq = 50; w.sigmaq = .8;
    n.b = .7; n.alpha = .75; n.muq = 30; n.sigmaq = .65;
    d.b = .9; d.alpha = .75; d.muq = 10; d.sigmaq = .5;
    featureParams = {w, n, d};
    clear w n d
end

rng(1); % for repeatibility
nBlocks = param.horizon / param.blockSize;
% features
if isfield(param, 'featureScenario')
    featureScenario = param.featureScenario
else
    for block = 1:nBlocks
        num = rand();
        thisFeatureIndex = find(num < featureCDF);
        featureScenario(block) = thisFeatureIndex(1);
    end
end
% actual inflows
for block = 1:nBlocks
    blockStart = (block - 1) * param.blockSize + 2;
    blockEnd = block * param.blockSize + 1;
    b = featureParams{featureScenario(block)}.b;
    alpha = featureParams{featureScenario(block)}.alpha;
    mu = featureParams{featureScenario(block)}.muq;
    sigma = featureParams{featureScenario(block)}.sigmaq;
    for i = blockStart:blockEnd
        % logN ARX(1) - Harms and Campbell, 1967
        param.inflows(i,1) = exp(log(mu) + b * (log(param.inflows(i-1,1)) - ...
            log(mu)) + sigma * randn(1) * sqrt(1 - alpha.^2));
        if param.inflows(i,1) < 0; param.inflows(i,1) = 0; end
    end
end

param.physicalScenario.featureScenario = featureScenario;
param.physicalScenario.features = features;
param.physicalScenario.featureProbability = featureProbability;