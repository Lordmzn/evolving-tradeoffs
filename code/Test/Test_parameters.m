function parameters = Test_parameters(parameters)

% add specific fields to parameters struct
parameters.blockSize = 50;
parameters.s0 = 50;
parameters.horizon = 50000;

% indicative number of steps in each negotiation
parameters.averageT = 50;

% update theta -> f(J) before each negotiation?
parameters.movingUtility = false;

% the policy functional shape
parameters.policy.evaluate = @Test_SOP;
parameters.thetaLimits = [0 200; 0 200];

% list the stepcosts here

% linear threshold overtopping
stepcost = [];
stepcost.epsilon = 10^-1;
stepcost.initialConcessionCoefficient = 0.1;
stepcost.concessionCoefficientLimits = [5e-03 .5];
% stepcost.concessionCoefficientLimits = [1e-01 .4];
stepcost.stepcost = @Test_stepcost_flooding;
stepcost.objective = @(stepcosts) mean(stepcosts(~isnan(stepcosts)));
parameters.agents.('flooding') = stepcost;

% squared irrigation deficit
stepcost = [];
stepcost.epsilon = 10^0;
stepcost.initialConcessionCoefficient = 0.1;
stepcost.concessionCoefficientLimits = [5e-03 2.5];
% stepcost.concessionCoefficientLimits = [5e-01 2];
stepcost.stepcost = @Test_stepcost_irrigation;
stepcost.objective = @(stepcosts) mean(stepcosts(~isnan(stepcosts)));
parameters.agents.('irrigation') = stepcost;

% outputting behavior
parameters.plotJvsTheta = 0; % 0 never; 1 only the 1st time; 2 always
parameters.verboseMediator = false;

% return the same struct
