function env = Test_environment()

% physical system
env.surface           = 1;    % reservoir surface [m2]
env.stepLength        = 1;    % conversion factor from m3/s -> m3/step

% objectives
env.irrigationDemand  = 50;   % water demand [m3/s]
env.floodingThreshold = 110;  % flooding threshold [m]