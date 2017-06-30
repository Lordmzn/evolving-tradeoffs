function result = Test_stepcost_flooding(nextState, otherVars)

persistent env
if isempty(env)
    env = Test_environment();
end

result = max(nextState - env.floodingThreshold, 0);
    