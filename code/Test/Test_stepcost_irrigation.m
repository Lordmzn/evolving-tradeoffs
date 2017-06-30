function result = Test_stepcost_irrigation(nextState, otherVars)

persistent env
if isempty(env)
    env = Test_environment();
end

result = max(env.irrigationDemand - otherVars(1), 0).^2;