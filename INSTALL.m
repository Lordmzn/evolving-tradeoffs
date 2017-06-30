disp('Installing directories...')
addpath(genpath('.'))
disp('Done!')

disp('Checking compiled code...')
if ~exist(['utils/paretofront.' mexext], 'file')
    cd utils
    mex paretofront.c
    cd ..
end
disp('Done!')