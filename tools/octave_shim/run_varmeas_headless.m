pkg load statistics
cd(getenv('SHIM')); addpath([getenv('SHIM') '/../..']); addpath([[getenv('SHIM') '/../..'] '/tests']); warning('off','all');
test_variable_numMeas;
