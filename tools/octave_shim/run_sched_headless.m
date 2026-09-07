pkg load statistics
cd(getenv('SHIM')); addpath('/home/user/Spoof-Detection'); warning('off', 'all');
tic; run_scheduler_test; fprintf('[elapsed %.1f s]\n', toc);
