pkg load statistics
cd(getenv('SHIM')); addpath([getenv('SHIM') '/../..']); warning('off','all');
runAllTests;
