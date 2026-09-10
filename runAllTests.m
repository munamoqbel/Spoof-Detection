function runAllTests()
%RUNALLTESTS  Pre-integration check: run every test in the recommended
%             order and summarise. Replaces the old test_all.m.
%
%   1. tests/test_monitors          engine math (unit + Monte-Carlo)
%   2. tests/test_variable_numMeas  varying satellite count, both kfMeas builders
%   3. tests/test_shadow_mode       guide Step 5: shadow mode, then attack
%   4. tests/test_feedback_split    host setKF split (GAIN) + drain in COAST
%   5. tests/test_error_handlers    non-PD solves, numMeas clamp, NaN input, anchor stamp
%   6. run_scheduler_test           100 Hz + 2 Hz architecture (T1-T6)
%
% Each test is run in its own workspace (they may 'clear'); its console
% output is captured and searched for its own PASS marker.

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'tests'));

tests = { ...
    'test_monitors',         '=== ALL MONITOR TESTS PASS ===';
    'test_variable_numMeas', '=== test_variable_numMeas PASS ===';
    'test_shadow_mode',      '=== test_shadow_mode PASS ===';
    'test_feedback_split',   '=== test_feedback_split PASS ===';
    'test_error_handlers',   '=== test_error_handlers PASS ===';
    'run_scheduler_test',    '=== ALL SCHEDULER TESTS PASS ==='};

fprintf('==============================================\n');
fprintf('  Spoof-detection pre-integration tests\n');
fprintf('==============================================\n');
nPass = 0;
for i = 1:size(tests, 1)
    fprintf('%-24s ... ', tests{i, 1});
    [txt, ok] = runOne(tests{i, 1}, tests{i, 2});
    if ok
        fprintf('PASS\n'); nPass = nPass + 1;
    else
        fprintf('FAIL\n%s\n', txt);
    end
end
fprintf('==============================================\n');
fprintf('  %d / %d passed\n', nPass, size(tests, 1));
if nPass == size(tests, 1)
    fprintf('  All tests passed. Wire per docs/HOST_2HZ_WIRING.m.\n');
else
    fprintf('  Fix failures before integrating.\n');
end
fprintf('==============================================\n');
close all
end

function [txt, ok] = runOne(name, marker)
txt = captureRun(name);
ok  = ~isempty(strfind(txt, marker)); %#ok<STREMP>
end

function txt = captureRun(name)
% own workspace: a test script's 'clear' wipes only this function's variables
try
    txt = evalc(name);
catch e
    txt = sprintf('ERROR: %s', e.message);
end
end
