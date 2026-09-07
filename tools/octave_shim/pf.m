function s = pf(ok)
% OCTAVE-ONLY shim: run_scheduler_test's local function, which Octave cannot see before use
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
