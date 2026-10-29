function test_gravityCheck()
%% test_gravityCheck.m
% tools/gravityCheck: the WGS84 normal gravity at the equator, at 45 degrees and at the pole against the tabulated
% values, the height term at 1000 m, the site taken from a run's truePos, and the error of a constant in mg.
% Run from anywhere:  octave --no-gui --quiet --eval test_gravityCheck
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'tools'));
ok = true;
chk = struct();
txt = evalc('r0 = gravityCheck([0 0]);');
chk.equator = abs(r0.gNormal - 9.7803253) < 2e-6;
txt = evalc('r45 = gravityCheck([pi/4 0], 9.80665);');
chk.lat45   = abs(r45.gNormal - 9.8061992) < 2e-6 && abs(r45.hostErrMg - (9.80665 - 9.8061992) / 9.80665e-3) < 0.01;
txt = evalc('r90 = gravityCheck([pi/2 0]);');
chk.pole    = abs(r90.gNormal - 9.8321849) < 2e-6;
txt = evalc('rh = gravityCheck([pi/4 1000]);');
chk.height  = abs((r45.gNormal - rh.gNormal) - (3.0877e-6 - 4.4e-9 * 0.5) * 1000 + 7.2e-13 * 1e6) < 1e-9;
run = struct('truePos', [linspace(0.60, 0.62, 100)', zeros(100, 1), linspace(100, 300, 100)']);
txt = evalc('rr = gravityCheck(run, 9.81);');
chk.site    = abs(rr.lat - 0.61) < 1e-9 && abs(rr.h - 200) < 1e-9 && rr.spreadMg > 0 && ...
              abs(rr.hostErrMg - (9.81 - rr.gNormal) / 9.80665e-3) < 1e-9 && ~isempty(strfind(txt, 'the host''s own value'));
chk.table   = numel(rr.candidateErrMg) == 5 && abs(rr.candidateErrMg(1) - (9.80665 - rr.gNormal) / 9.80665e-3) < 1e-9;
names = fieldnames(chk);
for i = 1:numel(names)
    if ~chk.(names{i}), fprintf('   FAILED sub-check: %s\n', names{i}); ok = false; end
end
fprintf('equator %.7f, 45 deg %.7f, pole %.7f m/s^2; 9.81 at the test site %+.2f mg\n', r0.gNormal, r45.gNormal, r90.gNormal, rr.hostErrMg);
if ok, fprintf('=== test_gravityCheck PASS ===\n'); else, fprintf('=== test_gravityCheck FAIL ===\n'); end
end
