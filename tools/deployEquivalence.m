function [same] = deployEquivalence()
%DEPLOYEQUIVALENCE  Check that deploy/ behaves exactly like the root copies.
%   Runs the harness scripted attack (recovery_nav_sim) and a shadow sweep
%   with varying numMeas, one alignment epoch and one NaN epoch, first with
%   the root files, then with deploy/ ahead of them on the path, and
%   compares every output field bitwise. Run from the repo root (MATLAB or
%   Octave with tools/octave_shim). Restores the path afterwards.

repo      = fileparts(fileparts(mfilename('fullpath')));
deployDir = fullfile(repo, 'deploy');

refRes = runScenario();

addpath(deployDir);
clear functions
depRes = runScenario();
rmpath(deployDir);
clear functions

same = isequaln(refRes, depRes);
if same
    fprintf('deploy/ and root copies: bitwise identical -> PASS\n');
else
    fprintf('deploy/ and root copies DIFFER -> FAIL\n');
    fprintf('  max |x_nav| diff %g\n', max(abs(refRes.x_nav(:) - depRes.x_nav(:))));
    for k = 1:numel(refRes.tels)
        if ~isequaln(refRes.tels{k}, depRes.tels{k})
            fprintf('  first telemetry difference at shadow epoch %d\n', k);
            break
        end
    end
end
end

function [res] = runScenario()
prmBoot = kujur_params();
n = prmBoot.n_states;
prm = kujur_params(double(CST_spfParam.WINDOW_LENGTH), CST_spfParam.CPI_THRESHOLD, ...
    CST_spfParam.K_FALSE_ALERT, CST_spfParam.K_MISSED_DETECTION);

% scripted attack with authority
[~, ~, H_all, ~, ~, Phi, Q, ~, z_all, V, P0] = generate_test_data(prmBoot, 42, 0.10, 0.5, [1 1 1]);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
out = recovery_nav_sim(z_all, H_all, V, propTel, P0, prm);
res.x_nav      = out.x_nav;
res.ev         = out.ev;
res.alarm_axis = out.alarm_axis;

% shadow sweep: varying numMeas, one alignment epoch, one NaN epoch
[~, ~, H_all, ~, ~, Phi, Q, scn, z_all, V, P0] = generate_test_data(prmBoot, 7, 0.05, 0.3, [1 0 1]);
propTel = STRUCT_SPF.setPropTel(Phi, Q);
N = scn.N_total;
m = size(z_all, 1);
kfX = zeros(n, 1);
kfP = P0;
tels = cell(N, 1);
for k = 1:N
    mk = m;
    if mod(k, 37) == 0
        mk = 2;
    end
    [kfX, kfP, y, ~, ~, xp, xpP] = kalman_update_step(kfX, kfP, z_all(:, k), H_all(:, :, k), propTel, V);
    kfMeas = STRUCT_SPF.kfMeasFromUpdate(y, H_all(:, :, k), V, mk, xp, xpP, kfX, kfP);
    if k == 150
        kfMeas.innovation(1) = NaN;
    end
    tels{k} = SPF_gate(kfMeas, propTel, k ~= 200, k == 1);
end
res.tels = tels;
end
