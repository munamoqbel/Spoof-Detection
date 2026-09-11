%******************************************************************************************
% DESCRIPTION:
% MATLAB-ONLY helper: screen and build SPF_gate with MATLAB Coder.
% Run from the repo root (or with the repo and the simulation's
% CST_gnssHybrid on the path). Produces a C static library plus the HTML
% code-generation report in <outDir>.
%
%   codegenSpoofMonitor              -> ./codegen_out
%   codegenSpoofMonitor('C:\tmp\spf')
%
% The entry-point argument types are taken from the STRUCT_SPF zero
% constructors, so the compiled layout is exactly the one the host builds
% with STRUCT_SPF.kfMeasFromUpdate / setPropTel:
%   kfMeas   fixed MAX_MEAS-padded arrays, numMeas uint8, 60-state vectors
%   propTel  accumPhi / accumQ [60 x 60]
%   navActive, resetRequest logical scalars
%
% ASSUMPTIONS AND LIMITATIONS:
% - Requires MATLAB Coder. Not for Octave.
% - STRUCT_SPF.kfMeasFromUpdate is compiled inside the host's 2 Hz
%   function; if the host passes exact-size (numMeas-row) arrays, declare
%   them there as bounded variable-size, e.g.
%   coder.typeof(zeros(51, 1), [51 1], [1 0])  (51 = CST_gnssHybrid.MAX_MEASURES).
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
function codegenSpoofMonitor(outDir)

if (nargin < 1)
    outDir = fullfile(pwd, 'codegen_out');
end
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    error('codegenSpoofMonitor:octave', 'MATLAB Coder is MATLAB-only.');
end

% 1. static readiness screen of the whole runtime call tree
coder.screener('SPF_gate');

% 2. entry-point argument types (fixed layouts)
kfMeasType  = coder.typeof(STRUCT_SPF.zeroKfMeas);
propTelType = coder.typeof(STRUCT_SPF.zeroPropTel);
navType     = coder.typeof(false);
resetType   = coder.typeof(false);

% 3. build a C library with the report
cfg = coder.config('lib');
cfg.TargetLang     = 'C';
cfg.GenerateReport = true;
cfg.LaunchReport   = false;

codegen('-config', cfg, 'SPF_gate', ...
    '-args', {kfMeasType, propTelType, navType, resetType}, '-d', outDir);

fprintf('SPF_gate built with MATLAB Coder into %s\n', outDir);

end
%------------------------------------------------------------------------
