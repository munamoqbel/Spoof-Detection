function octaveHeadless()
%OCTAVEHEADLESS  OCTAVE ONLY. Prepare a headless Octave session to run the
%   harness scripts: loads the statistics package, then materialises in a
%   TEMPORARY folder (never inside the repo, so nothing can shadow MATLAB):
%     * a constant-property stand-in for the CST_spfMode enumeration
%       (Octave has no enumeration classes),
%     * no-op plotting stubs (figure, plot, ...) so the scripts run headless,
%     * the scripts' local 'pf' helper (Octave cannot see script-local
%       functions).
%   The temp folder becomes the cwd (cwd has path precedence), then the repo
%   and its tests folder are added to the path.
%
%   Never call this from MATLAB.

if exist('OCTAVE_VERSION', 'builtin') == 0
    error('octaveHeadless:matlab', 'octaveHeadless is for GNU Octave only.');
end
pkg load statistics
warning('off', 'all');

here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', '..');

tmp = tempname(); mkdir(tmp);

% enumeration stand-in
copyfile(fullfile(here, 'CST_spfMode_octave.txt'), fullfile(tmp, 'CST_spfMode.m'));

% headless plotting stubs (files, so they survive 'clear functions')
stubs = {'figure','plot','subplot','semilogy','fill','yline','xline','ylabel', ...
         'xlabel','legend','grid','xlim','ylim','title','stairs','scatter', ...
         'yticks','yticklabels','hold','close'};
for k = 1:numel(stubs)
    fid = fopen(fullfile(tmp, [stubs{k} '.m']), 'w');
    fprintf(fid, 'function varargout = %s(varargin)\nvarargout = cell(1, nargout);\nend\n', stubs{k});
    fclose(fid);
end

% script-local helper used by the harness scripts
fid = fopen(fullfile(tmp, 'pf.m'), 'w');
fprintf(fid, 'function s = pf(ok)\nif ok, s = ''PASS''; else, s = ''FAIL''; end\nend\n');
fclose(fid);

cd(tmp);
addpath(repo, fullfile(repo, 'tests'));
end
