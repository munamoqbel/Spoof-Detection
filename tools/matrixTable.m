function matrixTable(rows)
%MATRIXTABLE  The injection-matrix table from the rows matrixRow returns.
%   matrixTable(rows)   rows a struct array (or cell) of matrixRow results, one line per run
%   Columns: the run, the scenario, the latch delay from tOn (s) and the monitor and axis that fired with its ratio,
%   the largest SS and CPI ratio from tOn to the hand-back and on the clean epochs, the nav error at the latch
%   (N E, m), the smallest CPI against the coast while the spoof was on (a gate that logs it; the dwell holds above 1),
%   the largest horizontal error from tOn to the hand-back, the error at tOff (N E), the hand-back delay
%   from tOff (s), the coast length (s), the vetoes, the false alarms. A dash is a number the run does not have.
%   MATLAB desktop tool, not codegen. Runs under Octave.
if iscell(rows), rows = [rows{:}]; end
fprintf('%-14s %-14s %7s %-8s %5s %7s %7s %7s %7s %7s %13s %7s %13s %8s %7s %4s %3s\n', 'run', 'scenario', 'latch', 'monitor', 'ratio', ...
    'SS in', 'CPI in', 'SS cln', 'CPI cln', 'cstCPI', 'err@latch N E', 'maxHor', 'err@off N E', 'handback', 'coast', 'veto', 'FA');
fprintf('%-14s %-14s %7s %-8s %5s %7s %7s %7s %7s %7s %13s %7s %13s %8s %7s %4s %3s\n', '', '', '[s]', '', '', '', '', '', '', 'min', '[m]', '[m]', '[m]', '[s]', '[s]', '', '');
for i = 1:numel(rows)
    r = rows(i);
    fprintf('%-14s %-14s %7s %-8s %5s %7s %7s %7s %7s %7s %13s %7s %13s %8s %7s %4d %3d\n', cut(r.name, 14), ...
        cut(sprintf('%s %g', r.kind, r.size), 14), num(r.latchDelay, '%.1f'), cut([r.monitor ' ' r.axis], 8), num(r.ratioAtLatch, '%.2f'), ...
        num(r.maxSsIn, '%.2f'), num(r.maxCpiIn, '%.2f'), num(r.maxSsClean, '%.2f'), num(r.maxCpiClean, '%.2f'), num(r.coastCpiMin, '%.2f'), ...
        pair(r.errAtLatch), num(r.maxHorizIn, '%.1f'), pair(r.errAtOff), num(r.handbackDelay, '%.1f'), num(r.coastLen, '%.1f'), ...
        r.nVeto, r.falseAlarms);
end
end

function [s] = num(x, fmt)
if isnan(x), s = '-'; else, s = sprintf(fmt, x); end
end

function [s] = pair(v)
if any(isnan(v(1:2))), s = '-'; else, s = sprintf('%+6.1f %+6.1f', v(1), v(2)); end
end

function [s] = cut(s, n)
if numel(s) > n, s = s(1:n); end
end
