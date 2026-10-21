%****************************************************************************
% DESCRIPTION:
% The spoof offsets at time t for one scenario of the matrix: a position
% step, a position ramp, a clock step or a clock ramp, switched on at
% tOn and off at tOff. The result goes to spoofInject once per epoch and
% per constellation. Off the interval every offset is zero, so the same
% call runs the clean part of the run too.
%
% INPUTS:
%   - t      host time of the epoch (s)
%   - scn    struct: kind ('step' | 'ramp' | 'clock' | 'clockRamp'),
%            tOn, tOff (s), size (m for step and clock, m/s for ramp and
%            clockRamp), dirNed [1 x 3] unit vector of the position
%            offset (unused by the clock kinds)
%
% OUTPUTS:
%   - inj    struct: dPosNed [1 x 3] (m), dVelNed [1 x 3] (m/s), dClk (m),
%            dClkRate (m/s)
%
% ASSUMPTIONS AND LIMITATIONS:
% A ramp starts at zero at tOn and grows at size m/s along dirNed, with
% the matching velocity, so the Doppler stays consistent with the
% position as a spoofer dragging a receiver would keep it. One scenario
% per run: a step on top of a ramp is two calls with two structs, added
% field by field by the caller.
%
% REQUIREMENT TRACEABILITY:
% The matrix of the gate tests: steps of 5, 10 and 20 m, ramps of 0.3, 1
% and 2 m/s and clock pushes, each against the latch time, the largest
% error before the latch and the hand-back.
%****************************************************************************
function [inj] = spoofScenario(t, scn)

dPosNed  = [0.0 0.0 0.0];
dVelNed  = [0.0 0.0 0.0];
dClk     = 0.0;
dClkRate = 0.0;

if (t >= scn.tOn) && (t < scn.tOff)
    switch scn.kind
        case 'step'
            dPosNed = scn.size * scn.dirNed;
        case 'ramp'
            dPosNed = scn.size * (t - scn.tOn) * scn.dirNed;
            dVelNed = scn.size * scn.dirNed;
        case 'clock'
            dClk = scn.size;
        case 'clockRamp'
            dClk     = scn.size * (t - scn.tOn);
            dClkRate = scn.size;
        otherwise
            error('spoofScenario:kind', 'kind must be step, ramp, clock or clockRamp; got %s', scn.kind);
    end
end % ELSE: outside the interval every offset stays zero

inj = struct('dPosNed', dPosNed, 'dVelNed', dVelNed, 'dClk', dClk, 'dClkRate', dClkRate);

end
