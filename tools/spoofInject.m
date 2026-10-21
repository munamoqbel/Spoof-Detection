%****************************************************************************
% DESCRIPTION:
% Inject a spoof into one constellation's slot arrays of the host, for
% one epoch: the ranges and range rates a receiver standing at the true
% position plus a position offset, moving with a velocity offset and
% carrying a clock offset and drift would measure. Call it in the host's
% 2 Hz function after the corrected ranges are formed and before the
% rows are built, once per constellation to spoof. The injection is exact
% in the geometry (no linearisation), so a step, a ramp or a clock push
% of the matrix is reproduced as a spoofer simulating that trajectory
% would. Unused slots (range 0 or not finite) are left as they are.
%
% INPUTS:
%   - geomRange      [S x 1] corrected pseudorange per slot (m), 0 unused
%   - geomRangeRate  [S x 1] range rate per slot (m/s), d(range)/dt
%   - satPos         [S x 3] satellite ECEF position per slot (m), the one
%                    the host pairs with the slot in h(x)
%   - satVel         [S x 3] satellite ECEF velocity per slot (m/s)
%   - posTrueEcef    [1 x 3] true receiver ECEF position (m)
%   - velTrueEcef    [1 x 3] true receiver ECEF velocity (m/s)
%   - latLon         [1 x 2] latitude, longitude (rad) of the receiver,
%                    the NED frame the offsets are given in
%   - inj            struct from spoofScenario: dPosNed [1 x 3] (m),
%                    dVelNed [1 x 3] (m/s), dClk (m), dClkRate (m/s)
%
% OUTPUTS:
%   - geomRange, geomRangeRate  the slot arrays with the spoof in them
%
% ASSUMPTIONS AND LIMITATIONS:
% The range rate convention is the host's: u' (vSat - vRec) plus the
% clock drift, u the unit line of sight from the receiver to the
% satellite. A slot whose satellite position is closer than 1000 km to
% the receiver is not a satellite and is skipped. The injection does not
% touch satPos, satVel, the validity flags or the row builder; a spoof
% that changes the satellite set is outside this function.
%
% REQUIREMENT TRACEABILITY:
% Matrix of the gate tests (steps, ramps, clock pushes on the clean
% twin): the injected offset is the quantity the SS test (Eq. 49) and the
% CPI (Eq. 17 to 35) must find, and spoofProfile reads it back from the
% ranges (its reference fit) for a check of the injection itself.
%****************************************************************************
function [geomRange, geomRangeRate] = spoofInject(geomRange, geomRangeRate, ...
    satPos, satVel, posTrueEcef, velTrueEcef, latLon, inj)

MIN_SAT_DISTANCE = 1.0e6;                                 % m: below it the slot holds no satellite

sinLat = sin(latLon(1));  cosLat = cos(latLon(1));
sinLon = sin(latLon(2));  cosLon = cos(latLon(2));
nedToEcef = [-sinLat * cosLon, -sinLon, -cosLat * cosLon; ...   % columns: N, E, D unit vectors in ECEF
             -sinLat * sinLon,  cosLon, -cosLat * sinLon; ...
              cosLat,           0.0,    -sinLat];
dPosEcef = (nedToEcef * inj.dPosNed(:))';
dVelEcef = (nedToEcef * inj.dVelNed(:))';
posSpoof = posTrueEcef + dPosEcef;

for slot = 1:numel(geomRange)
    losTrue  = satPos(slot, :) - posTrueEcef;
    losSpoof = satPos(slot, :) - posSpoof;
    rhoTrue  = sqrt(sum(losTrue.^2));
    rhoSpoof = sqrt(sum(losSpoof.^2));
    if isfinite(geomRange(slot)) && (geomRange(slot) ~= 0.0) && (rhoTrue > MIN_SAT_DISTANCE)
        uTrue  = losTrue  / rhoTrue;
        uSpoof = losSpoof / rhoSpoof;
        geomRange(slot)     = geomRange(slot) + (rhoSpoof - rhoTrue) + inj.dClk;
        geomRangeRate(slot) = geomRangeRate(slot) ...
            + (uSpoof - uTrue) * (satVel(slot, :) - velTrueEcef)' - uSpoof * dVelEcef' + inj.dClkRate;
    end % ELSE: an unused slot, or no satellite in it: left as it is
end

end
