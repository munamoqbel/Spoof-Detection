# Octave regression shim (NOT for MATLAB, NOT for the simulation)

Lets the harness and tests run headless under GNU Octave (tested 8.4 +
statistics package). Nothing in this folder may be added to a MATLAB path:
`CST_spfMode_octave.txt` is a constant-property stand-in for the
enumeration class, and `octaveHeadless.m` writes it plus no-op plotting
stubs (`figure`, `plot`, ...) into a temporary folder outside the repo for
the session only. Any copy of those stubs on a MATLAB path silently
suppresses all figures.

From the repo root:

    octave --no-gui --quiet --eval "addpath tools/octave_shim; octaveHeadless; runAllTests"
    octave --no-gui --quiet --eval "addpath tools/octave_shim; octaveHeadless; run_recovery"

Octave's `rng(42)`/`randn` stream differs from MATLAB's, so the random
geometry differs: expect the same design constants and the same qualitative
events, not bit-exact epochs. Bit-exact regression against `reference/*.txt`
needs MATLAB.
