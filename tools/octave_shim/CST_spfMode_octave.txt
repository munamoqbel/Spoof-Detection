% OCTAVE-ONLY SHIM (not part of the repo): Octave has no enumeration
% classes, so expose the same three constants as uint8 Constant properties.
classdef CST_spfMode
    properties (Constant = true)
        NOMINAL   = uint8(1);
        COAST     = uint8(2);
        PROBATION = uint8(3);
    end
end
