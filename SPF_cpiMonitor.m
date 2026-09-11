%******************************************************************************************
% DESCRIPTION:
% Cumulative Position-domain Innovation monitor, one complete window.
%
% INPUTS:
%   - innovationBuffer          [MAX_MEAS x N] gamma history of the window
%   - innovationCovBuffer       [MAX_MEAS x MAX_MEAS x N] S history
%   - obsMatrixBuffer           [MAX_MEAS x n x N] H history
%   - numMeasBuffer             [N x 1] valid measurement count per epoch
%                               (rows/cols beyond it are padding)
%   - axisIdx                   monitored position state index
%   Window length N and threshold T_N are CST_spfParam.WINDOW_LENGTH and
%   CST_spfParam.CPI_THRESHOLD; the buffers are N deep.
%
% OUTPUTS:
%   - cpiAlarm                  logical
%   - qStatistic                double Eq. 33 statistic
%   - xiHistory                 [N x 1] normalised projections (diagnostics)
%   - solveFault                logical, true if S of at least one epoch was
%                               not positive definite (that epoch counted
%                               xi = 0; host contract violation, e.g. a row
%                               with zero noise or a duplicated row)
%
% ASSUMPTIONS AND LIMITATIONS:
% An epoch with numMeas = 0 (no GNSS) contributes xi = 0 to the window.
% S is inverted through its Cholesky factor: S = Rc'*Rc, so with
% w_g = Rc'\gamma and w_f = Rc'\f,  f'S^-1 gamma = w_f'w_g  and
% f'S^-1 f = w_f'w_f  (no explicit inverse, no division by a pivot that can
% be zero). chol with two outputs never errors: p > 0 flags a non-PD S, and
% a relative pivot floor (CST_spfParam.PIVOT_REL_TOL) rejects an S that is
% singular up to rounding.
%
% REQUIREMENT TRACEABILITY:
% - PAPER MAPPING
%    Equations (Kujur 2024)
%    Eq. 17: gamma_u = f'*(S\gamma), f = H(:,axis_idx) = (H * u)
%    Eq. 20: sigma2 = f'*(S\f)
%    Eq. 29: xi = gamma_u / sqrt(sigma2) ~ N(0,1) under H0
%    Eq. 33: q = sum(xi.^2) ~ Gamma(N/2,2) under H0
%    Eq. 35: alarm iff q > cpi_threshold (T_N, precomputed offline)
%
%******************************************************************************************
%#codegen
function [cpiAlarm, qStatistic, xiHistory, solveFault] = SPF_cpiMonitor...
    (innovationBuffer, innovationCovBuffer, obsMatrixBuffer, numMeasBuffer, axisIdx)

% Define variables
windowLength = CST_spfParam.WINDOW_LENGTH;
cpiThreshold = CST_spfParam.CPI_THRESHOLD;
maxMeas      = uint8(size(innovationBuffer, 1));
xiHistory    = zeros(windowLength, 1);
qStatistic   = 0.0;
cpiAlarm     = false;
solveFault   = false;

for idx = 1:windowLength
    numMeas = min(numMeasBuffer(idx), maxMeas);       % defensive: never past the buffer
    xiNormalised = 0.0;

    if (numMeas > 0)
        % ------------------------------------------------------------------
        % Fixed-size formulation (MATLAB Coder without variable sizing):
        % work on the full MAX_MEAS layout. Rows/cols beyond numMeas are
        % zeroed and the padding diagonal of S is set to a positive value,
        % so S_full = blkdiag(S, padValue*I). Its Cholesky factor is
        % blkdiag(Rc, sqrt(padValue)*I), the whitened vectors are zero in
        % the padding, and every sum below equals the numMeas-row result.
        % ------------------------------------------------------------------
        innovation    = innovationBuffer(:, idx);                % gamma (padded). Eq.3
        innovationCov = innovationCovBuffer(:, :, idx);          % S (padded). Defined under Eq. 4
        projection    = obsMatrixBuffer(:, axisIdx, idx);        % f = H(:,axis) (padded). From Eq. 17

        padValue = max(max(diag(innovationCov)), 1.0);
        for rowIdx = (double(numMeas) + 1):double(maxMeas)
            innovation(rowIdx)    = 0.0;
            projection(rowIdx)    = 0.0;
            innovationCov(rowIdx, :) = 0.0;
            innovationCov(:, rowIdx) = 0.0;
            innovationCov(rowIdx, rowIdx) = padValue;
        end

        % ------------------------------------------------------------------------
        % Exception handler: S must be positive definite. chol never errors;
        % p > 0 means singular / indefinite / non-finite -> this epoch is
        % dropped (xi = 0) and the fault is reported.
        % ------------------------------------------------------------------------
        [cholFactor, cholFail] = chol(innovationCov);
        pivotOk = (cholFail == 0);
        if (pivotOk)
            % an exactly singular S can still factorise with a rounding-level
            % pivot: require every pivot to be above PIVOT_REL_TOL of the
            % largest diagonal (condition number below ~1e12). The padding
            % pivots equal padValue >= max(diag S) and never trip it.
            minPivot = min(diag(cholFactor)) ^ 2;
            pivotOk  = (minPivot > CST_spfParam.PIVOT_REL_TOL * max(diag(innovationCov)));
        end % ELSE is trivial
        if (pivotOk)
            wInnovation = cholFactor' \ innovation;      % Rc'^-1 gamma  (zero in the padding)
            wProjection = cholFactor' \ projection;      % Rc'^-1 f      (zero in the padding)

            gammaProjection  = 0.0;  % Eq. 17
            sigma2Projection = 0.0;  % Eq. 20
            for idxMeas = 1:double(maxMeas)
                gammaProjection  = gammaProjection + wProjection(idxMeas) * wInnovation(idxMeas);   % = w_f' * w_g (Eq.17); loop used for deterministic rounding
                sigma2Projection = sigma2Projection + wProjection(idxMeas) * wProjection(idxMeas);  % = w_f' * w_f (Eq.20); loop used for deterministic rounding
            end

            % Exception handler: axis unobservable this epoch (f = 0) -> xi = 0
            if (sigma2Projection > 0.0) && isfinite(sigma2Projection) && isfinite(gammaProjection)
                xiNormalised = gammaProjection / sqrt(sigma2Projection);   % Eq. 29
            end % ELSE is trivial
        else
            solveFault = true;
        end
    end % ELSE: no GNSS this epoch -> xi = 0

    xiHistory(idx) = xiNormalised;
    qStatistic     = qStatistic + (xiNormalised * xiNormalised);   % Eq. 33
end

if (qStatistic > cpiThreshold)    % alarm rule (section 2 under Eq. 4)
    cpiAlarm = true;
end % ELSE is trivial

end

%------------------------------------------------------------------------
