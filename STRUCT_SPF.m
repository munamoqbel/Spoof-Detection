%******************************************************************************************
% DESCRIPTION:
% Struct constructors for the spoofing-detection gate (definitions only).
% Every struct used by SPF_gate / SPF_protectedNav / SPF_monitorPool is built
% here so Coder sees one fixed layout per type. All vector/matrix fields are
% sized from CST_gnssHybrid.NO_STATES, CST_spfParam.MAX_MEAS and
% CST_spfParam.WINDOW_LENGTH.
%
% STATE CONVENTION (host error-state EKF, closed loop):
%   The gate never uses an absolute state. It accumulates the host KF's
%   UPDATE INCREMENTS  dx = postState - priorState  (= K*y), propagated with
%   the host's own interval matrices (propTel.accumPhi / accumQ). This is
%   invariant to how the host splits its estimate between the mechanization
%   feedback and the residual KF.states.
%
% ASSUMPTIONS AND LIMITATIONS:
%
% REQUIREMENT TRACEABILITY:
%
%******************************************************************************************
%#codegen
classdef STRUCT_SPF
    methods(Static)

        %% ---------------- inputs from the host ----------------
        function [kfMeas] = setKfMeas(innovation, innovationCov, obsMatrix, ...
                measNoiseCov, numMeas, priorState, postState, postCov)

            % Exception handler: never index past the fixed MAX_MEAS layout
            numMeasClamped = (double(numMeas) > double(CST_spfParam.MAX_MEAS));
            numMeasUsed    = uint8(min(double(numMeas), double(CST_spfParam.MAX_MEAS)));

            % Define structure (rows/cols 1:numMeas are valid)
            kfMeas = struct( ...
                'innovation',     innovation, ...      % y = z - h(xPrior)      [MAX_MEAS x 1]
                'innovationCov',  innovationCov, ...   % S = H P H' + R         [MAX_MEAS x MAX_MEAS]
                'obsMatrix',      obsMatrix, ...       % H(x)                   [MAX_MEAS x n]
                'measNoiseCov',   measNoiseCov, ...    % R                      [MAX_MEAS x MAX_MEAS]
                'numMeas',        numMeasUsed, ...     % valid rows this epoch (0 = no GNSS), <= MAX_MEAS
                'numMeasClamped', logical(numMeasClamped), ... % host passed more rows than MAX_MEAS
                'priorState',     priorState, ...      % x_bar into the update  [n x 1]
                'postState',      postState, ...       % x+ out of the update   [n x 1]
                'postCov',        postCov);            % P+                     [n x n]
        end

        function [kfMeas] = kfMeasFromUpdate(innovation, obsMatrix, measNoiseCov, ...
                numMeas, priorState, priorCov, postState, postCov)
            % Host convenience: kfUpdate exposes only y, H, R, numMeas and
            % x+, P+; the innovation covariance is formed HERE from the prior
            % covariance the host had going into the update (KF.covariance
            % as extrapolated at 100 Hz):  S = H P_bar H' + R.
            % Arrays may be exact-size (numMeas rows) or padded to MAX_MEAS.

            mMax = double(CST_spfParam.MAX_MEAS);
            n    = CST_gnssHybrid.NO_STATES;
            m    = min(double(numMeas), mMax);    % exception handler: rows beyond MAX_MEAS are dropped

            innovationCov = zeros(mMax, mMax);
            if (m > 0)
                H = obsMatrix(1:m, :);
                S = H * priorCov * H' + measNoiseCov(1:m, 1:m);
                innovationCov(1:m, 1:m) = (S + S') / 2;
            end % ELSE: no GNSS this epoch

            % pad the host arrays to the fixed layout
            innovationP = zeros(mMax, 1);  innovationP(1:m)     = innovation(1:m);
            obsMatrixP  = zeros(mMax, n);  obsMatrixP(1:m, :)   = obsMatrix(1:m, :);
            measNoiseP  = zeros(mMax);     measNoiseP(1:m, 1:m) = measNoiseCov(1:m, 1:m);

            kfMeas = STRUCT_SPF.setKfMeas(innovationP, innovationCov, obsMatrixP, ...
                measNoiseP, numMeas, priorState, postState, postCov);
        end

        function [kfMeas] = zeroKfMeas
            m = double(CST_spfParam.MAX_MEAS);
            n = CST_gnssHybrid.NO_STATES;
            kfMeas = STRUCT_SPF.setKfMeas(zeros(m, 1), zeros(m, m), zeros(m, n), ...
                zeros(m, m), 0, zeros(n, 1), zeros(n, 1), eye(n));
        end

        function [propTel] = setPropTel(accumPhi, accumQ)

            % Define structure: interval matrices accumulated on the 100 Hz
            % side since the previous GNSS epoch (spfAccumProp)
            propTel = struct( ...
                'accumPhi', accumPhi, ...
                'accumQ',   accumQ);
        end

        function [propTel] = zeroPropTel
            n = CST_gnssHybrid.NO_STATES;
            propTel = STRUCT_SPF.setPropTel(eye(n), zeros(n, n));
        end

        %% ---------------- monitor pool ----------------
        function [pool] = setMonitorPool(active, hadAlarm, windowAge, ...
                separation, coastCovariance, innovationBuffer, ...
                innovationCovBuffer, obsMatrixBuffer, numMeasBuffer)

            % Define structure
            pool = struct( ...
                'active',              logical(active), ...
                'hadAlarm',            logical(hadAlarm), ...
                'windowAge',           windowAge, ...
                'separation',          separation, ...        % d_w = KF - coast (increments) [n x N]
                'coastCovariance',     coastCovariance, ...   % P_C per window               [n x n x N]
                'innovationBuffer',    innovationBuffer, ...
                'innovationCovBuffer', innovationCovBuffer,...
                'obsMatrixBuffer',     obsMatrixBuffer, ...
                'numMeasBuffer',       numMeasBuffer);
        end

        function [pool] = zeroMonitorPool

            % Define variables (buffers sized to the measurement upper bound;
            % numMeasBuffer holds the valid row count per buffered epoch)
            windowLength = CST_spfParam.WINDOW_LENGTH;
            numStates = CST_gnssHybrid.NO_STATES;
            num_meas  = CST_spfParam.MAX_MEAS;

            active   = false(windowLength, 1);
            hadAlarm = false(windowLength, 1);
            windowAge = zeros(windowLength, 1);

            separation      = zeros(numStates, windowLength);
            coastCovariance = zeros(numStates, numStates, windowLength);

            innovationBuffer = zeros(num_meas, windowLength, windowLength);
            innovationCovBuffer = zeros(num_meas, num_meas, windowLength, ...
                windowLength);
            obsMatrixBuffer = zeros(num_meas, numStates, windowLength, ...
                windowLength);
            numMeasBuffer = zeros(windowLength, windowLength, 'uint8');

            % Set the output
            pool = STRUCT_SPF.setMonitorPool(active, hadAlarm, windowAge, ...
                separation, coastCovariance, innovationBuffer, ...
                innovationCovBuffer, obsMatrixBuffer, numMeasBuffer);
        end

        function [report] = setMonitorReport(alarmPerAxis, anyAlarm, ...
                maxProtectionLevel, ssAlarm, cpiAlarm, cleanCloseFound, ...
                cleanCloseSeparation, cleanCloseCovar, solveFault)

            % Define structure
            report = struct( ...
                'alarmPerAxis',         logical(alarmPerAxis), ...
                'anyAlarm',             logical(anyAlarm), ...
                'maxProtectionLevel',   maxProtectionLevel, ...
                'ssAlarm',              logical(ssAlarm), ...
                'cpiAlarm',             logical(cpiAlarm), ...
                'cleanCloseFound',      logical(cleanCloseFound), ...
                'cleanCloseSeparation', cleanCloseSeparation, ...
                'cleanCloseCovar',      cleanCloseCovar, ...
                'solveFault',           logical(solveFault));   % a CPI epoch had a non-PD S
        end

        function [report] = zeroMonitorReport

            % Init values
            axisAlarm = zeros(1, 3, 'logical');
            alarm = false;
            scalar = 0.0;
            zeroState = zeros(CST_gnssHybrid.NO_STATES, 1);
            zeroCovar = zeros(CST_gnssHybrid.NO_STATES, CST_gnssHybrid.NO_STATES);

            report = STRUCT_SPF.setMonitorReport(axisAlarm, alarm, ...
                scalar, alarm, alarm, alarm, zeroState, zeroCovar, alarm);
        end

        function [poolOut] = closeAllWindows(poolIn)

            % Set the output
            poolOut = poolIn;
            poolOut.active(:)    = false;
            poolOut.hadAlarm(:)  = false;
            poolOut.windowAge(:) = 0;

        end

        function [poolOut] = openWindow(poolIn, kfMeas)
            % Open a window on this epoch's post-update solution: the window's
            % coast is that solution propagated INS-only, so its separation
            % starts at zero and its covariance at P+ (paper E24).

            % Set the output
            poolOut = poolIn;

            % ---- find the first free slot ----
            slotFound = false;
            freeSlot = uint8(1);
            windowLength = CST_spfParam.WINDOW_LENGTH;

            for wIndx = 1:windowLength
                if (~slotFound) && (~poolIn.active(wIndx))
                    freeSlot = wIndx;
                    slotFound = true;
                end
            end

            % ---- initialise the window ----
            if (slotFound)
                poolOut.active(freeSlot)    = true;
                poolOut.hadAlarm(freeSlot)  = false;
                poolOut.windowAge(freeSlot) = 1;

                poolOut.separation(:, freeSlot)         = zeros(CST_gnssHybrid.NO_STATES, 1);
                poolOut.coastCovariance(:, :, freeSlot) = kfMeas.postCov;

                m = kfMeas.numMeas;
                poolOut.innovationBuffer(1:m, 1, freeSlot)          = kfMeas.innovation(1:m);
                poolOut.innovationCovBuffer(1:m, 1:m, 1, freeSlot)  = kfMeas.innovationCov(1:m, 1:m);
                poolOut.obsMatrixBuffer(1:m, :, 1, freeSlot)        = kfMeas.obsMatrix(1:m, :);
                poolOut.numMeasBuffer(1, freeSlot)                  = uint8(m);
            end

        end

        %% ---------------- SS result ----------------
        function [SPF_ssMonitor] = setSSmonitor(alarmPerAxis, anyAlarm, ...
                separation, sigmaSeparation, protectionLevel, ...
                maxProtectionLevel)

            % Define structure
            SPF_ssMonitor = struct( ...
                'alarmPerAxis',       logical(alarmPerAxis), ...
                'anyAlarm',           logical(anyAlarm), ...
                'separation',         separation, ...
                'sigmaSeparation',    sigmaSeparation, ...
                'protectionLevel',    protectionLevel, ...
                'maxProtectionLevel', maxProtectionLevel);
        end

        function [SPF_ssMonitor] = zeroSSmonitor

            % Init values
            axisAlarm = zeros(1, 3, 'logical');
            alarm = false;
            zeroVector = zeros(1, 3);
            scalar = 0.0;

            SPF_ssMonitor = STRUCT_SPF.setSSmonitor(axisAlarm, alarm, ...
                zeroVector, zeroVector, zeroVector, scalar);
        end

        %% ---------------- FSM state ----------------
        function [anchor] = setAnchor(valid, separation, covariance, epoch)
            % Certified-clean coast, kept live: separation = (host solution
            % now) - (anchor coast now), accumulated increments since the
            % anchor window opened; covariance = that coast's P_C.

            % Define structure
            anchor = struct( ...
                'valid',      logical(valid), ...
                'separation', separation, ...
                'covariance', covariance, ...
                'epoch',      uint32(epoch));
        end

        function [anchor] = zeroAnchor

            % Init values
            valid = false;
            separation = zeros(CST_gnssHybrid.NO_STATES, 1);
            covariance = zeros(CST_gnssHybrid.NO_STATES, CST_gnssHybrid.NO_STATES);
            epoch = uint32(0);

            anchor = STRUCT_SPF.setAnchor(valid, separation, covariance, epoch);

        end

        function [sys] = setSys(mode, coastCov, probSep, pool, anchor, ...
                dwellCount, probationCount, coastCount)

            % Define structure
            sys = struct( ...
                'mode',           mode, ...           % CST_spfMode
                'coastCov',       coastCov, ...       % covariance of the protected solution
                'probSep',        probSep, ...        % (active filter) - (coast), increments [n x 1]
                'pool',           pool, ...
                'anchor',         anchor, ...
                'dwellCount',     dwellCount, ...
                'probationCount', probationCount, ...
                'coastCount',     coastCount);
        end

        function [sys] = zeroSys

            % Init values
            mode = CST_spfMode.NOMINAL;
            n = CST_gnssHybrid.NO_STATES;
            coastCov = eye(n);
            probSep = zeros(n, 1);
            pool = STRUCT_SPF.zeroMonitorPool;
            anchor = STRUCT_SPF.zeroAnchor;
            dwellCount = 0.0;
            probationCount = 0.0;
            coastCount = 0.0;

            sys = STRUCT_SPF.setSys(mode, coastCov, probSep, pool, ...
                anchor, dwellCount, probationCount, coastCount);
        end

        function [altXCheck] = setAltXCheck(suspect, altDiff)

            % Define structure
            altXCheck = struct( ...
                'suspect', suspect, ...
                'altDiff', altDiff);

        end

        function [altXCheck] = zeroAltXCheck

            % Init values
            suspect = false;
            altDiff = 0.0;

            altXCheck = STRUCT_SPF.setAltXCheck(suspect, altDiff);
        end

        %% ---------------- telemetry / commands ----------------
        function [spoofTel] = setTel(info, kfCommand, nav)

            % Define structure
            spoofTel = struct( ...
                'info',      info, ...
                'kfCommand', kfCommand, ...
                'nav',       nav);
        end

        function [spoofTel] = zeroTel

            % Init values
            info = STRUCT_SPF.zeroInfo;
            kfCommand = STRUCT_SPF.zeroCommand;
            nav = STRUCT_SPF.zeroNav;

            % Define structure
            spoofTel = STRUCT_SPF.setTel(info, kfCommand, nav);
        end

        function [info] = setInfo(mode, ssAlarm, cpiAlarm, ...
                alarmPerAxis, maxProtectionLevel, qReval, revalComputed,...
                dwellCount, eventLatched, eventAnchorEpoch, eventProbationStarted,...
                eventProbationVetoed, eventHandback, anchorMissing, coastEpochs, ...
                inputFault, numMeasClamped, solveFault)

            % Define structure
            info = struct( ...
                'mode',                  mode, ...
                'ssAlarm',               logical(ssAlarm), ...
                'cpiAlarm',              logical(cpiAlarm), ...
                'alarmPerAxis',          logical(alarmPerAxis), ...
                'maxProtectionLevel',    maxProtectionLevel, ...
                'qReval',                qReval, ...
                'revalComputed',         logical(revalComputed), ...
                'dwellCount',            dwellCount, ...
                'eventLatched',          logical(eventLatched), ...
                'eventAnchorEpoch',      uint32(eventAnchorEpoch), ...
                'eventProbationStarted', logical(eventProbationStarted), ...
                'eventProbationVetoed',  logical(eventProbationVetoed), ...
                'eventHandback',         logical(eventHandback), ...
                'anchorMissing',         logical(anchorMissing), ...
                'coastEpochs',           coastEpochs, ...
                'inputFault',            logical(inputFault), ...     % non-finite input: gate reset this epoch
                'numMeasClamped',        logical(numMeasClamped), ... % host passed > MAX_MEAS rows (clamped)
                'solveFault',            logical(solveFault));        % S or residual covariance not PD this epoch
        end

        function [info] = zeroInfo

            % Init values
            mode                  = CST_spfMode.NOMINAL;
            ssAlarm               = false;
            cpiAlarm              = false;
            alarmPerAxis          = false(1, 3);
            maxProtectionLevel    = 0.0;
            qReval                = 0.0;
            revalComputed         = false;
            dwellCount            = 0;
            eventLatched          = false;
            eventAnchorEpoch      = uint32(0);
            eventProbationStarted = false;
            eventProbationVetoed  = false;
            eventHandback         = false;
            anchorMissing         = false;
            coastEpochs           = 0;
            inputFault            = false;
            numMeasClamped        = false;
            solveFault            = false;

            % Define structure
            info = STRUCT_SPF.setInfo(mode, ssAlarm, cpiAlarm, ...
                alarmPerAxis, maxProtectionLevel, qReval, revalComputed,...
                dwellCount, eventLatched, eventAnchorEpoch, eventProbationStarted,...
                eventProbationVetoed, eventHandback, anchorMissing, coastEpochs, ...
                inputFault, numMeasClamped, solveFault);
        end

        function [command] = setCommand(startTrial)
            % startTrial true = probation opens: the host must start its TRIAL
            % filter as a copy of the operational (coasting) KF, states and
            % covariance as the 100 Hz side extrapolated them.

            % Define structure
            command = struct( ...
                'startTrial', logical(startTrial));
        end

        function [command] = zeroCommand

            % Init values
            startTrial = false;

            % Define structure
            command = STRUCT_SPF.setCommand(startTrial);
        end

        function [nav] = setNav(applyCorrection, state, correction, covar, sigmaPosition)
            % applyCorrection true on LATCH and COMMIT epochs: (state, covar)
            % is a complete update result (x+, P+) for the OPERATIONAL KF,
            % expressed in the host's KF-state space. The host passes it to
            % its normal setKF exactly like any kfUpdate result.
            % correction = state - (this epoch's x+ of the active filter) is
            % the same information as an increment (diagnostics / harness).

            % Define structure
            nav = struct( ...
                'applyCorrection', logical(applyCorrection), ...
                'state',           state, ...         % [n x 1] x+ to hand to setKF
                'correction',      correction, ...    % [n x 1] increment form of the same
                'covar',           covar, ...         % [n x n] P+ / protected-solution covariance
                'sigmaPosition',   sigmaPosition);    % [3 x 1] 1-sigma of position
        end

        function [nav] = zeroNav

            % Init values
            n = CST_gnssHybrid.NO_STATES;
            nav = STRUCT_SPF.setNav(false, zeros(n, 1), zeros(n, 1), eye(n), zeros(3, 1));
        end

    end
end

%------------------------------------------------------------------------
