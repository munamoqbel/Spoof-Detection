%******************************************************************************************
% DESCRIPTION:
% Struct constructors for the spoofing-detection gate (one fixed layout per type, sized from
% CST_gnssHybrid.NO_STATES, CST_spfParam.MAX_MEAS and CST_spfParam.WINDOW_LENGTH).
% The gate works on update increments dx = postState - priorState (= K*y).
%******************************************************************************************
%#codegen
classdef STRUCT_SPF
    methods(Static)

        %% ---------------- inputs from the host ----------------
        function [kfMeas] = setKfMeas(innovation, innovationCov, obsMatrix, ...
                measNoiseCov, numMeas, priorState, postState, postCov)

            numMeasClamped = (double(numMeas) > double(CST_spfParam.MAX_MEAS));
            numMeasUsed    = uint8(min(double(numMeas), double(CST_spfParam.MAX_MEAS)));

            kfMeas = struct( ...
                'innovation',     innovation, ...      % y = z - h(xPrior)   [MAX_MEAS x 1]
                'innovationCov',  innovationCov, ...   % S = H P H' + R      [MAX_MEAS x MAX_MEAS]
                'obsMatrix',      obsMatrix, ...       % H                   [MAX_MEAS x n]
                'measNoiseCov',   measNoiseCov, ...    % R                   [MAX_MEAS x MAX_MEAS]
                'numMeas',        numMeasUsed, ...     % valid rows (0 = no GNSS)
                'numMeasClamped', logical(numMeasClamped), ...
                'priorState',     priorState, ...      % x_bar               [n x 1]
                'postState',      postState, ...       % x+                  [n x 1]
                'postCov',        postCov);            % P+                  [n x n]
        end

        function [kfMeas] = kfMeasFromUpdate(innovation, obsMatrix, measNoiseCov, ...
                numMeas, priorState, priorCov, postState, postCov)
            % pad the host arrays (exact-size or padded) to MAX_MEAS and form S = H P_bar H' + R

            mMax = double(CST_spfParam.MAX_MEAS);
            n    = CST_gnssHybrid.NO_STATES;
            m    = min(double(numMeas), mMax);

            innovationP = zeros(mMax, 1);
            obsMatrixP  = zeros(mMax, n);
            measNoiseP  = zeros(mMax, mMax);
            for rowIdx = 1:m                       % element loops: no variable-size slices for Coder
                innovationP(rowIdx) = innovation(rowIdx);
                for colIdx = 1:n
                    obsMatrixP(rowIdx, colIdx) = obsMatrix(rowIdx, colIdx);
                end
                for colIdx = 1:m
                    measNoiseP(rowIdx, colIdx) = measNoiseCov(rowIdx, colIdx);
                end
            end

            S = obsMatrixP * priorCov * obsMatrixP' + measNoiseP;
            innovationCov = (S + S') / 2;

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
            % interval matrices accumulated on the 100 Hz side since the previous epoch
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

            pool = struct( ...
                'active',              logical(active), ...
                'hadAlarm',            logical(hadAlarm), ...
                'windowAge',           windowAge, ...
                'separation',          separation, ...        % KF - coast per window [n x N]
                'coastCovariance',     coastCovariance, ...   % P_C per window        [n x n x N]
                'innovationBuffer',    innovationBuffer, ...
                'innovationCovBuffer', innovationCovBuffer,...
                'obsMatrixBuffer',     obsMatrixBuffer, ...
                'numMeasBuffer',       numMeasBuffer);
        end

        function [pool] = zeroMonitorPool

            windowLength = CST_spfParam.WINDOW_LENGTH;
            numStates = CST_gnssHybrid.NO_STATES;
            num_meas  = CST_spfParam.MAX_MEAS;

            active    = false(windowLength, 1);
            hadAlarm  = false(windowLength, 1);
            windowAge = zeros(windowLength, 1);

            separation      = zeros(numStates, windowLength);
            coastCovariance = zeros(numStates, numStates, windowLength);

            innovationBuffer    = zeros(num_meas, windowLength, windowLength);
            innovationCovBuffer = zeros(num_meas, num_meas, windowLength, windowLength);
            obsMatrixBuffer     = zeros(num_meas, numStates, windowLength, windowLength);
            numMeasBuffer       = zeros(windowLength, windowLength, 'uint8');

            pool = STRUCT_SPF.setMonitorPool(active, hadAlarm, windowAge, ...
                separation, coastCovariance, innovationBuffer, ...
                innovationCovBuffer, obsMatrixBuffer, numMeasBuffer);
        end

        function [report] = setMonitorReport(alarmPerAxis, anyAlarm, ...
                maxProtectionLevel, ssAlarm, cpiAlarm, cleanCloseFound, ...
                cleanCloseSeparation, cleanCloseCovar, solveFault, ssRatio, cpiRatio)

            report = struct( ...
                'alarmPerAxis',         logical(alarmPerAxis), ...
                'anyAlarm',             logical(anyAlarm), ...
                'maxProtectionLevel',   maxProtectionLevel, ...
                'ssAlarm',              logical(ssAlarm), ...
                'cpiAlarm',             logical(cpiAlarm), ...
                'cleanCloseFound',      logical(cleanCloseFound), ...
                'cleanCloseSeparation', cleanCloseSeparation, ...
                'cleanCloseCovar',      cleanCloseCovar, ...
                'solveFault',           logical(solveFault), ...
                'ssRatio',              ssRatio, ...     % [1x3] |d|/(k_FA sigma_SS), alarm at > 1
                'cpiRatio',             cpiRatio);       % [1x3] q/T_N, alarm at > 1
        end

        function [report] = zeroMonitorReport

            axisAlarm = zeros(1, 3, 'logical');
            alarm = false;
            scalar = 0.0;
            zeroState = zeros(CST_gnssHybrid.NO_STATES, 1);
            zeroCovar = zeros(CST_gnssHybrid.NO_STATES, CST_gnssHybrid.NO_STATES);

            report = STRUCT_SPF.setMonitorReport(axisAlarm, alarm, ...
                scalar, alarm, alarm, alarm, zeroState, zeroCovar, alarm, zeros(1, 3), zeros(1, 3));
        end

        function [poolOut] = closeAllWindows(poolIn)

            poolOut = poolIn;
            poolOut.active(:)    = false;
            poolOut.hadAlarm(:)  = false;
            poolOut.windowAge(:) = 0;

        end

        function [poolOut] = openWindow(poolIn, kfMeas)
            % open a window on this epoch's post-update solution (separation 0, covariance P+)

            poolOut = poolIn;

            slotFound = false;
            freeSlot = uint8(1);
            windowLength = CST_spfParam.WINDOW_LENGTH;

            for wIndx = 1:windowLength
                if (~slotFound) && (~poolIn.active(wIndx))
                    freeSlot = wIndx;
                    slotFound = true;
                end
            end

            if (slotFound)
                poolOut.active(freeSlot)    = true;
                poolOut.hadAlarm(freeSlot)  = false;
                poolOut.windowAge(freeSlot) = 1;

                poolOut.separation(:, freeSlot)         = zeros(CST_gnssHybrid.NO_STATES, 1);
                poolOut.coastCovariance(:, :, freeSlot) = kfMeas.postCov;

                poolOut.innovationBuffer(:, 1, freeSlot)       = kfMeas.innovation;
                poolOut.innovationCovBuffer(:, :, 1, freeSlot) = kfMeas.innovationCov;
                poolOut.obsMatrixBuffer(:, :, 1, freeSlot)     = kfMeas.obsMatrix;
                poolOut.numMeasBuffer(1, freeSlot)             = kfMeas.numMeas;
            end

        end

        %% ---------------- SS result ----------------
        function [SPF_ssMonitor] = setSSmonitor(alarmPerAxis, anyAlarm, ...
                separation, sigmaSeparation, protectionLevel, ...
                maxProtectionLevel)

            SPF_ssMonitor = struct( ...
                'alarmPerAxis',       logical(alarmPerAxis), ...
                'anyAlarm',           logical(anyAlarm), ...
                'separation',         separation, ...
                'sigmaSeparation',    sigmaSeparation, ...
                'protectionLevel',    protectionLevel, ...
                'maxProtectionLevel', maxProtectionLevel);
        end

        function [SPF_ssMonitor] = zeroSSmonitor

            axisAlarm = zeros(1, 3, 'logical');
            alarm = false;
            zeroVector = zeros(1, 3);
            scalar = 0.0;

            SPF_ssMonitor = STRUCT_SPF.setSSmonitor(axisAlarm, alarm, ...
                zeroVector, zeroVector, zeroVector, scalar);
        end

        %% ---------------- FSM state ----------------
        function [anchor] = setAnchor(valid, separation, covariance, epoch)
            % certified-clean coast, kept live: separation = host solution - anchor coast

            anchor = struct( ...
                'valid',      logical(valid), ...
                'separation', separation, ...
                'covariance', covariance, ...
                'epoch',      uint32(epoch));
        end

        function [anchor] = zeroAnchor

            valid = false;
            separation = zeros(CST_gnssHybrid.NO_STATES, 1);
            covariance = zeros(CST_gnssHybrid.NO_STATES, CST_gnssHybrid.NO_STATES);
            epoch = uint32(0);

            anchor = STRUCT_SPF.setAnchor(valid, separation, covariance, epoch);

        end

        function [sys] = setSys(mode, coastCov, probSep, pool, anchor, ...
                dwellCount, probationCount, coastCount, armed, armCount)

            sys = struct( ...
                'mode',           mode, ...           % CST_spfMode
                'coastCov',       coastCov, ...       % covariance of the protected solution
                'probSep',        probSep, ...        % (active filter) - (coast) [n x 1]
                'pool',           pool, ...
                'anchor',         anchor, ...
                'dwellCount',     dwellCount, ...
                'probationCount', probationCount, ...
                'coastCount',     coastCount, ...
                'armed',          logical(armed), ... % monitors active (after the warm-up)
                'armCount',       armCount);          % qualifying epochs so far
        end

        function [sys] = zeroSys

            mode = CST_spfMode.NOMINAL;
            n = CST_gnssHybrid.NO_STATES;
            coastCov = eye(n);
            probSep = zeros(n, 1);
            pool = STRUCT_SPF.zeroMonitorPool;
            anchor = STRUCT_SPF.zeroAnchor;
            dwellCount = 0.0;
            probationCount = 0.0;
            coastCount = 0.0;
            armed = false;
            armCount = 0.0;

            sys = STRUCT_SPF.setSys(mode, coastCov, probSep, pool, ...
                anchor, dwellCount, probationCount, coastCount, armed, armCount);
        end

        function [altXCheck] = setAltXCheck(suspect, altDiff)

            altXCheck = struct( ...
                'suspect', suspect, ...
                'altDiff', altDiff);

        end

        function [altXCheck] = zeroAltXCheck

            suspect = false;
            altDiff = 0.0;

            altXCheck = STRUCT_SPF.setAltXCheck(suspect, altDiff);
        end

        %% ---------------- telemetry / commands ----------------
        function [spoofTel] = setTel(info, kfCommand, nav)

            spoofTel = struct( ...
                'info',      info, ...
                'kfCommand', kfCommand, ...
                'nav',       nav);
        end

        function [spoofTel] = zeroTel

            info = STRUCT_SPF.zeroInfo;
            kfCommand = STRUCT_SPF.zeroCommand;
            nav = STRUCT_SPF.zeroNav;

            spoofTel = STRUCT_SPF.setTel(info, kfCommand, nav);
        end

        function [info] = setInfo(mode, ssAlarm, cpiAlarm, ...
                alarmPerAxis, maxProtectionLevel, qReval, revalComputed,...
                dwellCount, eventLatched, eventAnchorEpoch, eventProbationStarted,...
                eventProbationVetoed, eventHandback, anchorMissing, coastEpochs, ...
                inputFault, numMeasClamped, solveFault, ssRatio, cpiRatio, armed, ...
                resetInhibit, coastBudgetExceeded)

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
                'inputFault',            logical(inputFault), ...     % non-finite input this epoch
                'numMeasClamped',        logical(numMeasClamped), ... % host passed > MAX_MEAS rows
                'solveFault',            logical(solveFault), ...     % a covariance was unusable
                'ssRatio',               ssRatio, ...                 % [1x3] SS margin, alarm > 1
                'cpiRatio',              cpiRatio, ...                % [1x3] CPI margin, alarm > 1
                'armed',                 logical(armed), ...          % monitors active
                'resetInhibit',          logical(resetInhibit), ...   % host must not reset this epoch
                'coastBudgetExceeded',   logical(coastBudgetExceeded)); % coast longer than COAST_BUDGET_EPOCHS
        end

        function [info] = zeroInfo

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
            ssRatio               = zeros(1, 3);
            cpiRatio              = zeros(1, 3);
            armed                 = false;
            resetInhibit          = false;
            coastBudgetExceeded   = false;

            info = STRUCT_SPF.setInfo(mode, ssAlarm, cpiAlarm, ...
                alarmPerAxis, maxProtectionLevel, qReval, revalComputed,...
                dwellCount, eventLatched, eventAnchorEpoch, eventProbationStarted,...
                eventProbationVetoed, eventHandback, anchorMissing, coastEpochs, ...
                inputFault, numMeasClamped, solveFault, ssRatio, cpiRatio, armed, ...
                resetInhibit, coastBudgetExceeded);
        end

        function [command] = setCommand(startTrial)
            % startTrial: the host starts its trial filter as a copy of the coasting KF

            command = struct( ...
                'startTrial', logical(startTrial));
        end

        function [command] = zeroCommand

            startTrial = false;

            command = STRUCT_SPF.setCommand(startTrial);
        end

        function [nav] = setNav(applyCorrection, state, correction, covar, sigmaPosition)
            % applyCorrection (LATCH, COMMIT): (state, covar) go to the host's setKF as (x+, P+)

            nav = struct( ...
                'applyCorrection', logical(applyCorrection), ...
                'state',           state, ...         % [n x 1] x+
                'correction',      correction, ...    % [n x 1] same as an increment
                'covar',           covar, ...         % [n x n] P+
                'sigmaPosition',   sigmaPosition);    % [3 x 1] 1-sigma of position
        end

        function [nav] = zeroNav

            n = CST_gnssHybrid.NO_STATES;
            nav = STRUCT_SPF.setNav(false, zeros(n, 1), zeros(n, 1), eye(n), zeros(3, 1));
        end

    end
end
%------------------------------------------------------------------------
