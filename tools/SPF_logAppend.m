function spfLog = SPF_logAppend(spfLog, spfTel, spfMeas, nis)
%SPF_LOGAPPEND  Append one epoch of gate telemetry to the log from
%   SPF_logInit. Call right after SPF_gate in the 2 Hz function:
%
%     spfLog = SPF_logAppend(spfLog, spfTel, spfMeas, nis);
%
%   nis (optional, recommended): the host KF's normalised innovation
%   squared for this epoch, y' * S^-1 * y / numMeas, formed from the live
%   rows of spfMeas. It is the KF's own consistency check: near 1 when R
%   and P are honest, well above 1 when the measurements are noisier than
%   R says (jamming, multipath), which is what the monitors then see too.

k = spfLog.count + 1;
spfLog.count = k;
if k > spfLog.numEpochs, return; end          % log full: keep counting, stop storing

info = spfTel.info;
m = double(spfMeas.numMeas);
spfLog.mode(k)        = double(info.mode);
spfLog.numMeas(k)     = m;
spfLog.ssAlarm(k)     = info.ssAlarm;
spfLog.cpiAlarm(k)    = info.cpiAlarm;
spfLog.alarmAxis(k,:) = info.alarmPerAxis;
spfLog.ssRatio(k,:)   = info.ssRatio;
spfLog.cpiRatio(k,:)  = info.cpiRatio;
spfLog.pl(k)          = info.maxProtectionLevel;
spfLog.qReval(k)      = info.qReval;
spfLog.revalDone(k)   = info.revalComputed;
if m > 0
    spfLog.revalThr(k) = CST_spfParam.REVAL_THRESHOLD_TABLE(min(m, numel(CST_spfParam.REVAL_THRESHOLD_TABLE)));
end
spfLog.dwell(k)       = info.dwellCount;
spfLog.coastEpochs(k) = info.coastEpochs;
spfLog.inputFault(k)  = info.inputFault;
spfLog.solveFault(k)  = info.solveFault;
spfLog.clamped(k)     = info.numMeasClamped;
if nargin >= 4
    spfLog.nis(k) = nis;
elseif m > 0
    y = spfMeas.innovation(1:m); S = spfMeas.innovationCov(1:m, 1:m);
    spfLog.nis(k) = (y' * (S \ y)) / m;
end
if info.eventLatched,          spfLog.evLatch(end+1)     = k; end
if info.eventProbationStarted, spfLog.evProbStart(end+1) = k; end
if info.eventProbationVetoed,  spfLog.evVeto(end+1)      = k; end
if info.eventHandback,         spfLog.evHandback(end+1)  = k; end
if info.anchorMissing,         spfLog.evNoAnchor(end+1)  = k; end
end
