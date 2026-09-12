function spfLog = SPF_logInit(numEpochs)
%SPF_LOGINIT  Preallocate a per-epoch log of the gate telemetry (MATLAB
%   desktop use, not code generation). Fill it with SPF_logAppend once per
%   2 Hz epoch, then plot with SPF_plotLog.
%
%   spfLog = SPF_logInit(numEpochs)

spfLog.numEpochs   = numEpochs;
spfLog.count       = 0;
spfLog.mode        = zeros(numEpochs, 1);        % CST_spfMode as double (1 NOMINAL, 2 COAST, 3 PROBATION)
spfLog.numMeas     = zeros(numEpochs, 1);
spfLog.ssAlarm     = false(numEpochs, 1);
spfLog.cpiAlarm    = false(numEpochs, 1);
spfLog.alarmAxis   = false(numEpochs, 3);
spfLog.ssRatio     = zeros(numEpochs, 3);        % |d| / (k_FA sigma_SS), alarm at > 1
spfLog.cpiRatio    = zeros(numEpochs, 3);        % q / T_N of windows closed this epoch, alarm at > 1
spfLog.pl          = zeros(numEpochs, 1);        % max protection level [m]
spfLog.qReval      = zeros(numEpochs, 1);
spfLog.revalThr    = zeros(numEpochs, 1);        % chi-square threshold for this epoch's numMeas
spfLog.revalDone   = false(numEpochs, 1);
spfLog.dwell       = zeros(numEpochs, 1);
spfLog.coastEpochs = zeros(numEpochs, 1);
spfLog.nis         = zeros(numEpochs, 1);        % host KF consistency y' S^-1 y / numMeas (optional)
spfLog.inputFault  = false(numEpochs, 1);
spfLog.solveFault  = false(numEpochs, 1);
spfLog.clamped     = false(numEpochs, 1);
spfLog.evLatch     = [];                         % epochs of each event
spfLog.evProbStart = [];
spfLog.evVeto      = [];
spfLog.evHandback  = [];
spfLog.evNoAnchor  = [];
end
