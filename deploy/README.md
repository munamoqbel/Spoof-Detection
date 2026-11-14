# deploy/ — drop-in runtime set

The ten files here are the gate as it goes into the simulation: the same
code as the copies in the repository root, with a header of five sections
(description, inputs, outputs, assumptions and limitations, requirement
traceability) kept to 76 columns, and inline comments that say what a
line does and point to the paper's equation where there is one. Logic,
signatures, struct layouts and `%#codegen` are unchanged; a scripted
attack and a shadow sweep produce bitwise identical telemetry with either
set on the path. The root copies keep the long design notes.

| file | role |
|---|---|
| `SPF_gate.m` | persistent-state 2 Hz gate |
| `SPF_protectedNav.m` | NOMINAL / COAST / PROBATION FSM |
| `SPF_monitorPool.m` | overlapping-window bank |
| `SPF_ssMonitor.m`, `SPF_cpiMonitor.m` | per-window SS and CPI tests |
| `SPF_revalidation.m` | chi-square re-validation in COAST / PROBATION |
| `SPF_insCoast.m` | INS-only propagation helper |
| `STRUCT_SPF.m` | struct constructors, window open / close |
| `CST_spfParam.m` | constants (`ARM_EPOCHS` is the repository value 10; the host runs 120. `NOMINAL_BIAS` is 0 here; the host sets the nominal GNSS excursion it tolerates per axis, in metres) |
| `CST_spfMode.m` | mode enumeration |

Not included, the simulation provides its own: `CST_gnssHybrid`
(`NO_STATES`, `MAX_MEASURES`) and `matrixInv` (`[inv, invalid] = matrixInv(A)`).

The S and R the gate reads (`innovationCov`, `measNoiseCov`) carry the
nominal measurement noise on all the rows the host built. A host that
scales R up on large innovations or drops a satellite for some epochs
after a jump in its range (its spike filter) and so leaves rows out of its own update
hands the gate the nominal pair and its own pair beside it through
`STRUCT_SPF.setKfMeasBoth` or `kfMeasFromUpdateBoth` (`innovationCovHost`,
`measNoiseCovHost`, kept for the record, never read by a monitor); with
`setKfMeas` / `kfMeasFromUpdate` the two pairs are the same matrices.

The re-validation and the first probation window use the coast covariance
bounded by `COAST_SIGMA_MAX` (3 m) on the `COAST_BOUND_STATES`, the NED
position and the receiver clock bias (`STRUCT_SPF.boundCoastCov`; the
host sets the clock index as it sets `MONITORED_AXES`): a spoof larger than what the bound and the
range noise allow, about 24 m here, stays refused for as long as it lasts,
and clean GNSS is re-admitted after a real coast whose drift stays below
that, and the dwell also needs the CPI of the innovations against the
coast below T_N, so what the monitors would latch on is not re-admitted
(from about 14 m here); beyond that the gate stays in COAST until the
coast budget and the host's reset. The budget and the nav's own sigma use
the unbounded covariance.

The gate's own coast covariance: with `COAST_Q_DATASHEET` the gate
propagates every coast it keeps (windows, anchor, long coasts, the coast
in COAST) with a Q built from the IMU datasheet along the host's accumPhi
(`SPF_coastQ`: velocity from the VRW, attitude from the ARW, the bias
blocks from their in-run stability and walk, the position without noise,
the host's accumQ on the states the IMU does not drive), instead of the
host's tuning Q, which overstates an INS-only coast several times. The
host's P is read through `COAST_SIGMA_FLOOR_STATES / _VALUES`, raised on
the states its own consistency read shows it underestimates (the heading
of a wheeled vehicle without heading aiding), once per epoch, before a
coast starts from it and before the separation test subtracts it; at a
commit the host takes its own P back. `SPF_coastPropTel`, `SPF_coastQ`,
`SPF_coastStartCov` and `SPF_floorCov` are the four files; the flag, the
`IMU_*` values and the `COAST_*_STATES` in the host's layout are the
host's constants.

The return path: `RETURN_BANK` long INS-only coasts are kept beside the
anchor, the oldest spanning 80 to 120 s. At a latch and on every COAST
epoch the position the rows ask for is read against it: when the nav had
wandered more than `RETURN_WANDER` from that coast and the GNSS lands
within `RETURN_RADIUS` of it, the event is the end of a followed spoof,
not an onset (`eventReturn`): the nav is corrected back to the long coast,
velocity included, through `nav.applyCorrection` (in COAST too), and the
bound's floor becomes that coast's own sigma, so the clean GNSS
re-validates; the first probation window keeps the short bound.

Through the warm-up the monitors run in shadow: windows open and close as
in NOMINAL, the margins are logged, nothing latches and the alarm flags
stay down, and a shadow alarm restarts the arming count, so the gate arms
only after `ARM_EPOCHS` epochs on which it would not have alarmed.

Use: copy the folder contents next to the host code. Keep the root copies
off the simulation path; they are the annotated reference (design notes,
paper equation mapping, exception-handler rationale). Any change goes into
both sets; the tests in `runAllTests.m` run against the root copies, and
`tools/deployEquivalence.m` confirms the two sets still agree.
