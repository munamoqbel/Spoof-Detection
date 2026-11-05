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
scales R up on large innovations or leaves rows out of its own update
hands the gate the nominal pair and its own pair beside it through
`STRUCT_SPF.setKfMeasBoth` or `kfMeasFromUpdateBoth` (`innovationCovHost`,
`measNoiseCovHost`, kept for the record, never read by a monitor); with
`setKfMeas` / `kfMeasFromUpdate` the two pairs are the same matrices.

The re-validation and the first probation window use the coast covariance
bounded by `COAST_SIGMA_MAX` (3 m per monitored position axis,
`STRUCT_SPF.boundCoastCov`): a spoof larger than what the bound and the
range noise allow, about 24 m here, stays refused for as long as it lasts,
and clean GNSS is re-admitted after a real coast whose drift stays below
that, and the dwell also needs the CPI of the innovations against the
coast below T_N, so what the monitors would latch on is not re-admitted
(from about 14 m here); beyond that the gate stays in COAST until the
coast budget and the host's reset. The budget and the nav's own sigma use
the unbounded covariance.

Use: copy the folder contents next to the host code. Keep the root copies
off the simulation path; they are the annotated reference (design notes,
paper equation mapping, exception-handler rationale). Any change goes into
both sets; the tests in `runAllTests.m` run against the root copies, and
`tools/deployEquivalence.m` confirms the two sets still agree.
