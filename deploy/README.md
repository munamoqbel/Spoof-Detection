# deploy/ — drop-in runtime set

The ten files here are the gate as it goes into the simulation: the same
code as the copies in the repository root, with the header paragraphs and
inline commentary cut to a few lines each. Logic, signatures, struct
layouts and `%#codegen` are unchanged; a scripted attack and a shadow
sweep produce bitwise identical telemetry with either set on the path.

| file | role |
|---|---|
| `SPF_gate.m` | persistent-state 2 Hz gate |
| `SPF_protectedNav.m` | NOMINAL / COAST / PROBATION FSM |
| `SPF_monitorPool.m` | overlapping-window bank |
| `SPF_ssMonitor.m`, `SPF_cpiMonitor.m` | per-window SS and CPI tests |
| `SPF_revalidation.m` | chi-square re-validation in COAST / PROBATION |
| `SPF_insCoast.m` | INS-only propagation helper |
| `STRUCT_SPF.m` | struct constructors, window open / close |
| `CST_spfParam.m` | constants (`ARM_EPOCHS` is the repository value 10; the host runs 120) |
| `CST_spfMode.m` | mode enumeration |

Not included, the simulation provides its own: `CST_gnssHybrid`
(`NO_STATES`, `MAX_MEASURES`) and `matrixInv` (`[inv, invalid] = matrixInv(A)`).

Use: copy the folder contents next to the host code. Keep the root copies
off the simulation path; they are the annotated reference (design notes,
paper equation mapping, exception-handler rationale). Any change goes into
both sets; the tests in `runAllTests.m` run against the root copies, and
`tools/deployEquivalence.m` confirms the two sets still agree.
