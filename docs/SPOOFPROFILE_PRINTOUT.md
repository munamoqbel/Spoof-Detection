# Reading the spoofProfile printout

`prof = spoofProfile(spfMeasLog, tKf, loggedData, tRef, posRef, spfLog)` prints
its reading in the order below. Every block is answering one question; the
question is given first, then the lines of that block and how to read them.
Units are metres and seconds unless the header says m/s. "Epoch" is one 2 Hz
gate call; "record" is the spfMeas the gate received at that epoch; "log" is
the host's pairing log (the receiver's ranges per slot); "row" is one
measurement line the host built from the log.

## 1. What the log holds (`pairingLogUnpack:` lines)

    pairingLogUnpack: 759 epochs, 732 with GNSS; systems: 1 = gpsRange / navstarStatus (15 slots, median 5 used) ...

The epochs in the log, how many held any used range, and the systems with
the log fields they come from. A slot is used when its range is finite and
not zero and satValid is true. Then, only when they occur:

- `ranges with satValid false, not used`: the receiver gave a range the
  host flagged out.
- `corrValid false while ranges are present`: the receiver's corrections
  flag was down while it still gave ranges (per system, epochs).
- `WARNING ... sysIdLog count the rows per system differently`: the host's
  own system ids do not match the used slots; the values it holds are
  listed. The gate never reads sysIdLog, so this only affects the tool's
  per-system labels.

## 2. How the record's rows map to systems

    spoofProfile: 759 epochs; range rows per epoch median 12; an offset group is 2 rows or more, 10 m or more from the rest
       row systems from the pairing log: ... their count matched the record's range rows on 587 of 654 epochs

The tool gives each range row of a record the system of the used slots in
slot order. The count line says on how many epochs the number of used
slots equals the number of range rows; where they differ the host dropped
or added a row and the per-system labels on that epoch are unreliable. A
WARNING follows when that happens on more than one epoch in ten.

## 3. The two-group split (innovations)

    no offset group on 3 consecutive epochs: ...
    or: onset / end / size at onset / rows in the offset group / the offset rows belong to system / course of the offset

Per epoch the range innovations are sorted; when a gap wider than the
spread inside the groups splits them, the group farther from zero is the
offset group. This only works when clean and spoofed rows sit in the same
record and the spoof is one common offset. A position spoof gives every
satellite its own offset, so this block usually reads "no offset group" or
a spurious short one; read the fits instead.

## 4. Stretches (innovations)

    stretches with a system 10 m or more off (off epochs closer than 10 s belong to one stretch):
    1) 27.06 to 31.06 s (epochs 16 to 24), 9 off epochs, median innovation: system 3 -11.2 m on 9;
    2) 239.57 to 290.07 s (epochs 442 to 543), 16 off epochs, median innovation: system 1 -2.4 m on 16; ...   <- the tables below

An epoch is "off" when some system's mean range innovation is 10 m or
more. Off epochs closer than 10 s form one stretch. Each stretch gives its
span, its off epochs, and the median innovation per system with the
number of off epochs that system contributed. The strongest stretch gets
the tables. A short stretch at the start of a run is usually the filter's
own convergence, listed so that it does not pull the tables back there.

## 5. Mean range innovation per system over the strongest stretch

    t [s]  system 1 [m] (rows)  system 2 [m] (rows)  system 3 [m] (rows)
    239.57        +42.9 ( 5)        +248.1 ( 2)       -600.5 ( 1)
    240.07            - ( 0)            - ( 0)            - ( 0)

The first 8 epochs of the stretch in a row, then spread to its end. Each
cell is the mean range innovation of that system and, in brackets, the
number of range rows the host built for it. A dash with (0) is an epoch on
which the host built no row for that system: this is the host's row set
epoch by epoch. For a position spoof the mean depends on which satellites
are in the record, so it moves from epoch to epoch; the fit below is the
stable number.

    row sets in this stretch (102 epochs): no range rows on 86; systems 1 2 3 on 13 (...); systems 1 2 on 3 (...)

How often each combination of systems had rows in the stretch, with the
median innovation per system on those epochs.

## 6. The position and clock the rows ask for (innovations)

    system      N [m]    E [m]    D [m]   horiz   clock   rms  epochs
    system 1   +379.3  +1152.9    -0.7  1213.7   -71.5   0.0      13

A least-squares fix on each system's range rows of a record (4 rows or
more): the position step and clock the rows ask for relative to the host's
prediction. Medians over the stretch; "epochs" is how many records had a
fit. "rms" is what the rows leave after the fit: 0 means exactly
determined (4 rows) or noiseless data, a few metres means a consistent
position with code noise, tens of metres means the rows do not describe
one position. The sign follows the host's H convention, and the
comparison block (section 10) says whether it is the reverse of the
reference fit. The per-epoch table for the strongest system follows;
`<4` marks epochs with fewer than 4 rows of that system.

The note under every fit table: with all satellites above the horizon the
vertical and the clock (or the drift) point almost the same way in the
fit, so a pair of opposite signs is one small term seen twice. The
horizontal part and the rms are the firm numbers.

## 7. The Doppler rows (innovations)

    the Doppler rows: mean range-rate innovation per system over the stretch system 1 +0.03 m/s on 16; ...
    system     N [m/s]   E [m/s]   D [m/s]   horiz   drift    rms  epochs
    the same before the stretch and after it ...

The same reading on the range-rate rows: the velocity and clock drift they
ask for relative to the host's prediction, medians over the stretch, and
the medians over the minute before and after the stretch for comparison.
A static spoofed position reads near zero here; a spoofer dragging the
position without matching Doppler reads as a velocity in section 11 and
none here.

## 8. The reference-based profile (only with posRef)

    493 slot-epochs with a range more than 100 km from the slot's geometry dropped (first at 236.54 s, epoch 436): ...
    on 53 of these 53 epochs the record holds fewer range rows than the log has used slots: the host left a row out

Each used range minus the geometric range from the reference position to
the satellite in its slot. A difference beyond 100 km means the range and
the position in that slot are not of the same satellite or instant; such
slots are dropped and counted, and the second line says whether the host
also left those rows out.

    against the reference position: ... inter-system bias relative to system 2: system 1 -0.6 m; ...
    the clock read on 630 of 732 epochs with ranges (two systems within 5 m), interpolated across the others

What is left in the residual is the receiver clock, a constant bias per
system, and the spoof. The biases are run medians relative to the system
present most often. The clock is read on the epochs where two systems
agree within 5 m and interpolated across the others.

    stretches with a system 10 m or more off ...        (as section 4, on the offset against the reference)
    onset / end / size at onset ... on system 2 (median of the first 3 epochs), largest ...
    the other systems in this stretch: median offset ...
    on 53 of the 90 epochs of this stretch at least one range lies beyond 100 km ... (g in the timeline)
    on 90 of the 90 epochs the clock is interpolated (gaps up to 29.0 s): ...
    within-system spread of system 2 in this stretch: median 572.2 m (...)

Onset and end are the first and last off epoch of the strongest stretch.
The size is the median offset of the named system over its first 3
epochs, the largest its extreme. The garbage line says how many epochs of
the stretch hold a range beyond 100 km; the medians mix those epochs with
the consistent ones, so for a run like that read the fit tables, which
use only the consistent slots. The interpolation line says on how many
epochs no two systems agreed, so the clock there rests on its line. The
spread line is the key to the spoof's kind: at code noise all satellites
of the system carry the same offset (a clock-type spoof); hundreds of
metres means every satellite has its own offset (a position-type spoof)
or only some satellites are spoofed.

    t [s]  system 1 [m] (slots)  system 2 [m] (slots)  system 3 [m] (slots)   clock [m]  spread 2 [m]
    236.06        +708.9 ( 1)         +335.1 ( 3)         +433.7 ( 1)      -4.5i        429.6

The course: per system the median offset against the reference and the
number of usable slots, the receiver clock with `i` where it was
interpolated, and the within-system spread of the strongest system.

    row sets in this stretch (90 epochs): systems 2 on 53 (...); systems 1 2 3 on 37 (...)

As in section 5, on the log's usable slots instead of the record's rows.

## 9. The position and clock the ranges describe (reference)

    system      N [m]    E [m]    D [m]   horiz   clock   rms  epochs
    system 1   -381.2  -1158.4    -0.3  1219.5   +85.3   0.0      36
    per epoch for system 2: ...
    system 2 moves within 14.4 N 22.3 E 55.4 D m (10th to 90th percentile) on 29 fitted epochs

The same least-squares fix as section 6, on the log's usable slots of
each system, with the receiver clock and the system bias removed: the
position and clock the spoofed ranges describe, relative to the truth. The
"moves within" line is the spread of the fitted position over the stretch:
a few metres is a fixed false position, hundreds of metres a moving one.

    the fix before the stretch and after it, medians N E D clock [m] ...
                 before, 176.09 to 235.56 s                 after, 281.04 to 340.08 s
    system 1      -0.7     -7.2     +4.4     -3.0  (109)      +0.4     +2.3     +1.4     +0.5  (92)

The fix over the minute before the stretch and the minute after. A sky
that is clean again reads near zero after; a difference between before
and after is a shift in the receiver's output that the nav will follow at
the hand-back.

    the Doppler against the reference: ... the receiver clock drift read on 667 of 732 epochs ..., median +0.007 m/s; the rate offset per system over the stretch: ...
    system     N [m/s]   E [m/s]   D [m/s]   horiz   drift    rms  epochs
    the same before the stretch and after it ...

Section 7 against the reference: the log's range rates minus the
line-of-sight satellite velocity leave the receiver clock drift (read
where two systems agree within 0.25 m/s) and the velocity and drift each
system's rates describe. This needs the satellite velocity in satPosLog
and a stationary reference; otherwise the tool says so and skips it.

## 10. The two fits side by side

    system 1: the two fits agree within 6 m horizontal with the sign reversed on 13 off epochs (the position columns of the host's H carry +LOS ...)

Per system, the innovation fit (section 6) against the reference fit
(section 9) on the off epochs: the same, the same with the sign reversed
(the host's position error state is estimate minus truth), both in the
clock only, or apart by a stated distance (the host's prediction sat that
far from the reference, or the record's rows are not the log's slots).

## 11. The timeline

    226.07 rows  AAAAAAAAA. .......... .......A.. .A...A.... ..A..A....
           log   AAAAAAAAA. .......... AggAggAAgg AAggAAgggA gAAgAAgggA
           mode  NNNNNNNNNN NNNNNNNNNN NNNNNNNCCC CCCCCCCCCC CCCCCCCCCC

One character per epoch, 50 per line, from 10 s before the strongest
stretch to 10 s after it. `rows` is the record's row set: `.` no range
rows, `A` all systems, a digit that system alone, `x` some systems (listed
under the strip). `log` is the log's used slots in the same code, printed
only on lines where it differs from `rows`, with `g` for an epoch holding
a range beyond 100 km. `mode` is the gate mode from spfLog: `N` nominal,
`C` coast, `P` probation. Reading the three together shows what the
receiver gave, what the host built from it, and what the gate did.

    in this window the record holds no range row on 78 epochs where the log has used slots (the host built none), and rows on 0 epochs where the log has no used slot
    every range in the log was used in this window ... / the log holds ranges the host did not use (satValid false): ...
    the log's corrValid flag is false while ranges are present: ...
    gate events: latch at 239.57 s (epoch 442); probation start at 293.57 s (epoch 550); hand-back at 302.57 s (epoch 568)

The counts behind the strips, and the gate's events from spfLog.

## The result struct

`prof` carries every number printed: the per-epoch arrays (`perSys`,
`nPerSys`, `inj`, `nSlots`, `spread`, `clock`, `clockGap`, `fitInn`,
`fitRef`, `fitRateInn`, `fitRateRef`, `rateDrift`, `badEpochs`), the
stretches (`stretchesInn`, `stretchesRef`), the onset, end, size and
system of the strongest stretch (`injOnset`, `injEnd`, `injSize`,
`injMax`, `injSys`), the before and after medians (`fixBefore`,
`fixAfter`, `rateBefore`, `rateAfter`), `fitSign`, and the timeline counts
(`nRecNoRows`, `nLogNoSlots`). The fit arrays are [epochs x 6 x systems]
with the columns N, E, D, clock (or drift), rms, rows.
