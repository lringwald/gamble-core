# tests/

## One command for everything

    Rscript tests/run_all.R              # every test + HTML report and figure
    Rscript tests/run_all.R fast         # skip the two long suites
    Rscript tests/run_all.R no-report    # table only

Writes `output/report/test_suite.html` (self-contained, with the figure embedded) and
`output/report/test_suite.png`. The **mode is stamped on both artefacts**, so a green `fast` run can
never be mistaken for a green full run.

### Three statuses, not two

| status | meaning |
|---|---|
| `PASS` | assertions ran and all held |
| `FAIL` | non-zero exit, a `[FAIL]` line, **or no assertions at all** |
| `SKIP` | the file printed `SKIP:` because its inputs are legitimately absent |

A file that asserts **nothing** is a FAIL, not a pass. Exit code 0 with zero assertions means it did
not run — an unreadable file, an early return, a stub. Silent success is the failure mode this repo
keeps hitting, and a green runner that hides it is worse than no runner.

`SKIP` exists so a test whose fixtures are genuinely missing can say so out loud instead of being
scored as either. `projects/` and `output/` are git-ignored, so project-dependent tests
(`test_target_class_split.R`, which needs the BMLEH cascade built) skip on a fresh clone and run
wherever that build exists.

### Membership is by discovery

Every `tests/test_*.R` is found and run. The list inside `run_all.R` fixes only the **order** of the
ones that should go first; anything undeclared still runs and is reported as a note. This matters:
the list used to be the membership, hand-maintained, and three feature tests (`re_mean_shift`,
`joint_shrink`, `mundlak`) sat outside it and never ran at all.

The two long gates in `codes/` are appended unless `fast`. They are self-contained — they locate
*source* files, not fixtures — so they always run and never skip.

## The gates

    Rscript codes/test_suite_lu_pixel.R    # 37 checks — the gate for codes/mnlogit_rcpp_sym.R
    Rscript codes/test_suite_ls_count.R    # count-feature validation

`test_suite_lu_pixel.R` was called `master_test_suite_sym.R` until 2026-08-25; the rename left stale
instructions in several READMEs, since corrected.

> **That suite does NOT catch probability misallocation.** Its symmetric-HS checks only assert
> baseline-invariance, which a prior can satisfy while sending mass to the wrong class. For changes
> to the `symmetric_hs` path, also run the real-design root probe (see the symmetric-HS notes in the
> git log).

## The nested-framework tests

Synthetic validation with known-λ data-generating processes, so a failure means the code is wrong
rather than the data being hard:

    Rscript tests/test_nested_cut.R        # CUT / IV-imputation framework, λ recovery + CI coverage
    Rscript tests/test_nested_iv.R         # IV-nested framework, λ recovered and inside (0,1]
    Rscript tests/test_nested_cut_focal.R  # focal routing -> λ identification
