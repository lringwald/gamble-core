# tests/

Synthetic validation of the nested framework — known-λ data-generating processes, so a failure means
the code is wrong rather than the data being hard. Run from the repo root:

    Rscript tests/test_nested_cut.R        # CUT / IV-imputation framework, λ recovery + CI coverage
    Rscript tests/test_nested_iv.R         # IV-nested framework, λ recovered and inside (0,1]
    Rscript tests/test_nested_cut_focal.R  # focal routing -> λ identification

The sampler's own feature suite lives at `codes/master_test_suite_sym.R` (35 checks) and is the gate
for any change to `codes/mnlogit_rcpp_sym.R`:

    Rscript codes/master_test_suite_sym.R

NOTE: that suite does NOT catch probability misallocation — its symmetric-HS checks only assert
baseline-invariance, which a prior can satisfy while sending mass to the wrong class. For changes to
the symmetric_hs path, also run the real-design root probe (see the symmetric-HS notes in git log).
