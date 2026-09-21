#!/usr/bin/env Rscript
# =============================================================================
# test_nest_tree_auto.R — recursive nesting derived from class names
# =============================================================================
# The properties pinned here are the ones that make the output a usable nested-logit tree
# rather than a faithful transcription of the naming convention:
#   * no nest with a single child (its inclusive value is degenerate; lambda is unidentified)
#   * a level that splits nothing is skipped, not emitted as a one-child tier
#   * every class appears exactly once (coverage is what nest_tree_check enforces downstream)
# =============================================================================
source("codes/nest_trees.R")
np <- 0L; nf <- 0L
ok <- function(c_, m) { if (isTRUE(c_)) { np <<- np+1L; cat(sprintf("[PASS] %s\n", m)) }
                        else { nf <<- nf+1L; cat(sprintf("[FAIL] %s\n", m)) } }
leaves <- function(n) if (is.character(n)) n else unlist(lapply(n, leaves), use.names = FALSE)
depth  <- function(n, d = 1L) if (is.character(n)) d else max(vapply(n, depth, 0L, d + 1L))
nests  <- function(n) if (is.character(n)) character(0) else c(names(n)[vapply(n, is.list, TRUE)],
                                                               unlist(lapply(n, nests), use.names = FALSE))
one_child <- function(n) if (is.character(n)) FALSE else
  any(vapply(n, function(z) is.list(z) && length(z) == 1L, TRUE)) || any(vapply(n, one_child, TRUE))

cl <- c("Forest_primary","Forest_managed",
        "Crop_arable_wheat","Crop_arable_maize","Crop_arable_fodder_maize","Crop_arable_fodder_other",
        "Crop_perm_olives","Crop_perm_wine",
        "Urban","Water_inland")
t <- nest_tree_auto(cl)

ok(setequal(leaves(t), cl), "every class appears, none invented")
ok(!any(duplicated(leaves(t))), "no class appears twice")
ok(!one_child(t), "no nest has exactly one child")
ok("Crop" %in% names(t) && is.list(t$Crop), "a multi-member prefix becomes a nest")
ok(!("Urban" %in% nests(t)), "a single-member prefix stays a LEAF, not a one-child nest")
ok(!("Water" %in% nests(t)), "Water_inland alone does not create a Water nest")
ok(all(c("Crop_arable","Crop_perm") %in% names(t$Crop)), "the second token splits Crop (the level by_prefix misses)")
ok(depth(t) >= 3L, sprintf("recursion goes deeper than two levels (depth %d)", depth(t)))

# a level that discriminates nothing must be skipped rather than emitted
cl2 <- c("A_x_one","A_x_two","A_x_three")
t2 <- nest_tree_auto(cl2)
ok(!one_child(t2), "a shared middle token does not produce a chain of one-child nests")
ok(setequal(leaves(t2), cl2), "skipping a tier loses no class")

# classes shorter than the current depth sit beside deeper nests
cl3 <- c("B", "B_deep_one", "B_deep_two")
t3 <- nest_tree_auto(cl3)
ok(setequal(leaves(t3), cl3), "a short class coexists with its deeper siblings")

# max_depth caps the recursion
ok(depth(nest_tree_auto(cl, max_depth = 1L)) <= depth(t), "max_depth caps the tree")

# the AUTO scheme is reachable through the registry entry point
t4 <- build_nest_tree("AUTO", cl, check = FALSE)
ok(identical(leaves(t4), leaves(t)), "build_nest_tree('AUTO') routes to nest_tree_auto")

cat(sprintf("RESULT: %d/%d checks passed\n", np, np + nf))
if (nf > 0L) quit(status = 1)
