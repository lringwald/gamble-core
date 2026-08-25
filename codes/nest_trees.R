# =============================================================================
# nest_trees.R — flexible nesting-tree construction for ANY land-use classification
# =============================================================================
# A tree (consumed by nested_cut_fit / mnlogit_nested_iv) is a named list whose
# leaves are character vectors of fine-class names. This module builds such trees
# from either
#   (a) an explicit class -> ancestor-PATH map           -> nest_tree_from_paths()
#   (b) naming conventions (prefix tokens)                -> nest_tree_by_prefix()
# and validates coverage against a class set             -> nest_tree_check().
#
# NEST_TREES registry: scheme -> function(classes) -> tree, for the shipped
# branches (AGMIP / GLOBIOM / BIOCLIMA). Add a scheme by adding one entry.
# =============================================================================

# --- (a) assemble a tree from class -> ancestor path -------------------------
# paths: named list; names = fine classes; value = character vector of nest names
#        from the ROOT down to the leaf's parent (empty = a root-level singleton).
#   e.g. Wheat = c("Cropland","Arable"); Grassland = character(0)
nest_tree_from_paths <- function(paths) {
  stopifnot(is.list(paths), !is.null(names(paths)))
  build <- function(items) {                       # items: named list class -> remaining path
    term    <- names(items)[vapply(items, length, 0L) == 0L]         # terminal here
    nonterm <- items[vapply(items, length, 0L) > 0L]
    if (length(nonterm) == 0L) return(unname(term))                  # pure terminal nest -> char vec
    node <- list()
    toks <- vapply(nonterm, `[[`, character(1), 1L)                  # first token of each
    for (tk in unique(toks)) {
      sub <- nonterm[toks == tk]
      node[[tk]] <- build(lapply(sub, function(p) p[-1]))            # drop first token, recurse
    }
    for (cl in term) node[[cl]] <- cl                                # terminal-at-mixed-node -> singleton
    node
  }
  build(paths)
}

# --- (a2) assemble a tree from a MAPPING TABLE -------------------------------
# The authoritative route when the taxonomy is CURATED BY HAND in the mapping CSV
# rather than encoded in class names. Two columns per row:
#   class_col  the LEAF class      e.g. "Cropland_other"
#   nest_col   its ancestor path   e.g. "Cropland/permanent"; BLANK = root-level leaf
# This frees leaf names from having to encode their own nest -- the constraint that
# forces `nest_tree_by_prefix` to rename classes to restructure them.
#
# `classes` is the OBSERVED class set (the fitted columns of Y); the mapping is
# filtered to it, so unused mapping rows are harmless. Classes the mapping does not
# mention (e.g. the derived "no_choice", which is a collapse, not a LUM row) become
# root-level leaves with a message rather than an error -- they are legitimately
# outside the curated taxonomy.
nest_tree_from_mapping <- function(map_dt, classes,
                                   class_col = "GLOBIOM_subclass",
                                   nest_col  = "GLOBIOM_nest", sep = "/") {
  miss <- setdiff(c(class_col, nest_col), names(map_dt))
  if (length(miss)) stop("mapping lacks column(s): ", paste(miss, collapse = ", "),
                         " -- add them, or use a prefix scheme")
  cl <- make.names(trimws(as.character(map_dt[[class_col]])))
  nz <- trimws(as.character(map_dt[[nest_col]])); nz[is.na(nz)] <- ""
  keep <- cl %in% classes & cl != ""
  cl <- cl[keep]; nz <- nz[keep]

  # A leaf must not be assigned two different nests -- that is a data-entry error in the
  # CSV and would silently drop one assignment, so fail loudly with the offenders.
  u <- unique(data.frame(cl, nz, stringsAsFactors = FALSE))
  dup <- u$cl[duplicated(u$cl)]
  if (length(dup)) stop("class(es) mapped to CONFLICTING nests in ", class_col, "/", nest_col, ": ",
                        paste(sprintf("%s -> {%s}", unique(dup),
                              vapply(unique(dup), function(d) paste(u$nz[u$cl == d], collapse = " | "),
                                     character(1))), collapse = "; "))
  paths <- setNames(lapply(u$nz, function(s)
                    if (!nzchar(s)) character(0) else strsplit(s, sep, fixed = TRUE)[[1]]), u$cl)

  orphan <- setdiff(classes, names(paths))
  if (length(orphan)) {
    message(sprintf("nest_tree_from_mapping: %d class(es) absent from the mapping -> root leaves: %s",
                    length(orphan), paste(orphan, collapse = ", ")))
    paths <- c(paths, setNames(rep(list(character(0)), length(orphan)), orphan))
  }
  # A NEST NAME MUST NOT COLLIDE WITH A CLASS NAME. `nest_tree_from_paths` writes terminal classes into
  # the same named list as nest nodes, and terminals are written LAST -- so a class named e.g.
  # "no_choice" silently OVERWRITES a nest of that name, deleting every leaf inside it. ("no_choice" is
  # reserved: the driver rewrites NODATA and unresolved rows to it.) Fail loudly instead.
  clash <- intersect(unlist(paths, use.names = FALSE), classes)
  if (length(clash)) stop("nest name(s) collide with CLASS name(s): ", paste(clash, collapse = ", "),
      ". A nest and a leaf cannot share a name -- the leaf would silently overwrite the nest and drop ",
      "its members. Rename the nest (e.g. 'Waterbodies' instead of 'no_choice').")
  tr <- nest_tree_from_paths(paths[classes])
  # A curated nest is filtered to the OBSERVED classes, so a nest whose other members have zero area
  # can prune down to a single leaf. That nest is degenerate -- its inclusive value equals its own
  # utility, so its lambda is unidentified. Warn rather than fail: it is a data outcome, not a typo.
  .n_leaf <- function(n) if (is.character(n)) length(n) else sum(vapply(n, .n_leaf, numeric(1)))
  .scan <- function(node, prefix = "") for (nm in names(node)) {
    ch <- node[[nm]]
    if (is.list(ch) || length(ch) > 0L) {
      if (.n_leaf(ch) == 1L && is.list(node) && !identical(unname(ch), nm))
        warning(sprintf("curated nest '%s%s' has only ONE observed class (%s) -> degenerate, its lambda is unidentified",
                        prefix, nm, paste(unlist(ch), collapse = ",")), call. = FALSE, immediate. = TRUE)
      if (is.list(ch)) .scan(ch, paste0(prefix, nm, "/"))
    }
  }
  .scan(tr)
  tr
}

# --- (b) build a tree from naming prefixes -----------------------------------
# Groups classes by their first token (before `sep`). Single-member groups become
# root singletons. Groups named in `subsplit` are split a second level by their
# 2nd token. Works for schemes whose names encode Type[_Subtype][_level] (GLOBIOM,
# BIOCLIMA). `flat` forces the listed classes to stay root singletons.
nest_tree_by_prefix <- function(classes, sep = "_", subsplit = character(0), flat = character(0)) {
  tok  <- strsplit(classes, sep, fixed = TRUE)
  grp1 <- vapply(tok, `[`, character(1), 1L)
  n1   <- table(grp1)
  paths <- setNames(vector("list", length(classes)), classes)
  for (i in seq_along(classes)) {
    cl <- classes[i]
    if (cl %in% flat || n1[[grp1[i]]] == 1L) { paths[[cl]] <- character(0); next }   # singleton
    if (grp1[i] %in% subsplit && length(tok[[i]]) >= 2L)
      paths[[cl]] <- c(grp1[i], paste(grp1[i], tok[[i]][2], sep = sep))              # 2-level
    else
      paths[[cl]] <- grp1[i]                                                         # 1-level
  }
  nest_tree_from_paths(paths)
}

# --- validate + pretty-print --------------------------------------------------
.tree_leaves <- function(node) if (is.character(node)) node else unlist(lapply(node, .tree_leaves), use.names = FALSE)
nest_tree_check <- function(tree, classes, print = TRUE) {
  lv <- .tree_leaves(tree)
  miss <- setdiff(classes, lv); extra <- setdiff(lv, classes); dup <- lv[duplicated(lv)]
  ok <- length(miss) == 0 && length(extra) == 0 && length(dup) == 0
  if (print) {
    render <- function(node, name, ind) {
      if (is.character(node)) {
        if (length(node) == 1L) cat(sprintf("%s- %s\n", strrep("  ", ind), node))
        else cat(sprintf("%s+ %s  {%s}\n", strrep("  ", ind), name, paste(node, collapse = ", ")))
      } else {
        cat(sprintf("%s+ %s\n", strrep("  ", ind), name))
        for (nm in names(node)) render(node[[nm]], nm, ind + 1)
      }
    }
    for (nm in names(tree)) render(tree[[nm]], nm, 0)
    cat(sprintf("[coverage] %d leaves / %d classes | %s%s%s\n", length(lv), length(classes),
        if (length(miss)) paste0("MISSING: ", paste(miss, collapse = ","), " ") else "",
        if (length(extra)) paste0("EXTRA: ", paste(extra, collapse = ","), " ") else "",
        if (ok) "OK" else "FAIL"))
  }
  invisible(ok)
}

# --- shipped branch trees (scheme -> builder) --------------------------------
.agmip_tree <- function(classes) {
  arable <- c("Wheat","Barley","Maize","Rice","Other_cereals","Fresh_vegetables","Dry_pulses",
              "Potatoes","Sugar_beet","Sunflower","Soybeans","Rapeseed","Flax_cotton_hemp","Cropland_arable_other")
  perm   <- c("Grapes","Olives","Fruits","Nuts","Cropland_permanent_other")
  forest <- c("Forests_primary","Forests_managed")
  singles<- setdiff(classes, c(arable, perm, forest))          # Grassland, Natural_unmanaged, Built_up_area, no_choice
  paths <- c(
    setNames(rep(list(c("Cropland","Arable")),   length(arable)), arable),
    setNames(rep(list(c("Cropland","Permanent")),length(perm)),   perm),
    setNames(rep(list("Forests"),                length(forest)), forest),
    setNames(rep(list(character(0)),             length(singles)),singles)
  )
  nest_tree_from_paths(paths[intersect(names(paths), classes)])
}

# GLOBIOM_ORG: split the ORGANIC twins onto their own level under Cropland and Pasture.
# WHY. In the flat GLOBIOM tree those nests mix TWO dimensions at one level -- management intensity
# (LI/HI/other) crossed with organic status (the _O / _*O suffix) -- while Forests is a single
# management gradient. Measured on 25,000 pixels (2026-08-08): lambda_Forests 0.306 [0.227, 0.377]
# (healthy) but lambda_Cropland -0.131 [-0.245, -0.056] and lambda_Pasture 0.032 [-0.145, 0.200].
# A negative lambda is not RUM-consistent. Cropland also has the LOWEST IV~design R2 (0.70 vs 0.83
# for Forests), so collinearity does not explain it -- the pattern tracks WHICH NESTS MIX TWO
# FACTORS. A nest assumes its members share unobserved factors so the logsum summarises them; organic
# is a cross-cutting attribute, not a sub-choice at that level. This tree tests that reading by
# giving organic its own level: Cropland -> {conventional{...}, organic{...}}.
#
# RESULT 2026-08-08: **REFUTED -- this does NOT fix lambda.** Matched settings (6,000 px, M=5, 600
# iter, 1 chain), flat GLOBIOM vs GLOBIOM_ORG:
#     lambda_Cropland  -0.172 [-0.332, +0.019]   ->  -0.274 [-1.007, +0.804]
#     IV~design R2      0.70                     ->   0.82   (MORE collinear, not less)
# The new sub-nests are unremarkable (conventional 0.166 [-0.145,0.434], organic 0.122
# [-0.473,0.385], both spanning 0). So the "intensity x organic crammed into one level" reading was
# a coincidence of 3 observations, not the cause. Kept as a selectable scheme because the deeper
# structure may still be defensible on its own terms, but it is NOT the remedy for the negative
# lambda. Do not re-test this expecting a fix.
.globiom_org_tree <- function(classes) {
  base <- nest_tree_by_prefix(classes)
  is_org <- function(v) grepl("(_O|O)$", v)          # Cropland_LIO / _other_O, Pasture_HIO ...
  for (nm in intersect(c("Cropland", "Pasture"), names(base))) {
    leaves <- base[[nm]]
    if (!is.character(leaves)) next                   # already nested -> leave alone
    org <- leaves[is_org(leaves)]; conv <- setdiff(leaves, org)
    if (length(org) && length(conv))                  # only split when BOTH sides exist
      base[[nm]] <- list(conventional = conv, organic = org)
  }
  base
}

# GLOBIOM_CROP: Cropland -> {arable, permanent}. `Cropland_other` is NOT an intensity like LI/HI --
# it is PERMANENT cropland (LUM 3001-3005: orchards, vineyards, olive groves), i.e. a crop-TYPE axis.
# Under the plain prefix tree it competes head-to-head with the arable intensities, flattening two
# different dimensions into one sibling set. This restores the level WITHOUT renaming any class:
# `nest_tree_by_prefix(subsplit=)` would need an arable/permanent token in the NAME
# (Cropland_arable_LI ...), breaking every existing managed class name; an explicit grouper does not.
# Composes with the organic split (recurses into conventional/organic when org = TRUE).
.globiom_crop_tree <- function(classes, org = FALSE) {
  base <- if (org) .globiom_org_tree(classes) else nest_tree_by_prefix(classes)
  is_perm <- function(v) grepl("^Cropland_other(_O)?$", v)
  .split <- function(node) {
    if (!is.character(node)) return(lapply(node, .split))   # recurse into conventional/organic
    perm <- node[is_perm(node)]; ara <- setdiff(node, perm)
    # BOTH sides need >= 2 members. A one-alternative nest is degenerate: its inclusive value IS its
    # own utility, so that nest's lambda is unidentified. Under the organic split each side holds a
    # single permanent class, so this correctly leaves those levels flat instead of minting two of them.
    if (length(perm) >= 2L && length(ara) >= 2L) list(arable = ara, permanent = perm) else node
  }
  if (!is.null(base[["Cropland"]])) base[["Cropland"]] <- .split(base[["Cropland"]])
  base
}

NEST_TREES <- list(
  AGMIP    = .agmip_tree,                                                   # crops -> Arable/Permanent (explicit)
  GLOBIOM  = function(classes) nest_tree_by_prefix(classes),                # Forests/Cropland/Pasture by prefix
  GLOBIOM_ORG = .globiom_org_tree,                                          # + Cropland/Pasture -> {conventional,organic}
  GLOBIOM_CROP     = function(classes) .globiom_crop_tree(classes, org = FALSE),  # + Cropland->{arable,permanent}
  GLOBIOM_CROP_ORG = function(classes) .globiom_crop_tree(classes, org = TRUE),   # + both levels
  BIOCLIMA = function(classes) nest_tree_by_prefix(classes, subsplit = "Cropland")  # + Cropland->{arable,permanent}
)

# convenience: build (and check) the tree for a named scheme + its classes
build_nest_tree <- function(scheme, classes, check = TRUE) {
  b <- NEST_TREES[[toupper(scheme)]]
  if (is.null(b)) stop("Unknown scheme: ", scheme, " (have: ", paste(names(NEST_TREES), collapse=", "), ")")
  tree <- b(classes)
  if (check) nest_tree_check(tree, classes)
  tree
}
