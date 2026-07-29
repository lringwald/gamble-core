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

NEST_TREES <- list(
  AGMIP    = .agmip_tree,                                                   # crops -> Arable/Permanent (explicit)
  GLOBIOM  = function(classes) nest_tree_by_prefix(classes),                # Forests/Cropland/Pasture by prefix
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
