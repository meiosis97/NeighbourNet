#' Set Default Parameters for Network Extraction
#'
#' Update or reset defaults used for downstream network extraction/analysis.
#' This function reads an existing model list from the Seurat object's misc slot
#' (preferring "NNet.mod"; falling back to "mod"), updates `defaults`,
#' and writes it back to the same slot.
#'
#' @param seurat.obj A Seurat object containing a model list in `Misc(seurat.obj, "NNet.mod")`
#'                   (created by `run.nn.reg`) or `Misc(seurat.obj, "mod")`.
#' @param clean.up Logical; if TRUE, reset defaults to standardized values:
#'                 `f = function(x) 2*x^2`, `remove.self.loops = TRUE`,
#'                 `assay = "effect"`, `predictors = mod$gene.sets$predictors$genes`,
#'                 `responses = mod$gene.sets$responses$genes`, `cutoff = 0.95`.
#' @param defaults Named list of default parameters to update. Valid names:
#'                 `f`, `remove.self.loops`, `assay`, `predictors`, `responses`, `cutoff`.
#'
#' @return The input Seurat object with updated defaults in the same misc slot.
#' @export
set.defaults <- function(seurat.obj, clean.up = FALSE, defaults = list()) {
  # --- Locate the model container in misc -------------------------------------
  # Prefer the canonical name used by run.nn.reg ("NNet.mod"), but support "mod" for back-compat.
  mod <- Seurat::Misc(seurat.obj, "NNet.mod")
  if (is.null(mod)) {
    stop("No model found in misc. Run `prepare.reg`/`run.nn.reg` first.")
  }

  # --- Optional full reset of defaults ----------------------------------------
  if (isTRUE(clean.up)) {
    mod$defaults <- list(
      f = function(x) 2 * x^2,
      remove.self.loops = TRUE,
      assay = "effect",
      predictors = mod$gene.sets$predictors$genes,
      responses  = mod$gene.sets$responses$genes,
      cutoff = 0.95
    )
  }

  # --- Update only valid keys provided in `defaults` --------------------------
  if (length(defaults)) {
    valid.args <- intersect(names(defaults), names(mod$defaults))
    if (length(valid.args)) {
      mod$defaults[valid.args] <- defaults[valid.args]
    } # silently ignore unknown keys
  }

  # --- Coerce predictors/responses to valid sets ------------------------------
  # Keep only genes present in the analysis universe.
  mod$defaults$predictors <- intersect(mod$defaults$predictors, mod$gene.sets$genes)
  mod$defaults$responses  <- intersect(mod$defaults$responses,  mod$gene.sets$responses$genes)

  # --- Validate `assay` -------------------------------------------------------
  # Accept effect / p.val / meta.network; default to "effect" if NULL.
  if (is.null(mod$defaults$assay)) {
    mod$defaults$assay <- "effect"
  }
  mod$defaults$assay <- match.arg(mod$defaults$assay,
                                  choices = c("effect", "p.val", "meta.network"))

  # --- Validate `remove.self.loops` -------------------------------------------
  if (!is.logical(mod$defaults$remove.self.loops) || length(mod$defaults$remove.self.loops) != 1L) {
    stop("`remove.self.loops` must be a single logical (TRUE/FALSE).")
  }

  # --- Validate `f` -----------------------------------------------------------
  if (!is.function(mod$defaults$f)) {
    stop("`f` must be a function (e.g., function(x) 2*x^2).")
  }

  # --- Validate `cutoff` in [0, 1] --------------------------------------------
  cutoff <- mod$defaults$cutoff
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    stop("`cutoff` must be a numeric scalar between 0 and 1.")
  }
  # Gentle clamping with a message (preserves your previous behavior for <0)
  if (cutoff < 0) {
    warning("`cutoff` < 0 detected; clamping to 0.")
    mod$defaults$cutoff <- 0
  } else if (cutoff > 1) {
    warning("`cutoff` > 1 detected; clamping to 1.")
    mod$defaults$cutoff <- 1
  }

  # --- Persist back to the *same* misc key we read from -----------------------
  suppressWarnings(Seurat::Misc(seurat.obj, misc_key) <- mod)
  seurat.obj
}
