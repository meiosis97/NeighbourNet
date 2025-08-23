# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Set Default Parameters for Network Extraction
#'
#' This function updates or resets default settings used in downstream network
#' extraction and analysis within a Seurat object. Defaults include pruning
#' thresholds, response and predictor gene sets, co-expression measures, and
#' transformation functions for effect estimates.
#'
#' @param seurat.obj A \code{Seurat} object with a \code{NNet.mod} list stored in the \code{misc} slot.
#'                   This list is created by \code{\link{run.nn.reg}}.
#' @param clean.up A logical indicating whether to reset defaults to initial values
#'                 (effect-based assay, quadratic transformation, self-loop removal,
#'                 and pruning at cutoff 0.95). Default is \code{FALSE}.
#' @param defaults A named list of default parameters to update. Valid entries include:
#'   \itemize{
#'     \item{f} {A function mapping effect estimates to importance scores.
#'                      Default is \code{function(x) 2*x^2}.}
#'     \item{remove.self.loops} {Logical, whether to remove self-loops. Default: \code{TRUE}.}
#'     \item{assay} {Co-expression measure to use for downstream analysis.
#'                          Options: \code{"effect"}, \code{"p.val"}, or \code{"meta.network"}.}
#'     \item{predictors} {Character vector of predictor genes.}
#'     \item{responses} {Character vector of response genes.}
#'     \item{cutoff} {Numeric threshold for pruning, between 0 and 1. Default: \code{0.95}.}
#'   }
#'
#' @return A \code{Seurat} object with updated default settings stored in
#'         \code{misc$NNet.mod$defaults}.
#'
#' @details
#' The \code{defaults} list defines how networks are extracted and post-processed
#' in downstream analyses (e.g., pruning, TF activity inference, meta-network learning).
#' By setting \code{clean.up = TRUE}, all defaults are reset to standardized values.
#' Only valid arguments present in the current \code{NNet.mod$defaults} are updated.
#'
#' @examples
#' # Reset defaults to initial values
#' seurat.obj <- set.defaults(seurat.obj, clean.up = TRUE)
#'
#' # Update only the cutoff threshold and assay
#' seurat.obj <- set.defaults(
#'   seurat.obj,
#'   defaults = list(cutoff = 0.5, assay = "p.val")
#' )
#'
#' # Check updated defaults
#' Seurat::Misc(seurat.obj, "mod")$defaults
#'
#' @seealso \code{\link{run.nn.reg}}, \code{\link{prepare.reg}}
#'
#' @export
set.defaults <- function(seurat.obj,
                         clean.up = FALSE,
                         defaults = list()
                         ) {
  # Model extraction
  mod <- Seurat::Misc(seurat.obj, "NNet.mod")
  if (is.null(mod)) {
    stop("No model found in misc. Run prepare.reg / run.nn.reg first.")
  }

  # Reset defaults
  if (clean.up) {
    mod$defaults <- list(
      f = function(x) 2 * x^2,
      remove.self.loops = TRUE,
      assay = "effect",
      predictors = mod$gene.sets$predictors$genes,
      responses  = mod$gene.sets$responses$genes,
      cutoff = 0.95
    )
  }

  # Update only valid keys provided in `defaults`
  if (length(defaults)) {
    valid.args <- intersect(names(defaults), names(mod$defaults))
    if (length(valid.args)) {
      mod$defaults[valid.args] <- defaults[valid.args]
    }
  }

  # Coerce predictors/responses to valid sets
  mod$defaults$predictors <- intersect(mod$defaults$predictors, mod$gene.sets$genes)
  mod$defaults$responses  <- intersect(mod$defaults$responses,  mod$gene.sets$responses$genes)

  # Accept effect / p.val / meta.network; default to "effect" if NULL.
  if (is.null(mod$defaults$assay)) {
    mod$defaults$assay <- "effect"
  }
  mod$defaults$assay <- match.arg(mod$defaults$assay,
                                  choices = c("effect", "p.val", "meta.network"))

  # Validate `remove.self.loops`
  if (!is.logical(mod$defaults$remove.self.loops) || length(mod$defaults$remove.self.loops) != 1L) {
    stop("remove.self.loops must be a single logical (TRUE/FALSE).")
  }

  # Validate `f`
  if (!is.function(mod$defaults$f)) {
    stop("`f` must be a function (e.g., default: function(x) 2*x^2).")
  }

  # Validate `cutoff` in [0, 1]
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

  message("Default now is")
  str(mod$defaults)
  
  # Persist back to the same misc key we read from
  suppressWarnings(Seurat::Misc(seurat.obj, "NNet.mod") <- mod)
  seurat.obj
}
