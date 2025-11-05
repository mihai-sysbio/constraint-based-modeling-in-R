# Helper functions shared by the Gurobi solver integration --------------------

.gurobi_available <- function() {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop("The 'gurobi' package is required to use the Gurobi solver.", call. = FALSE)
  }
}

.gurobi_build_params <- function(method, tolerance) {
  params <- list(OutputFlag = 0L,
                 FeasibilityTol = tolerance)

  method <- tolower(method)
  params$Method <- switch(method,
                          primal = 0L,
                          dual = 1L,
                          barrier = 2L,
                          mip = NULL,
                          auto = NULL,
                          NULL)

  if (is.null(params$Method)) {
    params$Method <- NULL
  }

  params
}

.gurobi_constraint_entries <- function(row_mat, lb, ub, type) {
  senses <- character()
  rhs <- numeric()
  mats <- list()

  type <- toupper(type)
  lb <- ifelse(is.na(lb), -Inf, lb)
  ub <- ifelse(is.na(ub),  Inf, ub)

  add_entry <- function(s, r) {
    senses <<- c(senses, s)
    rhs <<- c(rhs, r)
    mats[[length(senses)]] <<- row_mat
  }

  switch(type,
         "E" = {
           if (is.finite(lb)) {
             add_entry("=", lb)
           } else if (is.finite(ub)) {
             add_entry("=", ub)
           }
         },
         "L" = {
           if (is.finite(lb)) {
             add_entry(">", lb)
           }
         },
         "U" = {
           if (is.finite(ub)) {
             add_entry("<", ub)
           }
         },
         "D" = {
           if (is.finite(lb)) {
             add_entry(">", lb)
           }
           if (is.finite(ub)) {
             add_entry("<", ub)
           }
         },
         "F" = {
           # free row, skip
         },
         {
           warning(sprintf("Unknown constraint type '%s', treated as equality.", type))
           if (is.finite(lb)) {
             add_entry("=", lb)
           } else if (is.finite(ub)) {
             add_entry("=", ub)
           }
         })

  list(sense = senses, rhs = rhs, mats = mats)
}

.gurobi_merge_rows <- function(rows, ncols) {
  if (length(rows) == 0L) {
    Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                         dims = c(0, ncols))
  } else {
    out <- rows[[1L]]
    if (length(rows) > 1L) {
      for (idx in seq.int(2L, length(rows))) {
        out <- Matrix::rBind(out, rows[[idx]])
      }
    }
    out
  }
}

.gurobi_modelsense <- function(direction) {
  ifelse(tolower(direction) == "max", "max", "min")
}

#' Structure of LPproblem_gurobi Class
#'
#' A class structure to link LP problems solved with the Gurobi optimizer.
#' Class is derived from \link{LPproblem}.
#'
#' @exportClass LPproblem_gurobi
setClass("LPproblem_gurobi",
         contains = "LPproblem")

.gurobi_validate_method <- function(method) {
  allowed <- c("auto", "primal", "dual", "barrier", "mip")
  method <- tolower(method)
  if (!method %in% allowed) {
    warning(sprintf("Method '%s' not available for Gurobi, using 'auto' instead.", method))
    return("auto")
  }
  method
}

#------------------------------------------------------------------------------#
# Function instances for the Gurobi solver                                     #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature = "LPproblem_gurobi",
          definition = function(.Object,
                                name,
                                method) {

            .gurobi_available()

            env <- new.env(parent = emptyenv())
            env$name <- name
            env$model <- NULL
            env$result <- NULL
            env$row_map <- list()
            env$status <- NULL
            env$params <- list()

            method <- .gurobi_validate_method(method)

            .Object@ptr <- env
            .Object@solver <- "gurobi"
            .Object@method <- method
            .Object@tol_bnd <- COBRAR_SETTINGS("TOLERANCE")

            return(.Object)
          }
)

#' @rdname loadLPprob-methods
#' @aliases loadLPprob,LPproblem_gurobi
setMethod("loadLPprob", signature(lp = "LPproblem_gurobi"),

          function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype, lpdir,
                   rub = NULL, ctype = NULL) {

            .gurobi_available()

            env <- lp@ptr

            if (is.null(rub)) {
              rub <- rlb
            }

            if (!inherits(mat, "dgCMatrix")) {
              mat <- as(mat, "dgCMatrix")
            }

            nCols <- as.integer(nCols)
            nRows <- as.integer(nRows)

            lb <- as.numeric(lb)
            ub <- as.numeric(ub)
            lb[is.na(lb)] <- -Inf
            ub[is.na(ub)] <- Inf

            obj <- as.numeric(obj)
            rlb <- as.numeric(rlb)
            rub <- as.numeric(rub)

            r_entries <- list()
            sense <- character()
            rhs <- numeric()
            row_map <- vector("list", length = nRows)

            for (i in seq_len(nRows)) {
              row_mat <- mat[i, , drop = FALSE]
              entries <- .gurobi_constraint_entries(row_mat, rlb[i], rub[i], rtype[i])

              if (length(entries$sense) > 0) {
                start_idx <- length(sense) + 1L
                sense <- c(sense, entries$sense)
                rhs <- c(rhs, entries$rhs)
                r_entries <- c(r_entries, entries$mats)
                end_idx <- length(sense)
                row_map[[i]] <- seq.int(start_idx, end_idx)
              } else {
                row_map[[i]] <- integer(0)
              }
            }

            model <- list()
            model$A <- .gurobi_merge_rows(r_entries, nCols)
            model$obj <- obj
            model$lb <- lb
            model$ub <- ub
            model$rhs <- rhs
            model$sense <- sense
            model$modelsense <- .gurobi_modelsense(lpdir)

            if (is.null(ctype)) {
              model$vtype <- rep("C", nCols)
            } else {
              model$vtype <- as.character(ctype)
            }

            env$model <- model
            env$row_map <- row_map
            env$result <- NULL
            env$status <- NULL
            env$params <- .gurobi_build_params(lp@method, lp@tol_bnd)

            invisible(lp)
          }
)

#' @rdname setObjDirection-methods
#' @aliases setObjDirection,LPproblem_gurobi
setMethod("setObjDirection", signature(lp = "LPproblem_gurobi"),
          function(lp, lpdir) {
            env <- lp@ptr
            env$model$modelsense <- .gurobi_modelsense(lpdir)
            invisible(lp)
          }
)

#' @rdname addCols-methods
#' @aliases addCols,LPproblem_gurobi
setMethod("addCols", signature(lp = "LPproblem_gurobi"),
          function(lp, ncols) {
            env <- lp@ptr
            ncols <- as.integer(ncols)
            if (ncols <= 0) {
              return(invisible(lp))
            }

            if (is.null(env$model)) {
              env$model <- list(A = Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                                         x = numeric(0),
                                                         dims = c(0, ncols)),
                                obj = numeric(ncols),
                                lb = rep(-Inf, ncols),
                                ub = rep(Inf, ncols),
                                rhs = numeric(0),
                                sense = character(0),
                                modelsense = "max",
                                vtype = rep("C", ncols))
            } else {
              add_mat <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                              dims = c(nrow(env$model$A), ncols))
              env$model$A <- Matrix::cBind(env$model$A, add_mat)
              env$model$obj <- c(env$model$obj, rep(0, ncols))
              env$model$lb <- c(env$model$lb, rep(-Inf, ncols))
              env$model$ub <- c(env$model$ub, rep(Inf, ncols))
              env$model$vtype <- c(env$model$vtype, rep("C", ncols))
            }

            invisible(lp)
          }
)

#' @rdname addRows-methods
#' @aliases addRows,LPproblem_gurobi
setMethod("addRows", signature(lp = "LPproblem_gurobi"),
          function(lp, nrows) {
            env <- lp@ptr
            nrows <- as.integer(nrows)
            if (nrows <= 0) {
              return(invisible(lp))
            }

            if (is.null(env$model)) {
              env$model <- list(A = Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                                         x = numeric(0),
                                                         dims = c(nrows, 0)),
                                obj = numeric(0),
                                lb = numeric(0),
                                ub = numeric(0),
                                rhs = rep(0, nrows),
                                sense = rep("=", nrows),
                                modelsense = "max",
                                vtype = character(0))
            } else {
              add_mat <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                              dims = c(nrows, ncol(env$model$A)))
              env$model$A <- Matrix::rBind(env$model$A, add_mat)
              env$model$rhs <- c(env$model$rhs, rep(0, nrows))
              env$model$sense <- c(env$model$sense, rep("=", nrows))
            }

            env$row_map <- c(env$row_map, replicate(nrows, integer(0), simplify = FALSE))

            invisible(lp)
          }
)

#' @rdname loadMatrix-methods
#' @aliases loadMatrix,LPproblem_gurobi
setMethod("loadMatrix", signature(lp = "LPproblem_gurobi"),
          function(lp, ne, ia, ja, ra) {
            stop("loadMatrix is not used for the Gurobi backend.")
          }
)

#' @rdname setColsBndsObjCoefs-methods
#' @aliases setColsBndsObjCoefs,LPproblem_gurobi
setMethod("setColsBndsObjCoefs", signature(lp = "LPproblem_gurobi"),
          function(lp, j, lb, ub, obj_coef, type = NULL) {
            env <- lp@ptr
            idx <- as.integer(j)
            env$model$lb[idx] <- as.numeric(lb)
            env$model$ub[idx] <- as.numeric(ub)
            env$model$obj[idx] <- as.numeric(obj_coef)
            env$model$lb[is.na(env$model$lb)] <- -Inf
            env$model$ub[is.na(env$model$ub)] <- Inf
            if (!is.null(type)) {
              env$model$vtype[idx] <- as.character(type)
            }
            invisible(lp)
          }
)

#' @rdname setColsKind-methods
#' @aliases setColsKind,LPproblem_gurobi
setMethod("setColsKind", signature(lp = "LPproblem_gurobi"),
          function(lp, j, kind) {
            env <- lp@ptr
            env$model$vtype[as.integer(j)] <- as.character(kind)
            invisible(lp)
          }
)

#' @rdname setRowsBnds-methods
#' @aliases setRowsBnds,LPproblem_gurobi
setMethod("setRowsBnds", signature(lp = "LPproblem_gurobi"),
          function(lp, i, lb, ub , type) {
            env <- lp@ptr

            inds <- as.integer(i)
            lb_vals <- rep(as.numeric(lb), length.out = length(inds))
            ub_vals <- rep(as.numeric(ub), length.out = length(inds))
            type_vals <- rep(toupper(as.character(type)), length.out = length(inds))

            for (pos in seq_along(inds)) {
              map <- env$row_map[[inds[pos]]]
              if (length(map) == 0L) {
                next
              }

              tval <- type_vals[pos]

              if (tval == "D") {
                if (length(map) >= 1L) {
                  env$model$sense[map[1L]] <- ">"
                  env$model$rhs[map[1L]] <- lb_vals[pos]
                }
                if (length(map) >= 2L) {
                  env$model$sense[map[2L]] <- "<"
                  env$model$rhs[map[2L]] <- ub_vals[pos]
                }
              } else if (tval == "L") {
                env$model$sense[map] <- ">"
                env$model$rhs[map] <- lb_vals[pos]
              } else if (tval == "U") {
                env$model$sense[map] <- "<"
                env$model$rhs[map] <- ub_vals[pos]
              } else if (tval == "E") {
                env$model$sense[map] <- "="
                val <- ifelse(is.na(lb_vals[pos]), ub_vals[pos], lb_vals[pos])
                env$model$rhs[map] <- val
              }
            }

            invisible(lp)
          }
)

#' @rdname solveLp-methods
#' @aliases solveLp,LPproblem_gurobi
setMethod("solveLp", signature(lp = "LPproblem_gurobi"),
          function(lp) {
            env <- lp@ptr

            params <- .gurobi_build_params(lp@method, lp@tol_bnd)
            env$params <- params

            result <- try(gurobi::gurobi(env$model, params = params), silent = TRUE)

            if (inherits(result, "try-error")) {
              env$result <- NULL
              env$status <- "ERROR"
              term <- conditionMessage(attr(result, "condition"))
              return(list(code = 1L,
                          term = term))
            }

            env$result <- result
            env$status <- result$status

            term <- switch(result$status,
                           OPTIMAL = "optimization process was successful",
                           INFEASIBLE = "problem has no feasible solution",
                           UNBOUNDED = "problem has unbounded solution",
                           INF_OR_UNBD = "no primal/dual feasible solution",
                           ITERATION_LIMIT = "iteration limit exceeded",
                           TIME_LIMIT = "time limit exceeded",
                           NUMERIC = "numerical difficulties encountered",
                           paste("solver returned status", result$status))

            code <- if (identical(result$status, "OPTIMAL")) 0L else 1L

            list(code = code,
                 term = term)
          }
)

#' @rdname getObjValue-methods
#' @aliases getObjValue,LPproblem_gurobi
setMethod("getObjValue", signature(lp = "LPproblem_gurobi"),
          function(lp) {
            env <- lp@ptr
            if (is.null(env$result) || is.null(env$result$objval)) {
              return(NA_real_)
            }
            env$result$objval
          }
)

#' @rdname getSolStat-methods
#' @aliases getSolStat,LPproblem_gurobi
setMethod("getSolStat", signature(lp = "LPproblem_gurobi"),
          function(lp) {
            env <- lp@ptr
            status <- env$status

            if (is.null(status)) {
              return(list(code = NA_integer_,
                          term = "solution status is unavailable"))
            }

            code <- switch(status,
                           OPTIMAL = 5L,
                           INFEASIBLE = 4L,
                           UNBOUNDED = 6L,
                           INF_OR_UNBD = 1L,
                           NA_integer_)

            term <- switch(status,
                           OPTIMAL = "solution is optimal",
                           INFEASIBLE = "problem has no feasible solution",
                           UNBOUNDED = "problem has unbounded solution",
                           INF_OR_UNBD = "solution is undefined",
                           paste("solver returned status", status))

            list(code = code,
                 term = term)
          }
)

#' @rdname getColsPrimal-methods
#' @aliases getColsPrimal,LPproblem_gurobi
setMethod("getColsPrimal", signature(lp = "LPproblem_gurobi"),
          function(lp) {
            env <- lp@ptr
            if (is.null(env$result) || is.null(env$result$x)) {
              return(rep(NA_real_, length(env$model$obj)))
            }
            env$result$x
          }
)

#' @rdname getRedCosts-methods
#' @aliases getRedCosts,LPproblem_gurobi
setMethod("getRedCosts", signature(lp = "LPproblem_gurobi"),
          function(lp) {
            env <- lp@ptr
            if (is.null(env$result) || is.null(env$result$rc)) {
              return(rep(NA_real_, length(env$model$obj)))
            }
            env$result$rc
          }
)

#' @rdname addSingleConstraint-methods
#' @aliases addSingleConstraint,LPproblem_gurobi
setMethod("addSingleConstraint", signature(lp = "LPproblem_gurobi"),
          function(lp, coeffs, lb, ub, type) {
            env <- lp@ptr

            coeffs <- as.numeric(coeffs)
            row_mat <- Matrix::sparseMatrix(i = rep(1L, length(coeffs)),
                                            j = seq_along(coeffs),
                                            x = coeffs,
                                            dims = c(1L, length(coeffs)))

            entries <- .gurobi_constraint_entries(row_mat, lb, ub, type)

            if (length(entries$sense) == 0L) {
              env$row_map <- c(env$row_map, list(integer(0)))
              return(invisible(lp))
            }

            merged <- .gurobi_merge_rows(entries$mats, ncol(env$model$A))
            env$model$A <- Matrix::rBind(env$model$A, merged)
            env$model$sense <- c(env$model$sense, entries$sense)
            env$model$rhs <- c(env$model$rhs, entries$rhs)

            start_idx <- length(env$model$sense) - length(entries$sense) + 1L
            end_idx <- length(env$model$sense)
            env$row_map <- c(env$row_map, list(seq.int(start_idx, end_idx)))

            invisible(lp)
          }
)

#' @rdname fvaJob-methods
#' @aliases fvaJob,LPproblem_gurobi
setMethod("fvaJob", signature(lp = "LPproblem_gurobi"),
          function(lp, ind) {
            env <- lp@ptr
            indices <- as.integer(ind)

            base_model <- env$model
            params <- .gurobi_build_params(lp@method, lp@tol_bnd)

            min_vals <- numeric(length(indices))
            max_vals <- numeric(length(indices))

            for (k in seq_along(indices)) {
              obj_vec <- numeric(length(base_model$obj))
              obj_vec[indices[k]] <- 1

              model_max <- base_model
              model_max$modelsense <- "max"
              model_max$obj <- obj_vec

              res_max <- try(gurobi::gurobi(model_max, params = params), silent = TRUE)
              if (inherits(res_max, "try-error") || !identical(res_max$status, "OPTIMAL")) {
                max_vals[k] <- NA_real_
              } else {
                max_vals[k] <- res_max$objval
              }

              model_min <- base_model
              model_min$modelsense <- "min"
              model_min$obj <- obj_vec

              res_min <- try(gurobi::gurobi(model_min, params = params), silent = TRUE)
              if (inherits(res_min, "try-error") || !identical(res_min$status, "OPTIMAL")) {
                min_vals[k] <- NA_real_
              } else {
                min_vals[k] <- res_min$objval
              }
            }

            data.frame(`min.flux` = min_vals,
                       `max.flux` = max_vals)
          }
)

#' @rdname deleteLP-methods
#' @aliases deleteLP,LPproblem_gurobi
setMethod("deleteLP", signature(lp = "LPproblem_gurobi"),
          function(lp) {
            env <- lp@ptr
            env$model <- NULL
            env$result <- NULL
            env$row_map <- list()
            env$status <- NULL
            env$params <- list()
            TRUE
          }
)

