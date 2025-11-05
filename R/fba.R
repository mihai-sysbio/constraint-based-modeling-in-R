#' Optimize a model with flux balance analysis
#'
#' Provides a user-facing name that mirrors
#' [cobrapy](https://opencobra.github.io/cobrapy/)'s
#' [`Model.optimize()`](https://opencobra.github.io/cobrapy/build/html/model.html#cobra.core.model.Model.optimize)
#' entry point while keeping the original function available for
#' backward compatibility.
#'
#' @param model Model of class \link{ModelOrg}
#'
#' @returns A \link{FluxPrediction-class} object with the predicted reaction
#' fluxes, reduced costs, objective value, and solver status information.
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- read_sbml_model(fpath)
#'
#' # aerobic growth
#' res_aero <- optimize_model(mod)
#' cat(" Growth rate:       ", res_aero@obj,"\n",
#'     "Acetate production:", res_aero@fluxes[mod@react_id == "EX_ac_e"],"\n")
#'
#' mod <- changeBounds(mod, react = "EX_o2_e", lb = 0) # before: -1000
#' res_anaero <- optimize_model(mod)
#' cat(" Growth rate:       ", res_anaero@obj,"\n",
#'     "Acetate production:", res_anaero@fluxes[mod@react_id == "EX_ac_e"],"\n")
#'
#' @family Flux prediction algorithms
#' @export
optimize_model <- function(model) {

  #----------------------------------------------------------------------------#
  # Initializing and defining LP problem                                       #
  #----------------------------------------------------------------------------#
  LPprob <- new(paste0("LPproblem_",COBRAR_SETTINGS("SOLVER")),
                name = paste0("LP_", model@mod_id),
                method = COBRAR_SETTINGS("METHOD"))

  loadLPprob(LPprob,
             nCols = react_num(model),
             nRows = met_num(model)+constraint_num(model),
             mat   = rbind(model@S, model@constraints@coeff),
             ub    = ifelse(abs(model@uppbnd)>COBRAR_SETTINGS("MAXIMUM"),
                            sign(model@uppbnd)*COBRAR_SETTINGS("MAXIMUM"),
                            model@uppbnd),
             lb    = ifelse(abs(model@lowbnd)>COBRAR_SETTINGS("MAXIMUM"),
                            sign(model@lowbnd)*COBRAR_SETTINGS("MAXIMUM"),
                            model@lowbnd),
             obj   = model@obj_coef,
             rlb   = c(rep(0, met_num(model)),
                       model@constraints@lb),
             rtype = c(rep("E", met_num(model)),
                       model@constraints@rtype),
             lpdir = substr(model@obj_dir,1,3),
             rub   = c(rep(NA, met_num(model)),
                       model@constraints@ub),
             ctype = NULL
  )

  #----------------------------------------------------------------------------#
  # Optimizing problem                                                         #
  #----------------------------------------------------------------------------#
  lp_ok   <- solveLp(LPprob)
  lp_stat <- getSolStat(LPprob)

  #----------------------------------------------------------------------------#
  # Retrieve predictions                                                       #
  #----------------------------------------------------------------------------#
  objRes <- getObjValue(LPprob)
  lp_fluxes <- getColsPrimal(LPprob)

  redCosts <- getRedCosts(LPprob)

  #----------------------------------------------------------------------------#
  # Delete LP-Problem and free associated memory                               #
  #----------------------------------------------------------------------------#
  deleteLP(LPprob)

  return(new("FluxPrediction",
             algorithm = "FBA",
             ok = lp_ok$code,
             ok_term = lp_ok$term,
             stat = lp_stat$code,
             stat_term = lp_stat$term,
             obj = objRes,
             obj_sec = NA_real_,
             fluxes = lp_fluxes,
             redCosts = redCosts))
}

#' @rdname optimize_model
#' @export
fba <- function(model) {
  .Deprecated("optimize_model", package = "cobrar")
  optimize_model(model)
}

