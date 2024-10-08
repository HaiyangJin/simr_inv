#' Estimate power by simulation.
#'
#' Perform a power analysis for a mixed model.
#'
#' @param fit a fitted model object (see \code{\link{doFit}}).
#' @param test specify the test to perform. By default, the first fixed effect in \code{fit} will be tested.
#'     (see: \link{tests}).
#' @param sim an object to simulate from. By default this is the same as \code{fit} (see \code{\link{doSim}}).
#' @param seed specify a random number generator seed, for reproducible results.
#' @param fitOpts extra arguments for \code{\link{doFit}}.
#' @param testOpts extra arguments for \code{\link{doTest}}.
#' @param simOpts extra arguments for \code{\link{doSim}}.
#' @param ... any additional arguments are passed on to \code{\link{simrOptions}}. Common options include:
#' \describe{
#'   \item{\code{nsim}:}{the number of simulations to run (default is \code{1000}).}
#'   \item{\code{alpha}:}{the significance level for the statistical test (default is \code{0.05}).}
#'   \item{\code{progress}:}{use progress bars during calculations (default is \code{TRUE}).}
#'   }
#' @examples
#' fm1 <- lmer(y ~ x + (1|g), data=simdata)
#' powerSim(fm1, nsim=10)
#'
#' @seealso \code{\link{print.powerSim}}, \code{\link{summary.powerSim}}, \code{\link{confint.powerSim}}
#' @export
powerSim <- function(

    fit,
    test = fixed(getDefaultXname(fit)),
    sim = fit,

    fitOpts = list(),
    testOpts = list(),
    simOpts = list(),

    seed,

    ...

    ) {

    opts <- simrOptions(...)
    on.exit(simrOptions(opts))

    nsim <- getSimrOption("nsim")
    alpha <- getSimrOption("alpha")
    nrow <- NA

    # START TIMING
    start <- proc.time()

    # setup
    if(!missing(seed)) set.seed(seed)

    # summarise the fitted models
    test <- wrapTest(test)
    #p <- maybe_laply(z, test, .text="Testing")

    # custom function to save more info from simulated fit #####################
    custom_inv_test <- function(fit, pval){
      # estimated marginal means
      emm_fit <- emmeans::emmeans(fit, ~ Congruency + Alignment + CA)
      
      # composite effect
      emm_cf <- emmeans::contrast(emm_fit,
                                  interaction = "pairwise",
                                  infer=TRUE) %>% 
        as.data.frame()
      pval$cf_p <- emm_cf$p.value
      pval$cf_lb <- emm_cf$asymp.LCL
      pval$cf_ub <- emm_cf$asymp.UCL
      
      # facilitation and interference
      emm_fi <- emmeans::contrast(emm_fit, 
                                  interaction = "pairwise", 
                                  by = "Congruency", 
                                  infer = TRUE, 
                                  adjust = "none")[1:2] %>% 
        as.data.frame()
      pval$fac_p <- emm_fi |> dplyr::filter(Congruency=="con") |> dplyr::pull(p.value)
      pval$fac_lb <- emm_fi |> dplyr::filter(Congruency=="con") |> dplyr::pull(asymp.LCL)
      pval$fac_ub <- emm_fi |> dplyr::filter(Congruency=="con") |> dplyr::pull(asymp.UCL)
      pval$int_p <- emm_fi |> dplyr::filter(Congruency=="inc") |> dplyr::pull(p.value)
      pval$int_lb <- emm_fi |> dplyr::filter(Congruency=="inc") |> dplyr::pull(asymp.LCL)
      pval$int_ub <- emm_fi |> dplyr::filter(Congruency=="inc") |> dplyr::pull(asymp.UCL)
      
      return(pval)
    }
    # custom function ends #####################################################
    
    f <- function() {

        # y <- doSim(sim, [opts])
        tag(y <- do.call(doSim, c(list(sim), simOpts)), tag="Simulating")

        # how many rows?
        ss <- fitOpts$subset
        nrow <<- length(if(is.null(ss)) y else y[ss])

        # fit <- doFit(y, fit, [opts])
        tag(z <- do.call(doFit, c(list(y, fit), fitOpts)), tag="Fitting")

        # doTest(fit, test, [opts])
        tag(pval <- do.call(doTest, c(list(z, test), testOpts)), tag="Testing")
        
        # custom added
        pval <- custom_inv_test(z, pval)

        return(pval)
    }

    p <- f()

    # END TIMING
    timing <- proc.time() - start

    # structure the return value
    rval <- list()

    rval $ x <- sum(p$value < alpha, na.rm=TRUE)
    rval $ n <- nsim

    #rval $ xname <- xname
    #rval $ effect <- fixef(sim)[xname] # can't guarantee this is available?

    rval $ text <- attr(test, "text")(fit, sim)
    rval $ description <- attr(test, "description")(fit, sim)

    rval $ pval <- p$value

    rval $ alpha <- alpha
    rval $ nrow <- nrow

    rval $ warnings <- p$warnings
    rval $ errors <- p$errors

    rval $ timing <- timing
    rval $ simrTag <- observedPowerWarning(sim)
    
    # custom added #############################################################
    rval $ cf_p <- p$cf_p
    rval $ cf_lb <- p$cf_lb
    rval $ cf_ub <- p$cf_ub
    rval $ fac_p <- p$fac_p
    rval $ fac_lb <- p$fac_lb
    rval $ fac_ub <- p$fac_ub
    rval $ int_p <- p$int_p
    rval $ int_lb <- p$int_lb
    rval $ int_ub <- p$int_ub
    # custom added ends ########################################################

    class(rval) <- "powerSim"

    .simrLastResult $ lastResult <- rval

    return(rval)
}

#' @export
plot.powerSim <- function(x, ...) stop("Not yet implemented.")
