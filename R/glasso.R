## Matt Galloway


#' @title Penalized precision matrix estimation
#' 
#' @description Penalized precision matrix estimation using the graphical lasso (glasso) algorithm.
#' Consider the case where \eqn{X_{1}, ..., X_{n}} are iid \eqn{N_{p}(\mu,
#' \Sigma)} and we are tasked with estimating the precision matrix,
#' denoted \eqn{\Omega \equiv \Sigma^{-1}}. This function solves the
#' following optimization problem:
#' \describe{
#' \item{Objective:}{
#' \eqn{\hat{\Omega}_{\lambda} = \arg\min_{\Omega \in S_{+}^{p}}
#' \left\{ Tr\left(S\Omega\right) - \log \det\left(\Omega \right) +
#' \lambda \left\| \Omega \right\|_{1} \right\}}}
#' }
#' where \eqn{\lambda > 0} and we define
#' \eqn{\left\|A \right\|_{1} = \sum_{i, j} \left| A_{ij} \right|}.
#' 
#' @details For details on the implementation of the 'glasso' function, see Tibshirani's website.
#' \url{http://statweb.stanford.edu/~tibs/glasso/}.
#' 
#' @param X option to provide a nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
#' @param nlam number of \code{lam} tuning parameters for penalty term generated from \code{lam.min.ratio} and \code{lam.max} (automatically generated). Defaults to 10.
#' @param lam.min.ratio smallest \code{lam} value provided as a fraction of \code{lam.max}. The function will automatically generate \code{nlam} tuning parameters from \code{lam.min.ratio*lam.max} to \code{lam.max} in log10 scale. \code{lam.max} is calculated to be the smallest \code{lam} such that all off-diagonal entries in \code{Omega} are equal to zero (\code{alpha} = 1). Defaults to 1e-2.
#' @param lam option to provide positive tuning parameters for penalty term. This will cause \code{nlam} and \code{lam.min.ratio} to be disregarded. If a vector of parameters is provided, they should be in increasing order. Defaults to NULL.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores must be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
#' @param tol convergence tolerance. Iterations will stop when the average absolute difference in parameter estimates in less than \code{tol} times multiple. Defaults to 1e-4.
#' @param maxit maximum number of iterations. Defaults to 1e4.
#' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged. This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param crit.cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
#' @param ... additional arguments to pass to \code{glasso}.
#' 
#' @return returns class object \code{CVglasso} which includes:
#' \item{Call}{function call.}
#' \item{Iterations}{number of iterations}
#' \item{Tuning}{optimal tuning parameters (lam and alpha).}
#' \item{Lambdas}{grid of lambda values for CV.}
#' \item{maxit}{maximum number of iterations for outer (blockwise) loop.}
#' \item{Omega}{estimated penalized precision matrix.}
#' \item{Sigma}{estimated covariance matrix from the penalized precision matrix (inverse of Omega).}
#' \item{Path}{array containing the solution path. Solutions will be ordered by ascending lambda values.}
#' \item{MIN.error}{minimum average cross validation error (cv.crit) for optimal parameters.}
#' \item{AVG.error}{average cross validation error (cv.crit) across all folds.}
#' \item{CV.error}{cross validation errors (cv.crit).}
#' 
#' @references
#' \itemize{
#' \item Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 'Sparse inverse covariance estimation with the graphical lasso.' \emph{Biostatistics} 9.3 (2008): 432-441.
#' \item Banerjee, Onureen, Ghauoui, Laurent El, and d'Aspremont, Alexandre. 2008. 'Model Selection through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data.' \emph{Journal of Machine Learning Research} 9: 485-516.
#' \item Tibshirani, Robert. 1996. 'Regression Shrinkage and Selection via the Lasso.' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}. JSTOR: 267-288.
#' \item Meinshausen, Nicolai and Buhlmann, Peter. 2006. 'High-Dimensional Graphs and Variable Selection with the Lasso.' \emph{The Annals of Statistics}. JSTOR: 1436-1462.
#' \item Witten, Daniela M, Friedman, Jerome H, and Simon, Noah. 2011. 'New Insights and Faster computations for the Graphical Lasso.' \emph{Journal of Computation and Graphical Statistics}. Taylor and Francis: 892-900.
#' \item Tibshirani, Robert, Bien, Jacob, Friedman, Jerome, Hastie, Trevor, Simon, Noah, Jonathan, Taylor, and Tibshirani, Ryan J. 'Strong Rules for Discarding Predictors in Lasso-Type Problems.' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}. Wiley Online Library 74 (2): 245-266.
#' \item Ghaoui, Laurent El, Viallon, Vivian, and Rabbani, Tarek. 2010. 'Safe Feature Elimination for the Lasso and Sparse Supervised Learning Problems.' \emph{arXiv preprint arXiv: 1009.4219}.
#' \item Osborne, Michael R, Presnell, Brett, and Turlach, Berwin A. 'On the Lasso and its Dual.' \emph{Journal of Computational and Graphical Statistics}. Taylor and Francis 9 (2): 319-337.
#' \item Rothman, Adam. 2017. 'STAT 8931 notes on an algorithm to compute the Lasso-penalized Gausssian likelihood precision matrix estimator.' 
#' }
#' 
#' @author Matt Galloway \email{gall0441@@umn.edu}
#' 
#' @seealso \code{\link{plot.CVglasso}}
#' 
#' @export
#' 
#' @examples
#' # generate data from a sparse matrix
#' # first compute covariance matrix
#' S = matrix(0.7, nrow = 5, ncol = 5)
#' for (i in 1:5){
#'  for (j in 1:5){
#'    S[i, j] = S[i, j]^abs(i - j)
#'  }
#' }
#'
#' # generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
#' Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' out = eigen(S, symmetric = TRUE)
#' S.sqrt = out$vectors %*% diag(out$values^0.5)
#' S.sqrt = S.sqrt %*% t(out$vectors)
#' X = Z %*% S.sqrt
#'
#' # lasso penalty CV
#' CVglasso(X)

# we define the CVglasso precision matrix estimation function
CVglasso = function(X = NULL, S = NULL, nlam = 10, lam.min.ratio = 0.01, 
    lam = NULL, diagonal = FALSE, path = FALSE, tol = 1e-04, maxit = 10000, 
    adjmaxit = NULL, K = 5, crit.cv = c("loglik", "AIC", "BIC"), 
    start = c("warm", "cold"), cores = 1, trace = c("progress", "print", 
        "none"), ...) {
    
    # checks
    if (is.null(X) && is.null(S)) {
        stop("Must provide entry for X or S!")
    }
    if (!all(lam > 0)) {
        stop("lam must be positive!")
    }
    if (!(all(c(tol, maxit, adjmaxit, K, cores) > 0))) {
        stop("Entry must be positive!")
    }
    if (!(all(sapply(c(tol, maxit, adjmaxit, K, cores, nlam, lam.min.ratio), 
        length) <= 1))) {
        stop("Entry must be single value!")
    }
    if (all(c(maxit, adjmaxit, K, cores)%%1 != 0)) {
        stop("Entry must be an integer!")
    }
    if (cores < 1) {
        stop("Number of cores must be positive!")
    }
    if (cores > 1 && path) {
        cat("\nParallelization not possible when producing solution path. Setting cores = 1...")
        cores = 1
    }
    K = ifelse(path, 1, K)
    if (cores > K) {
        cat("\nNumber of cores exceeds K... setting cores = K")
        cores = K
    }
    if (is.null(adjmaxit)) {
        adjmaxit = maxit
    }
    
    # match values
    crit.cv = match.arg(crit.cv)
    start = match.arg(start)
    trace = match.arg(trace)
    call = match.call()
    MIN.error = AVG.error = CV.error = NULL
    n = ifelse(is.null(X), nrow(S), nrow(X))
    
    # compute sample covariance matrix, if necessary
    if (is.null(S)) {
        S = (nrow(X) - 1)/nrow(X) * cov(X)
    }
    
    Sminus = S
    diag(Sminus) = 0
    
    # compute grid of lam values, if necessary
    if (is.null(lam)) {
        if (!((lam.min.ratio <= 1) && (lam.min.ratio > 0))) {
            cat("\nlam.min.ratio must be in (0, 1]... setting to 1e-2!")
            lam.min.ratio = 0.01
        }
        if (!((nlam > 0) && (nlam%%1 == 0))) {
            cat("\nnlam must be a positive integer... setting to 10!")
            nlam = 10
        }
        
        # calculate lam.max and lam.min
        lam.max = max(abs(Sminus))
        lam.min = lam.min.ratio * lam.max
        
        # calculate grid of lambda values
        lam = 10^seq(log10(lam.min), log10(lam.max), length = nlam)
        
    } else {
        
        # sort lambda values
        lam = sort(lam)
        
    }
    
    # perform cross validation, if necessary
    if ((length(lam) > 1) & (!is.null(X) || path)) {
        
        # run CV in parallel?
        if (cores > 1) {
            
            # execute CVP
            GLASSO = CVP(X = X, lam = lam, diagonal = diagonal, tol = tol, 
                maxit = maxit, adjmaxit = adjmaxit, K = K, crit.cv = crit.cv, 
                start = start, cores = cores, trace = trace, ...)
            MIN.error = GLASSO$min.error
            AVG.error = GLASSO$avg.error
            CV.error = GLASSO$cv.error
            
        } else {
            
            # execute CV_ADMMc
            if (is.null(X)) {
                X = matrix(0)
            }
            GLASSO = CV(X = X, S = S, lam = lam, diagonal = diagonal, 
                path = path, tol = tol, maxit = maxit, adjmaxit = adjmaxit, 
                K = K, crit.cv = crit.cv, start = start, trace = trace, 
                ...)
            MIN.error = GLASSO$min.error
            AVG.error = GLASSO$avg.error
            CV.error = GLASSO$cv.error
            Path = GLASSO$path
            
        }
        
        # print warning if lam on boundary
        if ((GLASSO$lam == lam[1]) && (length(lam) != 1) && !path) {
            cat("\nOptimal tuning parameter on boundary... consider providing a smaller lam value or decreasing lam.min.ratio!")
        }
        
        # specify initial estimate for Sigma
        if (diagonal) {
            
            # simply force init to be positive definite final diagonal
            # elements will be increased by lam
            init = S + GLASSO$lam
            
        } else {
            
            # provide estimate that is pd and dual feasible
            alpha = min(c(GLASSO$lam/max(abs(Sminus)), 1))
            init = (1 - alpha) * S
            diag(init) = diag(S)
            
        }
        
        # compute final estimate at best tuning parameters
        lam_ = GLASSO$lam
        GLASSO = glasso(s = S, rho = lam_, thr = tol, maxit = maxit, 
            penalize.diagonal = diagonal, start = "warm", w.init = init, 
            wi.init = diag(ncol(S)), trace = FALSE, ...)
        GLASSO$lam = lam_
        
        
    } else {
        
        # execute ADMM_sigmac
        if (length(lam) > 1) {
            stop("Must set specify X, set path = TRUE, or provide single value for lam.")
        }
        
        # specify initial estimate for Sigma
        if (diagonal) {
            
            # simply force init to be positive definite final diagonal
            # elements will be increased by lam
            init = S + lam
            
        } else {
            
            # provide estimate that is pd and dual feasible
            alpha = min(c(lam/max(abs(Sminus)), 1))
            init = (1 - alpha) * S
            diag(init) = diag(S)
            
        }
        
        GLASSO = glasso(s = S, rho = lam, thr = tol, maxit = maxit, 
            penalize.diagonal = diagonal, start = "warm", w.init = init, 
            wi.init = diag(ncol(S)), trace = FALSE, ...)
        GLASSO$lam = lam
        
    }
    
    
    # option to penalize diagonal
    if (diagonal) {
        C = 1
    } else {
        C = 1 - diag(ncol(S))
    }
    
    # compute penalized loglik
    loglik = (-n/2) * (sum(GLASSO$wi * S) - determinant(GLASSO$wi, 
        logarithm = TRUE)$modulus[1] + GLASSO$lam * sum(abs(C * GLASSO$wi)))
    
    
    # return values
    tuning = matrix(c(log10(GLASSO$lam), GLASSO$lam), ncol = 2)
    colnames(tuning) = c("log10(lam)", "lam")
    if (!path) {
        Path = NULL
    }
    
    returns = list(Call = call, Iterations = GLASSO$niter, Tuning = tuning, 
        Lambdas = lam, maxit = maxit, Omega = GLASSO$wi, Sigma = GLASSO$w, 
        Path = Path, Loglik = loglik, MIN.error = MIN.error, AVG.error = AVG.error, 
        CV.error = CV.error)
    
    class(returns) = "CVglasso"
    return(returns)
    
}






##-----------------------------------------------------------------------------------



#' @title Print CVglasso object
#' @description Prints CVglasso object and suppresses output if needed.
#' @param x class object CVglasso
#' @param ... additional arguments.
#' @keywords internal
#' @export
print.CVglasso = function(x, ...) {
    
    # print warning if maxit reached
    if (x$maxit <= x$Iterations) {
        print("Maximum iterations reached...!")
    }
    
    # print call
    cat("\n\nCall: ", paste(deparse(x$Call), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
    
    # print iterations
    cat("\nIterations:\n")
    print.default(x$Iterations, quote = FALSE)
    
    # print optimal tuning parameters
    cat("\nTuning parameter:\n")
    print.default(round(x$Tuning, 3), print.gap = 2L, quote = FALSE)
    
    # print loglik
    cat("\nLog-likelihood: ", paste(round(x$Loglik, 5), sep = "\n", 
        collapse = "\n"), "\n", sep = "")
    
    # print Omega if dim <= 10
    if (nrow(x$Omega) <= 10) {
        cat("\nOmega:\n")
        print.default(round(x$Omega, 5))
    } else {
        cat("\n(...output suppressed due to large dimension!)\n")
    }
    
}



##-----------------------------------------------------------------------------------




#' @title Plot CVglasso object
#' @description Produces a plot for the cross validation errors, if available.
#' @param x class object CVglasso
#' @param type produce either 'heatmap' or 'line' graph
#' @param footnote option to print footnote of optimal values. Defaults to TRUE.
#' @param ... additional arguments.
#' @export
#' @examples
#' # generate data from a sparse matrix
#' # first compute covariance matrix
#' S = matrix(0.7, nrow = 5, ncol = 5)
#' for (i in 1:5){
#'  for (j in 1:5){
#'    S[i, j] = S[i, j]^abs(i - j)
#'  }
#' }
#'
#' # generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
#' Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' out = eigen(S, symmetric = TRUE)
#' S.sqrt = out$vectors %*% diag(out$values^0.5)
#' S.sqrt = S.sqrt %*% t(out$vectors)
#' X = Z %*% S.sqrt
#' 
#' # produce line graph for CVglasso
#' plot(CVglasso(X))
#' 
#' # produce CV heat map for CVglasso
#' plot(CVglasso(X), type = 'heatmap')

plot.CVglasso = function(x, type = c("line", "heatmap"), footnote = TRUE, 
    ...) {
    
    # check
    type = match.arg(type)
    Means = NULL
    if (is.null(x$CV.error)) {
        stop("No cross validation errors to plot!")
    }
    
    if (type == "line") {
        
        # gather values to plot
        cv = cbind(expand.grid(lam = x$Lambdas, alpha = 0), Errors = as.data.frame.table(x$CV.error)$Freq)
        
        # produce line graph
        graph = ggplot(summarise(group_by(cv, lam), Means = mean(Errors)), 
            aes(log10(lam), Means)) + geom_jitter(width = 0.2, color = "navy blue") + 
            theme_minimal() + geom_line(color = "red") + labs(title = "Cross-Validation Errors", 
            y = "Error") + geom_vline(xintercept = x$Tuning[1], linetype = "dotted")
        
    } else {
        
        # augment values for heat map (helps visually)
        lam = x$Lambdas
        cv = expand.grid(lam = lam, alpha = 0)
        Errors = 1/(c(x$AVG.error) + abs(min(x$AVG.error)) + 1)
        cv = cbind(cv, Errors)
        
        # design color palette
        bluetowhite <- c("#000E29", "white")
        
        # produce ggplot heat map
        graph = ggplot(cv, aes(alpha, log10(lam))) + geom_raster(aes(fill = Errors)) + 
            scale_fill_gradientn(colours = colorRampPalette(bluetowhite)(2), 
                guide = "none") + theme_minimal() + labs(title = "Heatmap of Cross-Validation Errors") + 
            theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
                axis.ticks.x = element_blank())
        
    }
    
    if (footnote) {
        
        # produce with footnote
        graph + labs(caption = paste("**Optimal: log10(lam) = ", 
            round(x$Tuning[1], 3), sep = ""))
        
    } else {
        
        # produce without footnote
        graph
        
    }
}
