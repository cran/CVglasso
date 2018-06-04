## Matt Galloway


#' @title Parallel Cross Validation
#' @description Parallel implementation of cross validation.
#'
#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-2, 2, 0.2)}.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
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
#' @return returns list of returns which includes:
#' \item{lam}{optimal tuning parameter.}
#' \item{min.error}{minimum average cross validation error (cv.crit) for optimal parameters.}
#' \item{avg.error}{average cross validation error (cv.crit) across all folds.}
#' \item{cv.error}{cross validation errors (cv.crit).}
#' 
#' @keywords internal

# we define the CVP function
CVP = function(X = NULL, lam = 10^seq(-2, 2, 0.2), diagonal = FALSE, 
    tol = 1e-04, maxit = 10000, adjmaxit = NULL, K = 5, crit.cv = c("loglik", 
        "AIC", "BIC"), start = c("warm", "cold"), cores = 1, trace = c("progress", 
        "print", "none"), ...) {
    
    # match values
    crit.cv = match.arg(crit.cv)
    start = match.arg(start)
    trace = match.arg(trace)
    lam = sort(lam)
    
    # make cluster and register cluster
    num_cores = detectCores()
    if (cores > num_cores) {
        cat("\nOnly detected", paste(num_cores, "cores...", sep = " "))
    }
    if (cores > K) {
        cat("\nNumber of cores exceeds K... setting cores = K")
        cores = K
    }
    
    cluster = makeCluster(cores)
    registerDoParallel(cluster)
    
    # use cluster for each fold in CV
    n = nrow(X)
    ind = sample(n)
    k = NULL
    CV = foreach(k = 1:K, .packages = "CVglasso", .combine = "cbind", 
        .inorder = FALSE) %dopar% {
        
        # set progress bar
        if (trace == "progress") {
            progress = txtProgressBar(max = length(lam), style = 3)
        }
        
        # training set
        leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * n/K)]
        X.train = X[-leave.out, , drop = FALSE]
        X_bar = apply(X.train, 2, mean)
        X.train = scale(X.train, center = X_bar, scale = FALSE)
        
        # validation set
        X.valid = X[leave.out, , drop = FALSE]
        X.valid = scale(X.valid, center = X_bar, scale = FALSE)
        
        # sample covariances
        S.train = crossprod(X.train)/(dim(X.train)[1])
        S.valid = crossprod(X.valid)/(dim(X.valid)[1])
        
        # initial estimates
        init = S.train
        initOmega = diag(ncol(S.train))
        CV_error = array(0, length(lam))
        
        # initial sigma
        if (!diagonal) {
            
            # provide estimate that is pd and dual feasible
            Sminus = S.train
            diag(Sminus) = 0
            alpha = min(c(lam[1]/max(abs(Sminus)), 1))
            init = (1 - alpha) * S.train
            diag(init) = diag(S.train)
            
        }
        
        # loop over all tuning parameters
        for (i in 1:length(lam)) {
            
            # set temporary tuning parameter
            lam_ = lam[i]
            
            # update diagonal elements of init, if necessary
            if (diagonal) {
                diag(init) = diag(S.train) + lam_
            }
            
            # compute the penalized likelihood precision matrix estimator
            GLASSO = glasso(s = S.train, rho = lam_, thr = tol, maxit = maxit, 
                penalize.diagonal = diagonal, start = "warm", w.init = init, 
                wi.init = initOmega, trace = FALSE, ...)
            
            if (start == "warm") {
                
                # option to save initial values for warm starts
                init = GLASSO$w
                initOmega = GLASSO$wi
                maxit = adjmaxit
                
            }
            
            # compute the observed negative validation loglikelihood (close
            # enoug)
            CV_error[i] = (nrow(X)/2) * (sum(GLASSO$wi * S.valid) - 
                determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])
            
            # update for crit.cv, if necessary
            if (crit.cv == "AIC") {
                CV_error[i] = CV_error[i] + sum(GLASSO$wi != 0)
            }
            if (crit.cv == "BIC") {
                CV_error[i] = CV_error[i] + sum(GLASSO$wi != 0) * 
                  log(nrow(X))/2
            }
            
            # update progress bar
            if (trace == "progress") {
                setTxtProgressBar(progress, i)
                
                # if not quiet, then print progress lambda
            } else if (trace == "print") {
                cat("\nFinished lam = ", paste(lam_, sep = ""))
            }
        }
        
        # return CV errors
        return(CV_error)
        
    }
    
    # determine optimal tuning parameters
    AVG = apply(CV, 1, mean)
    best_lam = lam[which.min(AVG)]
    error = min(AVG)
    
    # stop cluster
    stopCluster(cluster)
    
    # return best lam and alpha values
    return(list(lam = best_lam, min.error = error, avg.error = AVG, 
        cv.error = CV))
    
}




