

#' Data transformation.
#' 
#' Transforms the data into the desired data format required by different 
#' estimation packages.
#' 
#' The \code{des} specified should be the full aggregate design. Thus if all 
#' participants responded to the same design, \code{des} will be a repetition of
#' that design matrix.
#' 
#' The responses in \code{y} should be succesive when there are multiple
#' respondents. There can be \code{n.sets} elements for each respondent with
#' discrete values indicating the chosen alternative for each set. Or there can
#' be \code{n.sets * n.alts} elements for each respondent with binary values
#' indicating for each alternative whether it was chosen or not. In the latter
#' case the \code{bin} argument should be \code{TRUE}.
#' 
#' \code{n.sets} indicates the number of sets each repondent responded to. It is
#' assumed that every responded responded to the same number of choice sets.
#' 
#' @param pkg Indicates the desired estimation package. Options are 
#'   \code{bayesm = \link[bayesm]{rhierMnlRwMixture}}, \code{ChoiceModelR = 
#'   \link[ChoiceModelR]{choicemodelr}}, \code{RSGHB = \link[RSGHB]{doHB}}, 
#'   \code{Mixed.Probit = \link[bayesm]{rbprobitGibbs}}, \code{mlogit = 
#'   \link[mlogit]{mlogit}}, and \code{Rchoice = \link[Rchoice]{Rchoice}}).
#' @param des A design matrix in which each row is a profile.
#' @inheritParams Modfed
#' @param y A numeric vector containing binary or discrete responses. See \code{bin}.
#' @param n.resp Numeric value indicating the number of respondents.
#' @param bin Logical value indicating whether the reponse vector contains 
#'   binary data (\code{TRUE}) or discrete data (\code{FALSE}). See \code{y}.
#' @param alt.names A character vector containing the names of the alternatives.
#'   The default = \code{NULL}
#' @return The data ready to be used by the specified package.
#' @examples 
#' idefix.data <- aggregate_design 
#' des <- as.matrix(idefix.data[, 3:8], ncol = 6)
#' y <- idefix.data[, 9]
#' bayesm.data <- Datatrans(pkg = "bayesm", des = des, y = y, 
#' n.alts = 2, n.sets = 8, n.resp = 7, bin = TRUE)
#' Mix.pro.data <- Datatrans(pkg = "Mixed.Probit", des = des, y = y,
#'  n.alts = 2, n.sets = 8, n.resp = 7, bin = TRUE)
#' mlogit.data <- Datatrans(pkg = "mlogit", des = des, y = y,
#'  n.alts = 2, n.sets = 8, n.resp = 7, bin = TRUE)
#' @export
Datatrans <- function(pkg, des, y, n.alts, n.sets, n.resp, bin, alt.names = NULL) {
  ## Errors 
  #des
  if(!is.matrix(des)){stop("'des' should be a matrix")}
  if(!isTRUE(all.equal((n.resp * n.sets * n.alts), nrow(des)))){
    stop("'n.resp' * 'n.sets' * 'n.alts' should equal nrow(des)")
  }
  # y
  if(!is.vector(y)){stop("'y' should be a vector")}
  #binary response data
  if(bin){
    if(!isTRUE(all.equal(length(y) %% nrow(des), 0))){
      stop("length(y) should equal nrow (des), when 'bin = TRUE'")
    }
    y <- BinDis(y = y, n.alts = n.alts, n.resp = n.resp)
    #discrete response data
  } else {
    y.bin <- DisBin(y, n.alts)
    if(!isTRUE(all.equal(nrow(y.bin), nrow(des)))){
      stop("'n.alts * nrow(y) should equal nrow(des), when 'bin = FALSE'")
    }
  }
  #pkg
  if(!pkg %in% c("bayesm", "ChoiceModelR", "RSGHB", "Mixed.Probit",
                 "mlogit", "Rchoice")){
    stop("'pkg' argument is not recognized")
  }
  #alt.names
  if(!is.null(alt.names)){
    if(!isTRUE(all.equal(length(alt.names), n.alts))){
      stop("length(alt.names) should equal 'n.alts'")
    }
  }
  rownames(des) <- NULL
  varnames <- colnames(des)
  n.par <- ncol(des)
  ###########################
  # Bayesm package 
  if(pkg == "bayesm") {
    # Initialize
    lgtdata <- NULL
    ni <- rep(n.sets, n.resp)
    csa <- n.sets * n.alts
    # For every respondent
    for (i in 1:n.resp) {
      # Obtain y
      ychoice <- NULL
      ybeg <- n.sets * (i - 1) + 1
      yend <- n.sets * i
      for (c in 1:n.sets) {
        ychoice[1:n.sets] <- y[ybeg:yend]
      }
      # Transform des into dataframe
      xmat <- NULL
      xbeg <- csa * (i - 1) + 1
      xend <- csa * i
      xmat[[i]] <- des[xbeg:xend, ]
      lgtdata[[i]] <- list(y = ychoice, X = xmat[[i]])
    }
    # The bayesmin function returns a list of 2
    bayesmdata <- list(p = n.alts, lgtdata = lgtdata)
    print("The dataset is ready to be used for bayesm package")
    return(bayesmdata)
  }
  # ChoiceModelR
  if(pkg == "ChoiceModelR") {
    # matrix y to 1 dim 
    y <- as.vector(y)
    y <- matrix(y, length(y))
    set <- rep(1:n.sets, each = n.alts, times = n.resp)
    id <- rep(1:n.resp, each = n.sets * n.alts)
    alt <- rep(c(1:n.alts), n.sets * n.resp)
    initialmat <- t(rbind(id, set, alt))
    xmat <- cbind(initialmat, des)
    # Create choice columns.
    newchoice <- y
    zeromat <- matrix(0, n.sets * n.resp, n.alts - 1)
    choicemat <- cbind(newchoice, zeromat)
    # This is the final y column representing choice.
    choicecol <- matrix(c(t(choicemat)))
    c.data <- cbind(xmat, choicecol)
    colnames(c.data)[ncol(c.data)] <- "Choice" 
    print("The dataset is ready to be used for ChoiceModelR package")
    return(c.data)
  }
  # RSGHB
  if (pkg == "RSGHB") {
    # matrix y to 1 dim 
    y <- as.vector(y)
    y <- matrix(y, length(y))
    id <- rep(1:n.resp, each = n.sets)
    ncs <- rep(1:n.sets, times = n.resp)
    rsghbinitialmat <- t(rbind(id, ncs))
    # Attribute matrix.
    rsghbattrmat <- NULL
    indset <- n.sets * n.resp
    for (cs in 1:indset) {
      beg <- n.alts * (cs - 1) + 1
      end <- n.alts * cs
      xtemp <- NULL
      for(col in 1:n.par) {
        xtemp <- cbind(xtemp, t(des [beg:end, col]))
      }
      rsghbattrmat <- rbind(rsghbattrmat, xtemp)
    }
    RSGHBchoice <- y
    RSGHBdata <- data.frame(cbind(rsghbinitialmat, rsghbattrmat, RSGHBchoice))
    colnames(RSGHBdata)[1] <- "ID"
    colnames(RSGHBdata)[2] <- "Choice Set"
    colnames(RSGHBdata)[ncol(RSGHBdata)] <- "Choice"
    if(!is.null(varnames)){
     colnames(RSGHBdata)[3:(ncol(RSGHBdata)-1)] <- paste(paste("alt", rep(1:n.alts, each = length(varnames)), sep =""), varnames, sep=".")
    }
    print("The dataset is ready to be used for RSGHB package")
    return(RSGHBdata)
  }
  # Mixed Probit estimation.
  if (pkg == "Mixed.Probit") {
    # matrix y to 1 dim 
    y <- as.vector(y)
    y <- matrix(y, length(y))
    
    ynum <- nrow(y)
    yind <- NULL
    for (i in 1:ynum) {
      zerotemp <- matrix(0, n.alts, 1)
      index <- y[i, ]
      zerotemp[index] <- 1
      yind <- rbind(yind, zerotemp)
    }
    y <- array(t(yind), dim = c(n.alts, n.sets, n.resp))
    des <- array(t(des), dim = c(n.par, n.alts, n.sets, n.resp))
    mix.data <- list(y = y, X = des, nlgt = n.resp, nset = n.sets, n.alts = n.alts, nbeta = n.par)
    
    print("The dataset is ready to be used for Mixed Probit Estimation")
    return(mix.data)
  } 
  # mlogit package
  if(pkg == "mlogit") {
    # matrix y to 1 dim 
    if(is.null(alt.names)){
      alt.names <- vector(length = n.alts)
      for (i in 1:n.alts){ alt.names[i] <- paste("alternative", i, sep=".") }
    }
    y.bin <- DisBin(y = y, n.alts = n.alts)
    Choice <- as.logical(y.bin)
    id.var <- rep(1:n.resp, each = n.alts * n.sets)
    mdes <- as.data.frame(cbind(id.var, alt.names, des, Choice))
    mdes$Choice <- as.logical(Choice)
    rownames(mdes) <- NULL
    alt.ctes <- dplyr::select(mdes, dplyr::ends_with(".cte"))
    cte.names <- colnames(alt.ctes)
    logitdata <- mlogit::mlogit.data(mdes, choice = "Choice", shape = "long", 
                             varying = "cte.names", alt.var = "alt.names",
                             id.var = "id.var")
    print("The dataset is ready to be used for mlogit package")
    return(logitdata) 
  }
  # Rchoice
  if(pkg == "Rchoice") {
    # matrix y to 1 dim 
    if(is.null(alt.names)){
      alt.names <- vector(length = n.alts)
      for (i in 1:n.alts){ alt.names[i] <- paste("alternative", i, sep=".") }
    }
    y.id <- DisBin(y = y, n.alts = n.alts)
    colnames(y.id) <- "Choice"
    id.var <- rep(1:n.resp, each = n.alts * n.sets)
    Rdes <- as.data.frame(cbind(id.var, alt.names, des, y.id))
    print("The dataset is ready to be used for Rchoice package")
    return(Rdes) 
  }
}

#' Response generation
#' 
#' Function to generate responses given parameter values and a design matrix,
#' assuming a MNL model.
#' @param par Numeric vector containing parameter values.
#' @inheritParams SeqMOD
#' @param bin Indicates whether the returned value should be a binary vector or 
#'   a discrete value which denotes the chosen alternative.
#' @return Numeric vector indicating the chosen alternatives.
#' @examples 
#' # design: 3 dummy coded attributes, each 3 levels. There are 8 choice sets.
#' des <- example_design
#' set.seed(123)
#' true_par <- rnorm(6)
#' RespondMNL(par=true_par, des = des, n.alts = 2)
#' @export
RespondMNL <- function(par, des, n.alts, bin = TRUE) {
  if (!is.matrix(des)){
    stop("'des' should be a matrix")
  }
  # Error par is not vector
  if (!is.vector(par)) {
    stop("'par' should be a vector")
  }
  # Error n.alts 
  if ((nrow(des) %% n.alts) != 0) {
    stop("number of rows in 'des' is not a multiple of 'n.alts'")
  }
  # Error par
  if (ncol(des) != length(par)) {
    stop("length of 'par' does not match the number of columns in 'des'")
  }
  # Probability
  group <- rep(seq(1, nrow(des) / n.alts, 1), each = n.alts)
  u <- des %*% diag(par)
  u <- .rowSums(u, m = nrow(des), n = length(par))
  p <- exp(u) / rep(rowsum(exp(u), group), each = n.alts)
  # Choice
  n.sets <- nrow(des) / n.alts
  draws <- (0:(n.sets-1)) + (stats::runif(n.sets))
  choice <- findInterval(x = draws, vec = c(0, cumsum(p)))
  Y <- rep(0, length(p))
  Y[choice] <- 1
  # Return
  ifelse(bin, return(Y), return(choice))
}

# Binary to discrete choice data
# 
# Transforms a numeric matrix with binary choice data for each respondent 
# (columns), to a matrix with discrete values representing the chosen 
# alternatives.
BinDis <- function(y, n.alts, n.resp) {
  # y matrix.
  if (!is.matrix(y)) {
    y <- matrix(y, ncol = n.resp)
  }
  # Create choice sets.
  f <- function(x) {
    split(x, ceiling(seq_along(x) / n.alts))
  }
  cs <- apply(y, 2, f)
  # Error n.alts 
  for (i in 1:ncol(y)){
    if((length(unique(lengths(cs[[i]]))) != 1L)){
      stop("length of 'Y' vector does match expected length based on 'n.alts'")
    }
  }
  # Index 1.
  Ones <- function(x) {
    xx <- (x == 1)
    ione <- which(xx, arr.ind = TRUE)
    if(length(ione) > 1) {
      stop("Multiple alternatives are chosen per choice set.")
    }
    # When no choice was chosen
    if(isTRUE(all.equal(length(ione), 0))){
      stop("set(s) detected where no alternative was chosen.")
    }
    return(ione)
  }
  yy <- list()
  for(c in 1:ncol(y)){
    yy[[c]] <- lapply(cs[[c]], Ones)
  }
  # Rbind.
  ynom <- lapply(yy, rbind)
  y.nom <- matrix(unlist(ynom), ncol = ncol(y), byrow = FALSE)
  return(y.nom)
}


# DIScrete to BINary choice data.
DisBin <- function(y, n.alts){
  y <- as.vector(y)
  y.bin <- matrix(y, length(y))
  ynum <- nrow(y.bin)
  y.bin <- NULL
  for (i in 1:ynum) {
    zerotemp <- matrix(0, n.alts, 1)
    index <- y[i]
    zerotemp[index] <- 1
    y.bin <- rbind(y.bin, zerotemp)
  }
  return(y.bin)
}

