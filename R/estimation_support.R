

#' Data transformation.
#'
#' Transforms the data into the neccesary format in order to use estimation functions from different packages. 
#' @param pkg Indicates the required package for estimation (\code{1 = 
#'   \link[bayesm]{rhierMnlRwMixture}}, \code{2 =
#'   \link[ChoiceModelR]{choicemodelr}}, \code{3 = \link[RSGHB]{doHB}} and
#'   \code{4 = \link[bayesm]{rbprobitGibbs}}).
#' @param des A design matrix in which each row is a profile. Can be generated 
#'   with \code{\link{Modfed}}
#' @inheritParams Modfed
#' @param y A numeric matrix. Each column represents a respondent, there are 
#'   \code{n.sets} rows (discrete choice data), or \code{n.sets * n.alts} rows 
#'   (binary data). see \code{bin}
#' @param n.resp Numeric value indicating the number of respondents.
#' @param n.par Numeric value indicating the number of model parameters that 
#'   needs to be estimated.
#' @param no.choice Logical value indicating whether a no choice response could 
#'   be observed. This would be a \code{0} for each alternative.
#' @param bin Logical value indicating whether the reponse matrix contains 
#'   binary data (\code{TRUE}) or not (\code{FALSE}).
#' @return The data ready to be used by the specified package.
#' @examples 
#' # 3 Attributes, 2 are dummy coded and 1 continuous.
#' cs <- Profiles(lvls = c(2, 3, 2), coding = c("D", "C", "D"), c.lvls = list(c(2,4,6)))
#' p <- c(0.8, 0.2, -0.3) # parameter vector
#' # Generate design
#' des <- Modfed(cand.set = cs, n.sets = 8, n.alts = 2, alt.cte = c(0,0), par.draws = p)$des
#' # Generate responses
#' y <- RespondMNL(par = p, des = des, n.alts = 2)
#' y <- matrix(y, 16)
#' #  data 
#' Datatrans(pkg = 4, des = des, y = y, n.alts = 2, n.sets = 8, n.resp = 1,
#'           n.par = 3, no.choice = FALSE, bin = TRUE)
#' @export
Datatrans <- function(pkg, des, y, n.alts, n.sets, n.resp, n.par, no.choice, bin) {
  rownames(des) <- NULL
  # Transform data if binary.
  if (bin) {
    y <- BinDis(y = y, n.alts = n.alts, no.choice = no.choice)
  }
  # Transform y and des into matrix form.
  y <- as.matrix(y)
  des <- as.matrix(des)
  # Bayesm package. 
  if(pkg == 1) {
    bayesmin <- function(des, y, n.alts, n.sets, n.resp) {
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
      return(bayesmdata)
    }
    print("The dataset is ready to be used for bayesm package")
    return(bayesmin(des, y, n.alts, n.sets, n.resp))
  }
  # ChoiceModelR.
  else if(pkg == 2) {
    Cmodrin <- function(des, y, n.alts, n.sets, n.resp) {
      set <- rep(1:n.sets, each = n.alts, times = n.resp)
      id <- rep(1:n.resp, each = n.sets * n.alts)
      alt <- rep(c(1:n.alts), n.sets * n.resp)
      initialmat <- t(rbind(id, set, alt))
      xmat <- cbind(initialmat, des)
      # Make choice columns.
      newchoice <- y
      zeromat <- matrix(0, n.sets * n.resp, n.alts - 1)
      choicemat <- cbind(newchoice, zeromat)
      # This is the final y column representing choice.
      choicecol <- matrix(c(t(choicemat)))
      choicemodelrdata <- cbind(xmat, choicecol)
      return(choicemodelrdata)
    }
    c.data <- Cmodrin(des, y, n.alts, n.sets, n.resp)
    print("The dataset is ready to be used for ChoiceModelR package")
    return(c.data)
  }
  # RSGHB.
  else if (pkg == 3) {
    Rsg <- function(des, y, n.alts, n.sets, n.resp, n.par) {
      n.par <- ncol(des)
      rsghbid <- rep(1:n.resp, each = n.sets)
      ncs <- rep(1:n.sets, times = n.resp)
      rsghbinitialmat <- t(rbind(rsghbid, ncs))
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
      colnames(RSGHBdata)[[1]] <- "ID"
      colnames(RSGHBdata)[[2]] <- "Choice Set"
      cy <- ncol(RSGHBdata)
      colnames(RSGHBdata)[[cy]] <- "Choice"
      return(RSGHBdata)
    }
    rsg.data <- Rsg(des, y, n.alts, n.sets, n.resp, n.par)
    print("The dataset is ready to be used for RSGHB package")
    return(rsg.data)
  }
  # Mixed Probit estimation.
  else if (pkg == 4) {
    mxpin <- function(des, y, n.alts, n.sets, n.resp, n.par) {
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
      return(Data = list(y = y, X = des, nlgt = n.resp, nset = n.sets, n.alts = n.alts, nbeta = n.par))
    }
    mix.data <- mxpin(des, y, n.alts, n.sets, n.resp, n.par)
    print("The dataset is ready to be used for Mixed Probit Estimation")
    return(mix.data)
  } else {
    return(print("please specify: 1-bayesm, 2-ChoiceModelR, 3-RSGHB, 4-Mixed Probit") )
  }
}

#' Response generation
#' 
#' Function to generate responses given parameter values and a design matrix, 
#' assuming a MNL model.
#' @param par Numeric vector containing parameter values.
#' @inheritParams SeqDB
#' @param bin Indicates whether the returned value should be a binary vector or 
#'   a discrete value which denotes the chosen alternative.
#' @return Numeric vector indicating the chosen alternatives.
#' @examples 
#' # 3 Attributes, 2 are dummy coded and 1 continuous.
#' cs <- Profiles(lvls = c(2, 3, 2), coding = c("D", "C", "D"), c.lvls = list(c(2,4,6)))
#' p <- c(0.8, 0.2, -0.3) # parameter vector
#' # Generate design
#' des <- Modfed(cand.set = cs, n.sets = 8, n.alts = 2, alt.cte = c(0,0), par.draws = p)$des
#' # Generate responses
#' y <- RespondMNL(par = p, des = des, n.alts = 2)
#' @export
RespondMNL <- function(par, des, n.alts, bin = TRUE) {
  # Error par is not vector
  if (!is.vector(par)) {
    stop('par should be a vector.')
  }
  # Error n.alts 
  if ((nrow(des) %% n.alts) != 0) {
    stop('number of rows in des is not a multiple of n.alts.')
  }
  # Error par
  if (ncol(des) != length(par)) {
    stop("length of par vector does not match the number of parameters in the design.")
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



