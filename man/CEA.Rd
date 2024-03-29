% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CEA.R
\name{CEA}
\alias{CEA}
\title{Coordinate Exchange algorithm for MNL models.}
\usage{
CEA(
  lvls,
  coding,
  c.lvls = NULL,
  n.sets,
  n.alts,
  par.draws,
  alt.cte = NULL,
  no.choice = FALSE,
  start.des = NULL,
  parallel = TRUE,
  max.iter = Inf,
  n.start = 12,
  best = TRUE
)
}
\arguments{
\item{lvls}{A numeric vector which contains for each attribute the number
of levels.}

\item{coding}{Type of coding that needs to be used for each attribute.}

\item{c.lvls}{A list containing numeric vectors with the attribute levels for
each continuous attribute. The default is \code{NULL}.}

\item{n.sets}{Numeric value indicating the number of choice sets.}

\item{n.alts}{Numeric value indicating the number of alternatives per choice
set.}

\item{par.draws}{A matrix or a list, depending on \code{alt.cte}.}

\item{alt.cte}{A binary vector indicating for each alternative whether an
alternative specific constant is desired. The default is \code{NULL}.}

\item{no.choice}{A logical value indicating whether a no choice alternative
should be added to each choice set. The default is \code{FALSE}.}

\item{start.des}{A list containing one or more matrices corresponding to
initial start design(s). The default is \code{NULL}.}

\item{parallel}{Logical value indicating whether computations should be done
over multiple cores. The default is \code{TRUE}.}

\item{max.iter}{A numeric value indicating the maximum number allowed
iterations. The default is \code{Inf}.}

\item{n.start}{A numeric value indicating the number of random start designs
to use. The default is 12.}

\item{best}{A logical value indicating whether only the best design should be
returned. The default is \code{TRUE}.}
}
\value{
If \code{best = TRUE} the design with the lowest D(B)-error is returned.
If \code{best = FALSE}, the results of all (provided) start designs are
returned. \item{design}{A
numeric matrix wich contains an efficient design.} \item{error}{Numeric
value indicating the D(B)-error of the design.} \item{inf.error}{Numeric
value indicating the percentage of draws for which the D-error was
\code{Inf}.} \item{probs}{Numeric matrix containing the probabilities of
each alternative in each choice set. If a sample matrix was provided in
\code{par.draws}, this is the average over all draws.}
}
\description{
The algorithm improves an initial start design by considering changes on an
attribute-by-attribute basis. By doing this, it tries to minimize the
D(B)-error based on a multinomial logit model. This routine is repeated for
multiple starting designs.
}
\details{
Each iteration will loop through all profiles from the initial design,
evaluating the change in D(B)-error for every level in each attribute.
The algorithm stops when an iteration occured without replacing a profile or
when \code{max.iter} is reached.

By specifying a numeric vector in \code{par.draws}, the D-error will be
calculated and the design will be optimised locally. By specifying a matrix,
in which each row is a draw from a multivariate distribution, the DB-error
will be calculated, and the design will be optimised globally. Whenever there
are alternative specific constants, \code{par.draws} should be a list
containing two matrices: The first matrix containing the parameter draws for
the alternative specific constant parameters. The second matrix containing
the draws for the rest of the parameters.

The DB-error is calculated by taking the mean over D-errors. It could be that
for some draws the design results in an infinite D-error. The percentage of
draws for which this was true for the final design can be found in the output
\code{inf.error}.

Alternative specific constants can be specified in \code{alt.cte}. The length
of this binary vector should equal \code{n.alts}, were \code{0} indicates the
absence of an alternative specific constant and \code{1} the opposite.

\code{start.des} is a list with one or several matrices corresponding to
initial start design(s). In each matrix each row is a profile. The number of rows
equals \code{n.sets * n.alts}, and the
number of columns equals the number of columns of the design matrix + the
number of non-zero elements in \code{alt.cte}. Consider that for a
categorical attribute with \emph{p} levels, there are \emph{p - 1} columns in the design
matrix, whereas for a continuous attribute there is only one column. If
\code{start.des = NULL}, \code{n.start} random initial designs will be
generated. If start designs are provided, \code{n.start} is ignored.

If \code{no.choice} is \code{TRUE}, in each choice set an alternative with
one alternative specific constant is added. The return value of the
D(B)-error is however based on the design without the no choice option.

When \code{parallel} is \code{TRUE}, \code{\link[parallel]{detectCores}} will
be used to decide upon the number of available cores. That number minus 1
cores will be used to search for efficient designs. The computation time will
decrease significantly when \code{parallel = TRUE}.
}
\examples{
\donttest{
# DB-efficient designs
# 3 Attributes, all dummy coded. 1 alternative specific constant = 7 parameters
mu <- c(1.2, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) # Prior parameter vector
v <- diag(length(mu)) # Prior variance.
set.seed(123) 
pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
p.d <- list(matrix(pd[,1], ncol = 1), pd[,2:7])
CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), par.draws = p.d,
n.alts = 2, n.sets = 8, parallel = FALSE, alt.cte = c(0, 1))

# DB-efficient design with categorical and continuous factors
# 2 categorical attributes with 4 and 2 levels (effect coded) and 1 
# continuous attribute (= 5 parameters)
mu <- c(0.5, 0.8, 0.2, 0.4, 0.3) 
v <- diag(length(mu)) # Prior variance.
set.seed(123) 
pd <- MASS::mvrnorm(n = 3, mu = mu, Sigma = v) # 10 draws.
CEA(lvls = c(4, 2, 3), coding = c("E", "E", "C"), par.draws = pd,
c.lvls = list(c(2, 4, 6)), n.alts = 2, n.sets = 6, parallel = FALSE)

# DB-efficient design with start design provided.  
# 3 Attributes with 3 levels, all dummy coded (= 6 parameters).
mu <- c(0.8, 0.2, -0.3, -0.2, 0.7, 0.4) 
v <- diag(length(mu)) # Prior variance.
sd <- list(example_design)
set.seed(123)
ps <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), par.draws = ps,
n.alts = 2, n.sets = 8, parallel = FALSE, start.des = sd)
}
}
