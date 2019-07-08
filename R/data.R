#' Discrete choice design.
#' 
#' This discrete choice design is generated using the \code{\link{Modfed}} 
#' function. There are 8 choice sets, each containig 2 alternatives (rows). The
#' alternatives consist of 3 attributes (time, price, comfort) with 3
#' levels each, all of which are dummy coded (columns).
#' 
#' @docType data
#' @usage data(example_design)
#' @format A matrix with 16 rows and 6 columns.
#' @keywords data
"example_design"


#' Discrete choice design.
#' 
#' This discrete choice design is generated using the \code{\link{Modfed}} 
#' function. There are 8 choice sets, each containig 3 alternatives (rows). The
#' alternatives consist of 3 attributes (time, price, comfort) with 3 
#' levels each, all of which are dummy coded (columns). The first two colums are
#' alternative specific constants for alternative 1 and 2.
#' 
#' @docType data
#' @usage data(example_design2)
#' @format A matrix with 24 rows and 8 columns.
#' @keywords data
"example_design2"


#' Discrete choice aggregate design.
#' 
#' The dataset contains fictional data for seven participants, each
#' responding to eight choice sets with two alternatives. Each alternative
#' consists of three attributes, and each attribute contains three levels, which 
#' are dummy coded.
#' 
#' @docType data
#' @usage data(aggregate_design)
#' @format A matrix with 112 rows and 9 variables
#' @keywords data
"aggregate_design"


#' Discrete choice design with no choice option.
#' 
#' This discrete choice design is generated using the \code{\link{Modfed}} 
#' function. There are 8 choice sets, each containig 3 alternatives (rows), of
#' which one is a no choice option. The no choice option consist of an
#' alternative specific constant and zero's for all other attribute levels. There
#' are three attributes (time, price, comfort) with 3 levels each, all of which
#' are dummy coded (columns).
#' 
#' @docType data
#' @usage data(nochoice_design)
#' @format A matrix with 24 rows and 7 variables
#' @keywords data
"nochoice_design"


