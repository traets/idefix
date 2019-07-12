
#' idefix: efficient designs for discrete choice experiments.
#'
#' Generates efficient designs for discrete choice experiments based on the
#' Multinomial Logit (MNL) model, and individually adapted designs for the Mixed
#' Multinomial Logit model. The (adaptive) designs can be presented on screen and
#' choice data can be gathered using a shiny application.
#'
#'
#' \itemize{
#' \item To generate efficient designs using the Modified Federov algorithm, please consult the \link[idefix]{Modfed} documentation.
#' \item To generate efficient designs using the Coordinate Exchange algorithm, 
#' please consult the \link[idefix]{CEA} documentation.
#' \item To generate adaptive designs using the Modified Fedorov algorithm, 
#' please consult the \link[idefix]{SeqMOD} documentation.
#' \item To generate adaptive designs using the Coordinate Exchange algorithm, 
#' please consult the \link[idefix]{SeqCEA} documentation.
#' \item To generate a discrete choice survey on screen, please consult the \link[idefix]{SurveyApp} documentation. 
#' }
#'
#' @useDynLib idefix
#' @importFrom Rcpp sourceCpp
"_PACKAGE"