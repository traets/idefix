
#' idefix: efficient designs for discrete choice experiments.
#'
#' Generates efficient designs for discrete choice experiments based on the
#' MNL model, and individually adapted designs for the mixed
#' multinomial logit model. The (adaptive) designs can be presented on screen and
#' choice data can be gathered using a shiny application.
#'
#'
#' \itemize{
#' \item To generate efficient designs using a modified federov algorithm, please consult the \link[idefix]{Modfed} documentation.
#' \item To generate adaptive designs, please consult the \link[idefix]{SeqDB} documentation.
#' \item To generate a discrete choice survey on screen, please consult the \link[idefix]{SurveyApp} documentation. 
#' }
#'
#' @useDynLib idefix
#' @importFrom Rcpp sourceCpp
"_PACKAGE"