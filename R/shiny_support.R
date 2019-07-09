#'Shiny application to generate a discrete choice survey.
#'
#'This function starts a shiny application which puts choice sets on screen and 
#'saves the responses. The complete choice design can be provided in advance, or
#'can be generated sequentially adaptively, or can be a combination of both.
#'
#'A pregenerated design can be specified in \code{des}. This should be a matrix 
#'in which each row is a profile. This can be generated with \code{\link{Modfed}}
#'or  \code{\link{CEA}}, but it is not necessary.
#'
#'If \code{n.total} = \code{nrow(des)} / \code{length(alts)}, the specified 
#'design will be put on screen, one set after the other, and the responses will 
#'be saved. If \code{n.total} > (\code{nrow(des)} / \code{length(alts)}), first 
#'the specified design will be shown and afterwards the remaining sets will be 
#'generated adaptively. If \code{des} = \code{NULL}, \code{n.total} sets will be
#'generated adaptively. See \code{\link{SeqMOD}} for more information on adaptive
#'choice sets.
#'
#'Whenever adaptive sets will be generated, \code{prior.mean}, 
#'\code{prior.covar}, \code{cand.set} and \code{n.draws}, should be specified. 
#'These arguments are necessary for the underlying importance sampling algorithm
#'to update the prior preference distribution. \code{lower} and \code{upper} can
#'be used to specify lower and upper truncation points. See
#'\code{\link{ImpsampMNL}} for more details.
#'
#'The names specified in \code{alts} will be used to label the choice 
#'alternatives. The names specified in \code{atts} will be used to name the 
#'attributes in the choice sets. The values of \code{lvl.names} will be used to 
#'create the values in the choice sets. See \code{\link{Decode}} for more 
#'details. 
#'
#'The text specified in \code{buttons.text} will be displayed above the buttons 
#'to indicate the preferred choice (for example: "indicate your preferred 
#'choice"). The text specified in \code{intro.text} will be displayed before the
#'choice sets. This will generally be a description of the survey and some 
#'instructions. The text specified in \code{end.text} will be displayed after 
#'the survey. This will generally be a thanking note and some further 
#'instructions.
#'
#'A no choice alternative is coded as an alternative with 1 alternative specific
#'constant and zero's for all other attribute levels. If a no choice alternative
#'is present in \code{des}, or is desired when generating adaptive choice sets,
#'\code{no.choice} should be specified. This should be done with an integer,
#'indicating which alternative is the no choice option. This alternative will
#'not be presented on screen, but the option to select "no choice" will be. The
#'\code{alt.cte} argument should be specified accordingly, namely with a
#'\code{1} on the location of the \code{no.choice} option. See examples for
#' illustration.
#'
#'When \code{parallel} is \code{TRUE}, \code{\link[parallel]{detectCores}} will
#'be used to decide upon the number of available cores. That number minus 1 
#'cores will be used to search for the optimal adaptive choice set. For small problems 
#'(6 parameters), \code{parallel = TRUE} can be slower. For larger problems the
#'computation time will decrease significantly.
#'
#'When \code{reduce = TRUE}, the set of all potential choice sets will be
#'reduced to choice sets that have a unique information matrix. If no
#'alternative specific constants are used, \code{reduce} should always be
#'\code{TRUE}. When alternative specific constants are used \code{reduce} can be
#'\code{TRUE} so that the algorithm will be faster, but the combinations of
#'constants and profiles will not be evaluated exhaustively.
#'
#'When \code{algorithm = "MOD"}, the \code{\link[idefix]{SeqMOD}} function is 
#'used to select the choice set that minimizes the DB-error. Otherwise, 
#'\code{algorithm = "CEA"} uses the \code{\link[idefix]{SeqCEA}} function. 
#'
#'@param alts A character vector containing the names of the alternatives.
#'@param atts A character vector containing the names of the attributes.
#'@param n.total A numeric value indicating the total number of choice sets.
#'@param buttons.text A string containing the text presented together with the 
#'  option buttons.
#'@param no.choice An integer indicating which alternative should be a no choice
#'  alternative. The default is \code{NULL}.
#'@param intro.text A string containing the text presented before the choice 
#'  survey.
#'@param end.text A string containing the text presented after the choice 
#'  survey.
#'@param data.dir A character string with the directory denoting where the data
#'  needs to be written. The default is NULL
#'@param algorithm Can be either "MOD" or "CEA" to use the Modified Fedorov or the Coordinate Exchange algorithm, respectively. The default is "MOD".
#'@inheritParams Decode
#'@inheritParams Modfed
#'@inheritParams Profiles
#'@inheritParams SeqMOD
#'@inheritParams SeqCEA
#'@inheritParams ImpsampMNL
#'@importFrom Rdpack reprompt
#'@references \insertRef{ju}{idefix}
#'@return After completing the survey, two text files can be found in 
#'  \code{data.dir}. The file with "num" in the filename is a matrix with the 
#'  numeric choice data. The coded design matrix ("par"), presented during the 
#'  survey, together with the observed responses ("resp") can be found here. 
#'  Rownames indicate the setnumbers. The file with "char" in the filename is a 
#'  matrix with character choice data. The labeled design matrix ("par"), 
#'  presented during the survey, together with the observed responses ("resp") 
#'  can be found here. See \code{\link{LoadData}} to load the data.
#' @examples 
#' \donttest{
#'#### Present choice design without adaptive sets (n.total = sets in des)
#'# example design 
#'data("example_design") # pregenerated design
#'xdes <- example_design
#'### settings of the design 
#'code <- c("D", "D", "D")
#'n.sets <- 8
#'# settings of the survey
#'alternatives <- c("Alternative A", "Alternative B")
#'attributes <- c("Price", "Time", "Comfort")
#'labels <- vector(mode="list", length(attributes))
#'labels[[1]] <- c("$10", "$5", "$1")
#'labels[[2]] <- c("20 min", "12 min", "3 min")
#'labels[[3]] <- c("bad", "average", "good")
#'i.text <- "Welcome, here are some instructions ... good luck!"
#'b.text <- "Please choose the alternative you prefer"
#'e.text <- "Thanks for taking the survey"
#'dataDir <- getwd()
#'# Display the survey 
#'SurveyApp (des = xdes, n.total = n.sets, alts = alternatives, 
#'           atts = attributes, lvl.names = labels, coding = code, 
#'           buttons.text = b.text, intro.text = i.text, end.text = e.text, 
#'           algorithm = "MOD")
#'
#' #### Present choice design with partly adaptive sets (n.total > sets in des)
#' # example design 
#' data("example_design") # pregenerated design
#' xdes <- example_design
#' ### settings of the design 
#' code <- c("D", "D", "D")
#' n.sets <- 12
#' # settings of the survey
#' alternatives <- c("Alternative A", "Alternative B")
#' attributes <- c("Price", "Time", "Comfort")
#' labels <- vector(mode="list", length(attributes))
#' labels[[1]] <- c("$10", "$5", "$1")
#' labels[[2]] <- c("20 min", "12 min", "3 min")
#' labels[[3]] <- c("bad", "average", "good")
#' i.text <- "Welcome, here are some instructions ... good luck!"
#' b.text <- "Please choose the alternative you prefer"
#' e.text <- "Thanks for taking the survey"
#' # setting for adaptive sets 
#' levels <- c(3, 3, 3)
#' cand <- Profiles(lvls = levels, coding = code)
#' p.mean <- c(0.3, 0.7, 0.3, 0.7, 0.3, 0.7)
#' p.var <- diag(length(p.mean))
#' dataDir <- getwd()
#' # Display the survey 
#' SurveyApp(des = xdes, n.total = n.sets, alts = alternatives, 
#'           atts = attributes, lvl.names = labels, coding = code, 
#'           buttons.text = b.text, intro.text = i.text, end.text = e.text, 
#'           prior.mean = p.mean, prior.covar = p.var, cand.set = cand, 
#'           n.draws = 50, algorithm = "MOD")
#' # If CEA algorithm is desired, cand.set argument is not needed
#' SurveyApp(des = xdes, n.total = n.sets, alts = alternatives, 
#'           atts = attributes, lvl.names = labels, coding = code, 
#'           buttons.text = b.text, intro.text = i.text, end.text = e.text, 
#'           prior.mean = p.mean, prior.covar = p.var, n.draws = 50, 
#'           algorithm = "CEA")
#'           
#'#### Choice design with only adaptive sets (des=NULL)
#'# setting for adaptive sets 
#'levels <- c(3, 3, 3)
#'p.mean <- c(0.3, 0.7, 0.3, 0.7, 0.3, 0.7)
#'low = c(-Inf, -Inf, -Inf, 0, 0, -Inf)
#'up = rep(Inf, length(p.mean))
#'p.var <- diag(length(p.mean)) 
#'code <- c("D", "D", "D")
#'cand <- Profiles(lvls = levels, coding = code)
#'n.sets <- 12
#'# settings of the survey
#'alternatives <- c("Alternative A", "Alternative B")
#'attributes <- c("Price", "Time", "Comfort")
#'labels <- vector(mode="list", length(attributes))
#'labels[[1]] <- c("$10", "$5", "$1")
#'labels[[2]] <- c("20 min", "12 min", "3 min")
#'labels[[3]] <- c("bad", "average", "good")
#'i.text <- "Welcome, here are some instructions ... good luck!"
#'b.text <- "Please choose the alternative you prefer"
#'e.text <- "Thanks for taking the survey"
#'dataDir <- getwd()
#'# Display the survey 
#'SurveyApp(des = NULL, n.total = n.sets, alts = alternatives,
#'           atts = attributes, lvl.names = labels, coding = code, 
#'           buttons.text = b.text, intro.text = i.text, end.text = e.text, 
#'           prior.mean = p.mean, prior.covar = p.var, cand.set = cand, 
#'           lower = low, upper = up, n.draws = 50, algorithm = "MOD")
#' # If CEA algorithm is desired, cand.set argument is not needed
#'SurveyApp(des = NULL, n.total = n.sets, alts = alternatives,
#'          atts = attributes, lvl.names = labels, coding = code, 
#'          buttons.text = b.text, intro.text = i.text, end.text = e.text, 
#'          prior.mean = p.mean, prior.covar = p.var, 
#'          lower = low, upper = up, n.draws = 50, algorithm = "CEA")
#'          
#'#### Present choice design with a no choice alternative.
#'# example design 
#'data("nochoice_design") # pregenerated design
#'xdes <- nochoice_design
#'### settings of the design 
#'code <- c("D", "D", "D")
#'n.sets <- 8
#'# settings of the survey
#'alternatives <- c("Alternative A", "Alternative B", "None")
#'attributes <- c("Price", "Time", "Comfort")
#'labels <- vector(mode = "list", length(attributes))
#'labels[[1]] <- c("$10", "$5", "$1")
#'labels[[2]] <- c("20 min", "12 min", "3 min")
#'labels[[3]] <- c("bad", "average", "good")
#'i.text <- "Welcome, here are some instructions ... good luck!"
#'b.text <- "Please choose the alternative you prefer"
#'e.text <- "Thanks for taking the survey"
#'
#'# Display the survey 
#'SurveyApp(des = xdes, n.total = n.sets, alts = alternatives, 
#'           atts = attributes, lvl.names = labels, coding = code, 
#'           buttons.text = b.text, intro.text = i.text, end.text = e.text,
#'           no.choice = 3, alt.cte = c(0, 0, 1), algorithm = "MOD")
#'}
#'@import shiny
#'@export
SurveyApp <- function(des = NULL, n.total, alts, atts, lvl.names, coding,
                       alt.cte = NULL, no.choice = NULL, algorithm = "MOD", 
                      n.cs = NULL,
                      buttons.text, intro.text, end.text, data.dir = NULL,
                       c.lvls = NULL, prior.mean = NULL,
                       prior.covar = NULL, cand.set = NULL, n.draws = NULL, 
                       lower = NULL, upper = NULL, parallel = TRUE, reduce = TRUE) {
  # Initialize 
  sdata <- vector(mode = "list")
  surveyData <- vector(mode = "list")
  y.bin <- vector("numeric")
  resp  <- vector("character")
  n.atts <- length(atts)
  n.alts <- length(alts)
  n.levels <- as.vector(unlist(lapply(lvl.names,length)))
  choice.sets <- matrix(data = NA, nrow = n.total * n.alts, ncol = n.atts)
  buttons <- NULL
  sn <- 0
  if (is.null(des)) {
    n.init <- 0
  } else {
    n.init <- nrow(des) / n.alts
    if (!isTRUE(all.equal(n.init, as.integer(n.init)))) {
      stop("the number of rows of 'des' are not a multiple of length(alts)")
    }
  }
  if (is.null(alt.cte) || all(alt.cte == 0)) {
    alt.cte <- rep(0, n.alts)
    n.cte <- 0
    cte.des <- NULL
  } else {
    # Error 
    if (length(alt.cte) != n.alts) {
      stop("length(alts) does not match length(alt.cte)")
    }
    if (!all(alt.cte %in% c(0, 1))) {
      stop("'alt.cte' should only contain 0's or 1's.")
    }
    if (!any(alt.cte == 0)) {
      stop("'alt.cte' should at least contain 1 zero")
    }
    n.cte <- sum(alt.cte)
    if (!is.null(des)) {
      cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.init)
      if (!isTRUE(all.equal(cte.des, matrix(des[ , 1:n.cte], ncol = n.cte)))) {
        stop("the first column(s) of 'des' are different from what is expected based on 'alt.cte'")
      }
    }
  }
  # Error handling
  if (!is.null(no.choice)) {
    if (!is.numeric(no.choice)) {
      stop("'no.choice' should be an integer indicating the no choice alternative.")
    }
      if (!no.choice %% 1 == 0) {
        stop("'no.choice' should be an integer")
      }
      if (any(isTRUE(no.choice > (n.alts + 0.2)), isTRUE(no.choice < 0.2))) {
        stop("'no.choice' does not indicate one of the alternatives")
      }
      if (!isTRUE(all.equal(alt.cte[no.choice], 1))) {
        stop("the location of the 'no.choice' option in the 'alt.cte' vector should correspond with 1")
      }
  }
  if (!is.null(data.dir)) {
    if (!dir.exists(data.dir)) {
      stop("Directory 'data.dir' does not exist")
    }
  }
  if (n.total > n.init) {
    if (is.null(lower)) {
      lower <- rep(-Inf, length(prior.mean))
    }
    if (is.null(upper)) {
      upper <- rep(Inf, length(prior.mean))
    }
    if (!any(c(isTRUE(all.equal(length(prior.mean), length(lower))), isTRUE(all.equal(length(prior.mean), length(upper)))))) {
      stop("length 'prior.mean' should equal 'upper' and 'lower'")
    }
    if (algorithm == "MOD") {
      if (any(c(is.null(prior.mean), is.null(prior.covar), is.null(cand.set), 
               is.null(n.draws)))) {
      stop("When n.total is larger than the number of sets in argument des, arguments prior.mean, prior.covar, cand.set and n.draws should be specified.")
      }
      if (length(prior.mean) != ncol(cand.set) + sum(alt.cte)) {
        stop("Number of parameters in prior.mean does not match with cand.set + alt.cte")
      }
    } else if (algorithm == "CEA") {
      if (any(c(is.null(prior.mean), is.null(prior.covar), is.null(n.draws)))) {
      stop("When n.total is larger than the number of sets in argument des, arguments prior.mean, prior.covar and n.draws should be specified.")
      }
    }
    
    if (!isTRUE(all.equal(length(prior.mean), ncol(prior.covar)))) {
      stop("length of 'prior.mean' differs from number of columns 'prior.covar'")
    }
  } else {
    if (!is.null(prior.mean)) {
      warning("'prior.mean' will be ignored, since there are no adaptive sets.")
    } 
    if (!is.null(prior.covar)) {
      warning("'prior.covar' will be ignored, since there are no adaptive sets.")
    }
    if (algorithm == "MOD" & !is.null(cand.set)) {
      warning("'cand.set' will be ignored, since there are no adaptive sets.")
    }
    if (!is.null(lower) || !is.null(upper)) {
      warning("'lower' and 'upper' bound will be ignored, since there are no adaptive sets.")
    }
    if (!is.null(n.draws)) {
      warning("'n.draws' will be ignored, since there are no adaptive sets.")
    }
  }
  if (is.null(des)) {
    if (algorithm == "MOD") {
      fulldes <- matrix(data = NA, nrow = (n.alts * n.total), ncol = ncol(cand.set))
    } else {
      fulldes <- matrix()
    }
  } else {
    bs <- seq(1, (nrow(des) - n.alts + 1), n.alts)
    es <- c((bs - 1), nrow(des))[-1] 
    rowcol <- Rcnames(n.sets = n.init, n.alts = n.alts, alt.cte = alt.cte, 
                      no.choice = FALSE)
    rownames(des) <- rowcol[[1]]
    if (is.null(colnames(des))) {
      colnames(des) <- c(rowcol[[2]], paste("par", 1:(ncol(des) - n.cte), 
                                            sep = "."))
    }
    fulldes <- des
    # Error handling
    if (length(bs) != n.init) {
      stop("The number of rows in 'des' is not a multiple of length(atts)")
    }
    if ("no.choice.cte" %in% colnames(des)) {
      if (is.null(no.choice)) {
        warning("no.choice.cte column name detected in 'des' while 'no.choice = NULL'")
      }
    }
  }
  if (!(algorithm %in% c("MOD","CEA"))) {
    stop("algorithm should be 'MOD' or 'CEA'")
  }
  
  shinyApp(
    ### User interface
    ui <- fluidPage(
      # Put setnr on screen
      column(8, align = 'center', textOutput("set.nr")),
      # Put design on screen
      column(8, align = 'center', tableOutput("choice.set")),
      # Put answer options on screen
      column(8, align = 'center', uiOutput('buttons')), 
      # put introtext on screen
      column(8, align = 'center', textOutput('intro')),
      # Put action button on screen
      column(8, align = "center", actionButton("OK", "OK")),
      # put end text on screen
      column(8, align = 'center', textOutput('end'))
    ),
    ### Server
    server <- function(input, output) {
      # Count set number
      observeEvent(input$OK, {
        sn <<- sn + 1
      })
      # Set selection function
      Select <- function() {
        if (sn <= n.total) {
          # for initial sets 
          if (sn <= n.init) {
            set <- des[bs[sn]:es[sn], ]
          } else {
            ## sample drawing for adaptive sets
            # if First set
            if (sn == 1) {
              # sample draws from prior
              s <- tmvtnorm::rtmvnorm(n = n.draws, mean = prior.mean, 
                                      sigma = prior.covar, lower = lower, 
                                      upper = upper)
              w <- rep(1, nrow(s)) / nrow(s)
              if (sum(alt.cte) > 0.2) {
                s <- list(as.matrix(s[ , 1:sum(alt.cte)], ncol = sum(alt.cte)), 
                          s[ , -c(1:sum(alt.cte))])
              }
              # From second set
            } else {
              # Sample draws from updated posterior
              sam <- ImpsampMNL(n.draws = n.draws, prior.mean = prior.mean, 
                                prior.covar = prior.covar,
                                des = fulldes, n.alts = n.alts, y = y.bin, 
                                alt.cte = alt.cte, lower = lower, upper = upper)
              s <- sam$sample
              w <- sam$weights
            }
            ## Selecting set
            if (algorithm == "MOD") {
              # Select new set based on Modfed
              setobj <- SeqMOD(des = des, cand.set = cand.set, n.alts = n.alts, 
                               par.draws = s, prior.covar = prior.covar, 
                               alt.cte = alt.cte, weights = w, 
                               no.choice = no.choice, parallel = parallel, 
                               reduce = reduce)
            } else if (algorithm == "CEA") {
              setobj <- SeqCEA(des = des, lvls = n.levels, coding = coding,
                               n.alts = n.alts, par.draws = s, 
                               prior.covar = prior.covar, alt.cte = alt.cte,
                               weights = w, no.choice = no.choice, 
                               parallel = parallel, reduce = reduce)
            }
            set <- setobj$set
            db  <- setobj$db
            
            ## Design storage
            if (sn == 1) {
              rowcol <- Rcnames(n.sets = 1, n.alts = n.alts, alt.cte = alt.cte, no.choice = FALSE)
              rownames(set) <- rownames(set, do.NULL = FALSE, prefix = paste(paste("set", sn , sep = ""), "alt", sep = "."))
              colnames(set) <- c(rowcol[[2]], paste("par", 1:(ncol(set) - n.cte), sep = "."))
              fulldes <<- set
            } else {
              rowcol <- Rcnames(n.sets = 1, n.alts = n.alts, alt.cte = alt.cte, no.choice = FALSE)
              rownames(set) <- rownames(set, do.NULL = FALSE, prefix = paste(paste("set", sn , sep = ""), "alt", sep = "."))
              colnames(set) <- c(rowcol[[2]], paste("par", 1:(ncol(set) - n.cte), sep = "."))
              fulldes <<- rbind(fulldes, set)
            }
          }
          # Transform coded set to attribute level character set.
          choice.set <- Decode(des = set, n.alts = n.alts, lvl.names = lvl.names, coding = coding, 
                               alt.cte = alt.cte, c.lvls = c.lvls, no.choice = no.choice)[[1]]
          choice.set <- t(choice.set[ , 1:n.atts])
          # Fill in attribute names and alternatives names
          colnames(choice.set) <- alts
          rownames(choice.set) <- atts
          # Store uncoded choice set
          if (sn == 1) {
            choice.sets <<- choice.set
          } else {
            choice.sets <<- rbind(choice.sets, choice.set)
          }
          #return design 
          if (!is.null(no.choice)) {
            no.choice.set <- choice.set[ ,-no.choice]
            return(no.choice.set)
          } else {
            return(choice.set)
          }
        }
      }
      #When action button is clicked
      observeEvent(input$OK, {
        # survey phase 
        if (sn <= n.total ) {
          # Plot new choice set
          output$choice.set <-  renderTable(Select(), rownames = TRUE)
        }
        # Store responses and design
        if (sn > 1 && sn <= (n.total + 1)) {
          resp  <<- c(resp, input$survey)
          y.bin <<- Charbin(resp = resp, alts = alts, n.alts = n.alts)
          sdata[["bin.responses"]] <- y.bin
          sdata[["responses"]] <- resp
          sdata[["desing"]] <- fulldes
          sdata[["survey"]] <- choice.sets
          surveyData <<- sdata 
        } 
        # end phase 
        if (sn > n.total) {
          #Don't show choice set
          output$choice.set <-  renderTable(NULL)
        }
      })
      #Output response options after first action button click
      output$buttons <- renderUI({
        # radiobuttons
        if (input$OK > 0 && input$OK <= n.total) {
          return(list(radioButtons("survey", buttons.text,
                                   alts , inline = T, selected = "None")))
        }
      })
      # set nr
      observeEvent(input$OK, {
        if (sn < n.total) {
          output$set.nr <- renderText(paste(c("choice set:", sn, "/", n.total)))
        } else {output$set.nr <- renderText(NULL)}
      })
      # Introtext
      output$intro <- renderText(intro.text)
      observeEvent(input$OK, {
        output$intro <- renderText(NULL)
      })
      # End of survey
      observeEvent(input$OK, {
        # Display end text 
        if (input$OK > n.total) {
          # Display end text 
          output$end <- renderText(end.text)
        }
        # Quit application 
        if (input$OK > (n.total + 1)) {
          # Write data to file
          if (!is.null(data.dir)) {
            saveData(data = surveyData, data.dir = data.dir, n.atts = n.atts)
          }
          # Stop application 
          stopApp()
        }
      })
    }
  )
}



#' Coded design to readable design.
#' 
#' Transforms a coded design matrix into a design containing character attribute
#' levels, ready to be used in a survey. The frequency of each attribute level 
#' in the design is also included in the output.
#' 
#' \code{des} A design matrix, this can also be a single choice set. See for
#' example the output of \link[idefix]{Modfed} or \link[idefix]{CEA}.
#' 
#' In \code{lvl.names}, the number of character vectors in the list should equal
#' the number of attributes in de choice set. The number of elements in each 
#' character vector should equal the number of levels for that attribute.
#' 
#' Valid arguments for \code{coding} are \code{C}, \code{D} and \code{E}. When 
#' using \code{C} the attribute will be treated as continuous and no coding will
#' be applied. All possible levels of that attribute should then be specified in
#' \code{c.lvls}. If \code{D} (dummy coding) is used 
#' \code{\link{contr.treatment}} will be applied to that attribute. The first 
#' attribute wil be used as reference level.  For \code{E} (effect coding) 
#' \code{\link{contr.sum}} is applied, in this case the last attribute level is 
#' used as reference level.
#' 
#' If \code{des} contains columns for alternative specific constants, 
#' \code{alt.cte} should be specified. In this case, the first column(s) (equal 
#' to the number of nonzero elements in \code{alt.cte}) will be removed from 
#' \code{des} before decoding the alternatives.
#'
#' @param des A numeric matrix which represents the design matrix. Each row is a
#'   profile.
#' @param lvl.names A list containing character vectors with the values of each 
#'   level of each attribute.
#' @param coding A character vector denoting the type of coding used for each 
#'   attribute. See also \code{\link{Profiles}}.
#' @param alt.cte A binary vector indicating for each alternative if an 
#'   alternative specific constant is present. The default is \code{NULL}.
#' @param no.choice An integer indicating the no choice alternative. The default
#'   is \code{NULL}.
#' @inheritParams Profiles
#' @inheritParams Modfed
#' @return \item{design}{A character matrix which represents the design.} 
#'   \item{lvl.balance}{A list containing the frequency of appearance of each 
#'   attribute level in the design.}
#' @examples 
#' \donttest{
#' # Example without continuous attributes.
#' design <- example_design 
#' coded <- c("D", "D", "D") # Coding.
#' # Levels as they should appear in survey. 
#' al <- list(
#'   c("$50", "$75", "$100"), # Levels attribute 1.
#'   c("2 min", "15 min", "30 min"), # Levels attribute 2.
#'   c("bad", "moderate", "good") # Levels attribute 3.
#' ) 
#' # Decode
#' Decode(des = design, n.alts = 2, lvl.names = al, coding = coded) 

#' 
#' # Example with alternative specific constants
#' design <- example_design2 
#' coded <- c("D", "D", "D") # Coding.
#' # Levels as they should appear in survey. 
#' al <- list(
#'   c("$50", "$75", "$100"), # Levels attribute 1.
#'   c("2 min", "15 min", "30 min"), # Levels attribute 2.
#'   c("bad", "moderate", "good") # Levels attribute 3.
#' ) 
#' # Decode
#' Decode(des = design, n.alts = 3, lvl.names = al, coding = coded, alt.cte = c(1, 1, 0)) 
#' }
#' @export
Decode <- function(des, n.alts, lvl.names, coding, alt.cte = NULL, c.lvls = NULL, no.choice = NULL) {
  n.sets <- nrow(des) / n.alts
  if(!isTRUE(all.equal(n.sets %% 1, 0))){
    stop("'n.alts' does not seem correct based on nrow(des)")
  }
  if(!is.list(lvl.names)){
    stop("'lvl.names' should be a list")
  }
  if(!is.matrix(des)){
    stop("'des' should be a matrix")
  }
  if(!is.null(alt.cte)) {
    if(!isTRUE(all.equal(length(alt.cte), n.alts))){
      stop("the length of 'alt.cte' does not match 'n.alts'")
    }
    contins <- which(alt.cte == 1)
    if( !length(contins) == 0) {
      des <- des[, -(1:length(contins))]
    }
  }
  if(!is.null(no.choice)){
    if(!is.numeric(no.choice)){
      stop("'no.choice' should be an integer indicating the no choice alternative.")
    }
    if(!no.choice %% 1 == 0){
      stop("'no.choice' should be an integer")
    }
    if(any(isTRUE(no.choice > (n.alts + 0.2)), isTRUE(no.choice < 0.2))){
      stop("'no.choice' does not indicate one of the alternatives")
    }
    if(!isTRUE(all.equal(alt.cte[no.choice], 1))){
      stop("the location of the 'no.choice' option in the 'alt.cte' vector should correspond with 1")
    }
    ncsek <- 1:nrow(des)
    ncsek <- ncsek[-seq(no.choice, (n.sets * no.choice), no.choice)] 
    #des <- des[-ncsek, ]
  } else {
    ncsek <- 1:nrow(des)
  }
  N.alts <- nrow(des) # Number of total alternatives.
  n.att <- length(lvl.names) # Number of attributes.
  conts <- which(coding == "C") # Continuous levels. 
  # Create vector where each element denotes the number of levels for each attribute.
  lvls <- numeric(n.att) 
  for (i in 1:n.att) { 
    lvls[i] <- length(lvl.names[[i]])
  }
  # Generate all possible profiles coded and uncoded
  dc <- Profiles(lvls = lvls, coding = coding, c.lvls = c.lvls)
  # Create uncoded grid. 
  d <- as.data.frame(expand.grid(lvl.names))
  # Create new matrix for choice des with attribute level names 
  m <- matrix(data = NA, nrow = N.alts, ncol = n.att)
  # Error handling
  if (ncol(des) != ncol(dc)) {
    stop("ncol(des) does not equal the expected number based on 'lvl.names' and 'alt.cte'.")
  }
  # For each alternative look for matching profile  
  for (i in ncsek) {
    # if coded choice des, look for match in coded version first, then take uncoded equivalent.
    lev.num <- d[as.numeric(which(apply(dc, 1, function(x) all(x == des[i, ])))), ]
    lev.num <- as.numeric(lev.num)
    # Error handling
    if (any(is.na(lev.num))) { 
      stop("alternatives detected that do not match with the type of 'coding' provided. 
           Make sure 'no.choice' is specified if applicable.")
    }
    # For each attribute fill in the attribute level name
    for (c in 1:n.att) {m[i, c] <- lvl.names[[c]][lev.num[c]]}
  }
  #frequency levels 
  ml <- list()
  ml <- split(m, rep(1:ncol(m), each = nrow(m)))
  lvl.balance <- lapply(ml, table)
  names(lvl.balance) <- paste("attribute", 1:length(coding), "")
  #col row names names 
  rownames(m) <- rownames(des)
  return(list("design" = as.data.frame(m), "lvl.balance" = lvl.balance))
}

# Character vector to binary vector.
# 
# Transforms a character vector with responses into a binary vector. Each
# alternative in each choice set wil be either 0 or 1. If the
# alternative was not chosen 0, and 1 if it was. The function can be used for example in a
# shiny application to transform the response vector received from
# \code{\link[shiny]{radioButtons}} into a numeric vector that can be used for
# estimation.
# 
# The \code{n.alts} denotes the number of alternatives a respondent could
# choose from, without counting a possible no choice option.
# 
# If \code{no.choice} is \code{TRUE} the first alternative specified in 
# \code{alts} will be treated as a no choice option. If the no choice option 
# was chosen all alternatives are zero for that choice set.
# @param resp String vector containing input responses
# @param alts String vector containing all possible alternatives. The order
#   should be the same as the order of the design matrix.
# @param n.alts The number of alternatives per choice set.
# @param no.choice Logical value indicating whether a no.choice option is 
#   provided or not. The default = \code{FALSE}.
# @return A binary response vector with length equal to \code{length(resp) *
#   length(n.alts)}.
# @examples 
# \donttest{
# # Observed Responses 
# resp <- c("alt1", "alt3", "alt2", "no.choice", "alt1") 
# # All possible alternatives 
# alts <- c("no.choice", "alt1", "alt2", "alt3")
# # 3 alternatives + no.choice 
# Charbin(resp = resp, alts = alts, n.alts = 3, no.choice = TRUE)
# }
Charbin <- function (resp, alts, n.alts, no.choice = FALSE) {
  # Error resp not in altsions
  if (!all(resp %in% alts)) {
    stop("1 or more responses do not match the possible response options.")
  }
  # Error altsions
  if (length(alts) != (n.alts + no.choice)) {
    stop("Number of response options is not correct")
  }
  map <- match(resp, alts)
  l <- list()
  for(i in 1:length(map)){
    l[[i]] <- rep(0, n.alts)
    if (no.choice) {
      l[[i]][map[i] - 1] <- 1
    } else {
      l[[i]][map[i]] <- 1
    }
  }
  v <- unlist(l)
  return(v)
}


# Function to save the data gathered by shiny app
saveData <- function(data, data.dir, n.atts) {
  # Data manipulation 
  d <- as.data.frame(cbind(data$desing, resp = data$bin.responses))
  unc_resp <- rep(data$responses, each = n.atts) 
  unc_setnr <- rep(1:length(data$responses), each = n.atts)
  unc_d <- cbind(set = unc_setnr, data$survey, resp = unc_resp) 
  # Create unique file names
  numname <- sprintf("%s_num_data.txt", as.integer(Sys.time()))
  charname <- sprintf("%s_char_data.txt", as.integer(Sys.time()))
  # Write files to data.dir
  utils::write.table(
    x = d,
    file = file.path(data.dir, numname), 
    row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA
  )
  utils::write.table(
    x = unc_d,
    file = file.path(data.dir, charname), 
    row.names = TRUE, quote = FALSE, sep = "\t", col.names = NA
  )
}

#' Load numeric choice data from directory
#' 
#' Reads all individual choice data files, created by \code{\link{SurveyApp}}
#' function, from a directory and concatenates 
#' those files into a single data file. Files containing either "num" or "char"
#' will be read, with num indicating numeric data and char indicating character
#' data. For more information, see output of \code{\link{SurveyApp}}.
#' 
#' @param data.dir A character string containing the directory to read from.
#' @param type Character vector containing either num or char. 
#' @return A data frame containg the full design and all the responses of the 
#'   combined data files that were found. Different files are indicated by an
#'   ID variable.
#' @export
LoadData <- function(data.dir, type) {
  # ErrorS
  if(!type %in% c("num", "char")){
    stop("'type' must be either num or char")
  }
  if (!dir.exists(data.dir)) {
    stop("Directory 'data.dir' does not exist")
  }
  error <- character(0)
  if(identical(list.files(data.dir, full.names = TRUE, pattern = type), error)){
    stop("No files of the specified 'type' in 'data.dir'")
  } 
  # Read all files into list
  files <- list.files(data.dir, full.names = TRUE, pattern = type)
  data <- lapply(files, utils::read.table, stringsAsFactors = FALSE, sep = '\t', header = T) 
  # check same col 
  ncols <- unlist(lapply(data, function(x) return(ncol(x))))
  if(!isTRUE(all.equal(min(ncols), max(ncols)))){
    stop("'data.dir'contains files with different number of columns")
  }
  # matrix
  id_rows <- sapply(data, nrow)
  id <- rep(1:length(id_rows), id_rows)
  data <- lapply(data, as.data.frame)
  # Bind together
  data <- do.call(rbind, data)
  data <- cbind("ID" = id, data)
  return(data)
}

