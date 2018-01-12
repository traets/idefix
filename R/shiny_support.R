#'Shiny application to generate a discrete choice survey.
#'
#'This function starts a shiny application which puts choice sets on screen and 
#'saves the responses. The complete choice design can be provided in advance, or
#'can be generated sequentially adaptively, or can be a combination of both.
#'
#'A pregenerated design can be specified in \code{des}. This should be a matrix 
#'in which each row is a profile. This can be generated with \code{Modfed}, but 
#'is not necesarry.
#'
#'If \code{n.total} = \code{nrow(des)} / \code{length(alts)}, the specified
#'design will be put on screen, one set after the other, and the responses will
#'be saved. If \code{n.total} > (\code{nrow(des)} / \code{length(alts)}), first
#'the specified design will be shown and afterwards the remaining sets will be 
#'generated adaptively. If \code{des} = \code{NULL}, \code{n.total} sets will be
#'generated adaptively.
#'
#'Whenever adaptive sets will be generated, \code{crit}, \code{prior.mean}, 
#'\code{prior.covar}, \code{cand.set} and \code{m}, should be specified.
#'
#'The names specified in \code{alts} will be used to label the choice 
#'alternatives. The names specified in \code{atts} will be used to name the 
#'attributes in the choice sets. The values of \code{lvl.names} will be used to 
#'create the values in the choice sets. See \code{\link{Decode}} for more 
#'details. The number of draws sampeled from the posterior preference 
#'distribution in the importance sampling algorithm used for adaptive sets can 
#'be specified with \code{m}, where the number is 2^\code{m}.
#'
#'The text specified in \code{buttons.text} will be displayed above the buttons 
#'to indicate the preferred choice (for example: "indicate your preferred 
#'choice"). The text specified in \code{intro.text} will be displayed before the
#'choice sets. This will generally be a description of the survey and some 
#'instructions. The text specified in \code{end.text} will be displayed after 
#'the survey. This will generally be a thanking note and some further 
#'instructions.
#'
#'
#'@param alts A character vector containing the names of the alternatives.
#'@param atts A character vector containing the names of the attributes.
#'@param n.total A numeric value indicating the total number of choice sets.
#'@param buttons.text A string containing the text presented together with the 
#'  option buttons.
#'@param intro.text A string containing the text presented before the choice 
#'  survey.
#'@param end.text A string containing the text presented after the choice 
#'  survey.
#'@param data.dir A character string with the directory denoting where the data
#'  needs to be written. The default is NULL
#'@param crit A string containing eihter KL or DB indicating the adaptive
#'  criterion to be used.
#'@inheritParams Decode
#'@inheritParams Modfed
#'@inheritParams Profiles
#'@inheritParams SeqKL
#'@inheritParams ImpsampMNL
#'@importFrom Rdpack reprompt
#'@references \insertRef{crabbe}{mnldes}
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
#'# NOTE that the data will be saved in the current working directory. 
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
#'           data.dir = dataDir)
#'# Data 
#'data_num <- LoadData(data.dir = dataDir, type  = "num")
#'data_char <- LoadData(data.dir = dataDir, type = "char")
#'
#'#### Present choice design with adaptive sets (n.total > sets in des)
#'# NOTE that the data will be saved in the current working directory. 
#'# example design 
#'data("example_design") # pregenerated design
#'xdes <- example_design
#'### settings of the design 
#'code <- c("D", "D", "D")
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
#'# setting for adaptive sets 
#'levels <- c(3, 3, 3)
#'cand <- Profiles(lvls = levels, coding = code)
#'p.mean <- c(0.3, 0.7, 0.3, 0.7, 0.3, 0.7)
#'p.var <- diag(length(p.mean))
#'dataDir <- getwd()
#'# Display the survey 
#'SurveyApp (des = NULL, n.total = n.sets, alts = alternatives, atts =
#'attributes, lvl.names = labels, coding = code, buttons.text = b.text,
#'intro.text = i.text, end.text = e.text, data.dir = dataDir, crit= "KL",
#'prior.mean = p.mean, prior.covar = p.var, cand.set = cand, m = 6)
#'# Data 
#'data_num <- LoadData(data.dir = dataDir, type = "num")
#'data_char <- LoadData(data.dir = dataDir, type = "char")
#'
#'#### Choice design with only adaptive sets (des=NULL)
#'# NOTE that the data will be saved in the current working directory. 
#'# setting for adaptive sets 
#'levels <- c(3, 3, 3)
#'p.mean <- c(0.3, 0.7, 0.3, 0.7, 0.3, 0.7)
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
#'SurveyApp (des = NULL, n.total = n.sets, alts = alternatives, 
#'           atts = attributes, lvl.names = labels, coding = code, 
#'           buttons.text = b.text, intro.text = i.text, end.text = e.text, data.dir = dataDir, 
#'           crit= "KL", prior.mean = p.mean, prior.covar = p.var, cand.set = cand, m = 6)
#'# Data 
#'data_num <- LoadData(data.dir = dataDir, type = "num")
#'data_char <- LoadData(data.dir = dataDir, type = "char")
#'}
#'@import shiny
#'@export
SurveyApp <- function(des = NULL, n.total, alts, atts, lvl.names, coding, 
                      buttons.text, intro.text, end.text, data.dir = NULL,
                      c.lvls = NULL, crit = NULL, alt.cte = NULL, prior.mean = NULL,
                      prior.covar = NULL, cand.set = NULL, m = NULL) {
  # Initialize 
  sdata <- vector(mode = "list")
  surveyData <- vector(mode = "list")
  y.bin <- vector("numeric")
  resp  <- vector("character")
  n.atts <- length(atts)
  n.alts <- length(alts)
  choice.sets <- matrix(data = NA, nrow = n.total * n.alts, ncol = n.atts)
  buttons <- NULL
  sn <- 0
  
  if (is.null(alt.cte)) {
    alt.cte <- rep(0, n.alts)
    cte.des <- NULL
  } else {
    # Error 
    if (!all(alt.cte %in% c(0,1))){
      stop("alt.cte should only contain 0s or 1s.")
    }
  }
  if (is.null(des)) {
    n.init <- 0
    fulldes <- matrix(data = NA, nrow = (n.alts * n.total), ncol = ncol(cand.set))
  } else {
    n.init <- nrow(des) / n.alts 
    bs <- seq(1, (nrow(des) - n.alts + 1), n.alts)
    es <- c((bs - 1), nrow(des))[-1] 
    if (sum(alt.cte) > 0) {
      cte.des <- Altspec(alt.cte = alt.cte, n.sets = (nrow(des) / n.alts))
      colnames(cte.des) <- paste(paste("alt", which(alt.cte == 1), sep = ""), ".cte", sep = "")
    }
    colnames(des) <- paste("par", 1 : ncol(des), sep = ".")
    fulldes <- cbind(cte.des, des)
    # Error handling
    if (length(bs) != n.init) {
      stop("The number of design rows does not match the number of alternatives times the number of sets.")
    }
  }
  # Error handling
  if (!is.null(data.dir)) {
    if (!dir.exists(data.dir)) {
      stop("Directory data.dir does not exist")
    }
  }
  if (n.total > n.init) {
      if (any(c(is.null(prior.mean), is.null(prior.covar), is.null(cand.set), is.null(m), is.null(crit)))) {
        stop("When n.total is larger than the number of sets in argument des, arguments crit, prior.mean, prior.covar, cand.set and m should be specified.")
      }
    if (length(prior.mean) != ncol(cand.set) + sum(alt.cte)) {
      stop("Number of parameters in prior.mean does not match with cand.set + alt.cte")
    }
  } else {
    if (!is.null(prior.mean)) {
      warning("prior.mean will be ignored, since there are no adaptive sets.")
    } 
    if (!is.null(prior.covar)) {
      warning("prior.covar will be ignored, since there are no adaptive sets.")
    }
    if (!is.null(cand.set)) {
      warning("cand.set will be ignored, since there are no adaptive sets.")
    }
    if (sum(alt.cte) > 0) {
      warning("alt.cte will be ignored, since there are no adaptive sets.")
    }
    if (!is.null(m)) {
      warning("m will be ignored, since there are no adaptive sets.")
    }
  }
  if (crit =="DB" && is.null(des)) {
    stop("In order to use the DB criterion, an initial design has to be provided.")
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
      Select <- function () {
        if (sn <= n.total) {
          # for initial sets 
          if (sn <= n.init) {
            set <- des[bs[sn] : es[sn], ]
          } else {
            ## sample drawing for adaptive sets
            # if First set
            if (sn == 1) {
              # sample draws from prior
              s <- MASS::mvrnorm(n = 2 ^ m, mu = prior.mean, Sigma  = prior.covar)
              w <- rep(1, nrow(s)) / nrow(s)
              # From second set
            } else {
              # Sample draws from updated posterior
              sam <- ImpsampMNL(prior.mean = prior.mean, prior.covar = prior.covar, des = fulldes, n.alts = n.alts, y = y.bin, m = m)
              s <- sam$sample
              w <- sam$weights
            }
            
            ## Selecting set
            if (crit == "KL") {
              # Select new set based on KL info
              set <- SeqKL(cand.set = cand.set, n.alts = n.alts, par.draws = s, alt.cte = alt.cte, weights = w)$set
              #delete alt.cte if necessary
              if (sum(alt.cte) > 0) {
                set <- set[ , -(1 : (sum(alt.cte)))]
              }
            } else if (crit == "DB") {
              # Select new set based on DB 
              setobj <- SeqDBApp(des = des, cand.set = cand.set, n.alts = n.alts, par.draws = s, prior.covar = prior.covar, alt.cte = alt.cte, w = w)
              set <- setobj$set
              db  <- setobj$db
            } else {
              stop("Argument crit should eihter be KL or DB.")
            }
            
            ## Design storage
            if (sn == 1) { 
              rownames(set) <- rownames(set, do.NULL = FALSE, prefix = paste(paste("set", sn , sep = ""), "alt", sep = "."))
              colnames(set) <- paste("par", 1:ncol(set), sep = ".")
              des <<- set
              # with alt.cte
              altset <- Altspec(alt.cte, n.sets = 1)
              if (sum(alt.cte) > 0) {
                colnames(altset) <- paste(paste("alt", which(alt.cte == 1), sep = ""), ".cte", sep = "")
              }
              fullset <- cbind(altset, set)
              fulldes <<- fullset
            } else {
              rownames(set) <- rownames(set, do.NULL = FALSE, prefix = paste(paste("set", sn , sep = ""), "alt", sep = "."))
              colnames(set) <- paste("par", 1:ncol(set), sep = ".")
              des <<- rbind(des, set)
              # with alt.cte
              altset <- Altspec(alt.cte, n.sets = 1)
              if (sum(alt.cte) > 0) {
                colnames(altset) <- paste(paste("alt", which(alt.cte == 1), sep = ""), ".cte", sep = "")
              }
              fullset <- cbind(altset, set)
              fulldes <<- rbind(fulldes, fullset)
            }
          }
          # Transform coded set to attribute level character set.
            choice.set <- Decode(set = set, lvl.names = lvl.names, coding = coding, c.lvls = c.lvls)
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
          return(choice.set)
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
        if (sn > 1 && sn <= (n.total +1)) {
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
          if (!is.null(data.dir)){
            saveData(data = surveyData, data.dir = data.dir, n.atts = n.atts)
          }
          # Stop application 
          stopApp()
        }
      })
    }
  )
}


#' Coded choice set to character choice set.
#' 
#' Transforms a coded choice set into a choice set containing character attribute
#' levels, ready to be used in a survey.
#' 
#' In \code{lvl.names}, the number of character vectors in the list should equal
#' the number of attributes in de choice set. The number of elements in each 
#' character vector should equal the number of levels for that attribute.
#' 
#' Valid arguments for \code{coding} are \code{C}, \code{D} and \code{E}. When
#' using \code{C} the attribute will be treated as continuous and no coding will
#' be applied. All possible levels should then be specified in \code{c.lvls}. If
#' \code{D} (dummy coding) is used \code{\link{contr.treatment}} will be applied
#' to that attribute. The first attribute wil be used as reference level.  For
#' \code{E} (effect coding) \code{\link{contr.sum}} is applied, in this case the
#' last attributelevel is used as reference level.
#' 
#' @param set A numeric matrix which represents a choice set. Each row is a
#'   profile.
#' @param lvl.names A list containing character vectors with the values of each
#'   level of each attribute.
#' @param coding A character vector denoting the type of coding used for each
#'   attribute. See also \code{\link{Profiles}}.
#' @param alt.cte A binary vector indicating for each alternative if an 
#'   alternative specific constant is present. The default is \code{NULL}. 
#' @inheritParams Profiles
#' @inheritParams Modfed
#' @return A character matrix which represents the choice set.
#' @examples 
#' # Example without continuous attributes.
#' l <- c(3, 4, 2) # 3 Attributes.
#' c <- c("D", "E", "D") # Coding.
#' # All profiles.
#' p <- Profiles(lvls = l, coding = c)
#' cs <- p[c(4, 8), ] # Choice set 
#' # Levels as they should appear in survey. 
#' al <- list(
#'  c("$50", "$75", "$100"), # Levels attribute 1.
#'  c("2 min", "15 min", "30 min", "50 min"), # Levels attribute 2.
#'  c("bad", "good") # Levels attribute 3.
#' ) 
#' # Decode
#' Decode(set = cs, lvl.names = al, coding = c, alt.cte = c(0, 0)) 
#'
#' # Example with continuous attribute.
#' l <- c(3, 4, 2) # 3 Attributes.
#' c <- c("D", "C", "D") # Coding.
#' cl <- list(c(50, 75, 80, 100))
#' # All profiles.
#' p <- Profiles(lvls = l, coding = c, c.lvls = cl)
#' cs <- p[c(4, 8), ] # Set. 
#' a <- c(1, 0) # Alternative specific constant. 
#' cs <- cbind(a, cs) # set with alt.cte
#' # Levels as they should appear in survey. 
#' al <- list(
#'   c("$50", "$75", "$100"), # Levels attribute 1.
#'   c("50 min", "75 min", "80 min", "100 min"), # Levels attribute 2.
#'   c("bad", "good") # Levels attribute 3.
#' ) 
#' # Decode
#' Decode(set = cs, lvl.names = al, coding = c, alt.cte = c(1, 0), c.lvls = cl) 
Decode <- function(set, lvl.names, coding, alt.cte = NULL, c.lvls = NULL) {
  
  if(!is.null(alt.cte)) {
    contins <- which(alt.cte == 1)
    if( !length(contins) == 0) {
      set <- set[, -length(contins)]
    }
  }
  
  n.alts <- nrow(set) # Number of alternatives.
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
  # Create new matrix for choice set with attribute level names 
  m <- matrix(data = NA, nrow = n.alts, ncol = n.att)
  # Error handling
  if (ncol(set) != ncol(dc)) {
    stop("Number of columns of the set does not match expected number based on the other arguments.")
  }
  # For each alternative look for matching profile  
  for (i in 1:n.alts) {
    # if coded choice set, look for match in coded version first, then take uncoded equivalent.
    lev.num <- d[as.numeric(which(apply(dc, 1, function(x) all(x == set[i, ])))), ]
    lev.num <- as.numeric(lev.num)
    # Error handling
    if (any(is.na(lev.num))) { 
      stop('The set does not match with the type of coding provided')
    }
    # For each attribute fill in the attribute level name
    for (c in 1:n.att) {
      m[i, c] <- lvl.names[[c]][lev.num[c]]
    }
  }
  return(m)
}


#' Character vector to binary vector.
#' 
#' Transforms a character vector with responses into a binary vector. Each
#' alternative in each choice set wil be either 0 or 1. If the
#' alternative was not chosen 0, and 1 if it was. The function can be used for example in a
#' shiny application to transform the response vector received from
#' \code{\link[shiny]{radioButtons}} into a numeric vector that can be used for
#' estimation.
#' 
#' The \code{n.alts} denotes the number of alternatives a respondent could
#' choose from, without counting a possible no choice option.
#' 
#' If \code{no.choice} is \code{TRUE} the first alternative specified in 
#' \code{alts} will be treated as a no choice option. If the no choice option 
#' was chosen all alternatives are zero for that choice set.
#' @param resp String vector containing input responses
#' @param alts String vector containing all possible alternatives. The order
#'   should be the same as the order of the design matrix.
#' @param n.alts The number of alternatives per choice set.
#' @param no.choice Logical value indicating whether a no.choice option is 
#'   provided or not. The default = \code{FALSE}.
#' @return A binary response vector with length equal to \code{length(resp) *
#'   length(n.alts)}.
#' @examples 
#' # Observed Responses 
#' resp <- c("alt1", "alt3", "alt2", "no.choice", "alt1") 
#' # All possible alternatives 
#' alts <- c("no.choice", "alt1", "alt2", "alt3")
#' # 3 alternatives + no.choice 
#' Charbin(resp = resp, alts = alts, n.alts = 3, no.choice = TRUE)
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


# Binary to discrete choice matrix.
# 
# Transforms a numeric matrix with binary choice data for each respondent 
# (columns), to a matrix with discrete values representing the chosen 
# alternatives.
# @param y NUmeric matrix containing the binary choice data. Each column is a 
#   different ID.
# @param n.alts Numeric value indicating the number of alternatives per choice 
#   set.
# @param no.choice Logical value indicating whether a no choice response could 
#   be observed. This would be a \code{0} for each alternative.
# @return A matrix with discrete values, indicating the chosen alternatives per
#   ID.
# @examples  
# # Binary response data, 2 participants
# y <- matrix(data = c(0,1,1,0,0,0,0,1), ncol = 2, byrow = FALSE)
# # no choice = TRUE 
# BinDis(y = y, n.alts = 2, no.choice = TRUE)
BinDis <- function(y, n.alts, no.choice) {
  # y matrix.
  if (!is.matrix(y)) {
    stop('y should be a matrix.')
  }
  # Error no.choice
  if(!is.logical(no.choice)) {
    stop('no.choice should be logical.')
  }
  # Create choice sets.
  f <- function(x) {
    split(x, ceiling(seq_along(x) / n.alts))
  }
  cs <- apply(y, 2, f)
  # Error n.alts 
  for (i in 1:ncol(y)){
    if((length(unique(lengths(cs[[i]]))) != 1L)){
      stop('length of Y vector does match expected length based on nr of alternatives.')
    }
  }
  # Index 1.
  Ones <- function(x) {
    xx <- (x == 1)
    ione <- which(xx, arr.ind = TRUE)
    if(length(ione) > 1) {
      stop('Multiple alternatives are chosen per choice set. The response data or the number of alternatives is probably incorrect.')
    }
    # When no choice was chosen
    if (length(ione) == 0) {
      if (!no.choice) {
        stop('no choice responses detected while no.choice = FALSE.')
      }
      ione <- 0
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

# Sequential DB function for shiny
# 
# Small changes in alt.cte argument in comparison with the \code{SeqDB}
# function. This way the function can be easily used in the SurveyApp function
SeqDBApp <- function(des, cand.set, n.alts, par.draws, prior.covar, alt.cte, reduce = TRUE, w = NULL) {
  # Initialize.
  n.sets <- nrow(des) / n.alts
  # If no w, equal w.
  if (is.null(w)) {
    w <- rep(1, nrow(par.draws))
  }
  # Create alternative specific design if necessay.
  if (!all(alt.cte == 0)) {
    cte.set <- Altspec(alt.cte = alt.cte, n.sets = 1)
    cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
    kcte <- ncol(cte.des)
    fdes <- cbind(cte.des, des)
  } else {
    cte.set <- NULL
    kcte <- 0
    fdes <- des
  }
  # Error handling cte.des
  if (ncol(cand.set) + kcte != ncol(par.draws)) {
    stop("dimension of par.draws does not match the dimension of alt.cte + cand.set.")
  }
  # Handling par.draws.
  if (!(is.matrix(par.draws))) {
    par.draws <- matrix(par.draws, nrow = 1)
  }
  # Error identifying model.
  if (n.sets < ncol(par.draws)) {
    stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
  }
  # Error par.draws
  if (ncol(fdes) != ncol(par.draws)) {
    stop("Numbers of parameters in par.draws does not match the number of parameters in the design.")
  }
  # Starting and initializing values.
  i.cov <- solve(prior.covar)
  d.start <- apply(par.draws, 1, Derr, des = fdes,  n.alts = n.alts)
  db.start <- mean(d.start, na.rm = TRUE)
  full.comb <- gtools::combinations(n = nrow(cand.set), r = n.alts, repeats.allowed = !reduce)
  n.par <- ncol(par.draws)
  # For each potential set, select best. 
  db.errors <- apply(full.comb, 1, DBerrS, cand.set, par.draws, fdes, n.alts, cte.set, i.cov, n.par, w)
  comb.nr <- as.numeric(full.comb[which.min(db.errors), ])
  set <- cand.set[comb.nr, ]
  row.names(set) <- NULL
  db <- min(db.errors)
  #return best set and db error design.
  return(list(set = set, db.error = db))
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
#' Reads all individual choice data files from a directory and concatenates 
#' those files into a single data file. Files containing either "num" or "char"
#' will be read, with num indicating numeric data and char indicating character
#' data. for more information see output of \code{\link{SurveyApp}}.
#' 
#' @param data.dir A character string containing the directory to read from.
#' @param type Character vector containing either num or char. 
#' @return A data frame containg the full design and all the responses of the 
#'   combined data files that where found. Different files are indicated by an
#'   ID variable.
#' @export
LoadData <- function(data.dir, type) {
  # ErrorS
  if(!type %in% c("num", "char")){
    stop("type must be either num or char")
  }
  if (!dir.exists(data.dir)) {
    stop("Directory data.dir does not exist")
  }
  error <- character(0)
  if(identical(list.files(data.dir, full.names = TRUE, pattern = type), error)){
    stop('No files of the specified type in data.dir')
  } 
  # Read all files into list
  files <- list.files(data.dir, full.names = TRUE, pattern = type)
  data <- lapply(files, utils::read.table, stringsAsFactors = FALSE, sep = '\t', header = T) 
  # matrix
  id_rows <- sapply(data, nrow)
  id <- rep(1:length(id_rows), id_rows)
  data <- lapply(data, as.data.frame)
  # Bind together
  data <- do.call(rbind, data)
  data <- cbind("ID" = id, data)
  return(data)
}

