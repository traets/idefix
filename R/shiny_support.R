#' Shiny application to generate a discrete choice survey.
#' 
#' This function starts a shiny application which puts choice sets on screen and
#' saves the responses. The complete choice design can be provided in advance,
#' or can be generated sequentially adaptively, or will be a combination of
#' both.
#' 
#' @param alts
#' @param atts 
#' @param n.total
#' @param lvl.names
#' @param buttons.text
#' @param intro.text
#' @param end.text
#' @inheritParams Decode
#' @inheritParams Modfed
#' @inheritParams Profiles
#' 
#' @export 
SurveyApp <- function(des = NULL, n.total, alts, atts, lvl.names, alt.cte, coding,  c.lvls = NULL, p.mean = NULL, p.cov = NULL, cand.set = NULL, m = 6, buttons.text, intro.text, end.text) {
  
  # Libraries
  require(shiny)
  
  # Initialize 
  survey <<- vector(mode = "list")
  y.bin <- vector("numeric")
  resp  <- vector("character")
  n.atts <- length(atts)
  n.alts <- length(alts)
  buttons <- NULL
  sn <- 0
  
  if (is.null(des)) {
    n.init <- 0
  } else {
    n.init <- nrow(des) / n.alts 
    bs <- seq(1, (nrow(des) - n.alts + 1), n.alts)
    es <- c((bs - 1), nrow(des))[-1] 
    # Error handling
    if (length(bs) != n.init) {
      stop("The number of design rows does not match the number of alternatives times the number of sets.")
    }
    if (n.total > n.init) {
      if (any(c(is.null(p.mean), is.null(p.cov), is.null(cand.set)))) {
        stop("When n.total is larger than the number of sets in des, arguments p.mean, p.cov, and cand.set should be specified")
      }
    }
  }
  
  shinyApp(
    ### User interface
    ui <- fluidPage(
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
          # Take set 
          # for initial sets 
          if (sn <= n.init) {
            set <- design[bs[sn] : es[sn], ]
          } else {
            # for adaptive sets
            # if First set
            if (sn == 1) {
              # Draw samples from prior
              s <- MASS::mvrnorm(n = 50, mu = p.mean, Sigma  = p.cov)
              w <- rep(1, nrow(s)) / nrow(s)
              # From second set
            } else {
              # Draw samples from updated posterior
              sam <- ImpsampMNL(prior.mean = p.mean, prior.covar = p.cov, des = des, n.alts = n.alts, y = y.bin, m = m)
              s <- sam$samples
              w <- sam$weights
            }
            # Select new set based on KL info
            set <- SeqKL(cand.set = cand.set, n.alts = n.alts, par.draws = s, alt.cte = alt.cte, weights = w)$set
            # Design storage
            if (sn == 1) { 
              rownames(set) <- rownames(set, do.NULL = FALSE, prefix = paste(paste("set", sn , sep = ""), "alt", sep = "."))
              des <<- set
            } else {
              rownames(set) <- rownames(set, do.NULL = FALSE, prefix = paste(paste("set", sn , sep = ""), "alt", sep = "."))
              des <<- rbind(des, set)
            }
          }
          # Transform coded set to attribute level character set.
          choice.set <- Decode(set = set, lvl.names = lvl.names, coding = coding, alt.cte = alt.cte, c.lvls = c.lvls)
          choice.set <- t(choice.set[ , 1:n.atts])
          # Fill in attribute names and alternatives names
          colnames(choice.set) <- alts
          rownames(choice.set) <- atts
          #store uncoded choice set
          if (sn == 1) {
            choice.sets <<- choice.set
          } else {
            choice.sets <<- rbind(choice.sets, choice.set)
          }
          #return choice set and design
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
          survey$bin.responses <<- y.bin
          survey$responses <<- resp
          survey$design <<- des
          survey$survey <<- choice.sets
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
      # Introtext
      output$intro <- renderText(intro.text)
      observeEvent(input$OK, {
        output$intro <- renderText("")
      })
      # End of survey
      observeEvent(input$OK, {
        # display end text 
        if (input$OK > n.total) {
          # display end text 
          output$end <- renderText(end.text)
        }
        # Quit application 
        if (input$OK > (n.total + 1)) {
          # End app
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
#' @export
Decode <- function(set, lvl.names, coding, alt.cte, c.lvls = NULL) {
  # Delete alt.cte's 
  contins <- which(alt.cte == 1)
  if( !length(contins) == 0) {
    set <- set[, -length(contins)]
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
#' @export
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




