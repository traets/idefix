

#' Transform coded choice set to Attribute level choice set
#'
#' Transforms a coded choice set into a choice set containing the attribute levels.
#' @param set A numeric matrix which represents a choice set. Each row is a profile.
#' @param lvl.names A list containing the values of each level of each attribute.
#' @param coding Type of coding used in the given set. See \code{\link[stats]{model.matrix}} 
#' @param intercept Logical argument indicating whether an intercept is included. The default is False.
#' @return A character matrix which represents the choice set.
#' @examples 
#' #choice set that is dummy coded
#' choice.set <- matrix(data = c(0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1), 
#' nrow=2, byrow = T) 
#' att.levels <- vector(mode="list", 3) #level names
#' att.levels[[1]] <- c("$50", "$75", "$100") #levels attribute 1
#' att.levels[[2]] <- c("2min", "15min", "30min") #levels attribute 2
#' att.levels[[3]] <- c("bad", "average", "good") #levels attribut 3
#' coding.type = "contr.treatment" #coding type is dummy coded 
#' #Transform
#' Present(set = choice.set, lvl.names = att.levels, coding = coding.type) 
#' @export
Present <- function(set, lvl.names, coding, intercept= FALSE) {
  #error handling
  codings.types<-c("contr.treatment", "contr.helmert", "contr.poly", "contr.sum", "none")
  if (!(coding %in% codings.types)) {
    stop("Coding argument is incorrect.")
  } 
  n.alts <- nrow(set) #number of alternatives
  n.att <- length(lvl.names) #number of attributes
  #create vector where each element denotes the number of levels for each attribute
  lvls <- numeric(n.att) 
  for (i in 1:n.att) { 
    lvls[i] <- length(lvl.names[[i]])
  }
  #generate all possible profiles coded and uncoded
  dc <- Profiles(lvls = lvls, coding = coding, intercept = intercept)
  d <- Profiles(lvls = lvls, coding = "none", intercept = intercept)
  #create new matrix for choice set with attribute level names 
  m <- matrix(data = NA, nrow = n.alts, ncol = n.att)
  #error handling
  if (ncol(set) != ncol(dc)) {
    stop("Number of columns of the set does not match expected number based on the other arguments.")
  }
  #for each alternative look for matching profile  
  for (i in 1:n.alts) {
    # if coded choice set, look for match in coded version first, then take uncoded equivalent.
    if (coding != "none") {
      lev.num <- d[as.numeric(which(apply(dc, 1, function(x) all(x == set[i, ])))), ]
      lev.num <- as.numeric(lev.num)
    } else {
      lev.num <- as.numeric(choice.set[i, ])
    }
    #error handling
    if (any(is.na(lev.num))) { 
      stop('The set does not match with the type of coding provided')
    }
    #for each attribute fill in the attribute level name
    for (c in 1:n.att) {
      m[i,c] <- lvl.names[[c]][lev.num[c]]
    }
  }
  return(m)
}



#' Transform responses
#'
#' Transforms character input responses to binary response vector.
#' @param resp String vector containing input responses
#' @param resp_options String vector containing all possible responses.
#' The response options should be specified in increasing order, starting with the neutral response (if included).
#' @param n_alts The number of alternatives per choice set.
#' @param neutral Logical value indicating whether a neutral option is provided or not. Default = TRUE.
#' @return A binary response vector.
#' @export
map_resp<-function(resp, resp_options, n_alts, neutral=T){

  map<-match(resp, resp_options)
  l<-list()

  for(i in 1:length(map)){

    if (neutral){
      l[[i]] <- rep(0, n_alts)
      l[[i]][map[i]-1]<-1
    }else{
      l[[i]] <- rep(0, n_alts)
      l[[i]][map[i]]<-1
    }
  }
  v<-unlist(l)
  return(v)
}


#' Load from dropbox
#'
#' Load design from dropbox map
#'@param inputDir A dropbox directory that contains the design.
#'@return the file in that directory (or concatenated files).
#'@export
loaddrop <- function(dir) {

  filesInfo <- drop_dir(dir)
  filePaths <- filesInfo$path
  data <- lapply(filePaths, drop_read_csv, stringsAsFactors = FALSE)
  # Concatenate all data together into one data frame
  data <- do.call(rbind, data)
  return(data)
}


#' Save to dropbox
#'
#' Save design and responses to dropbox map
#' @param des A design matrix in which each row is a profile.
#' @param Y A response vector.
#' @export
savedrop <- function(des, Y, dir, filename) {

  data<-cbind(des, Y)
  fileName <- sprintf("%s.csv", filename )
  filePath <- file.path(tempdir(), fileName)
  write.csv(data, filePath, row.names = FALSE, quote = TRUE)

  drop_upload(filePath, dest = dir)
}


#' Binary to discrete response matrix
#'
#' Transforms a matrix with binary choice data for each respondent (columns),
#' to a matrix with discrete values representing the choices.
#' @param y Matrix containing the binary choice data.
#' @param n_alts The number of alternatives per choice set.
#' @return A matrix with discrete values, indicating the choices.
#' @export
bindis<-function(y, n_alts){

  #divide into choice sets
  fun1<-function(x) split(x, ceiling(seq_along(x)/n_alts))
  YY<-apply(y,2,fun1)

  #warning 1
  for (i in 1:ncol(y)){
    if((length(unique(lengths(YY[[i]]))) == 1L) == FALSE){
      stop('length of Y vector does match expected length based on nr of alternatives')
    }
  }

  #index 1's
  fun<-function(x){xx<-(x==1); indexone <-which(xx, arr.ind = TRUE); if(length(indexone) > 1){stop('Multiple alternatives are chosen per choice set.
                                                                                                   The response data or the number of alternatives is probably incorrect.')}; return(indexone)}
  Y_Y<-list()
  for(r in 1: ncol(y)){
    Y_Y[[r]]<-lapply(YY[[r]], fun)
  }

  #rbind
  fun3<-function(x) {as.numeric(rbind(x))}
  ynom<-lapply(Y_Y, fun3)


  y_nom<-matrix(unlist(ynom), ncol=ncol(y), byrow= FALSE)

  return(y_nom)

  }





