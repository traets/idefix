#' data transformation.
#'
#' Transforms the data into the neccesary format in order to use the bayesm package, ChoiceModelR package or the RSGHB package.
#' @param pkg Indicates the required package or estimation (1=bayesm, 2=ChoiceModelR, 3=RSGHB, 4=Mixed Probit estimation (Bayesm))
#' @param des A design matrix in which each row is a profile.
#' @param Y A numeric matrix. Each column represents a respondent, there are n_sets rows (discrete choice data), or n_sets*n_alts rows (binary data).
#' @param n_alts Numeric value indicating the number of alternatives per choice set.
#' @param n_sets Numeric value indicating the number of choide sets.
#' @param n_resp Numeric value indicating the number of respondents.
#' @param n_beta Numeric value indicating the number of model parameters that needs to be estimated.
#' @param bin Logical value indicating whether the reponse matrix contains binary data (TRUE) or not (FALSE).
#' @return The data ready to be used by the specified package.
#' @export
datatrans<-function(pkg, des, Y, n_alts, n_sets, n_resp, n_beta, bin){


  #transform data if binary
  if(bin){Y<-bindis(y=Y, n_alts = n_alts)}

  if(pkg==1){

    library(bayesm)
    print("The dataset is ready to be used for bayesm package")

    bayesmin=function(des, Y, n_alts ,n_sets ,n_resp){

      #transform Y and des  into matrix form
      Y=t(t(Y))
      des =t(t(des))

      lgtdata=NULL

      ni=rep(n_sets,n_resp)
      csa=n_sets*n_alts

      for (i in 1:n_resp)
      {
        #obtain Y
        ychoice=NULL
        ybeg=n_sets*(i-1)+1
        yend=n_sets*i
        for(c in 1:n_sets){ychoice[1:n_sets]=Y[ybeg:yend]}

        #transform des into dataframe
        xmat=NULL
        xbeg=csa*(i-1)+1
        xend=csa*i
        xmat[[i]]=des[xbeg:xend, ]

        lgtdata[[i]]=list(y=ychoice, X=xmat[[i]])
      }


      #the bayesmin function returns a list of 2
      return(bayesmdata=list(p=n_alts, lgtdata=lgtdata))

    }

    return(bayesmin(des, Y, n_alts, n_sets, n_resp))


  }

  else if(pkg==2){

    library(ChoiceModelR)
    print("The dataset is ready to be used for ChoiceModelR package")

    choicemodelrin=function(des, Y, n_alts, n_sets, n_resp){

      #transform Y and des  into matrix form
      Y=t(t(Y))
      des =t(t(des))

      set=rep(1:n_sets, each = n_alts, times=n_resp)
      id=rep(1:n_resp, each=n_sets*n_alts)
      alt=rep(c(1:n_alts), n_sets*n_resp)

      initialmat=t(rbind(id, set, alt))

      xmat=cbind(initialmat, des)

      #make choice columns
      newchoice=Y
      zeromat=matrix(0, n_sets*n_resp, n_alts-1)
      choicemat=cbind(newchoice, zeromat)

      #This is the final Y column representing choice
      choicecol=matrix(c(t(choicemat)))

      return(choicemodelrdata=cbind(xmat, choicecol))

    }

    return(choicemodelrin(des, Y, p, n_sets, n_resp))

  }

  else if (pkg==3){

    library(RSGHB)
    print("The dataset is ready to be used for RSGHB package")


    rsghbin=function(des, Y, n_alts, n_sets, n_resp,n_beta){

      Y=t(t(Y))
      des =t(t(des ))
      n_beta=ncol(des)

      rsghbid=rep(1:n_resp, each=n_sets)
      ncs=rep(1:n_sets, times=n_resp)
      rsghbinitialmat=t(rbind(rsghbid, ncs))

      #attribute matrix
      rsghbattrmat=NULL
      indset=n_sets*n_resp
      for (cs in 1:indset){
        beg=n_alts*(cs-1)+1
        end=n_alts*cs

        xtemp=NULL
        for(col in 1:n_beta){
          xtemp=cbind(xtemp, t(des [beg:end, col]))
        }
        rsghbattrmat=rbind(rsghbattrmat, xtemp)
      }

      RSGHBchoice=Y

      RSGHBdata=data.frame(cbind(rsghbinitialmat, rsghbattrmat, RSGHBchoice))
      colnames(RSGHBdata)[[1]]="ID"
      colnames(RSGHBdata)[[2]]="Choice Set"
      cy=ncol(RSGHBdata)
      colnames(RSGHBdata)[[cy]]="Choice"

      return(RSGHBdata)

    }

    return(rsghbin(des, Y, p, n_sets, n_resp, n_beta))

  }

  else if (pkg==4){
    print("The dataset is ready to be used for Mixed Probit Estimation")
    library(bayesm)

    mxpin=function(des, Y, n_alts, n_sets, n_resp, n_beta){

      des =t(t(des))
      Y=t(t(Y))

      ynum=nrow(Y)
      yind=NULL
      for (i in 1:ynum){
        zerotemp=matrix(0, n_alts, 1)
        index=Y[i, ]
        zerotemp[index]=1
        yind=rbind(yind,zerotemp)
      }

      Y= array(t(yind), dim=c(n_alts, n_sets, n_resp))
      des= array(t(des), dim=c(n_beta, n_alts, n_sets, n_resp))

      return(Data=list(y=Y, X=des, nlgt=n_resp, nset=n_sets, n_alts=n_alts, nbeta=n_beta))

    }

    return(mxpin(des, Y, n_alts, n_sets, n_resp, n_beta))

  }

  else{

    return(print("please specify: 1-bayesm, 2-ChoiceModelR, 3-RSGHB, 4-Mixed Probit") )

  }

}

#' Response generation
#'
#' Function to generate responses given parameter values and a design matrix according to MNL model.
#' @param par Vector containing parameter values.
#' @param set A numeric matrix which represents a choice set. Each row is a profile.
#' @param bin Indicates whether the returned value should be a binary vector or a discrete value which denotes the chosen alternative.
#' @return Binary response vector or discrete value indicating the choosing alternative.
#' @export
respond<-function (par, set, n_alts, bin=TRUE){

  par<-as.matrix(par)
  d <- as.matrix(set)

  #prob
  U <- d %*% t(par)
  expU <- exp(U)
  p <- expU/sum(expU)

  #choice
  choice<-findInterval(x=runif(1), vec=c(0,cumsum(p)))

  Y<-rep(0,length(p))
  Y[choice] <- 1

  #return
  ifelse(bin, return(Y), return(choice))

}

#roxygen2::roxygenise()
