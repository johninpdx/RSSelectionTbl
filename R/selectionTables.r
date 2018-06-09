#_______________________________________________________________________
# R script selectionTables.r                                           #
# written by Tom A.B. Snijders                                         #
# March 25, 2018
# Customized and roxygenized by J.M. Light
# June 7, 2018
#_______________________________________________________________________

#_______________________________________________________________________
# These are functions for constructing and presenting selection tables #
# for the interpretation of results for network dynamics               #
# obtained with the RSiena program.                                    #
# Also consult the manual!                                             #
#                                                                      #
# FUNCTIONS:                                                           #
## selectionMatrix <- function(x,xd,name,vname,levls, levls.alt=levls) #
#_______________________________________________________________________
#                                                                      #
# which creates a matrix containing the selection table for            #
# siena data set xd, sienaFit or sienaMeta object x,                   #
# actor covariate vname (should be a character string),                #
# dependent variable name (also a character string),                   #
# levels for ego levls, levels for alter levls.alt.                    #
# Mostly levls.alt will be the same as levls,                          #
# in which case it does not have to be specified.                      #
#                                                                      #
## selectionTable.se                                                   #
#_______________________________________________________________________
# Computes standard error of a linear combination of elements          #
# of the selection table.                                              #
#                                                                      #
# selectionTableWithMax                                                #
#_______________________________________________________________________
# Creates a data frame of the selection table together with            #
# the maximum per ego of the selection function over alters.           #
#                                                                      #
# selectionTable.norm                                                  #
# Computes the location of the social norm and its standard error      #
# for the model discussed in Snijders & Lomi (2018).                   #
#                                                                      #
# selectionTable.plot                                                  #
#_______________________________________________________________________
# Constructs a plot of the selection table using ggplot2.              #
#                                                                      #
# If these functions are used for a sienaMeta object x                 #
# (as opposed to a sienaFit object),                                   #
# xd should be one of the individual data sets used for creating x,    #
# (a single-group instead of multi-group data set)                     #
# with an average ('representative') value                             #
# for the mean of variable <name>.                                     #
# NOTE: a 'data set' is an object of class 'siena'   (a Siena data     #
#       object)                                                        #
# _____________________________________________________________________

# ------------------------------- selectionTable ------------------------------
#'Creates a selection (or influence) table (matrix)
#'
#'@param x sienaFit: Results from a single group analysis
#'@param xd siena: the RSiena  Data Object or Data Set that
#'  was used to generate x
#'@param name character: the name of the endogenous variable
#'  of interest (a network or behavior variable)
#'@param vname character: the actor variable name (i.e. predictor
#'  variable)
#'@param levls numeric: number of ego levels in the variable given by vname
#'@param levls.alt numeric: number of alter levels in vname (usually = levls, in
#'  which case you do not need to include it in the call)
#'@silent logical: If set to TRUE, does not return feedback indicating which
#'  parameters were found, etc. Default FALSE.
#'@return matrix: the selection (or influence) table with alter values of
#'  vname across the columns, and ego values down the rows.
#'@export
selectionTable <- function(x, xd, name, vname, levls, levls.alt=levls){
  #'@import RSiena
  df <- selectionTable.basis(x, xd, name, vname, levls, levls.alt)$df
  df$ego <- as.character(df$vego)
  df$valter <- as.numeric(as.character(df$valter))
  df$select <- as.numeric(as.character(df$select))
  df
}

selectionMatrix <- function(x,xd,name,vname,levls, levls.alt=levls){
  # Creates a matrix containing the selection table for
  # siena data set xd, sienaFit object x,
  # actor covariate vname (should be a character string),
  # dependent variable name (also a character string),
  # levels for ego levls, levels for alter levls.alt.
  stab <- selectionTable(x,xd,name,vname,levls,levls.alt)
  mat <- matrix(as.numeric(as.character(stab$select)),
                length(unique(stab$vego)),
                length(unique(stab$valter)), byrow=TRUE)
  colnames(mat) <- levls.alt
  rownames(mat) <- levls
  mat
}

# ------------------------------- selectionTableWithMax ------------------------------
#'Creates a dataframe with the selection table, plus max values
#'across alter values
#'
#'@details  Creates a data frame including the selection table for
#' siena data set xd, sienaFit object x,
#' actor covariate vname,
#' dependent variable name,
#' levels for ego levls, levels for alter levls.alt,
#' to which is appended the values for the maximum across alter,
#' where alter ranges over levls.alt (a grid with 200 points between
#' minimum and maximum is used).
#'
#'@param x sienaFit: Results from a single group analysis
#'@param xd siena: the RSiena  Data Object or Data Set that
#'  was used to generate x
#'@param name character: the name of the endogenous variable
#'  of interest (a network or behavior variable)
#'@param vname character: the actor variable name (i.e. predictor
#'  variable)
#'@param levls numeric: number of ego levels in the variable given by vname
#'@param levls.alt numeric (=levls): number of alter levels in vname (usually = levls, in
#'  which case you do not need to include it in the call)
#'@discrete logical(TRUE): If set to TRUE, appropriately rounds integers
#'  to get means.
#'@return dataframe: MUST DOCUMENT THIS STRUCTURE WHEN YOU CAN FIGURE IT OUT!
#'  Apparently, there are two types of rows; if 'kind'=1, it's the selection
#'  table, and if 2, it's the max values.
#'@export
selectionTableWithMax <- function(x, xd, name, vname, levls, levls.alt=levls,
                                  discrete=TRUE){
#'@import RSiena
  dsign <- function(d){
    0.5*(ifelse(d > 0, 1, 0) + ifelse(d >= 0, 1, 0))
  }
  st <- selectionTable.basis(x, xd, name, vname, levls, levls.alt)
  df1 <- st$df
  df1$kind <- rep(1,dim(df1)[1])
  vtheta <- st$vtheta
  vmean <- st$vmean
  Delta <- st$Delta
  vsmean <- st$vsmean
  # f6 duplicates the information in coeffs in selectionTable.basis.
  # This duplication is undesirable as programming style, but I do this anyway.
  # When maintaining this script, f6 and coeffs in selectionTable.basis
  # must stay in line.
  f6 <- function(ve,va){
    contr <- vtheta[1]*(va-vmean) + vtheta[2]*(va-vmean)*(va-vmean) +
      vtheta[3]*(ve-vmean) +
      vtheta[4]*(ve-vmean)*(ve-vmean) + vtheta[5]*(ve-vmean)*(va-vmean) +
      vtheta[6]*(Delta-(abs(va-ve)/Delta)-vsmean) +
      vtheta[7]*(va-ve) + vtheta[8]*(va-ve)*(va-ve) +
      vtheta[9]*dsign(ve-va)
    names(contr) <- 'contribution'
    contr}
  # vtheta contains the parameter values in x
  K <- length(levls)
  KA <- length(levls.alt)
  valter <- rep(levls.alt,K)
  vego <- rep(levls,each=KA)
  fact <- 1:K
  ego <- factor(rep(fact,each=KA))
  # Now calculate the maximum.
  minv <- min(levls.alt)
  maxv <- max(levls.alt)
  gridv <- minv + (0:200)*((maxv-minv)/200)
  # vsmean is not important because it is an additive constant
  if (discrete){
    altm <- as.integer(round(minv)):as.integer(round(maxv))
    egom <- rep(levls[1], length(altm)) # dummy
    maxm <- sapply(altm, function(x){max(f6(x, gridv))})
  } else {
    egom <- rep(NA, K*4)
    altm <- rep(NA, K*4)
    maxm <- rep(NA, K*4)
    for (i in 1:K){
      for (j in 1:4){
        egom[4*(i-1) + j] <- levls[i]
        x1 <- ifelse(i <= 1, levls[1], 0.5*(levls[i-1]+levls[i]))
        x2 <- ifelse(i >= K, levls[K], 0.5*(levls[i+1]+levls[i]))
        dd <- x2-x1
        alterij <- ((j-1)/3)*x2 + (1 - (j-1)/3)*x1
        altm[4*(i-1) + j] <- alterij
        #			maxm[4*(i-1) + j]  <- max(f6(alterij, gridv))
        maxm[4*(i-1) + j]  <- max(f6(levls[i], gridv))
      }
    }
  }
  df1$ego <- as.character(df1$vego)
  df2 <- data.frame(ego=egom,vego=egom,valter=altm,select=maxm,kind=2)
  df <- rbind(df1,df2)
  df$ego <- as.character(df$vego)
  #	df$valter <- as.numeric(as.character(df$valter)) # superfluousc
  #	df$select <- as.numeric(as.character(df$select)) # superfluous
  df
}

# ------------------------------- selectionTable.se ------------------------------
#'Calculates standard errors for a linear combination
#' of elements of a selection table
#'
#'@details The linear combination for which the standard error is computed is
#' sum_{h,k} ww[h,k] * selectionTable[h,k].
#'
#'@param x sienaFit: Results from a single group analysis
#'@param xd siena: the RSiena  Data Object or Data Set that
#'  was used to generate x
#'@param name character: the name of the endogenous variable
#'  of interest (a network or behavior variable)
#'@param vname character: the actor variable name (i.e. predictor
#'  variable)
#'@param levls numeric: number of ego levels in the variable given by vname
#'@param ww matrix: same format as the selection table of interest, containing
#'  the coefficients of the linear combination
#'@param levls.alt numeric (=levls): number of alter levels in vname
#' (usually = levls, in which case you do not need to include it in the call)
#'@return matrix: the selection (or influence) table se's, with alter values
#'  of 'vname' across the columns, and ego values down the rows.
#'@export
selectionTable.se <- function(x, xd, name, vname, levls, ww, levls.alt=levls){
#'@import RSiena
  if (class(x) == "sienaBayesFit"){
    stop('This function does not work for sienaBayesFit objects.\n')
  }
  if (!all(dim(ww) == c(length(levls), length(levls.alt)))){
    stop('Dimension of ww should be the lengths of levls and levls.alt.\n')
  }
  cat("Requested cell weights (row = ego; col = alter):\n")
  print(cbind(which(ww != 0, arr.ind=TRUE),
              value=ww[which(ww != 0, arr.ind=TRUE)]))
  sb <- selectionTable.basis(x, xd, name, vname, levls, levls.alt, silent=TRUE)
  veff.r <- sb$veff[sb$veff != 0] # effects included in x
  cth <- x$covtheta[veff.r, veff.r] # their covariance matrix
  cth[is.na(cth)] <- 0 # if any effects were fixed, they have NA in covtheta
  lincomb <- sum(as.vector(t(ww)) * sb$df["select"])
  # the desired linear combination of cell values of the selection matrix
  wt <- colSums(as.vector(t(ww)) * sb$coeff)
  # Vector of weights for the linear combination of parameters
  names(wt) <- names(sb$veff)
  wt.r <- wt[names(veff.r)] # should be restricted to what is in the model
  cat("Parameter estimates and their resulting weights are \n")
  print(rbind('param'=sb$vtheta[sb$veff != 0], 'weight'=wt.r))
  cat("Linear combination of cells of selection matrix", round(lincomb,4),
      "\nStandard error\n")
  print(sqrt((wt.r %*% cth %*% wt.r)[1,1]))
}

# ------------------------------- selectionTable.norm--------------------------
#'Calculates the location of the social norm and its standard error
#' in the attraction function
#'
#'@details Assumes an attraction function with terms
#' 'altX','altSqX','egoX','egoSqX','diffSqX'. Note that the attraction function
#'  should not include the 'egoXaltX' effect. See Snijders & Lomi(c2018),
#'  "Beyond Homophily", for more info.
#'
#'@param x sienaFit: Results from a single group analysis
#'@param xd siena: the RSiena  Data Object or Data Set that
#'  was used to generate x
#'@param name character: the name of the endogenous variable
#'  of interest (a network or behavior variable)
#'@param vname character: the actor variable name (i.e. predictor
#'  variable)
#'@return list: [1] The calculated norm, [2] the SE of the norm
#'@export
selectionTable.norm <- function(x,xd,name,vname){
#'@import RSiena
  stab <- selectionTable.basis(x,xd,name,vname,1:2)
  # the 1:2 is arbitrary, not used
  vtheta <- stab$vtheta
  vmean <- stab$vmean
  veff <- stab$veff
  if (vtheta["egoXaltX"] != 0){
    cat('Warning: effect of egoXaltX is ',vtheta["egoXaltX"],
        ', not equal to 0; \n')
    cat(' this function does not apply to such a sienaFit object.\n')
  }
  # calculate gradient
  grad <- matrix(0, length(x$theta), 1)
  grad[veff["altX"],1] <- -1/(2*vtheta["altSqX"])
  grad[veff["altSqX"],1] <- vtheta["altX"]/(2*vtheta["altSqX"]*vtheta["altSqX"])
  covtheta <- x$covtheta
  covtheta[is.na(covtheta)] <- 0
  if (vtheta["altSqX"] > 0){
    cat('Warning: the coefficient of alter squared is positive.\n')
    cat('This means the extremum of the attraction function is a minimum,\n')
    cat('and cannot be interpreted as a social norm.\n')
    flush.console()
  }
  list(vnorm = vmean - (vtheta["altX"]/(2*vtheta["altSqX"])),
       se.vnorm = sqrt(t(grad) %*% covtheta %*% grad))
}

# ------------------------------- selectionTable.plot ------------------------------
#'Plots the selection table associated with a sienaFit object & other spex
#'
#'@details
#'
#'@param x sienaFit: Results from a single group analysis
#'@param xd siena: the RSiena  Data Object or Data Set that
#'  was used to generate x
#'@param name character: the name of the endogenous variable
#'  of interest (a network or behavior variable)
#'@param vname character: the actor variable name (i.e. predictor
#'  variable)
#'@param levls numeric: number of ego levels in the variable given by vname
#'@param levls.alt numeric (=levls): number of alter levels in vname
#' (usually = levls, in which case you do not need to include it in the call)
#'@param quad logical(=TRUE): If plotting a quadratic function, set to TRUE;
#'  For similarity effects, set to FALSE.
#'@param separation numeric(=0): Higher values give the curves more separation
#'  to improve readability. An advisable nonzero value is -0.01.
#'@param bw logical(=FALSE): If true, plot is black & white, otherwise color.
#'@param withMax logical (=FALSE): If TRUE, also plot the maximum over selection
#'  probability over values of alters as a dashed line.
#'@return A ggplot2 object
#'@export
selectionTable.plot <- function(x, xd, name, vname, levls, levls.alt=levls,
                                quad=TRUE, separation=0, bw=FALSE, withMax=FALSE){
#'@import RSiena
  if (withMax) {
    vselect00 <- selectionTableWithMax(x, xd, name, vname, levls, levls.alt)
    vselect <- subset(vselect00, kind==1)
    vselect.max <- subset(vselect00, kind==2)
  } else {
    vselect <- selectionTable(x, xd, name, vname, levls, levls.alt)
  }
  vr <- max(vselect$select) - min(vselect$select) # only for separation
  vselect$select <- vselect$select + separation*vr*as.numeric(factor(vselect$ego))
  if (bw) {
    sp <- ggplot(vselect, aes(valter, select, group=ego, linetype=ego))
  } else {
    sp <- ggplot(vselect, aes(valter, select, group=ego, colour=ego))
  }
  if (quad) {
    gs <- geom_smooth(size=1.2, span=3)
  } else {
    gs <- geom_line(size=1.2)
  }
  if (bw) {
    sp <- sp + geom_point() + gs + scale_linetype_manual(values =
                                                           c('solid',  'longdash','dashed', 'twodash', 'dotdash', 'dotted'))
  } else {
    sp <- sp + geom_point() + gs + scale_colour_hue()
  }
  ssp <- sp + scale_x_continuous(breaks=levls.alt) +
    theme(legend.key=element_blank())+
    labs(x=paste(vname),
         y=paste('selection function'),
         title=paste('Effect',vname,'on',name),
         colour=paste(vname,'\n ego\n', sep='')) +
    #				colour='age \n ego\n') +
    theme_grey(base_size=26, base_family="") +
    theme(legend.key.width=unit(1, "cm")) +
    theme(plot.title=element_text(hjust=0.5))
  if (withMax){
    ssp <- ssp + geom_point(data=vselect.max, size=2, shape=8)#, color='black')
    # alternatives
    #			sp <- sp + geom_smooth(data=vselect.max, size=1.2, span=3, linetype='dashed')
    #			sp <- sp + geom_line(data=vselect.max, size=2.5, linetype='solid')
  }
  ssp
}
#______________________________________________________________________________
## ---------------------      NOT EXPORTED    ---------------------------------
#______________________________________________________________________________
#----------------------- selectionTable.basis ---------------------------------
#'Makes calculations required by other functions ('basics'). Not exported.
#'
#'@param x sienaFit or sienaMeta: Results from SAOM
#'@param xd siena: the RSiena  Data Object or Data Set that
#'  was used to generate x (if class(x) is sienaMeta, this must be one
#'  of the Data Sets used to generate it, hopefully representative
#'  w/r/t vname, and it cannot be a multigroup data set)
#'@param name character: the name of the endogenous variable
#'  of interest (a network or behavior variable)
#'@param vname character: the actor variable name (i.e. predictor
#'  variable)
#'@param levls numeric: number of ego levels in the variable given by vname
#'@param levls.alt numeric: number of alter levels in vname (usually = levls, in
#'  which case you do not need to include it in the call)
#'@silent logical: If set to TRUE, does not return feedback indicating which
#'  parameters were found, etc. Default FALSE.
#'@return list: includes various calculations (predicted values for levels of
#'   vname, that kind of thing) needed for constructing the actual selection
#'   table.
#'@example
#' siena data set xd, sienaFit object x,
#' actor variable vname (should be a character string),
#' dependent variable name (also a character string),
#' levls a range of the actor variable
#'
#'vname <- '...'
#'name <- '...'
#'levls <- 1:5
#'
#'png(filename=paste("selectionTable_",vname,".png",sep=""), width=1000,height=800)
#'selectionTable.plot(x, xd, name, vname, levls)
#'graphics.off()
#'@export
selectionTable.basis <- function(x, xd, name, vname, levls, levls.alt=levls,
                                 silent=FALSE){
  #'@import RSiena
  dsign <- function(d){
    0.5*(ifelse(d > 0, 1, 0) + ifelse(d >= 0, 1, 0))
  }
  cat("Dependent network",name, "; actor variable",vname,".\n")
  depvar <- FALSE
  thevar <- xd$cCovars[[vname]] #Is vname a cCovar?
  if (is.null(thevar)) {thevar <- xd$vCovars[[vname]]} #Is vname a vCovar?
  if (is.null(thevar)) {
    # Is vname a dep var?
    depvar <- TRUE
    thevar <- xd$depvars[[vname]]
  }
  if (is.null(thevar)){stop(paste('There is no actor variable <',vname,'>.'))}
  if (depvar){# Then the mean is not stored as an attribute, so calc it
    means <- colMeans(thevar, na.rm=TRUE)
    vmean <- mean(means)
  } else {
    vmean <- attr(thevar, 'mean')
  }
  vsmean <- attr(thevar, 'simMean')
  # perhaps attr(thevar, 'simMean') fails
  # if there are more than one dependent network?
  # is this attribute then a vector of length >= 2?
  # Note that vsmean is used only if the model includes avSim or totSim.
  Delta  <- attr(thevar, 'range')
  replace0 <- function(k){ifelse(length(k)==0,0,k)}
  efNames <- c('altX','altSqX','egoX','egoSqX','egoXaltX','simX','diffX','diffSqX','higher')
  veff <- sapply(efNames, function(s)
  {replace0(which(x$requestedEffects$name == name &
                    x$requestedEffects$interaction1 == vname &
                    x$requestedEffects$include &
                    x$requestedEffects$shortName==s))})
  # veff gives the indicators of the effects in x
  taketheta <- function(k){ifelse(k==0,0,x$theta[k])}
  vtheta <- sapply(veff,taketheta)
  if (!silent){
    cat('Parameters found are\n')
    print(vtheta[vtheta != 0.0])
    cat('\n')
  }
  flush.console()
  # vtheta contains the parameter values in x
  K <- length(levls)
  KA <- length(levls.alt)
  valter <- rep(levls.alt,K)
  vego <- rep(levls,each=KA)
  fact <- 1:K
  ego <- factor(rep(fact,each=KA))
  # coeffs duplicates the information in f6 in selectionTableWithMax.
  # This duplication is undesirable, but I do this anyway.
  # When maintaining this script, coeffs and f6 in selectionTable.WithMax
  # must stay in line.
  coeffs <- matrix(NA, K*KA, length(efNames))
  coeffs[,1] <- (valter - vmean)
  coeffs[,2] <- (valter - vmean)*(valter - vmean)
  coeffs[,3] <- (vego - vmean)
  coeffs[,4] <- (vego - vmean)*(vego - vmean)
  coeffs[,5] <- (vego - vmean)*(valter - vmean)
  coeffs[,6] <- (1-(abs(valter - vego)/Delta)-vsmean)
  coeffs[,7] <- (valter - vego)
  coeffs[,8] <- (valter - vego)*(valter - vego)
  coeffs[,9] <- dsign(vego - valter)
  select <- coeffs %*% vtheta
  df <- data.frame(ego,vego,valter,select)
  list(df=df, veff=veff, vtheta=vtheta, vmean=vmean, vsmean=vsmean,
       Delta=Delta, coeffs=coeffs)
}
