####
# Model fitting
####

getMdl <- function(x, delta6=TRUE, modelName='E', nComp=2) {
  if (delta6) {
    w6 <- mean(x == 6)
    x <- x[x > 6]
  } else {
    w6 <- NA
  }
  list(w6=w6, mClu=densityMclust(x, G=nComp, modelNames=modelName))
}

getEligibleMdl <- function(x, delta6=TRUE, modelName='E', nComp=2, sigmaE) {
  mdl <- getMdl(x, delta6, modelName, nComp)
  # If sigma_error > min(sigma_comps) for the best model, prohibit this number of components:
  while (any(mdl$mClu$parameters$variance$sigmasq <= (sigmaE ^ 2))) {
    nComp <- setdiff(nComp, mdl$mClu$G)
    stopifnot(length(nComp) > 0)
    mdl <- getMdl(x, delta6, modelName, nComp)
  }
  return(mdl)
}


getParams <- function(mdl) {
  return(mdl$mClu$parameters)
}


getWComb <- function(mdl, type) {
  wComb <- mdl$mClu$parameters$pro
  if (!is.na(mdl$w6)) {
    wComb = (1 - mdl$w6) * wComb
    if (type == 'phenotype') {
      wComb[1] = mdl$w6 + wComb[1]
    } else if (type == 'component') {
      wComb <- c(mdl$w6, wComb)
    }
  }
  return(wComb)
}



####
# Printing
####


getSummaryTable <- function(mdl) {
  p <- getParams(mdl)
  smmry <- data.frame(mean=p$mean, sd=sqrt(p$variance$sigmasq), weight=p$pro)
  if (!is.na(mdl$w6)) {
    smmry$weight <- smmry$weight * (1 - mdl$w6)
    smmry <- rbind(data.frame(mean=6, sd=0, weight=mdl$w6), smmry)
  }
  smmry <- data.frame(component=1:nrow(smmry), smmry)
  colnames(smmry)[2:3] <- paste0(colnames(smmry)[2:3], '/mm')
  return(smmry)
}


printSummary <- function(mdl) {
  smmry <- getSummaryTable(mdl)
  cat('\n\n')
  print(kable(smmry, row.names=FALSE, digits=2))
  cat('\n\n')
}


####
# Calculations
####


getDensity <- function(mdl, type) {
  xx <- seq(from=5.5, to=40.5, length.out=1000)  # points at which density is computed and plotted
  densComp <- predict(mdl$mClu, xx, what='cdens')
  for (i in 1:mdl$mClu$G) {
    densComp[, i] = mdl$mClu$parameters$pro[i] * densComp[, i]  # scale densities s.t. they sum to 1
  }
  if (!is.na(mdl$w6)) {
    densComp = densComp * (1 - mdl$w6)
    idx6 <- 5.5 < xx & xx <= 6.5
    if (type == 'phenotype') {
      densComp[idx6, 1] = densComp[idx6, 1] + mdl$w6
    } else if (type == 'component') {
      densComp <- cbind(idx6 * mdl$w6, densComp)
    }
  }
  dens <- data.frame(diameter=xx, mix=rowSums(densComp), densComp)
  names(dens) <- getMdlColNames(mdl, type)
  return(dens)
}


jointProb <- function(mdl, sigma, xx, y) {
  stopifnot(length(y) == 1)
  
  pX <- predict(mdl$mClu, xx, what='dens')
  if (!is.na(mdl$w6)) {
    pX <- (1 - mdl$w6) * pX
  }
  pEps <- dnorm(y - xx, mean=0, sd=sigma)
  pXY <- pX * pEps
  
  if (!is.na(mdl$w6) & 5.5 < y & y <= 6.5) {
    pXY <- pXY + mdl$w6 * (5.5 < xx & xx <= 6.5)
  }
  
  return(pXY)
}



####
# Plotting
####


masterPlot <- function(x, mdl, cbp, sigma) {
  pltHist <- plotHistCount(x, cbp)
  pltHistAndFit <- plotHistAndFit(x, mdl, cbp, type='component')
  pltDens <- plotDensity(x, mdl, cbp, type='component')
  # pltCdf <- plotCdf(x, mdl, cbp, type='component')
  pltQQ <- myQQPlot(x, mdl)
  if (!is.na(cbp$S)) {
    tmp <- plotDeconvolution(x, mdl, cbp=cbp, sigma=c(sigma))
    pltDeconv <- tmp$plot
    uncertain <- tmp$uncertain
    errorProbs <- tmp$errorProbs
    pltThresh <- tmp$plotThresh
    dUncertainPlot <- tmp$dUncertainPlot
  }
  
  # Arrange figure for this report:
  leg <- getLegend(pltDens)
  plt <- arrangeGrob(arrangeGrob(pltHist + theme(legend.position="none"), pltDens + theme(legend.position="none"), pltQQ, heights=rep(1, 3)), arrangeGrob(pltDeconv), ncol=2)
  # plt <- arrangeGrob(arrangeGrob(pltHist + theme(legend.position="none"), leg, pltDens + theme(legend.position="none"), pltQQ, heights=c(1, 0.3, rep(1, 2))), arrangeGrob(pltDeconv), ncol=2)
  
  # Arrange figure for manuscript 1:
  pltMS1 <- arrangeGrob(pltHistAndFit + theme(legend.position="none"), pltQQ, pltThresh)
  
  # Arrange figure for manuscript 2:
  pltMS2 <- pltHist + theme(legend.position="none") + geom_rect(data=dUncertainPlot, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", alpha=0.3, inherit.aes=FALSE)
  # pltMS2 <- arrangeGrob(arrangeGrob(pltHist + theme(legend.position="none") + geom_rect(data=dUncertainPlot, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", alpha=0.3, inherit.aes=FALSE), pltThresh))
  
  
  return(list(plot=plt, uncertain=uncertain, errorProbs=errorProbs, plotMS1=pltMS1, plotMS2=pltMS2))
}


plotHistCount <- function(x, cbp) {
  return(ggplot(data.frame(diameter=x)) + geom_vline(xintercept=c(cbp$R, cbp$S) - 0.5, linetype=2, na.rm=TRUE) + geom_bar(aes(x=diameter), width=0.9) + xlab('observed diameter / mm') + theme(legend.position="top") + coord_cartesian(xlim=c(5.5, 40.5)))
}

plotHistDens <- function(x, cbp) {
  return(ggplot(data.frame(diameter=x)) + geom_vline(xintercept=c(cbp$R, cbp$S) - 0.5, linetype=2, na.rm=TRUE) + geom_bar(aes(x=diameter, y=(..count..)/sum(..count..)), width=0.9) + xlab('observed diameter / mm') + ylab('density') + theme(legend.position="top") + coord_cartesian(xlim=c(5.5, 40.5)))
}



plotHistAndFit <- function(x, mdl, cbp, type) {
  stopifnot(type %in% c('phenotype', 'component'))
  stopifnot(type != 'phenotype' | mdl$mClu$G == 2)
  
  dens <- getDensity(mdl, type)
  dens <- dens[, c('diameter', 'mix')]
  
  colsComp <- getColsComp(ncol(dens) - 2)
  return(plotHistDens(x, cbp) + geom_line(data=dens, aes(x=diameter, y=mix)))
}



plotDensity <- function(x, mdl, cbp, type) {
  stopifnot(type %in% c('phenotype', 'component'))
  stopifnot(type != 'phenotype' | mdl$mClu$G == 2)
  
  dens <- getDensity(mdl, type)
  dDens <- melt(dens, id.vars='diameter', variable.name=type, value.name='density')
  
  colsComp <- getColsComp(ncol(dens) - 2)
  return(plotHistDens(x, cbp) + geom_line(data=dDens, aes_string(x="diameter", y="density", color=type), size=1) + scale_colour_manual(values=colsComp))
  # return(ggplot(data.frame(diameter=x)) + geom_vline(xintercept=c(cbp$R, cbp$S) - 0.5, linetype=2, na.rm=TRUE) + geom_bar(aes(x=diameter, y=(..count..)/sum(..count..)), width=0.9) + geom_line(data=dDens, aes_string(x="diameter", y="density", color=type), size=1) + scale_colour_manual(values=colsComp) + xlab('observed diameter / mm') + ylab('density') + theme(legend.position="top"))
}



plotCdf <- function(x, mdl, cbp, type) {
  stopifnot(type %in% c('phenotype', 'component'))
  stopifnot(type != 'phenotype' | mdl$mClu$G == 2)
  
  p <- getParams(mdl)
  sd <- sqrt(ifelse(rep(mdl$mClu$modelName == 'E', mdl$mClu$G), rep(p$variance$sigmasq, mdl$mClu$G), p$variance$sigmasq))  # standard deviations of the different components
  xx <- seq(from=5.5, to=40.5, length.out=1000)  # points at which cdf is computed and plotted
  cdf <- sapply(1:mdl$mClu$G, function(i) pnorm(xx, mean=p$mean[i], sd=sd[i]))
  if (!is.na(mdl$w6)) {
    if (type == 'phenotype') {
      cdf[, 1] = apply(cbind(xx > 6, cdf[, 1]), 1, weighted.mean, w=c(mdl$w6, (1 - mdl$w6) * p$pro[1]))
    } else if (type == 'component') {
      cdf <- cbind((xx > 6) * 1, cdf)
    }
  }
  cdfMix <- apply(cdf, 1, weighted.mean, w=getWComb(mdl, type))
  tmp <- data.frame(diameter=xx, mix=cdfMix, cdf)
  names(tmp) <- getMdlColNames(mdl, type)
  dCdf <- melt(tmp, id.vars='diameter', variable.name=type, value.name='cdf')
  
  colsComp <- getColsComp(nlevels(dCdf[, 2]) - 1)
  return(ggplot(data.frame(diameter=x), aes(x=diameter)) + geom_vline(xintercept=c(cbp$R, cbp$S) - 0.5, linetype=2, na.rm=TRUE) + stat_ecdf(colour='grey35', size=1) + geom_line(data=dCdf, aes_string(x="diameter", y="cdf", color=type), size=1) + scale_colour_manual(values=colsComp) + xlab('observed diameter / mm') + ylab('cdf') + theme(legend.position ="top") + coord_cartesian(xlim=c(5.5, 40.5)))
}



plotDeconvolution <- function(x, mdl, cbp, sigma=c(0, 1)) {
  stopifnot(min(sigma) >= 0)
  sigma <- sigma[min(mdl$mClu$parameters$variance$sigmasq) > sigma ^ 2]
  stopifnot(mdl$mClu$parameters$variance$sigmasq > max(sigma ^ 2))
  
  #############
  # Top panel #
  #############
  
  # Compute density of true diameter:
  dDens <- data.frame()
  for (s in sigma) {
    thisMdl <- mdl
    thisMdl$mClu$parameters$variance$sigmasq <- mdl$mClu$parameters$variance$sigmasq - s ^ 2
    
    dens <- getDensity(thisMdl, type='component')
    dDens <- rbind(dDens, data.frame(diameter=dens$diameter, sigma=s, density=dens$mix))
  }
  dDens$sigma <- factor(dDens$sigma)
  
  # Plot:
  if (length(sigma) == 1) {
    cols <- 1
  } else {
    cols <- brewer.pal(length(sigma) + 1, 'Blues')[-1]
  }
  pltDens <- ggplot(data.frame(diameter=unname(x))) + geom_vline(xintercept=c(cbp$R, cbp$S) - 0.5, linetype=2, na.rm=TRUE) + geom_bar(aes(x=diameter, y=(..count..)/sum(..count..))) + geom_line(data=dDens, aes(x=diameter, y=density, color=sigma)) + scale_colour_manual(labels=c('true'), values=cols) + xlab('diameter / mm') + ylab('density') + guides(colour=guide_legend(title="")) + theme(legend.position="top")
  
  ####################
  # Bottom panel old #
  ####################
  
  #   # Compute p(X=>t|Y=y):
  #   t <- cbp$S - 0.5
  #   yy <- seq(from=5.5, by=0.1, to=40.5)
  #   pY <- predict(mdl$mClu, yy, what='dens')
  #   if (!is.na(mdl$w6)) {
  #     pY = mdl$w6 * (5.5 < yy & yy <= 6.5) + (1 - mdl$w6) * pY
  #   }
  #   bw <- min(0.01, sigma[sigma > 0] / 10)
  #   xx <- seq(from=t, by=bw, to=45)
  #   dP <- data.frame()
  #   for (s in sigma) {
  #     p <- rep(NA, length(yy))
  #     if (s == 0) {
  #       p <- t <= yy
  #     } else {
  #       thisMdl <- mdl
  #       thisMdl$mClu$parameters$variance$sigmasq <- mdl$mClu$parameters$variance$sigmasq - s ^ 2
  #       
  #       for (i in 1:length(yy)) {
  #         pXY <- jointProb(thisMdl, s, xx, yy[i])
  #         p[i] = sum(pXY) / pY[i]
  #       }
  #       p <- p * bw
  #     }
  #     dP <- rbind(dP, data.frame(diameter=yy, sigma=s, p=p))
  #   }
  #   dP$sigma <- factor(dP$sigma)
  #   
  #   # Zone of methodological uncertainty for largest sigma:
  #   rangeUncertain <- range(dP$diameter[dP$sigma == max(sigma) & 0.01 < dP$p & dP$p < 0.99])
  #   uncertainOut <- c(width=diff(rangeUncertain), mean(rangeUncertain[1] <= x & x <= rangeUncertain[2]), zmuMin=ceiling(rangeUncertain[1]), zmuMax=floor(rangeUncertain[2]))
  #   dUncertainPlot <- data.frame(xmin=rangeUncertain[1], xmax=rangeUncertain[2], ymin=-Inf, ymax=Inf)
  #   
  #   # Plot:
  #   pltThresh <- ggplot(dP, aes(x=diameter, y=p, color=sigma)) + geom_rect(data=dUncertainPlot, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=2, alpha=0.2, inherit.aes = FALSE) + geom_vline(xintercept=c(cbp$R, cbp$S) - 0.5, linetype=2, na.rm=TRUE) + geom_line() + scale_colour_manual(values=cols) + xlab('observed diameter / mm') + ylab(paste('P(S | observed diameter)')) + guides(colour=guide_legend(title=expression(sigma[epsilon]/"mm"))) + theme(legend.position="none") # + scale_y_continuous(breaks=seq(0, 1, 0.1)) + scale_x_continuous(breaks=seq(6, 40, 2))
  
  
  ####################
  # Bottom panel new #
  ####################
  
  # Compute p(X=>t|Y=y):
  yy <- seq(from=5.5, by=0.1, to=40.5)
  pY <- predict(mdl$mClu, yy, what='dens')
  if (!is.na(mdl$w6)) {
    pY = mdl$w6 * (5.5 < yy & yy <= 6.5) + (1 - mdl$w6) * pY
  }
  bw <- min(0.01, sigma[sigma > 0] / 10)
  xxR <- seq(from=1, by=bw, to=(cbp$R - 0.5))
  xxI <- seq(from=(cbp$R - 0.5), by=bw, to=(cbp$S - 0.5))
  xxS <- seq(from=(cbp$S - 0.5), by=bw, to=45)
  dP <- data.frame()
  for (s in sigma) {
    p <- rep(NA, length(yy))
    if (s == 0) {
      exit('code to be written')
      p <- t <= yy
    } else {
      thisMdl <- mdl
      thisMdl$mClu$parameters$variance$sigmasq <- mdl$mClu$parameters$variance$sigmasq - s ^ 2
      
      for (i in 1:length(yy)) {
        if (yy[i] < cbp$R - 0.5) {
          pXY <- jointProb(thisMdl, s, xxS, yy[i])
        } else if (yy[i] >= cbp$S - 0.5) {
          pXY <- jointProb(thisMdl, s, xxR, yy[i])
        } else {
          pXY <- 0
        }
        p[i] = sum(pXY) / pY[i]
      }
      p <- 1 - p * bw
    }
    dP <- rbind(dP, data.frame(diameter=yy, sigma=s, p=p))
  }
  dP$sigma <- factor(dP$sigma)
  
  # Zone of methodological uncertainty for largest sigma:
  rangeUncertain <- range(c(cbp$R - 0.5, cbp$S - 0.5, dP$diameter[dP$sigma == max(sigma) & dP$p < 0.99]))
  zmuMin <- ceiling(rangeUncertain[1])
  zmuMax <- ceiling(rangeUncertain[2])
  uncertainOut <- c(width=diff(rangeUncertain), mean(zmuMin <= x & x < zmuMax), zmuMin=zmuMin, zmuMax=zmuMax)
  dUncertainPlot <- data.frame(xmin=rangeUncertain[1], xmax=rangeUncertain[2], ymin=-Inf, ymax=Inf)
  
  # Plot:
  pltThresh <- ggplot(dP, aes(x=diameter, y=p, color=sigma)) + geom_rect(data=dUncertainPlot, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='grey35', alpha=0.3, inherit.aes = FALSE) + geom_vline(xintercept=c(cbp$R, cbp$S) - 0.5, linetype=2, na.rm=TRUE) + geom_line() + scale_colour_manual(values=cols) + xlab('observed diameter / mm') + ylab("forecast probability") + guides(colour=guide_legend(title=expression(sigma[epsilon]/"mm"))) + theme(legend.position="none") + coord_cartesian(ylim=c(0, 1)) # + scale_y_continuous(breaks=seq(0, 1, 0.1)) + scale_x_continuous(breaks=seq(6, 40, 2))
  
  
  ############################
  # Bottom panel alternative #
  ############################
  
  #   # Compute p(X=x|Y=>t):
  #   t <- cbp$S - 0.5
  #   
  #   xx <- seq(from=5.5, by=0.1, to=40.5)
  #   
  #   dP <- data.frame()
  #   for (s in sigma) {
  #     p <- rep(NA, length(xx))
  #     
  #     thisMdl <- mdl
  #     thisMdl$mClu$parameters$variance$sigmasq <- mdl$mClu$parameters$variance$sigmasq - s ^ 2
  #     
  #     pX <- predict(thisMdl$mClu, xx, what='dens')
  #     if (!is.na(thisMdl$w6)) {
  #       pX <- pX * (1 - thisMdl$w6)
  #     }
  #     pEps <- 1 - pnorm(t, mean=xx, sd=s)
  #     stopifnot(t > 6)
  #     pS <- 1 - cdfMclust(thisMdl$mClu, t)$y
  #     if (!is.na(thisMdl$w6)) {
  #       pS <- pS * (1 - thisMdl$w6)
  #     }
  #     p = pX * pEps / pS
  #     
  #     dP <- rbind(dP, data.frame(diameter=xx, sigma=s, p=p))
  #   }
  #   dP$sigma <- factor(dP$sigma)
  #   
  #   # Plot:
  #   pltThresh2 <- ggplot(dP, aes(x=diameter, y=p, color=sigma)) + geom_vline(xintercept=c(cbp$R, cbp$S) - 0.5, linetype=2, na.rm=TRUE) + geom_line() + scale_colour_manual(values=cols) + xlab('true diameter / mm') + ylab(paste('P(true diameter | pred=S)')) + guides(colour=guide_legend(title=expression(sigma[epsilon]/"mm"))) + theme(legend.position="none") # + scale_y_continuous(breaks=seq(0, 1, 0.1)) + scale_x_continuous(breaks=seq(6, 40, 2))
  
  
  ###################################
  # Print summary and combine plots #
  ###################################
  
  cat(sep='', '\n\nOfficial CBPs: ', cbp$R, ' mm and ', cbp$S, ' mm. ZMU: [', rangeUncertain[1], ' mm, ', rangeUncertain[2], ' mm].\n\n')
  
  # Combined figure:
  leg <- getLegend(pltDens)
  if (length(sigma) == 1) {
    pltDeconv <- arrangeGrob(pltDens + theme(legend.position="none"), pltThresh + theme(legend.position="none"), heights=rep(1, 2))
  } else {
    pltDeconv <- arrangeGrob(leg, pltDens + theme(legend.position="none"), pltThresh + theme(legend.position="none"), heights=c(0.2, rep(1, 2)))
  }
  # pltDeconv <- arrangeGrob(leg, pltDens + theme(legend.position="none"), pltThresh + theme(legend.position="none"), pltThresh2 + theme(legend.position="none"), heights=c(0.2, rep(1, 3)))
  
  
  #######################################
  # Probabilities of missclassification #
  #######################################
  
  s <- max(sigma)
  thisMdl <- mdl
  thisMdl$mClu$parameters$variance$sigmasq <- mdl$mClu$parameters$variance$sigmasq - s ^ 2
  
  # Compute probability of very major error (prediction=S, truth=R):
  #   bw <- 0.05
  #   xx <- seq(from=0, by=bw, to=(cbp$R - 0.5))
  #   yy <- seq(from=(cbp$S - 0.5), by=bw, to=45)
  #   pVME <- 0
  #   for (y in yy) {
  #     pXY <- jointProb(thisMdl, s, xx, y)
  #     pVME <- pVME + sum(pXY)
  #   }
  #   pVME <- pVME * bw ^ 2
  #   cat(sep='', '\n\nMethod 1: The estimated probability of a very major error is ', pVME, ' for sigma_E=', s, ' mm.\n\n')
  #   
  #   pX <- predict(thisMdl$mClu, xx, what='dens')
  #   if (!is.na(thisMdl$w6)) {
  #     pX <- pX * (1 - thisMdl$w6)
  #   }
  #   pEps <- 1 - pnorm(cbp$S - 0.5, mean=xx, sd=s)
  #   pVME <- sum(pX * pEps) * bw
  #   cat(sep='', '\n\nMethod 2: The estimated probability of a very major error is ', pVME, ' for sigma_E=', s, ' mm.\n\n')
  #   
  #   integrand <- function(x) {
  #     pX <- predict(thisMdl$mClu, x, what='dens')
  #     if (!is.na(thisMdl$w6)) {
  #       pX <- pX * (1 - thisMdl$w6)
  #     }
  #     pEps <- 1 - pnorm(cbp$S - 0.5, mean=x, sd=s)
  #     return(pX * pEps)
  #   }
  #   pVME <- integrate(integrand, -Inf, cbp$R - 0.5)
  #   cat(sep='', '\n\nMethod 3: The estimated probability of a very major error is ', pVME$value, ' ± ', pVME$abs.error, ' for sigma_E=', s, ' mm.\n\n')
  #   pVME <- pVME$value
  
  # More systematic approach:
  getIntegrand <- function(errBndry, errtype) {
    stopifnot(errtype %in% c('<', '>'))  # '<' / '>': error E has to be smaller / greater than errBndry
    integrand <- function(x) {
      pX <- predict(thisMdl$mClu, x, what='dens')
      if (!is.na(thisMdl$w6)) {
        pX <- pX * (1 - thisMdl$w6)
      }
      pEps <- pnorm(errBndry, mean=x, sd=s)
      if (errtype == '>') {
        pEps <- 1 - pEps
      }
      return(pX * pEps)
    }
    return(integrand)
  }
  
  errorProbs <- data.frame()
  
  for (case in c('official CBPs', 'CBP_S increased by 2mm', 'I for ZMU')) {
    for (errorType in c('major', 'very major')) {
      if (errorType == 'major') {
        if (case == 'I for ZMU') {
          xBndry <- c(max(cbp$S - 0.5, zmuMax - 0.5), Inf)
        } else if (case == 'CBP_S increased by 2mm') {
          xBndry <- c(cbp$S + 2 - 0.5, Inf)
        } else {
          xBndry <- c(cbp$S - 0.5, Inf)
        }
        epsDir <- '<'
        if (case == 'I for ZMU') {
          epsBndry <- min(cbp$R - 0.5, zmuMin - 0.5)
        } else {
          epsBndry <- cbp$R - 0.5
        }
      } else if (errorType == 'very major') {
        if (case == 'I for ZMU') {
          xBndry <- c(-Inf, min(cbp$R - 0.5, zmuMin - 0.5))
        } else {
          xBndry <- c(-Inf, cbp$R - 0.5)
        }
        epsDir <- '>'
        if (case == 'official CBPs') {
          epsBndry <- cbp$S - 0.5
        } else if (case == 'CBP_S increased by 2mm') {
          epsBndry <- cbp$S + 2 - 0.5
        } else if (case == 'I for ZMU') {
          epsBndry <- max(cbp$S - 0.5, zmuMax - 0.5)
        }
      } else {
        stopifnot(2)
      }
      pE <- integrate(getIntegrand(epsBndry, epsDir), xBndry[1], xBndry[2], rel.tol=1e-14, abs.tol=1e-14)
      # cat(sep='', '\n\n', pE$value, ' ± ', pE$abs.error, '.\n\n')
      errorProbs <- rbind(errorProbs, data.frame(case=case, errorType=errorType, errorProbability=pE$value, lower=(pE$value - pE$abs.error), upper=(pE$value + pE$abs.error)))
    }
  }
  
  printErrorProbs(errorProbs)
  
  return(list(plot=pltDeconv, uncertain=uncertainOut, errorProbs=errorProbs, plotThresh=pltThresh, dUncertainPlot=dUncertainPlot))
}


myQQPlot <- function(x, mdl) {
  p <- (1:length(x)) / length(x)
  if (!is.na(mdl$w6)) {
    pJump <- (1 - mdl$w6) * cdfMclust(mdl$mClu, 6)$y
    theorQuant <- rep(6, length(p))
    for (i in 1:2) {
      if (i == 1) {
        idx <- p < pJump
      } else if (i == 2) {
        idx <- (pJump + mdl$w6) < p
      }
      if (any(idx)) {
        theorQuant[idx] = quantileMclust(mdl$mClu, (p[idx] - mdl$w6 * (i == 2)) / (1 - mdl$w6))
      }
    }
  } else {
    theorQuant <- quantileMclust(mdl$mClu, p)
  }
  dQuantiles <- data.frame(sample=sort(x), theoretical=theorQuant)
  return(ggplot(dQuantiles, aes(x=theoretical, y=sample)) + geom_point(size=0.1) + geom_abline(colour='grey35') + xlab('theoretical\nquantile / mm') + ylab('sample\nquantile / mm') + coord_fixed())
}



getColsComp <- function(nComp) {
  #   if (nComp == 1) {
  #     colsComp <- brewer.pal(3, 'Set2')[2]
  #   } else if (nComp == 2) {
  #     colsComp <- brewer.pal(3, 'Set2')[2:3]
  #   } else {
  #     colsComp <- brewer.pal(nComp, 'Set2')
  #   }
  colsComp <- brewer.pal(max(nComp, 3), 'Set2')[1:nComp]
  colsComp <- c(1, colsComp)
}



getLegend<-function(myggplot){
  # This function returns the legend of a ggplot-objet.
  # Function taken from www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page-r-software-and-data-visualization
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

####
# Varia
####

getMdlColNames <- function(mdl, type) {
  if (type == 'phenotype') {
    nms <- c('diameter', 'both', 'non-wt', 'wt')
  } else if (type == 'component') {
    nms <- c('diameter', 'mix', 1:(mdl$mClu$G + !is.na(mdl$w6)))
  }
  return(nms)
}



formatErrorProbs <- function(errorProbs) {
  tmp <- format(errorProbs[, c('case', 'errorType', 'errorProbability')], scientific=TRUE, digits=2)
  tmp$errorType <- paste(tmp$errorType, 'error')
  tmp <- spread(tmp, errorType, errorProbability)[c(3, 1, 2), ]
  colnames(tmp)[1] <- 'basis for prediction'
  return(tmp)
}

printErrorProbs <- function(errorProbs) {
  tmp <- formatErrorProbs(errorProbs)
  cat('\n\n')
  print(kable(tmp, row.names=FALSE, caption="Probabilities of misclassification errors."))
  cat('\n\n')
}




