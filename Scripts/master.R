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
  stopifnot(sigmaE > 0)
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



getWComb <- function(mdl) {
  wComb <- mdl$mClu$parameters$pro
  if (!is.na(mdl$w6)) {
    wComb <- c(mdl$w6, (1 - mdl$w6) * wComb)
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



####
# Calculations
####



getDensity <- function(mdl, xx) {
  densComp <- predict(mdl$mClu, xx, what='cdens')
  for (i in 1:mdl$mClu$G) {
    densComp[, i] = mdl$mClu$parameters$pro[i] * densComp[, i]  # scale densities s.t. they sum to 1
  }
  if (!is.na(mdl$w6)) {
    idx6 <- 5.5 < xx & xx <= 6.5
    densComp <- cbind(idx6 * mdl$w6, densComp * (1 - mdl$w6))
  }
  dens <- data.frame(diameter=xx, mix=rowSums(densComp), densComp)
  names(dens) <- getMdlColNames(mdl)
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



getForcastProbability <- function(mdl, sigma, cbp, yy) {
  # Computes forecast probability, i.e. the probability tht no major or very major error occurs:
  pY <- predict(mdl$mClu, yy, what='dens')
  if (!is.na(mdl$w6)) {
    pY = mdl$w6 * (5.5 < yy & yy <= 6.5) + (1 - mdl$w6) * pY
  }
  bw <- min(0.01, sigma / 10)
  xxR <- seq(from=1, by=bw, to=(cbp$R - 0.5))
  xxI <- seq(from=(cbp$R - 0.5), by=bw, to=(cbp$S - 0.5))
  xxS <- seq(from=(cbp$S - 0.5), by=bw, to=45)
  p <- rep(NA, length(yy))
  
  thisMdl <- mdl
  thisMdl$mClu$parameters$variance$sigmasq <- mdl$mClu$parameters$variance$sigmasq - sigma ^ 2
  
  for (i in 1:length(yy)) {
    if (yy[i] < cbp$R - 0.5) {
      pXY <- jointProb(thisMdl, sigma, xxS, yy[i])
    } else if (yy[i] >= cbp$S - 0.5) {
      pXY <- jointProb(thisMdl, sigma, xxR, yy[i])
    } else {
      pXY <- 0
    }
    p[i] = sum(pXY) / pY[i]
  }
  p <- 1 - p * bw
  data.frame(diameter=yy, p=p)
}

####
# Plotting
####


masterPlot <- function(x, mdl, cbp, sigma) {
  pltHist <- plotHistCount(x, cbp)
  pltHistAndFit <- plotHistAndFit(x, mdl, cbp)
  pltDens <- plotDensity(x, mdl, cbp)
  pltQQ <- myQQPlot(x, mdl)
  if (!is.na(cbp$S)) {
    tmp <- plotDeconvolution(x, mdl, cbp=cbp, sigma=c(sigma))
    pltDeconv <- tmp$plot
    uncertain <- tmp$uncertain
    errorProbs <- tmp$errorProbs
    pltThresh <- tmp$plotThresh
    dUncertainPlot <- tmp$dUncertainPlot
  }
  
  # Arrange figure for report:
  leg <- getLegend(pltDens)
  plt <- arrangeGrob(arrangeGrob(pltHist + theme(legend.position="none"), pltDens + theme(legend.position="none"), pltQQ, heights=rep(1, 3)), arrangeGrob(pltDeconv), ncol=2)
  
  # Arrange figure for manuscript 1:
  pltMS1 <- arrangeGrob(pltHistAndFit + theme(legend.position="none"), pltQQ, pltThresh)
  
  # Arrange figure for manuscript 2:
  pltMS2 <- pltHist + theme(legend.position="none") + geom_rect(data=dUncertainPlot, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", alpha=0.3, inherit.aes=FALSE)
  
  
  return(list(plot=plt, uncertain=uncertain, errorProbs=errorProbs, plotMS1=pltMS1, plotMS2=pltMS2))
}



plotHistCount <- function(x, cbp) {
  postprocessPlot(ggplot(data.frame(diameter=x)) + geom_bar(aes(x=diameter), width=0.9), cbp)
}



plotHistDens <- function(x, cbp) {
  postprocessPlot(ggplot(data.frame(diameter=x)) + geom_bar(aes(x=diameter, y=(..count..)/sum(..count..)), width=0.9) + ylab('density'), cbp)
}



plotHistAndFit <- function(x, mdl, cbp) {
  dens <- getDensity(mdl, seq(from=5.5, to=40.5, length.out=1000))
  dens <- dens[, c('diameter', 'mix')]
  
  colsComp <- getColsComp(ncol(dens) - 2)
  plotHistDens(x, cbp) + geom_line(data=dens, aes(x=diameter, y=mix))
}



plotHistAndTrue <- function(x, mdl, cbp, sigma=1) {
  stopifnot(sigma >= 0)
  stopifnot(mdl$mClu$parameters$variance$sigmasq > sigma ^ 2)

  # Compute density of *true* diameter:
  thisMdl <- mdl
  thisMdl$mClu$parameters$variance$sigmasq <- mdl$mClu$parameters$variance$sigmasq - sigma ^ 2
  dDens <- getDensity(thisMdl, seq(from=5.5, to=40.5, length.out=1000)) %>% select(diameter, mix) %>% rename(density=mix)
  addCBPsToPlot(ggplot(data.frame(diameter=unname(x))) + geom_bar(aes(x=diameter, y=(..count..)/sum(..count..))) + geom_line(data=dDens, aes(x=diameter, y=density)) + xlab('diameter / mm') + ylab('density'), cbp)
}



plotDensity <- function(x, mdl, cbp) {
  dens <- getDensity(mdl, seq(from=5.5, to=40.5, length.out=1000))
  dDens <- melt(dens, id.vars='diameter', variable.name='component', value.name='density')
  
  colsComp <- getColsComp(ncol(dens) - 2)
  plotHistDens(x, cbp) + geom_line(data=dDens, aes(x=diameter, y=density, color=component)) + scale_colour_manual(values=colsComp)
}



plotCdf <- function(x, mdl, cbp) {
  p <- getParams(mdl)
  sd <- sqrt(ifelse(rep(mdl$mClu$modelName == 'E', mdl$mClu$G), rep(p$variance$sigmasq, mdl$mClu$G), p$variance$sigmasq))  # standard deviations of the different components
  xx <- seq(from=5.5, to=40.5, length.out=1000)  # points at which cdf is computed and plotted
  cdf <- sapply(1:mdl$mClu$G, function(i) pnorm(xx, mean=p$mean[i], sd=sd[i]))
  if (!is.na(mdl$w6)) {
    cdf <- cbind((xx > 6) * 1, cdf)
  }
  cdfMix <- apply(cdf, 1, weighted.mean, w=getWComb(mdl))
  tmp <- data.frame(diameter=xx, mix=cdfMix, cdf)
  names(tmp) <- getMdlColNames(mdl)
  dCdf <- melt(tmp, id.vars='diameter', variable.name='component', value.name='cdf')
  
  colsComp <- getColsComp(nlevels(dCdf[, 2]) - 1)
  postprocessPlot(ggplot(data.frame(diameter=x), aes(x=diameter)) + stat_ecdf(colour='grey35', size=1) + geom_line(data=dCdf, aes(x=diameter, y=cdf, color=component), size=1) + scale_colour_manual(values=colsComp) + ylab('cdf'))
}



postprocessPlot <- function(plt, cbp) {
  addCBPsToPlot(plt +
             xlab('observed diameter / mm') +
             theme(legend.position="top") +
             coord_cartesian(xlim=c(5.5, 40.5)), cbp)
}



addCBPsToPlot <- function(plt, cbp) {
  plt + geom_vline(xintercept=c(cbp$R, cbp$S) - 0.5, linetype=2, na.rm=TRUE)
}



plotDeconvolution <- function(x, mdl, cbp, sigma=1) {
  stopifnot(sigma >= 0)
  stopifnot(mdl$mClu$parameters$variance$sigmasq > sigma ^ 2)
  
  #############
  # Top panel #
  #############
  
  pltDens <- plotHistAndTrue(x, mdl, cbp, sigma)
  
  ####################
  # Bottom panel new #
  ####################
  
  dP <- getForcastProbability(mdl, sigma, cbp, yy=seq(from=5.5, by=0.1, to=40.5))
  
  # Zone of methodological uncertainty:
  rangeUncertain <- range(c(cbp$R - 0.5, cbp$S - 0.5, dP$diameter[dP$p < 0.99]))
  zmuMin <- ceiling(rangeUncertain[1])
  zmuMax <- ceiling(rangeUncertain[2])
  uncertainOut <- c(width=diff(rangeUncertain), mean(zmuMin <= x & x < zmuMax), zmuMin=zmuMin, zmuMax=zmuMax)
  dUncertainPlot <- data.frame(xmin=rangeUncertain[1], xmax=rangeUncertain[2], ymin=-Inf, ymax=Inf)
  
  # Plot:
  pltThresh <- addCBPsToPlot(ggplot(dP, aes(x=diameter, y=p)) + geom_rect(data=dUncertainPlot, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='grey35', alpha=0.3, inherit.aes=FALSE) + geom_line() + xlab('observed diameter / mm') + ylab("forecast probability") + coord_cartesian(ylim=c(0, 1)), cbp)

  
  
  ###################################
  # Print summary and combine plots #
  ###################################
  
  cat(sep='', '\n\nOfficial CBPs: ', cbp$R, ' mm and ', cbp$S, ' mm. ZMU: [', rangeUncertain[1], ' mm, ', rangeUncertain[2], ' mm].\n\n')
  
  # Combined figure:
  pltDeconv <- arrangeGrob(pltDens + theme(legend.position="none"), pltThresh + theme(legend.position="none"), heights=rep(1, 2))
  
  
  #######################################
  # Probabilities of missclassification #
  #######################################
  
  thisMdl <- mdl
  thisMdl$mClu$parameters$variance$sigmasq <- mdl$mClu$parameters$variance$sigmasq - sigma ^ 2
  
  # More systematic approach:
  getIntegrand <- function(errBndry, errtype) {
    stopifnot(errtype %in% c('<', '>'))  # '<' / '>': error E has to be smaller / greater than errBndry
    integrand <- function(x) {
      pX <- predict(thisMdl$mClu, x, what='dens')
      if (!is.na(thisMdl$w6)) {
        pX <- pX * (1 - thisMdl$w6)
      }
      pEps <- pnorm(errBndry, mean=x, sd=sigma)
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



getMdlColNames <- function(mdl) {
  c('diameter', 'mix', 1:(mdl$mClu$G + !is.na(mdl$w6)))
}


