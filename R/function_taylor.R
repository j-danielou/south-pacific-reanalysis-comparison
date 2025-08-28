taylor.diagram = function (ref, model, label, add = FALSE, 
                           col = "red", pch = 19, labpos=3,
                           pos.cor = TRUE, xlab = "", 
                           ylab = "", reflabel="Observed",
                           main = "", 
                           show.gamma = TRUE, ngamma = 5, 
                           gamma.col = 8, sd.arcs = 0, 
                           ref.sd = TRUE, 
                           grad.corr.lines = c(2, 4, 6, 8, 9, 9.5, 9.9)/10, 
                           pcex = 0.7, tcex=0.7, normalize = TRUE, 
                           mar = c(5, 5, 3, 3), test = "sst", ...) 
{
  
  grad.corr.full = seq(0,1, by=0.2)
  R = cor(ref, model, use = "pairwise")
  sd.r = sd(ref)
  sd.f = sd(model)
  if (normalize) {
    sd.f = sd.f/sd.r
    sd.r = 1
  }
  #maxsd = ifelse(test == TRUE, 11, 1.5 * max(sd.f, sd.r))
  maxsd = ifelse(test == "sss", 11, 2)
  oldpar = par("mar", "xpd", "xaxs", "yaxs")
  if (!add) {
    if (pos.cor) {
      if (nchar(ylab) == 0) 
        ylab = "Standard deviation"
      par(mar = mar)
      plot(0, xlim = c(0, maxsd), ylim = c(0, maxsd), xaxs = "i", 
           yaxs = "i", axes = FALSE, main = main, xlab = xlab, 
           ylab = ylab, type = "n", ...)
      if (grad.corr.lines[1]!=0) {
        for (gcl in grad.corr.lines) {
          lines(c(0, maxsd*gcl), c(0, maxsd * sqrt(1 - gcl^2)), lty = 3)
        }
      }
      segments(c(0, 0), c(0, 0), c(0, maxsd), c(maxsd,0))
      axis.ticks = pretty(c(0, maxsd))
      axis.ticks = axis.ticks[axis.ticks <= maxsd]
      axis(1, at = axis.ticks)
      axis(2, at = axis.ticks, las = 1)
      if (sd.arcs[1]) {
        if (length(sd.arcs) == 1) sd.arcs = axis.ticks
        for (sdarc in sd.arcs) {
          xcurve = cos(seq(0, pi/2, by = 0.03))*sdarc
          ycurve = sin(seq(0, pi/2, by = 0.03))*sdarc
          lines(xcurve, ycurve, col = "blue", lty = 3)
        }
      }
      if (show.gamma[1]!=0) {
        if(length(show.gamma) > 1) {
          gamma = show.gamma
        } else {
          gamma = pretty(c(0, maxsd), n = ngamma)[-1]
        }
        if(gamma[length(gamma)] > maxsd) {
          gamma = gamma[-length(gamma)]
        } 
        labelpos = seq(45, 70, length.out = length(gamma))
        for (gindex in seq_along(gamma)) {
          xcurve = cos(seq(0, pi, by = 0.03))*gamma[gindex] + 
            sd.r
          endcurve = which(xcurve < 0)
          endcurve = ifelse(length(endcurve), min(endcurve) - 
                              1, 105)
          ycurve = sin(seq(0, pi, by = 0.03)) * gamma[gindex]
          maxcurve = xcurve * xcurve + ycurve * ycurve
          startcurve = which(maxcurve > maxsd * maxsd)
          startcurve = ifelse(length(startcurve), max(startcurve) + 
                                1, 0)
          lines(xcurve[startcurve:endcurve], ycurve[startcurve:endcurve], 
                col = gamma.col)
          boxed.labels(xcurve[labelpos[gindex]], ycurve[labelpos[gindex]], 
                       gamma[gindex], border = FALSE)
        }
      }
      xcurve = cos(seq(0, pi/2, by = 0.01)) * maxsd
      ycurve = sin(seq(0, pi/2, by = 0.01)) * maxsd
      lines(xcurve, ycurve)
      bigtickangles = acos(seq(0.1, 0.9, by = 0.1))
      medtickangles = acos(seq(0.05, 0.95, by = 0.1))
      smltickangles = acos(seq(0.91, 0.99, by = 0.01))
      segments(cos(bigtickangles) * maxsd, sin(bigtickangles) * 
                 maxsd, cos(bigtickangles) * 0.97 * maxsd, sin(bigtickangles) * 
                 0.97 * maxsd)
      par(xpd = TRUE)
      if (ref.sd) {
        xcurve = cos(seq(0, pi/2, by = 0.01)) * sd.r
        ycurve = sin(seq(0, pi/2, by = 0.01)) * sd.r
        lines(xcurve, ycurve)
      }
      points(sd.r, 0)
      text(cos(c(bigtickangles, acos(c(0.95, 0.99)))) * 
             1.05 * maxsd, sin(c(bigtickangles, acos(c(0.95, 
                                                       0.99)))) * 1.05 * maxsd, c(seq(0.1, 0.9, by = 0.1), 
                                                                                  0.95, 0.99))
      text(maxsd*0.8, maxsd*0.8, srt=-45, "Correlation")
      segments(cos(medtickangles) * maxsd, sin(medtickangles) * 
                 maxsd, cos(medtickangles) * 0.98 * maxsd, sin(medtickangles) * 
                 0.98 * maxsd)
      segments(cos(smltickangles) * maxsd, sin(smltickangles) * 
                 maxsd, cos(smltickangles) * 0.99 * maxsd, sin(smltickangles) * 
                 0.99 * maxsd)
    }
    else {
      x = ref
      y = model
      R = cor(x, y, use = "pairwise.complete.obs")
      E = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
      xprime = x - mean(x, na.rm = TRUE)
      yprime = y - mean(y, na.rm = TRUE)
      sumofsquares = (xprime - yprime)^2
      Eprime = sqrt(sum(sumofsquares)/length(complete.cases(x)))
      E2 = E^2 + Eprime^2
      if (add == FALSE) {
        maxray = 1.5 * max(sd.f, sd.r)
        plot(c(-maxray, maxray), c(0, maxray), type = "n", 
             asp = 1, bty = "n", xaxt = "n", yaxt = "n", 
             xlab = xlab, ylab = ylab, main = main)
        discrete = seq(180, 0, by = -1)
        listepoints = NULL
        for (i in discrete) {
          listepoints = cbind(listepoints, maxray * 
                                cos(i * pi/180), maxray * sin(i * pi/180))
        }
        listepoints = matrix(listepoints, 2, length(listepoints)/2)
        listepoints = t(listepoints)
        lines(listepoints[, 1], listepoints[, 2])
        lines(c(-maxray, maxray), c(0, 0))
        lines(c(0, 0), c(0, maxray))
        for (i in grad.corr.lines) {
          lines(c(0, maxray * i), c(0, maxray * sqrt(1 - 
                                                       i^2)), lty = 3)
          lines(c(0, -maxray * i), c(0, maxray * sqrt(1 - 
                                                        i^2)), lty = 3)
        }
        for (i in grad.corr.full) {
          text(1.05 * maxray * i, 1.05 * maxray * sqrt(1 - 
                                                         i^2), i, cex = 0.6)
          text(-1.05 * maxray * i, 1.05 * maxray * sqrt(1 - 
                                                          i^2), -i, cex = 0.6)
        }
        seq.sd = seq.int(0, 2 * maxray, by = (maxray/10))[-1]
        for (i in seq.sd) {
          xcircle = sd.r + (cos(discrete * pi/180) * 
                              i)
          ycircle = sin(discrete * pi/180) * i
          for (j in 1:length(xcircle)) {
            if ((xcircle[j]^2 + ycircle[j]^2) < (maxray^2)) {
              points(xcircle[j], ycircle[j], col = "darkgreen", 
                     pch = ".")
              if (j == 10) 
                text(xcircle[j], ycircle[j], signif(i, 
                                                    2), cex = 0.5, col = "darkgreen")
            }
          }
        }
        seq.sd = seq.int(0, maxray, length.out = 5)
        for (i in seq.sd) {
          xcircle = (cos(discrete * pi/180) * i)
          ycircle = sin(discrete * pi/180) * i
          if (i) 
            lines(xcircle, ycircle, lty = 3, col = "blue")
          text(min(xcircle), -0.03 * maxray, signif(i, 
                                                    2), cex = 0.5, col = "blue")
          text(max(xcircle), -0.03 * maxray, signif(i, 
                                                    2), cex = 0.5, col = "blue")
        }
        text(0, -0.08 * maxray, "Standard Deviation", 
             cex = 0.7, col = "blue")
        text(0, -0.12 * maxray, "Centered RMS Difference", 
             cex = 0.7, col = "darkgreen")
        points(sd.r, 0, pch = 19, bg = "darkgreen", cex = 0.7)
      }
      S = (2 * (1 + R))/(sd.f + (1/sd.f))^2
    }
  }
  crmsd = sqrt(sd.r^2 + sd.f^2 - 2 * sd.r * sd.f * R)
  points(sd.f * R, sd.f * sin(acos(R)), pch = pch, col = col, 
         cex = pcex)
  text(sd.f * R, sd.f * sin(acos(R)),  labels=label, 
       cex = tcex, pos=labpos, col=col) 
  if (!add) text(sd.r, 0,  labels=reflabel, cex = tcex, adj=c(-0.1,1.1)) 
  invisible(oldpar)
  return(list(sd.f = sd.f, 
              R = R, 
              crmsd = crmsd))
}


boxed.labels = function (x, y = NA, labels, bg = ifelse(match(par("bg"), "transparent", 
                                                              0), "white", par("bg")), border = TRUE, xpad = 1.2, ypad = 1.2, 
                         srt = 0, cex = 1, adj = 0.5, ...) 
{
  oldcex = par("cex")
  par(cex = cex)
  if (is.na(y) && is.list(x)) {
    y = unlist(x[[2]])
    x = unlist(x[[1]])
  }
  box.adj = adj + (xpad - 1) * cex * (0.5 - adj)
  if (srt == 90 || srt == 270) {
    bheights = strwidth(labels)
    theights = bheights * (1 - box.adj)
    bheights = bheights * box.adj
    lwidths = rwidths = strheight(labels) * 0.5
  }
  else {
    lwidths = strwidth(labels)
    rwidths = lwidths * (1 - box.adj)
    lwidths = lwidths * box.adj
    bheights = theights = strheight(labels) * 0.5
  }
  args = list(x = x, y = y, labels = labels, srt = srt, adj = adj, 
              col = ifelse(colSums(col2rgb(bg) * c(1, 1.4, 0.6)) < 
                             350, "white", "black"))
  args = modifyList(args, list(...))
  rect(x - lwidths * xpad, y - bheights * ypad, x + rwidths * 
         xpad, y + theights * ypad, col = bg, border = border)
  do.call(text, args)
  par(cex = oldcex)
}


