# This is a small script to plot the QTL mapping results.
library(qtl)

# SCRIPT PARAMETERS
# -----------------
threshold <- 0.05  # Significance threshold.
ymax      <- 9     # Height of vertical axis.

# Load the QTL mapping results.
load("../results/gwscan.rr.Rdata")
n <- length(phenotypes)

# PLOT QTL MAPPING RESULTS
# ------------------------
x11(width = 19.5,height = 5.5,title = "QTL mapping results")
par(ps = 14,family = "HersheySans",font.axis = 2,font.lab = 2,
    font.main = 2,cex.main = 1)
par(mai = c(0.6,0.6,0.3,0.3))
layout(matrix(1:(3*n),nrow = 3,ncol = n))

# Repeat for each phenotype.
for (i in 1:n) {
  phenotype <- phenotypes[i]

  # PLOT RESULTS FROM F2 DATA
  # -------------------------
  # Plot the QTL mapping results from qtl using the F2 cross.
  plot(gwscan$F2.qtl,lodcolumn = i,incl.markers = FALSE,lwd = 4,
       bandcol = "powderblue",col = "dodgerblue",ylim = c(0,ymax),
       gap = 0,xlab = "",ylab = "LOD",font.lab = 2,font.main = 2,
       xaxp = c(1,10,1),main = paste(phenotype,"(F2 only)"))
  
  # Plot the QTL mapping results from QTLRel using the F2 cross.
  plot(gwscan$F2.rel,lodcolumn = i,incl.markers = FALSE,lwd = 1,
       bandcol = "powderblue",col = "darkblue",ylim = c(0,ymax),
       gap = 0,add = TRUE)

  # Add the significance threshold to the plot.
  add.threshold(gwscan$F2.qtl,perms = perms$F2,alpha = threshold,gap = 0,
                lodcolumn = i,col = "darkorange",lty = "dotted")

  # PLOT RESULTS FROM F34 DATA
  # --------------------------
  # Plot the QTL mapping results from QTLRel using the F34 cross.
  plot(gwscan$F34,lodcolumn = i,incl.markers = FALSE,lwd = 2,
       bandcol = "powderblue",col = "darkblue",ylim = c(0,ymax),
       gap = 0,xlab = "",ylab = "LOD",font.lab = 2,
       main = paste(phenotype,"(F34 only)"))

  # Add the significance threshold to the plot.
  add.threshold(gwscan$F34,perms = perms$F34,alpha = threshold,gap = 0,
                lodcolumn = i,col = "darkorange",lty = "dotted")

  # PLOT RESULTS FROM F2 + F34 DATA
  # -------------------------------
  # Plot the QTL mapping results from QTLRel using the combined F2 and
  # F34 samples.
  plot(gwscan$combined,lodcolumn = i,incl.markers = FALSE,lwd = 2,
       bandcol = "powderblue",col = "darkblue",ylim = c(0,ymax),
       gap = 0,xlab = "",ylab = "LOD",font.lab = 2,
       main = paste(phenotype,"(F2 + F34)"))  

  # Add the significance threshold to the plot.
  add.threshold(gwscan$combined,perms = perms$combined,alpha = threshold,
                gap = 0,lodcolumn = i,col = "darkorange",lty = "dotted")
}

