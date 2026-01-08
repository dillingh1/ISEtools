#' @title Plot ISE parameter values
#' @description Plots histograms of ISE parameter values a, b, c, sigma, and LOD (alpha, beta or S/N) for the equation y = a + b log(x + c) + error, where error follows a normal distribution with mean 0 and standard deviation sigma.
#' @param x ISE description (e.g. object of class ISEdescription)
#' @param breaks Approximate number of bins for histograms, defaults to 20	
#' @param ... Other arguments to be passed through to plotting (histogram) functions
#' @author Peter Dillingham, \email{peter.dillingham@@otago.ac.nz}
#' @seealso \code{\link{describeISE}}
plot.ISEdescription = function(x, breaks=20, ...) {
###
# plots histograms of ISE parameters and LOD if the sample coda was kept in describeISE (keep.coda=T)
###
	ISEdescription = x
	coda.exists = exists('ahat.coda', where = ISEdescription)
	if (!coda.exists) {
		message('No relevant plot.  Use keep.coda = TRUE in describeISE')
		return()
	}
	if (coda.exists) {
		if(ISEdescription$LOD.info$type == "alpha, beta") {
			ISE.LOD.label = bquote("log LOD"[paste(alpha,"=",.(ISEdescription$LOD.info$alpha), ",", beta,"=",.(ISEdescription$LOD.info$beta))])
		} else if(ISEdescription$LOD.info$type == "S/N") {
			ISE.LOD.label = bquote("log LOD"[paste("S/N=",.(ISEdescription$LOD.info$SN))])
		}
		R = length(ISEdescription$ahat)
		
		# Plot parameter estimates, one plot per ISE
		# Return mfrow values on exit
		oldpar <- graphics::par(no.readonly = TRUE) 
		on.exit(graphics::par(oldpar), add = TRUE)
		graphics::par(mfrow=c(3,2))

		if (R == 1) {
		  graphics::plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
		  graphics::text(0.5, 0.5, "ISE #1" )
		  graphics::hist(ISEdescription$ahat.coda, main="", xlab = expression(a), breaks=breaks, ...)
		  graphics::hist(ISEdescription$bhat.coda, main="", xlab = expression(b), breaks=breaks, ...)
		  graphics::hist(log10(ISEdescription$chat.coda), main="", xlab = expression(paste(log[10], c)), breaks=breaks, ...)
		  graphics::hist(ISEdescription$sigmahat.coda, main="", xlab = expression(sigma), breaks=breaks, ...)
		  graphics::hist(pmin(ISEdescription$LOD.coda, 0.01), main="", xlab = ISE.LOD.label, breaks=breaks, ...)
			
		}
		if (R > 1) {
			for (i in 1:R) {
			  graphics::plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
			  graphics::text(0.5, 0.5, paste("ISE #", i, sep="") )
			  graphics::hist(ISEdescription$ahat.coda[,i], main="", xlab = expression(a), breaks=breaks, ...)
			  graphics::hist(ISEdescription$bhat.coda[,i], main="", xlab = expression(b), breaks=breaks, ...)
			  graphics::hist(log10(ISEdescription$chat.coda[,i]), main="", xlab = expression(paste(log[10], c)), breaks=breaks, ...)
			  graphics::hist(ISEdescription$sigmahat.coda[,i], main="", xlab = expression(sigma), breaks=breaks, ...)
			  graphics::hist(pmin(ISEdescription$LOD.coda[,i], 0.01), main="", xlab = ISE.LOD.label, breaks=breaks, ...)
			}
		}
	}
}
