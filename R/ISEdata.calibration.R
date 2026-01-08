ISEdata.calibration <-
function(filename.calibration, calibration.only) {
##
# A sub-function called by loadISEdata for when there is only calibration data (i.e. the goal is to characterise an ISE, not estimate
#   activity for unknowns
##
    # List of required column names for the calibration data and the two experimental data formats
    req_calib = c("ISEID","log10x","emf")
	
	# Load the calibration files
    data.calib = utils::read.delim(filename.calibration, header=TRUE, sep="\t",
                               quote="\"", dec=".", fill=TRUE, comment.char="",
                               stringsAsFactors = FALSE, check.names = FALSE)
     
    if (!all(req_calib %in% names(data.calib))) {
      stop("Calibration file must contain columns: ", paste(req_calib, collapse=", "))
    }

	### Format data from the calibration file
	N = nrow(data.calib)
	
	# Confirm that ISE IDs were entered correctly (consecutive, starting at 1)
	ise.tmp = data.calib$ISEID
	if (!is.numeric(ise.tmp) || anyNA(ise.tmp) || any(ise.tmp %% 1 != 0)) stop("ISEID must be integer.")
	if (!setequal(sort(unique(ise.tmp)), seq_len(max(ise.tmp)))) stop("ISEID must be consecutive starting at 1.")
	R = max(data.calib$ISEID)
	
	ISEID = data.calib$ISEID
	log10x = data.calib$log10x
	emf = data.calib$emf
	
	if (R == 1) { data.out = list(N=N, R=R, log10x = log10x, emf=emf, 
		stdadd = NA, calibration.only = calibration.only, data.calib = data.calib, data.exp = NA) }
	if (R > 1) { data.out = list(N=N, R=R, ISEID=ISEID, log10x = log10x, emf=emf, 
		stdadd = NA, calibration.only = calibration.only, data.calib = data.calib, data.exp = NA) }

	return(data.out)
}
