#' @title Annotate spectra (MS/MS or pseudo-MS/MS) stored in a MSnExp object
#'
#' @description
#' This function annotates fragmentation spectra generated in LC-MS experiments
#' or pseudo-MS/MS spectra from AIF data or ISF from MS1 scans which are loaded
#' into an MSnExp data object from the MSnbase package.
#' The spectra are generally acquired in LC-MS or MS imaging experiments. 
#'
#' @author Goncalo Graca (Imperial College London)
#'
#' @param MSnExpObj A MSn experiment data ("MSnExp"), usually obtained after
#' importing spectra in .mgf format.
#' @param libs Fragment libraries to use. Specify one of default libraries
#' provided as data object with the package (\code{LipidPos}, \code{LipidNeg}, 
#' \code{MetabolitesPos} and \code{MetabolitesNeg}) or the full path to 
#' user-defined libraries.
#' @param RTs Optional data.frame with Lipid/metabolites classes Retention
#' Times in seconds.
#' @param checkIsotope Whether or not to check the isotope type;
#' default is set to FALSE as the precursor ions will likely be the 
#' monoisotopic mass of the isotopic series.
#' @param tolerance Tolerance in ppm for the candidate search.
#' @param maxMZdiff Maximum m/z difference between candidate fragments and
#' pseudo-MS/MS or AIF ions in Da.
#' @param matchWeight weight of the fragment matches to the final score;
#' value between 0 and 1; the remaining fraction of the weight comes from the
#' candidate m/z error.
#' @return For each spectra in the MSnExp object the function will return 
#' a list containing: a data frame with with rank 1 annotations 
#' (\code{global}), the the date and time of annotation, a data frame with the 
#' annotation options. For each feature the following lists are returned: 
#' ranked annotations for each feature (\code{rankedResults}), 
#' the corresponding ranked matched spectra (\code{rankedSpectra}), 
#' the pseudo-MS/MS spectra (\code{pseudoMSMS}), in-source spectra 
#' (\code{inSourceSpectra}) and AIF spectrum (\code{AIFspectra}).
#' @examples
#' # read the spectra to annotate from an .mgf file
#' mgfFile <- system.file("extdata", "pseudoMSMS.mgf", 
#' package="MetaboAnnotatoR")
#' MSnExpObj <- MSnbase::readMgfData(mgfFile)
#' # Read the default lipid positive libraries
#' data("LipidPos")
#' # Run the annotation procedure
#' annotations <- annotateSpectra(MSnExpObj, libs="LipidPos", RTs="none", 
#' checkIsotope=FALSE)
#' @export
annotateSpectra <- function(MSnExpObj,
                            libs="LipidPos", RTs="none",
                            checkIsotope=TRUE, tolerance=25,
                            maxMZdiff=0.01, matchWeight=0.5){
	
    ## organise spectra related data from the MSnExpObj -----------------------
	# extract mz and intensity values from the spec object
    # create spectrum objects with mz and intensity
    s <- lapply(seq_along(spectra(MSnExpObj)), 
            function(x) data.frame(mz=mz(MSnExpObj)[[x]], 
                                   into=intensity(MSnExpObj)[[x]]))

    # organise spectra data to specData list
    specData <- list(precursor=precursorMz(MSnExpObj), 
                    spectra=s,
                    rt=rtime(MSnExpObj))
    
    ## Targets list, i.e. list of precursors and RT (0 for imaging MS spectra)
    targets <- data.frame(feature.mz=specData$precursor, 
                            feature.rt=specData$rt)
    
    ## Initialize results object------------------------
    results <- initializeResultsSpec(targets)

    ## RTs intervals specified? --------------------------------------------
    if(RTs == "none") {
        message("No RT information provided...")
    } else if(is.data.frame(RTs)){
        message("Using user provided RT information...")
    } else stop("RTs must be a data.frame")

    ## load libraries ---------------------------------------------------------
    if(libs == "LipidPos") {
        if(!exists("LipidPos")) {
            stop("LipidPos not found, please use data(LipidPos)")
        } else libraries <- LipidPos
    } else if(libs == "LipidNeg") {
        if(!exists("LipidNeg")) {
            stop("LipidNeg not found, please use data(LipidNeg)")
        } else libraries <- LipidNeg
    } else if(libs == "MetabolitesPos") {
        if(!exists("MetabolitesPos")) {
            stop("MetabolitesPos not found, please use data(MetabolitesPos)")
        } else libraries <- MetabolitesPos
    } else if(libs == "MetabolitesNeg") {
        if(!exists("MetabolitesNeg")) {
            stop("MetabolitesNeg not found, please use data(MetabolitesNeg)")
        } else libraries <- MetabolitesNeg
    } else libraries <- loadLibs(libs)
    
    libfiles <- libraries$libfiles
    lib <- libraries$lib

    ## process each feature from the targets table-----------------------------
    for(i in seq_len(nrow(targets))){
        progNote <- paste("... Processing feature", i, "of", nrow(targets), 
                            "...")
        message(progNote)
        fmz <- targets[i,1]
        frt <- targets[i,2]

        # read MS/MS spectra --------------------------------------------------
        pseudoSpec <- specData$spectra[[i]]
        highCESpec <- pseudoSpec
        # store spectra in the annotation results object
        if(!is.null(pseudoSpec)) {
            results$querySpec[[i]] <- pseudoSpec
            } else results$querySpec[[i]] <- NA

        # Isotope check --------------------
        if(!checkIsotope) iso <- 0 else {
            iso <- checkIsotope(fmz, frt, pseudoSpec)
        }

        # Search Libraries --------------------
        if(is.null(pseudoSpec)) { next
        } else {
            message("Searching candidates...")
            candidates <- searchLib(lib, libfiles, fmz-iso*1.0034, frt,
                                    tolerance=tolerance, RTs, pseudoSpec)
        }
        #Compare fragments between Library candidates and high-collision-energy
        #pseudo-MS/MS spectra --------------------
        if(is.null(pseudoSpec) & 
            is.null(highCESpec) | 
            length(unlist(candidates)) == 0){
            result <- NULL
        } else {
            message("Matching fragments to pseudo-MS/MS...")
            output <- mapply(compFrag, candidates, 
                                lapply(as.numeric(names(candidates)),
                                function(x) lib[[x]]),
                                MoreArgs=list(fmz, frt, iso, highCESpec, 
                                                pseudoSpec,
                                                maxMZdiff=maxMZdiff,
                                                matchWeight=matchWeight), 
                                SIMPLIFY=FALSE)
            result <- do.call(rbind, lapply(output, "[[", 1))
            specMatch <- unlist(lapply(output, "[[", 2), recursive = FALSE)
            specMatch <- specMatch[!(specMatch) == "NULL"]
        }
        # Score ranking --------------------
        if(is.null(result)) {
            rankedResult <- targets[i,c(1,2)]
            rankedResult[c("metabolite", "feature.type", "ion.type", "isotope",
                            "mz.metabolite", "matched.mz", "mz.error", 
                            "pseudoMSMS", "fraction", "score")] <- NA
            rankedSpec <- NULL

        } else {
            output <- rankScore(result, specMatch)
            rankedSpec <- output$rankedSpecMatch
            rankedResult <- output$rankedResult
            # type of ion isotope
            rankedResult$isotope <- paste("M+", iso, sep="")
            # pseudoMSMS flag
            rankedResult$pseudoMSMS <- "TRUE"
            results$rankedResult[[i]] <- rankedResult
            results$rankedSpectra[[i]] <- rankedSpec
        }
        ## Store highest rank annotation in global results------
        results$global[i,] <- storeAnnotations(global=results$global[i,], 
                                rankedResult[1,])
    }
    # store options
    df <- data.frame(dataType="MSnExp",
                    polarity=unique(polarity(MSnExpObj)),
                    libraries=libs,
                    RTs=RTs,
                    checkIsotope=checkIsotope, matchWeight=matchWeight,
                    tolerance=paste(tolerance, "ppm"),
                    maxMZdiff=paste(maxMZdiff, "Da"), row.names="parameter")
    results$options <- as.data.frame(t(df))

    message('Job done!')

    return(results)
}


## Helper functions------------------------------------------------------------

## Initialization of results folder and global results table
initializeResultsSpec <- function(targets){
    # create objects to store results and metadata
    Date <- Sys.Date()
    Time <- format(Sys.time(), "%X")
    rankedResult <- list()
    rankedSpectra <- list()
    querySpec <- list()
    options <- NULL
    # create table to store global results
    global <- targets
    global[,c("metabolite", "feature.type", "ion.type", "isotope",
                "mz.metabolite", "matched.mz", "mz.error",
                "pseudoMSMS", "fraction", "score")] <- NA
    # return global results table and results path as list
    results <- list(global=global,
                    Date=Date, 
                    Time=Time,
                    options=options,
                    rankedResult=rankedResult, 
                    rankedSpectra=rankedSpectra, 
                    querySpec=querySpec)
    return(results)
}
