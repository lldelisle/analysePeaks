#' Create a GRList with a reference GR and a sample GR annotated with GRanges provided in arguments
#'
#'@param e name of the experiment usually sample1____sample2(____sample3 etc...)
#'@param ovF a list which contains the filtered overlaps between the samples. Should at least contains a item named e.
#'@param allExpe a GRangeList with at least each sample in `e`
#'@param grLScores a GRangeList with annotations for which you would like to know the best score overlapping the GRanges of experiments
#'@param grLDistance a GRangeList with annotations for which you would like to know the distance to the closest item
#'@param nameOfRef among the samples which one is the reference, the others will be merged as "replicates"
#'@param useSummitPMFlanking logical value to specify if you prefer to use the full region of the peak or only the summit +/- the `flankingSize`
#'@param flankingSize integer value used to extend the summit if `useSummitPMFlanking` is TRUE
#'@return A list with:
#'`stringSet` a string with all samples other than the reference separated by comma
#'`myGRs` a GRangeList of 2 elements which are
#'       1. The peaks which are shared by all samples which are considered as "replicates"
#'       2. The peaks of the reference
#'`nameOfRef` the name of the reference used
#'`namesOfGrLScores` the names of grLScores
#'`namesOfGrLDistance` the names of grLDistance
#'`useSummitPMFlanking` the value of useSummitPMFlanking which was used
#'`flankingSize` the value of flankingSize which was used
#'@details The GRanges of myGRs will have in mcols a new information which is
#'bestNameOfTheGRScore for each GR in grLScores and
#'distanceToTheNearestNameOfTheGR for each GR in grLDistance
#'The first item of myGRs will also have inUniqueRef which is the index of the merged item of the Reference.
#'The second item of myGRs will also have inNoneOfTheSet
#'@importFrom GenomicRanges mcols start end flank findOverlaps distanceToNearest
#'@export
createMyGRs <- function(e, ovF, allExpe,
                        grLScores = NULL, grLDistance = NULL,
                        nameOfRef = "ChIP",
                        useSummitPMFlanking = T, flankingSize = 150){
  if ( ! e %in% names(ovF)){
    stop(e, "is not part of the names of ovF")
  }
  # df contains the filtered overlaps
  df <- ovF[[e]]
  samplesToCheck <- colnames(df)
  if (!nameOfRef %in% samplesToCheck){
    stop(nameOfRef, " is not part of ", e)
  }
  # stringSet contains the name of other samples which are not the Ref
  stringSet <- paste(setdiff(samplesToCheck, nameOfRef), collapse = ",")
  df$nbNA <- apply(df, 1, function(v){
    sum(is.na(v))
  })
  df$RefisNA <- as.numeric(is.na(df[, nameOfRef]))
  # name1 is the name of the first non Ref sample
  name1 <- setdiff(samplesToCheck, nameOfRef)[1]
  # We are now selecting the indices of name1 items
  # which overlap with all other samples which are not the Ref
  i_fullOverlap <- df[, name1][(df$nbNA - df$RefisNA) == 0]
  # grSetRep contains the subset of the GRange of name1 which
  # Overlap with all other samples which are not the Ref
  grSetRep <- allExpe[[name1]][i_fullOverlap]
  # We store in this GRanges the index of the overlapped item in Ref
  # (in the metadata "inUniqueRef")
  grSetRep$inUniqueRef <- df[match(i_fullOverlap, df[, name1]), nameOfRef]
  # We now change the coordiantes of the GRange and keep
  # the summit +/- the flankingSize
  grMySummitsExtended <- grSetRep
  if (useSummitPMFlanking){
    grMySummits <- grSetRep
    GenomicRanges::start(grMySummits) <- GenomicRanges::start(grSetRep) +
      grSetRep$relativeSummit
    GenomicRanges::end(grMySummits) <- GenomicRanges::start(grMySummits)
    grMySummitsExtended <- GenomicRanges::flank(grMySummits,
                                                width = flankingSize,
                                                both = T)
  }
  grRef <- allExpe[[nameOfRef]]
  # We annotate the reference GRanges to
  # Specify when it is a totally specific peak
  grRef$inNoneOfTheSet <- TRUE
  grRef$inNoneOfTheSet[na.omit(df[df$nbNA != (length(samplesToCheck) - 1),
                                  nameOfRef])] <- F
  grRefSummitsExtended <- grRef
  if (useSummitPMFlanking){
    grRefSummits <- grRef
    GenomicRanges::start(grRefSummits) <- GenomicRanges::start(grRef) +
      grRef$relativeSummit
    GenomicRanges::end(grRefSummits) <- GenomicRanges::start(grRefSummits)
    grRefSummitsExtended <- GenomicRanges::flank(grRefSummits,
                                                 width = flankingSize,
                                                 both = T)
  }
  myGRs <- list(grMySummitsExtended, grRefSummitsExtended)
  # We will annotate the GRanges with the features which are in parameters
  for (iGR in 1:2){
    for (jGR in 1:length(grLScores)){
      annotScoreOverlap <-
        as.data.frame(GenomicRanges::findOverlaps(myGRs[[iGR]],
                                                  grLScores[[jGR]]))
      annotScoreOverlapByScore <-
        aggregate(list(
          score = grLScores[[jGR]]$score[annotScoreOverlap$subjectHits]
        ), by = list(peak = annotScoreOverlap$queryHits), FUN = max)
      nameJ <- names(grLScores)[jGR]
      GenomicRanges::mcols(myGRs[[iGR]])[, paste0("best", nameJ, "Score")] <- 0
      GenomicRanges::mcols(myGRs[[iGR]])[annotScoreOverlapByScore$peak,
                                         paste0("best", nameJ, "Score")] <-
        annotScoreOverlapByScore$score
    }
    for (jGR in 1:length(grLDistance)){
      nameJ <- names(grLDistance)[jGR]
      GenomicRanges::mcols(myGRs[[iGR]])[,
                                         paste0("distanceToNearest",
                                                nameJ)] <-
        as.data.frame(
          GenomicRanges::distanceToNearest(myGRs[[iGR]],
                                           grLDistance[[jGR]]))$distance
    }
  }
  return(list("stringSet" = stringSet, "myGRs" = myGRs,
              "nameOfRef" = nameOfRef, "namesOfGrLScores" = names(grLScores),
              "namesOfGrLDistance" = names(grLDistance),
              "useSummitPMFlanking" = useSummitPMFlanking,
              "flankingSize" = flankingSize))
}


#' Give category names from thresholds
#'
#' @param thresholdValues a numeric vector sorted in decreasing order which define the borders of categories
#' @param lastCateName a string which give the name of the category when it is below the last value of `thresholdValues`
#' @param unit a string which will be pasted after the threshold values
#' @importFrom utils head
#' @export
cateNamesFromThreshold <- function(thresholdValues, lastCateName, unit) {
  return(c(paste0("above", thresholdValues[1], unit),
           paste0(thresholdValues[-1], "-",
                  utils::head(thresholdValues, length(thresholdValues) - 1),
                  unit),
           lastCateName))
}

#' Add annotations to the list to be able to plot categories
#'
#' @param myGRAndAttributes Should be the output of \link[analysePeaks]{createMyGRs}
#' @param thresholdsForgrLScores a list with the threshold to use to make categories for the features where score is important
#' @param thresholdsForgrLDistance a list with the threshold to use to make categories for the features where distance is important
#' @return a list like myGRAndAttributes but with more attributes.
#' @importFrom GenomicRanges mcols
#' @export
annotateWithCate <- function(myGRAndAttributes, thresholdsForgrLScores,
                             thresholdsForgrLDistance){
  allCates <- list()
  for (jGR in 1:length(thresholdsForgrLScores)){
    nameJ <- names(thresholdsForgrLScores)[jGR]
    thresholdValues <- sort(thresholdsForgrLScores[[jGR]], decreasing = T)
    if (median(thresholdValues) > 1e3){
      cateNames <- analysePeaks::cateNamesFromThreshold(
        thresholdValues = thresholdValues / 1e3,
        lastCateName = paste0("no", nameJ),
        unit = "k")
    } else {
      cateNames <- analysePeaks::cateNamesFromThreshold(
        thresholdValues = thresholdValues,
        lastCateName = paste0("no", nameJ),
        unit = "")
    }
    whatString <- paste("score of", nameJ)
    if (myGRAndAttributes[["useSummitPMFlanking"]]){
      whatString <- paste(whatString,
                          "in summit +/-",
                          myGRAndAttributes[["flankingSize"]],
                          "bp")
    }
    allCates[[length(allCates) + 1]] <- list(
      nameOfOriCol = paste0("best", nameJ, "Score"),
      nameOfCol = paste0("best", nameJ, "ScoreCate"),
      cateNames = cateNames,
      thresholdValues = thresholdValues,
      what = whatString)
  }
  for (jGR in 1:length(thresholdsForgrLDistance)){
    nameJ <- names(thresholdsForgrLDistance)[jGR]
    thresholdValues <- sort(thresholdsForgrLDistance[[jGR]], decreasing = T)
    if (median(thresholdValues) > 1e3){
      cateNames <- analysePeaks::cateNamesFromThreshold(
        thresholdValues = thresholdValues / 1e3,
        lastCateName = paste0("at", nameJ),
        unit = "kb")
    } else {
      cateNames <- analysePeaks::cateNamesFromThreshold(
        thresholdValues = thresholdValues,
        lastCateName = paste0("at", nameJ),
        unit = "bp")
    }
    whatString <- paste("distance to the closest", nameJ)
    allCates[[length(allCates) + 1]] <- list(
      nameOfOriCol = paste0("distanceToNearest", nameJ),
      nameOfCol = paste0("distanceToNearest", nameJ, "Cate"),
      cateNames = cateNames,
      thresholdValues = thresholdValues,
      what = whatString)
  }
  myGRs <- myGRAndAttributes[["myGRs"]]
  for (i in 1:length(allCates)){
    nameOfOriCol <- allCates[[i]][["nameOfOriCol"]]
    nameOfCol <- allCates[[i]][["nameOfCol"]]
    cateNames <- allCates[[i]][["cateNames"]]
    thresholdValues <- allCates[[i]][["thresholdValues"]]
    whatString <- allCates[[i]][["what"]]
    for (iGR in 1:2){
      if (! nameOfOriCol %in% names(GenomicRanges::mcols(myGRs[[iGR]]))){
        stop(nameOfOriCol, " is not part of the attributes.")
      }
      GenomicRanges::mcols(myGRs[[iGR]])[, nameOfCol] <- cateNames[1]
      for (i in 1:length(thresholdValues)){
        GenomicRanges::mcols(myGRs[[iGR]])[
          GenomicRanges::mcols(myGRs[[iGR]])[, nameOfOriCol] <=
            thresholdValues[i],
          nameOfCol] <- cateNames[i + 1]
      }
    }
  }
  myGRAndAttributes[["myGRs"]] <- myGRs
  return(c(myGRAndAttributes, list(allCates = allCates)))
}


#' Extract common and specific peaks in narrowPeak format
#'
#' @param myGRAndAttributes Should be the output of \link[analysePeaks]{createMyGRs}
#' @param outputDirectory a path to a folder where will be written the narrowPeaks files.
#' @param useOriginalNarrowPeak logical to specify if instead of the coordinates in myGRAndAttributes which are by default the summit +/-150bp. You want the original coordinates as in the input (default is F).
#' @return Does not return anything but write 4 files into the `outputDirectory`:
#' Ref_in_at_least_one_Set.narrowPeak
#' Ref_in_not_a_single_Set.narrowPeak
#' share_by_all_Set_no_Ref.narrowPeak
#' shared_by_all_Set_and_Ref.narrowPeak
#' @importFrom GenomicRanges start end
#' @importFrom usefulLDfunctionsGR narrowPeakDFFromGR
#' @export
writeNarrowPeaksFromMyGRsAndAnnot <- function(myGRAndAttributes,
                                              outputDirectory,
                                              useOriginalNarrowPeak = F){
  if (!"myGRs" %in% names(myGRAndAttributes)){
    stop("myGRs is not part of the input myGRAndAttributes.")
  }
  myGRs <- myGRAndAttributes[["myGRs"]]
  dir.create(outputDirectory, showWarnings = F, recursive = T)
  grsToExport <- list()
  outputFileNames <- list()
  # First: shared by all Set and Ref
  grsToExport[[1]] <- subset(myGRs[[1]], !is.na(inUniqueRef))
  outputFileNames[[1]] <- file.path(outputDirectory,
                                    "shared_by_all_Set_and_Ref.narrowPeak")
  # Second: specific to Set:
  grsToExport[[2]] <- subset(myGRs[[1]], is.na(inUniqueRef))
  outputFileNames[[2]] <- file.path(outputDirectory,
                                    "share_by_all_Set_no_Ref.narrowPeak")
  # Third: Ref common with at least one Set
  grsToExport[[3]] <- subset(myGRs[[2]], !inNoneOfTheSet)
  outputFileNames[[3]] <- file.path(outputDirectory,
                                    "Ref_in_at_least_one_Set.narrowPeak")
  # Fourth: Ref specific
  grsToExport[[4]] <- subset(myGRs[[2]], inNoneOfTheSet)
  outputFileNames[[4]] <- file.path(outputDirectory,
                                    "Ref_in_not_a_single_Set.narrowPeak")

  for (i in 1:length(grsToExport)){
    if (useOriginalNarrowPeak){
      myMat <- matrix(unlist(strsplit(names(grsToExport[[i]]), ":|-|_")),
                      ncol = 4, byrow = T)
      GenomicRanges::start(grsToExport[[i]]) <- as.numeric(myMat[, 2])
      GenomicRanges::end(grsToExport[[i]]) <- as.numeric(myMat[, 3])
    }
    dfNarrowPeak <- usefulLDfunctionsGR::narrowPeakDFFromGR(grsToExport[[i]])
    dfNarrowPeak <- dfNarrowPeak[order(dfNarrowPeak$seqnames,
                                       dfNarrowPeak$start), ]
    cat("track type=narrowPeak\n", file = outputFileNames[[i]])
    write.table(dfNarrowPeak,
                file = outputFileNames[[i]],
                row.names = F, col.names = F, quote = F, sep = "\t",
                append = T)
  }
}

#' Export myGRs with all annotations in tab-delimited text files.
#'
#' @param myGRAndAttributes Should be the output of \link[analysePeaks]{createMyGRs} or \link[analysePeaks]{annotateWithCate}
#' @param outputDirectory a path to a folder where will be written the text files.
#' @param useOriginalNarrowPeak logical to specify if instead of the coordinates in myGRAndAttributes which are by default the summit +/-150bp. You want the original coordinates as in the input (default is F).
#' @return Does not return anything but write 2 files into the `outputDirectory`:
#' Ref_table.txt
#' Set_table.txt
#' @details The coordinates are 1-based closed.
#' @importFrom GenomicRanges start end
#' @export
writeTablesFromMyGRsAndAnnot <- function(myGRAndAttributes,
                                              outputDirectory,
                                              useOriginalNarrowPeak = F){
  if (!"myGRs" %in% names(myGRAndAttributes)){
    stop("myGRs is not part of the input myGRAndAttributes.")
  }
  myGRs <- myGRAndAttributes[["myGRs"]]
  dir.create(outputDirectory, showWarnings = F, recursive = T)
  outputFileNames <- list()
  outputFileNames[[1]] <- file.path(outputDirectory,
                                    "Set_table.txt")
  outputFileNames[[2]] <- file.path(outputDirectory,
                                    "Ref_table.txt")
  for (i in 1:2){
    if (useOriginalNarrowPeak){
      myMat <- matrix(unlist(strsplit(names(myGRs[[i]]), ":|-|_")),
                      ncol = 4, byrow = T)
      GenomicRanges::start(myGRs[[i]]) <- as.numeric(myMat[, 2])
      GenomicRanges::end(myGRs[[i]]) <- as.numeric(myMat[, 3])
    }
    df <- base::as.data.frame(myGRs[[i]])
    write.table(df, file = outputFileNames[[i]],
                row.names = F, col.names = T, quote = F, sep = "\t")
  }
}
