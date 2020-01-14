#' Plot a barplot from a named vector
#'
#' @param data a named vector where values are number and names are category1(which can have dot).category2(no dot possible here)
#' @param colorsForGraphs a named vector with the colors which should be used in the barplot (names should correspond to category2)
#' @param orderForCategory1 a vector with the names of the category1 in the order it should be plotted
#' @param legend a logical to precise if the legend should be added (default is TRUE)
#' @param putLegendToTheTopRightBorder a logical to precise if the legend should be until the extremity of the page (default is TRUE)
#' @param args.legend list of arguments for the legend (default is topright inset=(-0.2, 0), bg= "white")
#' @param las a numeric value between 0 and 3 to specify the style of axis label (default is 2)
#' @param ... other arguments for barplot
#' @details Will plot a barplot with one bar per category1
#' @return invisible x values for the barplot
#' @importFrom reshape cast
#' @importFrom graphics par barplot
#' @export
#' @examples
#'topTs <- c(5, 10, 15, 20) * 1e3
#'colorsForGraphs<-grDevices::rainbow(2+length(topTs))
#'names(colorsForGraphs)<-c(paste0("overlapTop",sort(topTs/1e3),"k"),"overlap","specific")
#'v <- c(15809, rep(5000, 4), c(9316, 4762, 4327, 3831, 5051, 36638))
#'names(v) <- paste(c(rep("ChIP1", 5), rep("ChIP2", 6)),
#'                  c(rep(paste0("overlap", c("", "Top5k", "Top10k", "Top15k", "Top20k")), 2),
#'                    "specific"),
#'                  sep = ".")
#'barplotFromNamedVector(v, colorsForGraphs = colorsForGraphs)
barplotFromNamedVector <- function(data, colorsForGraphs = NULL,
                                   orderForCategory1 = NULL,
                                   legend = T,
                                   putLegendToTheTopRightBorder = T,
                                   args.legend = list(x = "topright",
                                                      inset = c(-0.2, 0),
                                                      bg = "white"),
                                   las = 2, ...){
  # Split the names of data:
  # first row: everything before the last dot (category1)
  # and second row: after the last dot (category2)
  namesCate <- sapply(names(data), function(s){
    v <- strsplit(s, "\\.")[[1]]
    n <- length(v)
    return(c(paste(v[1:(n - 1)], collapse = "."), v[n]))
  })
  # Put these names with the values in a dataframe
  temp.df <- data.frame(value = data, t(namesCate))
  # Reshape the dataframe to have a matrix
  mat.tmp <- reshape::cast(temp.df, X1 ~ X2)
  # Remove the first column which is the category1 and get the transpose
  mat <- t(mat.tmp[, 2:ncol(mat.tmp)])
  # Replace NA by 0 to be able to do the barplot
  mat[is.na(mat)] <- 0
  # Put right names
  rownames(mat) <- colnames(mat.tmp)[-1]
  colnames(mat) <- mat.tmp$X1
  if (!is.null(colorsForGraphs)){
    # If the colors are provided,
    # use them to select the rownames and order them.
    newOrder <- intersect(names(colorsForGraphs), rownames(mat))
    mat <- mat[newOrder, ]
    myColor <- colorsForGraphs[rownames(mat)]
  } else {
    myColor <- NULL
  }
  if (! is.null(orderForCategory1)) {
    # If the order is provided, use it
    mat <- mat[, intersect(orderForCategory1, colnames(mat))]
  }
  maxChar <- max(sapply(colnames(mat), nchar))
  graphics::par(mar = c(4 + maxChar / 3, 4, 4, 6), xpd = TRUE)
  if (putLegendToTheTopRightBorder){
    # I get the position of the plotting region
    v <- graphics::par()$plt
    # I evaluate how much I could shift the legend to the right:
    possibleShift <- (1 - v[2]) / (v[2] - v[1])
    if (is.null(args.legend)){
      args.legend <- list(x = "topright",
                          inset = c(-possibleShift, 0))
    } else {
      args.legend[["x"]] <- "topright"
      args.legend[["inset"]] <- c(-possibleShift, 0)
    }
  }
  b <- graphics::barplot(mat, col = myColor, legend = legend,
                         args.legend = args.legend, las = las, ...)
  return(invisible(b))
}

#' Plot different plots to illustrate a comparison between 2 overlaps
#'
#' @param name1 the name of the first experiment to use
#' @param name2 the name of the second experiment to use
#' @param ovF the list with at least the filtered overlap or the dataframe with cluster info for name1 and name2
#' @param allExpe a GRanges list with at least name1 and name2 items with a dataframe with a column called score
#' @param step a integer which will be used to bin the items of name1 or name2 into categories (default is 5000).
#' @param fontsize base fontsize for the heatmaps (default is 10)
#' @return Will do a lot of plots but do not return anything.
#' @details `ovF` can be obtained by putting into a list the result of \link[usefulLDfunctionsGR]{filterByScore},
#' The name of the item in the list should be name1____name2
#' @importFrom pheatmap pheatmap
#' @importFrom stats cor.test na.omit
#' @importFrom grDevices rainbow rgb
#' @importFrom graphics plot
#' @importFrom GenomicRanges mcols
#' @importFrom usefulLDfunctions subsetByNamesOrIndices
#' @export
plotPairwiseComparison <- function(name1, name2, ovF, allExpe,
                                   step = 5000, fontsize = 10){
  sortedNames <- sort(c(name1, name2))
  name1 <- sortedNames[1]
  name2 <- sortedNames[2]
  # e is the name of the experiment
  e <- paste(c(name1, name2), collapse = "____")
  # smallExpe is a subset of allExpe
  smallExpe <- usefulLDfunctions::subsetByNamesOrIndices(allExpe,
                                                         c(name1, name2))
  # temp.df contains the filtered overlap between name1 and name2
  temp.df <- ovF[[e]]
  if (is.null(temp.df)){
    stop(paste(e, "is not part of ovF."))
  }
  # If ovF was obtained by cluster, we need to apply some modifications:
  if (class(temp.df[, 1]) == "list"){
    temp.df <- as.data.frame(apply(temp.df, 2, function(l){
      sapply(l, function(v){
        ifelse(test = (length(v) == 0),
               yes = NA,
               no = min(v))
      })
    }))
    iSample1 <- sort(stats::na.omit(temp.df[, 1]))
    iSample2 <- sort(stats::na.omit(temp.df[, 2]))
    smallExpe[[1]] <- smallExpe[[1]][iSample1]
    smallExpe[[2]] <- smallExpe[[2]][iSample2]
    temp.df[!is.na(temp.df[, 1]), 1] <- rank(stats::na.omit(temp.df[, 1]))
    temp.df[!is.na(temp.df[, 2]), 2] <- rank(stats::na.omit(temp.df[, 2]))
  }
  maxValue <- max(temp.df, na.rm = T)
  steps <- seq(step, maxValue + step, step)
  cate <- rep(paste0(steps / 1e3, "k"), each = step)
  cate <- factor(cate, levels = c(unique(cate), "nonOverlap"))
  # To be able to make plots more easily
  # The exact item index will be converted
  # to a factor which contains its "binned" category.
  temp.dfCate <- temp.df
  for (i in 1:2){
    temp.dfCate[, i] <- cate[temp.dfCate[, i]]
    temp.dfCate[is.na(temp.dfCate[, i]), i] <- "nonOverlap"
  }
  # We will now summary the number of time each duo appears:
  t <- table(temp.dfCate)
  t <- t[rowSums(t) > 0, colSums(t) > 0]
  # Another statistics is the overlapping proportion taking
  # only the first category,
  # The two first etc...
  propOverlap <- sapply(1:(min(nrow(t), ncol(t)) - 2),
                        function(i){
                          sum(t[1:i, 1:i] / (step * i))
                        })
  # We plot it.
  # First keeping y extremities to 0 and 1
  graphics::plot( (1:length(propOverlap)) * step,
                  propOverlap, type = "b",
                  xlab = "top peaks considered",
                  ylab = "Overlapping proportion",
                  ylim = c(0, 1),
                  main = "Overlapping proportion using top peaks")
  # Then letting free the y axis.
  graphics::plot( (1:length(propOverlap)) * step, propOverlap, type = "b",
                  xlab = "top peaks considered",
                  ylab = "Overlapping proportion",
                  main = "Overlapping proportion using top peaks")
  # We now plot the proportion of each category of name1
  # for each category of name2
  graphics::plot(t(t), color = c(grDevices::rainbow( (nrow(t) - 1)), "white"),
                 main = paste0("correlation between the ranking of\n",
                               name1, " and ", name2))
  # We now plot the proportion of each category of name2
  # for each category of name1
  graphics::plot(t, color = c(grDevices::rainbow( (ncol(t) - 1)), "white"),
                 main = paste0("correlation between the ranking of\n",
                               name1, " and ", name2))
  # We will now plot the same but without proportion
  # just with colors with pheatmap
  tWithNames <- t
  rownames(tWithNames) <- paste0(name1, "_", rownames(tWithNames))
  colnames(tWithNames) <- paste0(name2, "_", colnames(tWithNames))
  pheatmap::pheatmap(t(tWithNames), cluster_cols = F, cluster_rows = F,
                     display_numbers = T, number_format = "%d",
                     main = paste0("correlation between the ranking of\n",
                                   name1, " and ", name2),
                     fontsize = fontsize)
  pheatmap::pheatmap(tWithNames, cluster_cols = F, cluster_rows = F,
                     display_numbers = T, number_format = "%d",
                     main = paste0("correlation between the ranking of\n",
                                   name1, " and ", name2),
                     fontsize = fontsize)
  # We will now plot the correlation without binning
  temp.dfNoNA <- temp.df
  for (c in colnames(temp.dfNoNA)){
    # When there is no overlap we arbitrary put
    # the rank to 1.2 x the maximum rank.
    temp.dfNoNA[is.na(temp.dfNoNA[, c]), c] <- max(temp.dfNoNA[, c],
                                                   na.rm = T) * 1.2
  }
  corTest <- stats::cor.test(temp.dfNoNA[, 1], temp.dfNoNA[, 2])
  graphics::plot(temp.dfNoNA[, c(name1, name2)], pch = 16,
                 col = grDevices::rgb(0, 0, 0, 0.03),
                 main = "rank correlation",
                 sub = paste0("cor=", format(corTest$estimate, digits = 2)))
  graphics::plot(temp.dfNoNA[, c(name2, name1)], pch = 16,
                 col = grDevices::rgb(0, 0, 0, 0.03),
                 main = "rank correlation",
                 sub = paste0("cor=", format(corTest$estimate, digits = 2)))
  # We will now correlate the scores
  temp.dfScore <- temp.df
  temp.dfScore[, ] <- 0
  for (c in colnames(temp.dfNoNA)){
    temp.dfScore[!is.na(temp.df[, c]), c] <-
      smallExpe[[c]]$score[temp.df[!is.na(temp.df[, c]), c]]
  }
  corTestS <- stats::cor.test(temp.dfScore[, 1], temp.dfScore[, 2])
  graphics::plot(temp.dfScore[, c(name1, name2)], pch = 16,
                 col = grDevices::rgb(0, 0, 0, 0.03),
                 main = "score correlation",
                 sub = paste0("cor=", format(corTestS$estimate, digits = 2)))
  graphics::plot(temp.dfScore[, c(name2, name1)], pch = 16,
                 col = grDevices::rgb(0, 0, 0, 0.03),
                 main = "score correlation",
                 sub = paste0("cor=", format(corTestS$estimate, digits = 2)))
  # And finally the signalValue if available
  temp.dfSValue <- temp.df
  temp.dfSValue[, ] <- 0
  allSValue <- TRUE
  for (c in colnames(temp.dfNoNA)){
    if ("signalValue" %in% colnames(GenomicRanges::mcols(allExpe[[c]]))){
      temp.dfSValue[!is.na(temp.df[, c]), c] <-
        smallExpe[[c]]$signalValue[temp.df[!is.na(temp.df[, c]), c]]
    } else {
      allSValue <- FALSE
    }
  }
  if (allSValue){
    corTestS <- stats::cor.test(temp.dfSValue[, 1], temp.dfSValue[, 2])
    graphics::plot(temp.dfSValue[, c(name1, name2)], pch = 16,
                   col = grDevices::rgb(0, 0, 0, 0.03),
                   main = "signal Value correlation",
                   sub = paste0("cor=", format(corTestS$estimate, digits = 2)))
    graphics::plot(temp.dfSValue[, c(name2, name1)], pch = 16,
                   col = grDevices::rgb(0, 0, 0, 0.03),
                   main = "signal Value correlation",
                   sub = paste0("cor=", format(corTestS$estimate, digits = 2)))
  }

}


#' Plot histograms for the scores or the distance of features in the Set and in the Reference
#' As well as when it is in Set and in the Set but not in the Reference and vice-versa
#'
#' @param myGRAndAttributes Should be the output of \link[analysePeaks]{createMyGRs} or \link[analysePeaks]{createMyGRsUsingSimpleOverlap}
#' @return Plot histograms but do not return anything
#' @importFrom GenomicRanges mcols
#' @importFrom graphics hist legend
#' @importFrom grDevices rgb
#' @importFrom stats quantile median
#' @export
plotClassicalHistogramForMyGRs <- function(myGRAndAttributes){
  myGRs <- myGRAndAttributes[["myGRs"]]
  # I deduce from the namesOfGrLScores and namesOfGrLDistance
  # The interests
  myInterests <- list()
  for (nameJ in myGRAndAttributes[["namesOfGrLScores"]]){
    whatString <- paste("score of", nameJ)
    if (myGRAndAttributes[["useSummitPMFlanking"]]){
      whatString <- paste(whatString,
                          "in summit +/-",
                          myGRAndAttributes[["flankingSize"]],
                          "bp")
    }
    myInterests[[paste0("s_", nameJ)]] <- list(
      what = whatString,
      columnName = paste0("best", nameJ, "Score")
    )
  }
  for (nameJ in myGRAndAttributes[["namesOfGrLDistance"]]){
    whatString <- paste("distance to the closest", nameJ)
    myInterests[[paste0("d_", nameJ)]] <- list(
      what = whatString,
      columnName = paste0("distanceToNearest", nameJ)
    )
  }
  # For each interest we do histograms
  for (i in 1:length(myInterests)){
    what <- myInterests[[i]][["what"]]
    myCol <- myInterests[[i]][["columnName"]]
    allData <- unlist(sapply(myGRs, function(gr){
      GenomicRanges::mcols(gr)[, myCol]
    }))
    # Usually I do not want to plot all values, only like a boxplot would do:
    xmax <- min(max(allData), stats::quantile(allData, probs = 0.75) +
                  1.5 * diff(stats::quantile(allData, probs = c(0.25, 0.75))))
    # But sometimes the data is full of 0
    if (xmax != max(allData)){
      if (stats::median(allData) == 0){
        xmax <- min(max(allData), stats::quantile(allData[allData > 0], probs = 0.75) +
                      1.5 * diff(stats::quantile(allData[allData > 0],
                                          probs = c(0.25, 0.75))))
      }
      # And sometimes it is still very low
      if (xmax == min(allData)){
        xmax <- max(allData)
      }
    }
    h <- graphics::hist(allData[allData <= xmax], plot = F)
    breakD <- h$breaks[2] - h$breaks[1]
    # I will now plot the histograms of the values
    for (iGR in 1:2){
      for (freqV in c(T, F)){
        if (freqV){
          alpha <- 1
          title <- ""
          firstColor <- grDevices::rgb(0, 0, 1, alpha)
          firstLegend <- "shared"
        } else {
          alpha <- 0.5
          title <- "Frequency of "
          firstColor <- grDevices::rgb(1, 0, 0, alpha)
          firstLegend <- "all"
        }
        graphics::hist(GenomicRanges::mcols(myGRs[[iGR]])[, myCol],
                       freq = freqV,
                       breaks = seq(0, max(allData) + breakD, breakD),
                       xlim = c(0, xmax),
                       col = firstColor,
                       main = paste0(title,
                                     what, "\n in ", c("Set", "Ref")[iGR]),
                       xlab = "",
                       sub = paste0("set:", myGRAndAttributes[["stringSet"]],
                                    " ref:", myGRAndAttributes[["nameOfRef"]]))
        if (iGR == 1){
          if ("inUniqueRef" %in% colnames(GenomicRanges::mcols(myGRs[[1]]))){
            mySubsetLogic <-
              is.na(GenomicRanges::mcols(myGRs[[1]])[, "inUniqueRef"])
          } else {
            mySubsetLogic <- GenomicRanges::mcols(myGRs[[1]])[,
                                                              "inNoneOfTheRef"]
          }
        } else {
          mySubsetLogic <- GenomicRanges::mcols(myGRs[[2]])[, "inNoneOfTheSet"]
        }
        graphics::hist(GenomicRanges::mcols(myGRs[[iGR]])[mySubsetLogic,
                                                          myCol],
                       breaks = seq(0, max(allData) + breakD, breakD),
                       freq = freqV,
                       col = grDevices::rgb(0, 1, 0, alpha),
                       add = T)
        graphics::legend(x = "topright",
                         legend = c(firstLegend,
                                    "specific"),
                         fill = c(firstColor,
                                  grDevices::rgb(0, 1, 0, alpha)),
                         bg = "white")
      }
    }
  }
}


#' Plot barplot for categories of a specific feature for a Set and a Ref
#'
#' @param myGRs a GRangeList of 2 elements which are the set and the ref and should have as meta column `nameOfColWithCate`
#' @param nameOfColWithCate a string which gives the column to display
#' @param cateNames the names of the categories, in the good order.
#' @param stringSet the name of the set
#' @param what the name of the measure
#' @param nameOfRef the name of the reference
#' @param only5categoriesBP logical whether to plot only the barplots with 5 categories (default is FALSE).
#' @return Plot but do not return anything
#' @importFrom graphics par barplot
#' @importFrom grDevices rainbow
#' @importFrom GenomicRanges mcols
#' @export
plotAllBarPlotForCategoriesFromMyGR <- function(myGRs, nameOfColWithCate,
                                                cateNames, stringSet,
                                                what, nameOfRef,
                                                only5categoriesBP = FALSE){
  if (length(myGRs) < 2){
    stop("Wrong myGRs\n")
  }
  if (!"inUniqueRef" %in% colnames(GenomicRanges::mcols(myGRs[[1]]))){
    if (!"inNoneOfTheRef" %in% colnames(GenomicRanges::mcols(myGRs[[1]]))){
      stop("myGRs[[1]] does not contains inUniqueRef nor inNoneOfTheRef\n")
    } else {
      mySubsetLogic <- GenomicRanges::mcols(myGRs[[1]])[, "inNoneOfTheRef"]
      sharedLabel <- "sharedWithAtLeastOneRef"
      t2Labels <- c("shared\nby all Set\nand at least\none Ref",
                    "shared\nby all Set\nnone of the Ref",
                    "shared\nby all Set\nall",
                    "shared\nby all Ref\nin not a single\nof the set",
                    "shared\nby all Ref\nall")
    }
  } else {
    mySubsetLogic <- is.na(myGRs[[1]]$inUniqueRef)
    sharedLabel <- "sharedWithRef"
    t2Labels <- c("shared\nby all Set\nand Ref",
                  "shared\nby all Set\n no Ref",
                  "all shared\nby all Set",
                  "Ref\nin not a single\nof the set",
                  "all Ref")
  }
  if (!"inNoneOfTheSet" %in% colnames(GenomicRanges::mcols(myGRs[[2]]))){
    stop("myGRs[[2]] does not contains inNoneOfTheSet\n")
  }
  for (i in 1:2){
    if (!nameOfColWithCate %in% colnames(GenomicRanges::mcols(myGRs[[i]]))){
      stop(nameOfColWithCate, " is not in myGRs[[", i, "]]\n")
    }
  }
  # We create a table with the values of the column in the Set
  # And the other dimension is the specificity of the interval/peak
  # As we put them as factors all categories will be displayed
  # even if they are empty
  t <- table(factor(GenomicRanges::mcols(myGRs[[1]])[, nameOfColWithCate],
                    levels = rev(cateNames)),
             mySubsetLogic)
  colnames(t) <- c(sharedLabel, "specific")
  maxChar <- max(sapply(colnames(t), nchar))
  graphics::par(mar = c(4 + maxChar / 3, 4, 4, 6), xpd = TRUE)
  # I get the position of the plotting region
  v <- graphics::par()$plt
  # I evaluate how much I could shift the legend to the right:
  possibleShift <- (1 - v[2]) / (v[2] - v[1])
  if ( !only5categoriesBP ){
    # The table is displayed as barplot
    graphics::barplot(t, legend = T,
                      args.legend = list(x = "topright",
                                         inset = c(-possibleShift, 0),
                                         bg = "white"),
                      main = paste0(what, "\n for peaks of the Set\n"),
                      col = grDevices::rainbow(length(cateNames)),
                      sub = paste0("set:", stringSet, " ref:", nameOfRef))
    # Also as proportion
    graphics::barplot(prop.table(t, 2), legend = T,
                      args.legend = list(x = "topright",
                                         inset = c(-possibleShift, 0),
                                         bg = "white"),
                      main = paste0("Proportion of peaks in category for \n",
                                    what, "\n for peaks of the Set"),
                      col = grDevices::rainbow(length(cateNames)),
                      sub = paste0("set:", stringSet, " ref:", nameOfRef))
  }
  t2 <- t
  # We do another table without splitting for the specificity
  temp.t <- table(factor(GenomicRanges::mcols(myGRs[[1]])[, nameOfColWithCate],
                         levels = rev(cateNames)))
  t2 <- cbind(t2, temp.t)
  # We now add the info for the ref
  for (NinNoneOfTheSet in c(F, NA)){
    temp.t <- table(
      factor(GenomicRanges::mcols(subset(
        myGRs[[2]],
        !inNoneOfTheSet %in% NinNoneOfTheSet))[, nameOfColWithCate],
        levels = rev(cateNames)))
    t2 <- cbind(t2, temp.t)
  }
  colnames(t2) <- t2Labels
  maxChar <- 15
  # We plot everything
  graphics::par(mar = c(4 + maxChar / 3, 4, 4, 6), xpd = TRUE)
  # I get the position of the plotting region
  v <- graphics::par()$plt
  # I evaluate how much I could shift the legend to the right:
  possibleShift <- (1 - v[2]) / (v[2] - v[1])
  graphics::barplot(t2, legend = T, main = paste0(what, "\nset:", stringSet,
                                                  "\nref:", nameOfRef),
                    col = grDevices::rainbow(length(cateNames)),
                    args.legend = list(x = "topright",
                                       inset = c(- possibleShift, 0),
                                       bg = "white"),
                    las = 2)
  # Also as proportion
  graphics::barplot(prop.table(t2, 2), legend = T,
                    main = paste0("Proportion of peaks in category for \n",
                                  what, "\nset:", stringSet,
                                  "\nref:", nameOfRef),
                    col = grDevices::rainbow(length(cateNames)),
                    args.legend = list(x = "topright",
                                       inset = c(-possibleShift, 0),
                                       bg = "white"),
                    las = 2)
}

#' Plot pheatmaps for categories of 2 specific features for a Set and a Ref
#'
#' @param myGRs a GRangeList of 2 elements which are the set and the ref and should have 2 meta columns `nameOfColWithCate1` and `nameOfColWithCate2`
#' @param nameOfColWithCate1 a string which gives the column to display
#' @param cateNames1 the names of the categories, in the good order, for the first feature
#' @param nameOfColWithCate2 a string which gives the column to display for the first feature
#' @param what1 the name of the first measure
#' @param cateNames2 the names of the categories, in the good order, for the second feature
#' @param stringSet the name of the set for the second feature
#' @param what2 the name of the second measure
#' @param nameOfRef the name of the reference
#' @param fontsize base fontsize for the heatmaps (default is 10)
#' @param display_numbers logical whether to display the numbers in the pheatmaps (default is TRUE)
#' @param whichPheatmaps a vector containing the pheatmaps to plot between 1 to 5 with:
#' 1=shared by all Set and at least one Ref or shared by all Set and Ref
#' 2=shared by all Set none of the Ref or shared by all Set no Ref
#' 3=shared by all Set: all or all shared by all Set
#' 4=shared by all Ref in not a single of the set or Ref in not a single of the set
#' 5=shared by all Ref: all or all Ref
#' (default is 1:5)
#' @param plotProportion logical whether to display the pheatmap 1, 2, and 4 should be plot as proportion of pheatmap 3 and 5 (default is FALSE)
#' @return Plot but do not return anything
#' @importFrom pheatmap pheatmap
#' @export
plotAllPheatmapsFor2CategoriesFromMyGR <- function(myGRs, nameOfColWithCate1,
                                                   cateNames1, what1,
                                                   nameOfColWithCate2,
                                                   cateNames2,
                                                   what2,
                                                   stringSet, nameOfRef,
                                                   fontsize = 10,
                                                   display_numbers = T,
                                                   whichPheatmaps = 1:5,
                                                   plotProportion = F){
  if (length(myGRs) < 2){
    stop("Wrong myGRs\n")
  }
  if (any(whichPheatmaps > 5) | any(whichPheatmaps < 1)){
    stop("Invalid whichPheatmap:", whichPheatmaps, " should be between 1 and 5.")
  }
  if (!"inUniqueRef" %in% colnames(GenomicRanges::mcols(myGRs[[1]]))){
    if (!"inNoneOfTheRef" %in% colnames(GenomicRanges::mcols(myGRs[[1]]))){
      stop("myGRs[[1]] does not contains inUniqueRef nor inNoneOfTheRef\n")
    } else {
      mySubsetLogic <- GenomicRanges::mcols(myGRs[[1]])[, "inNoneOfTheRef"]
      labels <- c("shared by all Set and at least one Ref",
                  "shared by all Set none of the Ref",
                  "shared by all Set: all",
                  "shared by all Ref in not a single of the set",
                  "shared by all Ref: all")
    }
  } else {
    mySubsetLogic <- is.na(myGRs[[1]]$inUniqueRef)
    labels <- c("shared by all Set and Ref",
                "shared by all Set no Ref",
                "all shared by all Set",
                "Ref in not a single of the set",
                "all Ref")
  }
  if (!"inNoneOfTheSet" %in%
      colnames(GenomicRanges::mcols(myGRs[[2]]))){
    stop("myGRs[[2]] does not contains inNoneOfTheSet\n")
  }
  for (i in 1:2){
    if (!nameOfColWithCate1 %in%
        colnames(GenomicRanges::mcols(myGRs[[i]]))){
      stop(nameOfColWithCate1, " is not in myGRs[[", i, "]]\n")
    }
    if (!nameOfColWithCate2 %in%
        colnames(GenomicRanges::mcols(myGRs[[i]]))){
      stop(nameOfColWithCate2, " is not in myGRs[[", i, "]]\n")
    }
  }
  inputsGR <- list()
  inputsGRWhat <- list()
  # We want to do 5 pheatmaps:
  # First: shared by all Set and Ref
  inputsGR[[1]] <- subset(myGRs[[1]], ! mySubsetLogic)
  inputsGRWhat[[1]] <- labels[1]
  # Second: specific to Set:
  inputsGR[[2]] <- subset(myGRs[[1]], mySubsetLogic)
  inputsGRWhat[[2]] <- labels[2]
  # Third: all Set
  inputsGR[[3]] <- myGRs[[1]]
  inputsGRWhat[[3]] <- labels[3]
  # Fourth: Ref specific
  inputsGR[[4]] <- subset(myGRs[[2]], inNoneOfTheSet)
  inputsGRWhat[[4]] <- labels[4]
  # Fifth: Ref
  inputsGR[[5]] <- myGRs[[2]]
  inputsGRWhat[[5]] <- labels[5]

  for (i in whichPheatmaps){
    # We make the table using the 2 categories
    t1 <- table(factor(GenomicRanges::mcols(inputsGR[[i]])[,
                                                           nameOfColWithCate1],
                       levels = rev(cateNames1)),
                factor(GenomicRanges::mcols(inputsGR[[i]])[,
                                                           nameOfColWithCate2],
                       levels = rev(cateNames2)))
    if(plotProportion && ! i %in% c(3, 5)){
      if(i < 3){
        t2 <- table(factor(GenomicRanges::mcols(inputsGR[[3]])[,
                                                               nameOfColWithCate1],
                           levels = rev(cateNames1)),
                    factor(GenomicRanges::mcols(inputsGR[[3]])[,
                                                               nameOfColWithCate2],
                           levels = rev(cateNames2)))
      } else {
        t2 <- table(factor(GenomicRanges::mcols(inputsGR[[5]])[,
                                                               nameOfColWithCate1],
                           levels = rev(cateNames1)),
                    factor(GenomicRanges::mcols(inputsGR[[5]])[,
                                                               nameOfColWithCate2],
                           levels = rev(cateNames2)))

      }
      line2 <- paste0(round(sum(t1) / sum(t2) * 100), "% of peaks")
      t1 <- t1 / t2
      t1[is.nan(t1)] <- 0
      number_format <- "%.2f"
    } else {
      line2 <- paste0(sum(t1), " peaks")
      number_format <- "%d"
    }
    # We plot it
    pheatmap::pheatmap(t1, cluster_rows = F, cluster_cols = F,
                       main = paste0(inputsGRWhat[[i]], "\n", line2, "\n", what1, " vs ", what2,
                                     "\nset:", stringSet, "\nref:", nameOfRef),
                       display_numbers = T,
                       number_format = number_format,
                       fontsize = fontsize)
  }
}



#' Plot barplots and pheatmaps for the scores or the distance of features in the Set and in the Reference
#' As well as when it is in Set and in the Set but not in the Reference and vice-versa
#'
#' @param myGRAndAttributes Should be the output of \link[analysePeaks]{annotateWithCate}
#' @param fontsize base fontsize for the heatmaps (default is 10)
#' @param plotBarPlots logical whether to plot the barplots (default is TRUE)
#' @param only5categoriesBP logical whether to plot only the barplots with 5 categories (default is FALSE).
#' @param plotPheatmaps logical whether to plot the pheatmaps (default is TRUE)
#' @param whichPheatmaps a vector containing the pheatmaps to plot between 1 to 5 with:
#' 1=shared by all Set and at least one Ref or shared by all Set and Ref
#' 2=shared by all Set none of the Ref or shared by all Set no Ref
#' 3=shared by all Set: all or all shared by all Set
#' 4=shared by all Ref in not a single of the set or Ref in not a single of the set
#' 5=shared by all Ref: all or all Ref
#' (default is 1:5)
#' @param plotProportion logical whether to display the pheatmap 1, 2, and 4 should be plot as proportion of pheatmap 3 and 5 (default is FALSE)
#' @param allCates a vector of string with the categories to plot, if NULL all categories in `myGRAndAttributes` are used (default is NULL)
#' @param display_numbers logical whether to display the numbers in the pheatmaps (default is TRUE)
#' @return Plot barplots and pheatmaps but do not return anything
#' @importFrom GenomicRanges mcols
#' @export
plotCateComparisonSetAndRef <- function(myGRAndAttributes, fontsize = 10,
                                        plotBarPlots = TRUE,
                                        only5categoriesBP = FALSE,
                                        plotPheatmaps = TRUE,
                                        display_numbers = TRUE,
                                        whichPheatmaps = 1:5,
                                        allCates = NULL,
                                        plotProportion = FALSE){
  if (is.null(allCates)){
    allCates <- myGRAndAttributes[["allCates"]]
  }
  myGRs <- myGRAndAttributes[["myGRs"]]
  if (plotBarPlots){
    for (i in 1:length(allCates)){
      plotAllBarPlotForCategoriesFromMyGR(myGRs = myGRs,
                                          nameOfColWithCate =
                                            allCates[[i]][["nameOfCol"]],
                                          cateNames =
                                            allCates[[i]][["cateNames"]],
                                          stringSet =
                                            myGRAndAttributes[["stringSet"]],
                                          what = allCates[[i]][["what"]],
                                          nameOfRef =
                                            myGRAndAttributes[["nameOfRef"]],
                                          only5categoriesBP = only5categoriesBP)

    }
  }
  if (plotPheatmaps){
    if (length(allCates) > 1){
      for (i in 1:length(allCates)){
        for (j in setdiff(1:length(allCates), i)){
          plotAllPheatmapsFor2CategoriesFromMyGR(myGRs,
                                                 allCates[[i]][["nameOfCol"]],
                                                 allCates[[i]][["cateNames"]],
                                                 allCates[[i]][["what"]],
                                                 allCates[[j]][["nameOfCol"]],
                                                 allCates[[j]][["cateNames"]],
                                                 allCates[[j]][["what"]],
                                                 myGRAndAttributes[["stringSet"]],
                                                 myGRAndAttributes[["nameOfRef"]],
                                                 fontsize = fontsize,
                                                 display_numbers = display_numbers,
                                                 whichPheatmaps = whichPheatmaps,
                                                 plotProportion = plotProportion)
        }
      }
    }
  }
}
