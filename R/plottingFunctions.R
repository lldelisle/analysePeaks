#' Plot a barplot from a named vector
#'
#' @param data a named vector where values are number and names are category1(which can have dot).category2(no dot possible here)
#' @param colorsForGraphs a named vector with the colors which should be used in the barplot (names should correspond to category2)
#' @param orderForCategory1 a vector with the names of the category1 in the order it should be plotted
#' @param legend a logical to precise if the legend should be added (default is TRUE)
#' @param args.legend list of arguments for the leged (default is topright inset=(-0.2, 0), bg= "white")
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
                                   args.legend = list(x = "topright",
                                                      inset = c(-0.2, 0),
                                                      bg = "white"),
                                   las = 2, ...){
  # require package reshape
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
  b <- graphics::barplot(mat, col = myColor, legend = legend,
                         args.legend = args.legend, las = las, ...)
  return(invisible(b))
}

#' Plot different plots to illustrate a comparison between 2 overlaps
#'
#' @param name1 the name of the first experiment to use
#' @param name2 the name of the second experiment to use
#' @param ovF the list with at least the filtered overlap for name1 and name2
#' @param allExpe a GRanges list with at least name1 and name2 items with a dataframe with a column called score
#' @param step a integer which will be used to bin the items of name1 or name2 into categories (default is 5000).
#' @param fontsize base fontsize for the heatmaps (default is 10)
#' @return Will do a lot of plots but do not return anything.
#' @details `ovF` can be obtained by putting into a list the result of \link[usefulLDfunctionsGR]{filterByScore},
#' The name of the item in the list should be name1____name2
#' @importFrom pheatmap pheatmap
#' @import stats
#' @importFrom grDevices rainbow rgb
#' @importFrom graphics plot
#' @export
plotPairwiseComparison <- function(name1, name2, ovF, allExpe, step = 5000, fontsize = 10){
  # e is the name of the experiment
  e <- paste(sort(c(name1, name2)), collapse = "____")
  # temp.df contains the filtered overlap between name1 and name2
  temp.df <- ovF[[e]]
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
  corTest <- cor.test(temp.dfNoNA[, 1], temp.dfNoNA[, 2])
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
      allExpe[[c]]$score[temp.df[!is.na(temp.df[, c]), c]]
  }
  corTestS <- cor.test(temp.dfScore[, 1], temp.dfScore[, 2])
  graphics::plot(temp.dfScore[, c(name1, name2)], pch = 16,
                 col = grDevices::rgb(0, 0, 0, 0.03),
                 main = "score correlation",
                 sub = paste0("cor=", format(corTestS$estimate, digits = 2)))
  graphics::plot(temp.dfScore[, c(name2, name1)], pch = 16,
                 col = grDevices::rgb(0, 0, 0, 0.03),
                 main = "score correlation",
                 sub = paste0("cor=", format(corTestS$estimate, digits = 2)))
}


#' Plot histograms for the scores or the distance of features in the Set and in the Reference
#' As well as when it is in Set and in the Set but not in the Reference and vice-versa
#'
#' @param myGRAndAttributes Should be the output of \link[analysePeaks]{createMyGRs}
#' @return Plot histograms but do not return anything
#' @importFrom GenomicRanges mcols
#' @importFrom graphics hist legend
#' @importFrom grDevices rgb
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
    xmax <- min(max(allData), quantile(allData, probs = 0.75) +
                  1.5 * diff(quantile(allData, probs = c(0.25, 0.75))))
    # But sometimes the data is full of 0
    if (xmax != max(allData)){
      if (median(allData) == 0){
        xmax <- min(max(allData), quantile(allData[allData > 0], probs = 0.75) +
                      1.5 * diff(quantile(allData[allData > 0], probs = c(0.25, 0.75))))
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
        graphics::hist(GenomicRanges::mcols(myGRs[[iGR]])[, myCol],
                       freq = freqV,
                       breaks = seq(0, max(allData) + breakD, breakD),
                       xlim = c(0, xmax), col = grDevices::rgb(1, 0, 0, 0.5),
                       main = paste0(c("Frequency of\n",
                                       "")[1 + as.numeric(freqV)],
                                     what, "\n in ", c("Set", "Ref")[iGR]),
                       xlab = "",
                       sub = paste0("set:", myGRAndAttributes[["stringSet"]],
                                    " ref:", myGRAndAttributes[["nameOfRef"]]))
        if (iGR == 1){
          mySubsetLogic <-
            is.na(GenomicRanges::mcols(myGRs[[1]])[, "inUniqueRef"])
        } else {
          mySubsetLogic <- GenomicRanges::mcols(myGRs[[2]])[, "inNoneOfTheSet"]
        }
        graphics::hist(GenomicRanges::mcols(myGRs[[iGR]])[mySubsetLogic,
                                                          myCol],
                       breaks = seq(0, max(allData) + breakD, breakD),
                       freq = freqV,
                       col = grDevices::rgb(0, 1, 0, 0.5), add = T)
        graphics::legend(x = "topright", legend = c("all", "specific"),
                         fill = c(grDevices::rgb(1, 0, 0, 0.5),
                                  grDevices::rgb(0, 1, 0, 0.5)),
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
#' @return Plot but do not return anything
#' @importFrom graphics par barplot
#' @importFrom grDevices rainbow
#' @export
plotAllBarPlotForCategoriesFromMyGR <- function(myGRs, nameOfColWithCate,
                                                cateNames, stringSet,
                                                what, nameOfRef){
  if (length(myGRs) < 2){
    stop("Wrong myGRs\n")
  }
  if (!"inUniqueRef" %in% colnames(GenomicRanges::mcols(myGRs[[1]]))){
    stop("myGRs[[1]] does not contains inUniqueRef\n")
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
             is.na(myGRs[[1]]$inUniqueRef))
  colnames(t) <- c("sharedWithRef", "specific")
  maxChar <- max(sapply(colnames(t), nchar))
  graphics::par(mar = c(4 + maxChar / 3, 4, 4, 6), xpd = TRUE)
  # The table is displayed as barplot
  graphics::barplot(t, legend = T,
                    args.legend = list(x = "topright", inset = c(-0.2, 0),
                                       bg = "white"),
                    main = paste0(what, "\n for peaks of the Set\n"),
                    col = grDevices::rainbow(length(cateNames)),
                    sub = paste0("set:", stringSet, " ref:", nameOfRef))
  # Also as proportion
  graphics::barplot(prop.table(t, 2), legend = T,
                    args.legend = list(x = "topright", inset = c(-0.2, 0),
                                       bg = "white"),
                    main = paste0("Proportion of peaks in category for \n",
                                  what, "\n for peaks of the Set"),
                    col = grDevices::rainbow(length(cateNames)),
                    sub = paste0("set:", stringSet, " ref:", nameOfRef))
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
  colnames(t2) <- c("shared\nby all Set\nand Ref",
                    "shared\nby all Set\n no Ref",
                    "all shared\nby all Set",
                    "Ref\nin not a single\nof the set",
                    "all Ref")
  maxChar <- 15
  # We plot everything
  graphics::par(mar = c(4 + maxChar / 3, 4, 4, 6), xpd = TRUE)
  graphics::barplot(t2, legend = T, main = paste0(what, "\nset:", stringSet,
                                                  "\nref:", nameOfRef),
                    col = grDevices::rainbow(length(cateNames)),
                    args.legend = list(x = "topright", inset = c(-0.2, 0),
                                       bg = "white"),
                    las = 2)
  # Also as proportion
  graphics::par(mar = c(4 + maxChar / 3, 4, 5, 6), xpd = TRUE)
  graphics::barplot(prop.table(t2, 2), legend = T,
                    main = paste0("Proportion of peaks in category for \n",
                                  what, "\nset:", stringSet,
                                  "\nref:", nameOfRef),
                    col = grDevices::rainbow(length(cateNames)),
                    args.legend = list(x = "topright", inset = c(-0.2, 0),
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
#' @return Plot but do not return anything
#' @importFrom pheatmap pheatmap
#' @export
plotAllPheatmapsFor2CategoriesFromMyGR <- function(myGRs, nameOfColWithCate1,
                                                   cateNames1, what1,
                                                   nameOfColWithCate2,
                                                   cateNames2,
                                                   what2,
                                                   stringSet, nameOfRef,
                                                   fontsize = 10){
  if (length(myGRs) < 2){
    stop("Wrong myGRs\n")
  }
  if (!"inUniqueRef" %in%
      colnames(GenomicRanges::mcols(myGRs[[1]]))){
    stop("myGRs[[1]] does not contains inUniqueRef\n")
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
  inputsGR[[1]] <- subset(myGRs[[1]], !is.na(inUniqueRef))
  inputsGRWhat[[1]] <- "shared by all Set and Ref"
  # Second: specific to Set:
  inputsGR[[2]] <- subset(myGRs[[1]], is.na(inUniqueRef))
  inputsGRWhat[[2]] <- "shared by all Set no Ref"
  # Third: all Set
  inputsGR[[3]] <- myGRs[[1]]
  inputsGRWhat[[3]] <- "all shared by all Set"
  # Fourth: Ref specific
  inputsGR[[4]] <- subset(myGRs[[2]], inNoneOfTheSet)
  inputsGRWhat[[4]] <- "Ref in not a single of the set"
  # Fifth: Ref
  inputsGR[[5]] <- myGRs[[2]]
  inputsGRWhat[[5]] <- "all Ref"

  for (i in 1:length(inputsGR)){
    # We make the table using the 2 categories
    t1 <- table(factor(GenomicRanges::mcols(inputsGR[[i]])[,
                                                           nameOfColWithCate1],
                       levels = rev(cateNames1)),
                factor(GenomicRanges::mcols(inputsGR[[i]])[,
                                                           nameOfColWithCate2],
                       levels = rev(cateNames2)))
    # We plot it
    pheatmap::pheatmap(t1, cluster_rows = F, cluster_cols = F,
                       main = paste0(inputsGRWhat[[i]], "\n", sum(t1),
                                     " peaks\n", what1, " vs ", what2,
                                     "\nset:", stringSet, "\nref:", nameOfRef),
                       fontsize = fontsize)
  }
}



#' Plot barplots and pheatmaps for the scores or the distance of features in the Set and in the Reference
#' As well as when it is in Set and in the Set but not in the Reference and vice-versa
#'
#' @param myGRAndAttributes Should be the output of \link[analysePeaks]{annotateWithCate}
#' @param fontsize base fontsize for the heatmaps (default is 10)
#' @return Plot barplots and pheatmaps but do not return anything
#' @importFrom GenomicRanges mcols
#' @export
plotCateComparisonSetAndRef <- function(myGRAndAttributes, fontsize = 10){
  allCates <- myGRAndAttributes[["allCates"]]
  myGRs <- myGRAndAttributes[["myGRs"]]
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
                                          myGRAndAttributes[["nameOfRef"]])

  }
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
                                             fontsize = fontsize)
    }
  }
}
