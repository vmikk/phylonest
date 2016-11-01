# Permutation test for wapqe
#' @title Permutation test for the apportionment of quadratic entropy (\code{\link{wapqe}}).
#'
#' @param df Dataframe or matrix with sites as rows and species as columns. Entries are abundances of species within sites.
#' @param dis Dissimilarity among species (NULL or class 'dist').
#' @param structures Data frame that contains the name of the group (row) of an level (column) to which the site belongs. Sites in structures should be in the same order as in df. Default is NULL.
#' @param formula Quadratic entropy formula ("QE", "EDI"). "QE" is default.
#' @param wopt Site weighting type ("even", "speciesab"). Default is "even".
#' @param level Level to test. Provide a number between 1 and 1+the number of columns in structures. The number is discarded if the parameter 'structures' is set to NULL.
#' @param nrep The number of permutations.
#' @param alter Alternative hypothesis type ("greater" (default), "less" or "two-sided").
#' @param tol A tolerance threshold (a value less than tol is considered equal to zero).
#'
#' @details
#' Level:
#' If structures is different from NULL then 1 means test for differences among sites, within the levels of the first factor given in parameter 'structures' (column 1), 2 means test for differences among levels of the first factor given in parameter ‘structures’ (column 1) but within levels of the second factor given in parameter 'structures' (column 2) (if available), etc.
#'
#' @return A list of class 'randtest', see \code{\link{randtest}}.
#' @author Sandrine Pavoine, Eric Marcon, Carlo Ricotta.
#' @references Pavoine, S., Marcon, E. and Ricotta, C. (2016), ‘Equivalent numbers’ for species, phylogenetic or functional diversity in a nested hierarchy of multiple scales. Methods Ecol Evol, 7: 1152–1163. DOI:10.1111/2041-210X.12591
#' @seealso \code{\link{wapqe}}, \code{\link{randtest}}.
#'
#' @examples
#' data(macroloire)
#' # Taxonomic dissimilarities among species:
#' dTaxo <- dist.taxo(macroloire$taxo)^2/2
#' dTaxo <- dTaxo/max(dTaxo)
#' # Size-based dissimilarities among species
#' dSize <- dist.prop(macroloire$traits[ ,1:4], method = 2)
#' # Dissimilarities among species in terms of feeding categories
#' dFeed <- dist.prop(macroloire$traits[ ,5:11], method = 2)
#' # Dissimilarities among species in terms of both size and feeding categories
#' dSF <- (dSize+dFeed)/2
#' 
#' # Table with sites as rows (stations), species as columns and abundances as entries
#' ab <- as.data.frame(t(macroloire$fau))
#' # Table with sites as rows and one column only. Entries indicate the geological region associated with each site
#' stru <- macroloire$envir["Morphoregion"]
#' 
#' # Tests for dissimilarities among sites within regions:
#' 
#' aGS1 <- randtestapqe(ab, , stru, formula = "QE", level=1, nrep=999)
#' aGS1
#' plot(aGS1)
#' aTaxo1 <- randtestapqe(ab, dTaxo, stru, formula = "QE", level=1, nrep=999)
#' aTaxo1
#' plot(aTaxo1)
#' aSize1 <- randtestapqe(ab, dSize, stru, formula = "QE", level=1, nrep=999)
#' aSize1
#' plot(aSize1)
#' aFeed1 <- randtestapqe(ab, dFeed, stru, formula = "QE", level=1, nrep=999)
#' aFeed1
#' plot(aFeed1)
#' aSF1 <- randtestapqe(ab, dSF, stru, formula = "QE", level=1, nrep=999)
#' aSF1
#' plot(aSF1)
#' 
#' # Tests for dissimilarities among regions:
#' 
#' aGS2 <- randtestapqe(ab, , stru, formula = "QE", level=2, nrep=999)
#' aGS2
#' plot(aGS2)
#' aTaxo2 <- randtestapqe(ab, dTaxo, stru, formula = "QE", level=2, nrep=999)
#' aTaxo2
#' plot(aTaxo2)
#' aSize2 <- randtestapqe(ab, dSize, stru, formula = "QE", level=2, nrep=999)
#' aSize2
#' plot(aSize2)
#' aFeed2 <- randtestapqe(ab, dFeed, stru, formula = "QE", level=2, nrep=999)
#' aFeed2
#' plot(aFeed2)
#' aSF2 <- randtestapqe(ab, dSF, stru, formula = "QE", level=2, nrep=999)
#' aSF2
#' plot(aSF2)
#' 
#' @export

randtestapqe <- function(df, dis = NULL, structures = NULL, formula = c("QE", "EDI"), 
                         wopt = c("even", "speciesab"), level = 1, nrep = 99,
                         alter = c("greater", "less", "two-sided"), tol = 0.00000001) {
  dfold <- df
  df <- df[rowSums(df) > 0, ]
  ncomm <- nrow(df)
  if (!is.null(structures)) {
    if (!inherits(structures, "data.frame")) 
      stop("structures should be a data frame or NULL")
    if (!nrow(structures) == nrow(dfold)) 
      stop("incorrect number of rows in structures")
    structures <- structures[rowSums(dfold) > 0, , drop = FALSE]
    structures <- as.data.frame(apply(structures, 2, factor))
    if (!is.null(rownames(structures)) & !is.null(rownames(df))) {
      e <- sum(abs(match(rownames(df), rownames(structures)) - (1:ncomm)))
      if (e > 0.00000001) 
        warning("be careful that rownames in df should be in the same order as rownames in structures")
    }
    if (ncol(structures) > 1) 
      .checknested(structures)
  }
  dfbrut <- df
  P <- as.data.frame(sweep(df, 1, rowSums(df), "/"))
  if (wopt[1] == "speciesab") {
    w <- rowSums(df)/sum(df)
  } else if (wopt[1] == "even") {
    if (is.null(structures)) 
      w <- rep(1/nrow(df), nrow(df)) else {
        nc <- ncol(structures)
        fun <- function(i) {
          x <- table(structures[, i], structures[, i - 1])
          x[x > 0] <- 1
          x <- rowSums(x)
          v <- x[structures[, i]]
          v <- 1/v
          return(v)
        }
        if (ncol(structures) == 1) {
          firstw <- table(structures[, 1])
          w <- 1/firstw[structures[, 1]]/length(levels(structures[, 1]))
        } else {
          listw <- lapply(2:nc, fun)
          firstw <- table(structures[, 1])
          firstw <- 1/firstw[structures[, 1]]
          finalw <- 1/length(levels(structures[, ncol(structures)]))
          forw <- cbind.data.frame(as.vector(firstw), as.vector(listw), as.vector(rep(finalw, nrow(structures))))
          w <- apply(forw, 1, prod)
        }
      }
    df <- P * w
  } else if (is.numeric(wopt) & length(wopt) == nrow(df) & sum(wopt) > tol) {
    if (!is.null(names(w)) & all(rownames(df) %in% w)) 
      w <- w[rownames(df)]
    w <- w/sum(w)
    if (any(w <= tol)) {
      warnings("sites with weights of zero in w have been removed")
      df <- df[w > tol, ]
      structures <- structures[w > tol, ]
      w <- w[w > tol]
      w <- w/sum(w)
    }
    df <- P * w
  } else stop("incorrect definition of wopt")
  if (is.null(dis)) {
    dis <- as.dist((matrix(1, ncol(df), ncol(df)) - diag(rep(1, ncol(df)))))
    attributes(dis)$Labels <- colnames(df)
    formula <- "QE"
  }
  if (!inherits(dis, "dist")) 
    stop("dis must be of class dist")
  if (!formula[1] %in% c("QE", "EDI")) 
    stop("formula can be either QE or EDI")
  if (any(!colnames(df) %in% attributes(dis)$Labels)) 
    stop("column names in df are missing in dis") else {
      d <- as.matrix(dis)[colnames(df), colnames(df)]
      if (formula[1] == "EDI") {
        op <- options()$warn
        options(warn = -1)
        if (!is.euclid(as.dist(d))) 
          stop("dis should be Euclidean")
        options(warn = op)
        d <- d^2/2  # Euclidean Diversity Index
      } else {
        op <- options()$warn
        options(warn = -1)
        if (!is.euclid(sqrt(as.dist(d)))) 
          stop("dis should be squared Euclidean")
        options(warn = op)
      }
    }
  d <- as.dist(d)
  alter <- alter[1]
  if (is.null(structures)) {
    fun <- function(i) {
      dfperm <- as.data.frame(sapply(dfbrut, sample))
      if (any(rowSums(dfperm) < tol)) 
        return(NA)
      if (wopt[1] != "speciesab") {
        dfperm <- as.data.frame(sweep(dfperm, 1, rowSums(dfperm), "/"))
        dfperm <- dfperm * w
      }
      op <- options()$warn
      options(warn = -1)
      a <- apqe(as.data.frame(t(dfperm)), dis = sqrt(2 * d), NULL)$results
      options(warn = op)
      res <- a[1, ]/a[2, ]
      return(res)
    }
    ressim <- sapply(1:nrep, fun)
    op <- options()$warn
    options(warn = -1)
    a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), NULL)$results
    options(warn = op)
    obs <- a[1, ]/a[2, ]
    sim <- ressim[!is.na(ressim)]
    res <- as.randtest(obs = obs, sim = sim, alter = alter)
    res$call <- match.call()
  } else if (level == 1) {
    aggr.permut <- function(x) {
      listval <- split(1:ncomm, as.factor(structures[, 1]))
      posiori <- as.vector(unlist(listval))
      fun0 <- function(v) {
        if (length(v) == 1) 
          return(v) else return(sample(v))
      }
      listval2 <- lapply(listval, fun0)
      posiori2 <- as.vector(unlist(listval2))
      x[posiori] <- x[posiori2]
      return(x)
    }
    fun <- function(i) {
      dfperm <- sapply(dfbrut, aggr.permut)
      rownames(dfperm) <- rownames(df)
      dfperm <- as.data.frame(dfperm)
      if (any(rowSums(dfperm) < tol)) 
        return(NA) else {
          if (wopt[1] != "speciesab") {
            dfperm <- as.data.frame(sweep(dfperm, 1, rowSums(dfperm), "/"))
            dfperm <- dfperm * w
          }
          op <- options()$warn
          options(warn = -1)
          a <- apqe(as.data.frame(t(dfperm)), dis = sqrt(2 * d), structures)$results
          options(warn = op)
          res <- a[nrow(a) - 2, ]/a[nrow(a) - 1, ]
          return(res)
        }
    }
    ressim <- sapply(1:nrep, fun)
    op <- options()$warn
    options(warn = -1)
    a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures)$results
    options(warn = op)
    obs <- a[nrow(a) - 2, ]/a[nrow(a) - 1, ]
    sim <- ressim[!is.na(ressim)]
    res <- as.randtest(obs = obs, sim = sim, alter = alter)
    res$call <- match.call()
  } else if ((level - 1) == ncol(structures) & level == 2) {
    fun <- function(i) {
      e <- sample(ncomm)
      strusim <- structures[e, , drop = FALSE]
      rownames(strusim) <- rownames(structures)
      op <- options()$warn
      options(warn = -1)
      a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures = strusim)$results
      options(warn = op)
      res <- a[nrow(a) - 1 - level, ]/a[nrow(a) - level, ]
      return(res)
    }
    ressim <- sapply(1:nrep, fun)
    op <- options()$warn
    options(warn = -1)
    a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures)$results
    options(warn = op)
    obs <- a[nrow(a) - 1 - level, ]/a[nrow(a) - level, ]
    res <- as.randtest(obs = obs, sim = ressim, alter = alter)
    res$call <- match.call()
  } else if ((level - 1) == ncol(structures)) {
    strulev <- structures[!duplicated(structures[level - 2]), level - 1]
    names(strulev) <- unique(structures[level - 2])
    fun <- function(i) {
      strusim <- structures
      strulevperm <- sample(strulev)
      names(strulevperm) <- names(strulev)
      strusim[, level - 1] <- strulevperm[strusim[, level - 2]]
      rownames(strusim) <- rownames(df)
      op <- options()$warn
      options(warn = -1)
      a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures = strusim)$results
      options(warn = op)
      res <- a[nrow(a) - 1 - level, ]/a[nrow(a) - level, ]
      return(res)
    }
    ressim <- sapply(1:nrep, fun)
    op <- options()$warn
    options(warn = -1)
    a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures)$results
    options(warn = op)
    obs <- a[nrow(a) - 1 - level, ]/a[nrow(a) - level, ]
    res <- as.randtest(obs = obs, sim = ressim, alter = alter)
    res$call <- match.call()
  } else if (level == 2) {
    strulev <- as.character(structures[, level - 1])
    names(strulev) <- paste("c", 1:ncomm)
    strulevsup <- structures[, level]
    listbase <- split(strulev, strulevsup)
    fun <- function(i) {
      fun0 <- function(x) {
        strulevperm <- sample(x)
        names(strulevperm) <- names(x)
        return(strulevperm)
      }
      listbase2 <- lapply(listbase, fun0)
      names(listbase2) <- NULL
      listbase2 <- unlist(listbase2)
      strusim <- structures
      strusim[, level - 1] <- as.factor(listbase2[names(strulev)])
      rownames(strusim) <- rownames(df)
      op <- options()$warn
      options(warn = -1)
      a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures = strusim)$results
      options(warn = op)
      res <- a[nrow(a) - 1 - level, ]/a[nrow(a) - level, ]
      return(res)
    }
    ressim <- sapply(1:nrep, fun)
    op <- options()$warn
    options(warn = -1)
    a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures)$results
    options(warn = op)
    obs <- a[nrow(a) - 1 - level, ]/a[nrow(a) - level, ]
    res <- as.randtest(obs = obs, sim = ressim, alter = alter)
    res$call <- match.call()
  } else {
    if ((level - 1) > ncol(structures)) 
      stop("level should be between 1 and ", ncol(structures) + 1)
    strulev <- as.character(structures[!duplicated(structures[level - 2]), level - 1])
    names(strulev) <- unique(structures[level - 2])
    strulevsup <- structures[!duplicated(structures[level - 2]), level]
    listbase <- split(strulev, strulevsup)
    fun <- function(i) {
      fun0 <- function(x) {
        strulevperm <- sample(x)
        names(strulevperm) <- names(x)
        return(strulevperm)
      }
      listbase2 <- lapply(listbase, fun0)
      names(listbase2) <- NULL
      listbase2 <- unlist(listbase2)
      strusim <- structures
      strusim[, level - 1] <- as.factor(listbase2[strusim[, level - 2]])
      rownames(strusim) <- rownames(df)
      op <- options()$warn
      options(warn = -1)
      a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures = strusim)$results
      options(warn = op)
      res <- a[nrow(a) - 1 - level, ]/a[nrow(a) - level, ]
      return(res)
    }
    ressim <- sapply(1:nrep, fun)
    op <- options()$warn
    options(warn = -1)
    a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures)$results
    options(warn = op)
    obs <- a[nrow(a) - 1 - level, ]/a[nrow(a) - level, ]
    res <- as.randtest(obs = obs, sim = ressim, alter = alter)
    res$call <- match.call()
  }
  return(res)
}