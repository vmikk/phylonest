# Permutation test for the first type of diversity decomposition
#' @title Permutation test for the first type of diversity decomposition (\code{\link{EqRS}}).
#'
#' @param df Dataframe or matrix with sites as rows and species as columns. Entries are abundances of species within sites.
#' @param dis Dissimilarity among species (NULL or class 'dist').
#' @param structures Data frame that contains the name of the group (row) of an level (column) to which the site belongs. Sites in structures should be in the same order as in df. Default is NULL.
#' @param formula Quadratic entropy formula ("QE", "EDI"). "QE" is default.
#' @param option Rescaling type ("eq", "normed1" or "normed2").
#' @param level Level to test. Provide a number between 1 and 1+the number of columns in structures. The number is discarded if the parameter 'structures' is set to NULL.
#' @param nrep The number of permutations.
#' @param alter Alternative hypothesis type ("greater" (default), "less" or "two-sided").
#' @param tol A tolerance threshold (a value less than tol is considered equal to zero).
#'
#' @details
#'
#' @return A list of class 'randtest', see \code{\link{randtest}}.
#' @author Sandrine Pavoine, Eric Marcon, Carlo Ricotta.
#' @references Pavoine, S., Marcon, E. and Ricotta, C. (2016), ‘Equivalent numbers’ for species, phylogenetic or functional diversity in a nested hierarchy of multiple scales. Methods Ecol Evol, 7: 1152–1163. DOI:10.1111/2041-210X.12591
#' @seealso \code{\link{EqRS}}, \code{\link{randtest}}.
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
#' 
#' @export

randtestEqRS <- function(df, dis = NULL, structures = NULL, formula = c("QE", "EDI"), 
                         option = c("normed1", "normed2", "eq"), level = 1, nrep = 99,
                         alter = c("greater", "less", "two-sided"), tol = 0.00000001) {
  if (option[1] == "eq") 
    indexk <- 0 else indexk <- 2
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
    alter <- alter[1]
    if (is.null(structures)) {
      fun <- function(i) {
        dfperm <- as.data.frame(sapply(df, sample))
        if (any(rowSums(dfperm) < tol)) 
          return(NA) else {
            res <- EqRS(dfperm, dis = dis, NULL, formula = formula, option = option, tol = tol)[1, ]
            return(res)
          }
      }
      ressim <- sapply(1:nrep, fun)
      obs <- EqRS(df, dis = dis, NULL, formula = formula, option = option, tol = tol)[1, ]
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
        dfperm <- sapply(df, aggr.permut)
        rownames(dfperm) <- rownames(df)
        dfperm <- as.data.frame(dfperm)
        if (any(rowSums(dfperm) < tol)) 
          return(NA) else {
            res <- EqRS(dfperm, dis = dis, structures, formula = formula, option = option, tol = tol)
            res <- res[nrow(res) - 2 + indexk, ]
            return(res)
          }
      }
      ressim <- sapply(1:nrep, fun)
      obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
      obs <- obs[nrow(obs) - 2 + indexk, ]
      sim <- ressim[!is.na(ressim)]
      res <- as.randtest(obs = obs, sim = sim, alter = alter)
      res$call <- match.call()
    } else if ((level - 1) == ncol(structures) & level == 2) {
      fun <- function(i) {
        e <- sample(ncomm)
        strusim <- structures[e, , drop = FALSE]
        rownames(strusim) <- rownames(df)
        res <- EqRS(df, dis = dis, structures = strusim, formula = formula, option = option, tol = tol)
        res <- res[1, ]
        return(res)
      }
      ressim <- sapply(1:nrep, fun)
      obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
      obs <- obs[1, ]
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
        res <- EqRS(df, dis = dis, structures = strusim, formula = formula, option = option, tol = tol)
        res <- res[1, ]
        return(res)
      }
      ressim <- sapply(1:nrep, fun)
      obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
      obs <- obs[1, ]
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
        res <- EqRS(df, dis = dis, structures = strusim, formula = formula, option = option, tol = tol)
        res <- res[nrow(res) - 3 + indexk, ]
        return(res)
      }
      ressim <- sapply(1:nrep, fun)
      obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
      obs <- obs[nrow(obs) - 3 + indexk, ]
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
        res <- EqRS(df, dis = dis, structures = strusim, formula = formula, option = option, tol = tol)
        res <- res[nrow(res) - 2 + indexk - level + 1, ]
        return(res)
      }
      ressim <- sapply(1:nrep, fun)
      obs <- EqRS(df, dis = dis, structures, formula = formula, option = option, tol = tol)
      obs <- obs[nrow(obs) - 2 + indexk - level + 1, ]
      res <- as.randtest(obs = obs, sim = ressim, alter = alter)
      res$call <- match.call()
    }
    return(res)
}