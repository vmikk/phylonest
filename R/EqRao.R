# Third decomposition of diversity
#' @title Third decomposition of diversity introduced in Pavoine et al. 2016.
#' @description Use when the interest is in beta diversity expressed in terms of pairwise dissimilarities among sites (within regions) and among regions and/or sites within a region have different weights and regions have different weights (due, e.g., to uneven sampling pressures, uneven size of sites or regions).
#'
#' @param df Dataframe or matrix with sites as rows and species as columns. Entries are abundances of species within sites.
#' @param dis Dissimilarity among species (NULL or class 'dist').
#' @param structures Data frame that contains the name of the group (row) of an level (column) to which the site belongs. Sites in structures should be in the same order as in df. Default is NULL.
#' @param option Rescaling type ("eq", "normed1" or "normed2").
#' @param formula Quadratic entropy formula ("QE", "EDI"). "QE" is default.
#' @param wopt Site weighting type ("even", "speciesab"). Default is "even".
#' @param tol A tolerance threshold (a value less than tol is considered equal to zero).
#' @param metmean Mean type - "arithmetic" or "harmonic" (default).
#'
#' @details
#' Rescaling types:
#' - "eq" - the diversity components are given in terms of equivalent number of species, sites, regions etc.
#' - "normed1" - the normed components of diversity will be returned with formula (1 – 1 / E) / (1 - 1 / Emax)
#' - "normed2" - the normed components of diversity will be returned with formula (E – 1) / (Emax - 1).
#' For Eα and Eγ, Emax=S (the number of species in the data set).
#'
#' Formula type:
#' 
#'
#'
#' Site weighting type:
#' If wopt = "speciesab", then the sites will be weighted by their sum of species’ abundances.
#' If wopt ="even", the sites will be evenly weighted within the factors defined by the parameter 'structures'.
#'
#' For the associated permutation test see \code{\link{randtestEqRao}}.
#'
#' @return A data frame with each component of the selected diversity decomposition.
#' @author Sandrine Pavoine, Eric Marcon, Carlo Ricotta.
#' @references Pavoine, S., Marcon, E. and Ricotta, C. (2016), ‘Equivalent numbers’ for species, phylogenetic or functional diversity in a nested hierarchy of multiple scales. Methods Ecol Evol. DOI:10.1111/2041-210X.12591
#' @seealso \code{\link{randtestEqRao}}, \code{\link{EqRS}}, \code{\link{EqRSintra}}.
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
#' EqRao(ab, , stru, option="eq")
#' EqRao(ab, dTaxo, stru, formula = "QE", option="eq")
#' EqRao(ab, dSize, stru, formula = "QE", option="eq")
#' EqRao(ab, dFeed, stru, formula = "QE", option="eq")
#' EqRao(ab, dSF, stru, formula = "QE", option="eq")
#' 
#' EqRao(ab, , stru, option="normed2")
#' EqRao(ab, dTaxo, stru, formula = "QE", option="normed2")
#' EqRao(ab, dSize, stru, formula = "QE", option="normed2")
#' EqRao(ab, dFeed, stru, formula = "QE", option="normed2")
#' EqRao(ab, dSF, stru, formula = "QE", option="normed2")
#'
#' @export

EqRao <- function(df, dis = NULL, structures = NULL, option = c("eq", "normed1", "normed2"), 
                  formula = c("QE", "EDI"), wopt = c("even", "speciesab"), tol = 0.00000001,
                  metmean = c("harmonic", "arithmetic")) {
  metmean <- metmean[1]
  if (!option[1] %in% c("eq", "normed1", "normed2")) 
    stop("unavailable option, please modify; option can be eq, normed1 or normed2")
  if (!(is.data.frame(df) | is.matrix(df))) 
    stop("df must be a matrix or a data frame")
  if (!(is.data.frame(df))) 
    df <- as.data.frame(df)
  if (is.null(colnames(df)) | is.null(rownames(df))) 
    stop("df must have both row names and column names")
  if (length(colnames(df)[colSums(df) > 0]) == 1) 
    stop("df must contain at least two species with positive occurrence in at least one site")
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
      if (any(d > 1)) 
        d <- d/max(d)
    }
  nsp <- ncol(df)
  if (any(rowSums(df) < tol)) 
    warnings("empty sites have been removed")
  if (sum(df) < tol) 
    stop("df must contain nonnegative values and the sum of its values must be positive")
  dfold <- df
  df <- df[rowSums(df) > tol, ]
  ncomm <- nrow(df)
  if (nrow(df) == 1) 
    stop("df must contain at least two non-empty sites")
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
          forw <- cbind.data.frame(firstw, listw, rep(finalw, nrow(structures)))
          w <- apply(forw, 1, prod)
        }
      }
    df <- P * w
  } else if (is.numeric(wopt) & length(wopt) == nrow(df) & sum(wopt) > tol) {
    if (!is.null(names(wopt)) & all(rownames(df) %in% wopt)) 
      w <- wopt[rownames(df)]
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
  ncomm <- nrow(df)
  if (!is.null(structures)) {
    if (length(levels(factor(structures[, 1]))) == 1) 
      stop("All sites belong to a unique level in the first column of structures, remove this first column in structures")
    if (length(levels(factor(structures[, 1]))) == nrow(df)) 
      stop("Each site belongs to a distinct level in the first column of structures, this first column is useless, remove it and re-run")
  }

  discintra <- function(samples) {
    Ng <- as.matrix(samples)
    lesnoms <- colnames(samples)
    sumcol <- colSums(Ng)
    Lg <- t(t(Ng)/sumcol)
    colnames(Lg) <- lesnoms
    Pg <- as.matrix(apply(Ng, 2, sum)/sum(Ng))
    rownames(Pg) <- lesnoms
    deltag <- as.matrix(apply(Lg, 2, function(x) t(x) %*% d %*% x))
    ug <- matrix(1, ncol(Lg), 1)
    numdg2 <- t(Lg) %*% d %*% Lg - 1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
    numdg2[numdg2 < tol] <- 0
    dendg2 <- 1 - 1/2 * (deltag %*% t(ug) + ug %*% t(deltag))
    dendg2[dendg2 < tol] <- 1
    dg2 <- numdg2/dendg2
    dg2 <- dg2 - diag(diag(dg2))
    colnames(dg2) <- lesnoms
    rownames(dg2) <- lesnoms
    return(dg2)
  }
  
  if (is.null(structures)) {
    div <- sapply(as.data.frame(t(P)), .diversity, d=d)
    if (metmean == "arithmetic") 
      alphaEq <- sum(w * (1/(1 - div))) else alphaEq <- 1/(1 - sum(w * div))
      beta <- t(w) %*% discintra(as.data.frame(t(P))) %*% w
      beta <- 1/(1 - beta)
      gamma <- alphaEq * beta
      if (option[1] == "eq") {
        res <- cbind.data.frame(c(beta, alphaEq, gamma))
        rownames(res) <- c("beta", "alpha", "gamma")
        colnames(res) <- "Equivalent numbers"
      } else if (option[1] == "normed2") {
        beta <- (beta - 1)/(ncomm - 1)
        res <- cbind.data.frame(beta)
        rownames(res) <- c("beta")
        colnames(res) <- "Normed contributions to diversity"
      } else if (option[1] == "normed1") {
        beta <- (1 - 1/beta)/(1 - 1/ncomm)
        res <- cbind.data.frame(beta)
        rownames(res) <- c("beta")
        colnames(res) <- "Normed contributions to diversity"
      }
  } else {
    div <- sapply(as.data.frame(t(P)), .diversity, d=d)
    if (metmean == "arithmetic") 
      alphaEq <- sum(w * (1/(1 - div))) else alphaEq <- 1/(1 - sum(w * div))
      if (!(is.data.frame(structures))) 
        stop("structures should be either NULL or a data frame")
      if (any(!rownames(df) %in% rownames(structures))) 
        stop("row names in df are missing in structures")
      structures <- structures[rownames(df), , drop = FALSE]
      nc <- ncol(structures)
      if (nc > 1) {
        namrow <- rownames(structures)
        structures <- as.data.frame(sapply(structures, factor))
        rownames(structures) <- namrow
      } else {
        namfac <- colnames(structures)
        namrow <- rownames(structures)
        structures <- as.data.frame(as.factor(structures[, 1]))
        colnames(structures) <- namfac
        rownames(structures) <- namrow
      }
      if (nc > 1) {
        for (i in 1:(nc - 1)) {
          tab <- table(structures[, c(i, i + 1)])
          tab2 <- apply(tab, 1, function(x) sum(x != 0))
          if (any(tab2 != 1)) {
            stop(paste("non hierarchical design in structures, column", i, "in", i + 1))
          }
          if (length(levels(structures[, i])) == length(levels(structures[, i + 1]))) 
            stop("In structures, column ", i, " and ", i + 1, " are equivalent, remove one of them and rerun the function")
        }
      } else {
        if (length(levels(structures)) == nrow(df)) 
          stop("Leave structures as a NULL object because there are as many different levels in structures than there are rows in df; structures is uninformative")
        if (length(levels(structures)) == 1) 
          stop("Leave structures as a NULL object because there is only one level in structures; structures is uninformative")
      }
      W <- lapply(structures, function(x) tapply(w, x, sum))
      df1 <- sapply(df, function(x) tapply(x, structures[, nc], sum))
      P1 <- sweep(df1, 1, rowSums(df1), "/")
      P1 <- as.data.frame(t(P1))
      beta1 <- t(W[[nc]]) %*% discintra(P1) %*% W[[nc]]
      beta1Eq <- 1/(1 - beta1)
      funbeta <- function(i) {
        if (i == 1) {
          listdf <- split(df, structures[, i])
          fun <- function(x) {
            if (nrow(x) == 1) {
              if (option[1] == "eq") 
                return(1) else return(0)
            } else {
              x <- as.data.frame(t(x))
              Px <- as.data.frame(sweep(x, 2, colSums(x), "/"))
              wx <- colSums(x)/sum(x)
              betax <- t(wx) %*% discintra(Px) %*% wx
              betaEqx <- 1/(1 - betax)
              if (option[1] == "eq") 
                return(betaEqx) else if (option[1] == "normed1") {
                  betaNx <- (1 - 1/betaEqx)/(1 - 1/ncol(x))
                  return(betaNx)
                } else if (option[1] == "normed2") {
                  betaNx <- (betaEqx - 1)/(ncol(x) - 1)
                  return(betaNx)
                }
            }
          }
          resinterm <- unlist(lapply(listdf, fun))
          if (metmean == "arithmetic") 
            return(sum(W[[i]] * resinterm)) else return(1/sum(W[[i]]/resinterm))
        } else {
          dfi <- as.data.frame(sapply(df, function(x) tapply(x, structures[i - 1], sum)))
          fac <- structures[!duplicated(structures[, i - 1]), i]
          names(fac) <- unique(structures[, i - 1])
          fac <- fac[rownames(dfi)]
          listdf <- split(dfi, fac)
          fun <- function(x) {
            if (nrow(x) == 1) {
              if (option[1] == "eq") 
                return(1) else return(0)
            } else {
              x <- as.data.frame(t(x))
              Px <- as.data.frame(sweep(x, 2, colSums(x), "/"))
              wx <- colSums(x)/sum(x)
              betax <- t(wx) %*% discintra(Px) %*% wx
              betaEqx <- 1/(1 - betax)
              if (option[1] == "eq") 
                return(betaEqx) else if (option[1] == "normed1") {
                  betaNx <- (1 - 1/betaEqx)/(1 - 1/ncol(x))
                  return(betaNx)
                } else if (option[1] == "normed2") {
                  betaNx <- (betaEqx - 1)/(ncol(x) - 1)
                  return(betaNx)
                }
            }
          }
          resinterm <- unlist(lapply(listdf, fun))
          if (metmean == "arithmetic") 
            return(sum(W[[i]] * resinterm)) else return(1/sum(W[[i]]/resinterm))
        }
      }
      resbeta <- sapply(nc:1, funbeta)
      if (option[1] == "eq") {
        gammaEq <- beta1Eq * prod(resbeta) * alphaEq
        res <- cbind.data.frame(c(beta1Eq, resbeta, alphaEq, gammaEq))
        inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep = ""))
        intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep = ""))
        intrainter <- paste(inter[-(nc + 1)], intra[-1])
        rownames(res) <- c(inter[nc + 1], rev(intrainter), intra[1], "Gamma")
        colnames(res) <- "Equivalent numbers"
      } else if (option[1] == "normed1") {
        beta1N <- (beta1Eq - 1)/(length(levels(structures[, nc])) - 1)
        res <- cbind.data.frame(c(beta1N, resbeta))
        inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep = ""))
        intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep = ""))
        intrainter <- paste(inter[-(nc + 1)], intra[-1])
        rownames(res) <- c(inter[nc + 1], rev(intrainter))
        colnames(res) <- "Normed contributions to diversity"
      } else if (option[1] == "normed2") {
        beta1N <- (beta1Eq - 1)/(length(levels(structures[, nc])) - 1)
        res <- cbind.data.frame(c(beta1N, resbeta))
        inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep = ""))
        intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep = ""))
        intrainter <- paste(inter[-(nc + 1)], intra[-1])
        rownames(res) <- c(inter[nc + 1], rev(intrainter))
        colnames(res) <- "Normed contributions to diversity"
      }
  }
  return(res)
}
