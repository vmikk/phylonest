# Second decomposition of diversity
#' @title Second decomposition of diversity introduced in Pavoine et al. 2016.
#' @description Use when the interest is in multiple-site (within regions) and multiple-region beta diversity, sites within a region can be given equal weights and regions can be given equal weights but the sampling design is uneven (different numbers of sites within regions).
#'
#' @param df Dataframe or matrix with sites as rows and species as columns. Entries are abundances of species within sites.
#' @param dis Dissimilarity among species (NULL or class 'dist').
#' @param structures Data frame that contains the name of the group (row) of an level (column) to which the site belongs. Sites in structures should be in the same order as in df. Default is NULL.
#' @param option Rescaling type ("eq", "normed1" or "normed2").
#' @param formula Quadratic entropy formula ("QE", "EDI"). "QE" is default.
#' @param tol A tolerance threshold (a value less than tol is considered equal to zero).
#' @param metmean Mean type - "arithmetic" or "harmonic" (default).
#'
#' @details
#' For the associated permutation test see \code{\link{randtestEqRSintra}}.
#'
#' @return A data frame with each component of the selected diversity decomposition.
#' @author Sandrine Pavoine, Eric Marcon, Carlo Ricotta.
#' @references Pavoine, S., Marcon, E. and Ricotta, C. (2016), ‘Equivalent numbers’ for species, phylogenetic or functional diversity in a nested hierarchy of multiple scales. Methods Ecol Evol. DOI:10.1111/2041-210X.12591
#' @seealso \code{\link{randtestEqRSintra}}, \code{\link{EqRS}}, \code{\link{EqRao}}.
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
#' EqRSintra(ab, , stru, option="eq")
#' EqRSintra(ab, dTaxo, stru, formula = "QE", option="eq")
#' EqRSintra(ab, dSize, stru, formula = "QE", option="eq")
#' EqRSintra(ab, dFeed, stru, formula = "QE", option="eq")
#' EqRSintra(ab, dSF, stru, formula = "QE", option="eq")
#' 
#' EqRSintra(ab, , stru, option="normed2")
#' EqRSintra(ab, dTaxo, stru, formula = "QE", option="normed2")
#' EqRSintra(ab, dSize, stru, formula = "QE", option="normed2")
#' EqRSintra(ab, dFeed, stru, formula = "QE", option="normed2")
#' EqRSintra(ab, dSF, stru, formula = "QE", option="normed2")
#' 
#' @export

EqRSintra <- function(df, dis = NULL, structures = NULL, option = c("eq", "normed1", "normed2"), 
                      formula = c("QE", "EDI"), tol = 0.00000001, metmean = c("harmonic", "arithmetic")) {
  metmean <- metmean[1]
  if (!option[1] %in% c("eq", "normed1", "normed2")) 
    stop("unavailable option, please modify; option can be eq, normed1, or normed2")
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
        forw <- cbind.data.frame(firstw, listw, rep(finalw, nrow(structures)))[,-1]  # error was here
        w <- apply(forw, 1, prod)
      }
    }
  df <- P * w
  if (!is.null(structures)) {
    if (length(levels(factor(structures[, 1]))) == 1) 
      stop("All sites belong to a unique level in the first column of structures, remove this first column in structures")
    if (length(levels(factor(structures[, 1]))) == nrow(df)) 
      stop("Each site belongs to a distinct level in the first column of structures, this first column is useless, remove it and re-run")
  }

  op <- options()$warn
  options(warn = -1)
  a <- apqe(as.data.frame(t(df)), sqrt(2 * as.dist(d)), structures)
  options(warn = op)
  dfprop <- sweep(df, 1, rowSums(df), "/")
  divsites <- sapply(as.data.frame(t(dfprop)), .diversity, d=d)
  if (metmean == "arithmetic") {
    alphaEQ <- sum(w * (1/(1 - divsites)))
  } else {
    alphaEQ <- 1/(1 - sum(w * divsites))
  }
  nlev <- nrow(a$results)
  beta1 <- (1 - sum(a$results[-c(1, nlev), 1]))/(1 - a$results[nlev, 1])
  if (nlev == 3) {
    gamma <- beta1 * alphaEQ
    if (option[1] == "eq") {
      res <- cbind.data.frame(c(beta1, alphaEQ, gamma))
      rownames(res) <- c("Inter-sites", "Intra-sites", "Gamma")
      colnames(res) <- "Equivalent numbers"
    } else if (option[1] == "normed1") {
      beta1 <- (1 - 1/beta1)/(1 - 1/ncomm)
      res <- cbind.data.frame(beta1)
      rownames(res) <- c("Inter-sites")
      colnames(res) <- "Normed inter-site diversity"
    } else if (option[1] == "normed2") {
      beta1 <- (beta1 - 1)/(ncomm - 1)
      res <- cbind.data.frame(beta1)
      rownames(res) <- c("Inter-sites")
      colnames(res) <- "Normed inter-site diversity"
    }
    return(res)
  } else {
    nc <- ncol(structures)
    funnc <- function(i) {
      dfilist <- split(df, as.factor(structures[, i]))
      if (i > 1) {
        strulist <- split(structures, as.factor(structures[, i]))
        strulist <- lapply(strulist, function(x) x[, -(i:nc), drop = FALSE])
      }
      wi <- tapply(w, as.factor(structures[, i]), sum)
      fun <- function(j) {
        tab <- dfilist[[j]]
        if (i > 1) 
          thestru <- strulist[[j]] else thestru <- NULL
          op <- options()$warn
          options(warn = -1)
          ai <- apqe(as.data.frame(t(tab)), sqrt(2 * as.dist(d)), thestru)
          nlev <- nrow(ai$results)
          numi <- ai$results[nlev, 1] - ai$results[1, 1]
          gammai <- ai$results[nlev, 1]
          options(warn = op)
          ratioi <- (1 - numi)/(1 - gammai)
          if (option[1] == "eq") {
            return(ratioi)
          } else if (option[1] == "normed1") {
            if (!is.null(thestru)) 
              len <- length(levels(factor(thestru[, ncol(thestru)]))) else len <- nrow(tab)
              if ((ratioi - 1) < 0.0000000001) 
                res <- 0 else res <- (1 - 1/ratioi)/(1 - 1/len)
                return(res)
          } else if (option[1] == "normed2") {
            if (!is.null(thestru)) 
              len <- length(levels(factor(thestru[, ncol(thestru)]))) else len <- nrow(tab)
              if ((ratioi - 1) < 0.0000000001) 
                res <- 0 else res <- (ratioi - 1)/(len - 1)
                return(res)
          }
      }
      resi <- unlist(lapply(1:length(dfilist), fun))
      if (metmean == "arithmetic") 
        res <- sum(wi * resi) else res <- 1/sum(wi * (1/resi))
      return(res)
    }
    resbeta <- sapply(nc:1, funnc)
    if (option[1] == "eq") {
      gamma <- beta1 * prod(resbeta) * alphaEQ
      res <- cbind.data.frame(c(beta1, resbeta, alphaEQ, gamma))
      inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep = ""))
      intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep = ""))
      intrainter <- paste(inter[-(nc + 1)], intra[-1])
      rownames(res) <- c(inter[nc + 1], rev(intrainter), intra[1], "Gamma")
      colnames(res) <- "Equivalent numbers"
    } else if (option[1] == "normed1") {
      beta1N <- (1 - 1/beta1)/(1 - 1/length(levels(structures[, nc])))
      res <- cbind.data.frame(c(beta1N, resbeta))
      inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep = ""))
      intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep = ""))
      intrainter <- paste(inter[-(nc + 1)], intra[-1])
      rownames(res) <- c(inter[nc + 1], rev(intrainter))
      colnames(res) <- "Normed contributions to beta diversity"
    } else if (option[1] == "normed2") {
      beta1N <- (beta1 - 1)/(length(levels(structures[, nc])) - 1)
      res <- cbind.data.frame(c(beta1N, resbeta))
      inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep = ""))
      intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep = ""))
      intrainter <- paste(inter[-(nc + 1)], intra[-1])
      rownames(res) <- c(inter[nc + 1], rev(intrainter))
      colnames(res) <- "Normed contributions to beta diversity"
    }
    return(res)
  }
}
