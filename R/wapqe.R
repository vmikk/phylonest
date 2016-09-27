# Apportionment of quadratic entropy
#' @title Apportionment of quadratic entropy.
#'
#' @param df Dataframe or matrix with sites as rows and species as columns. Entries are abundances of species within sites.
#' @param dis Dissimilarity among species (NULL or class 'dist').
#' @param structures Data frame that contains the name of the group (row) of an level (column) to which the site belongs. Sites in structures should be in the same order as in df. Default is NULL.
#' @param formula Quadratic entropy formula ("QE", "EDI"). "QE" is default.
#' @param wopt Site weighting type ("even", "speciesab"). Default is "even".
#' @param tol A tolerance threshold (a value less than tol is considered equal to zero).
#'
#' @details
#' For the associated permutation test see \code{\link{randtestapqe}}.
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
#' wapqe(ab, , stru, formula = "QE")
#' wapqe(ab, dTaxo, stru, formula = "QE")
#' wapqe(ab, dSize, stru, formula = "QE")
#' wapqe(ab, dFeed, stru, formula = "QE")
#' wapqe(ab, dSF, stru, formula = "QE")
#' 
#' @export

wapqe <- function(df, dis = NULL, structures = NULL, formula = c("QE", "EDI"),
                  wopt = c("even", "speciesab"), tol = 0.00000001) {
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
  a <- apqe(as.data.frame(t(df)), dis = sqrt(2 * d), structures = structures)
  return(a$results)
}