# First type of diversity decomposition
#' @title First decomposition of diversity introduced in Pavoine et al. 2016.
#' @description Use when the interest is in multiple-site and multiple-region beta diversity, the sampling design is even (same number of sites within regions), sites within a region can be given equal weights and regions can be given equal weights.
#'
#' @param df Dataframe or matrix with sites as rows and species as columns. Entries are abundances of species within sites.
#' @param dis Dissimilarity among species (NULL or class 'dist').
#' @param structures Data frame that contains the name of the group (row) of an level (column) to which the site belongs. Sites in structures should be in the same order as in df. Default is NULL.
#' @param option Rescaling type ("eq", "normed1" or "normed2").
#' @param formula Quadratic entropy formula ("QE", "EDI"). "QE" is default.
#' @param tol A tolerance threshold (a value less than tol is considered equal to zero).
#'
#' @details
#' Rescaling types (argument 'option'):
#' \itemize{
#' \item \strong{"eq"} - the diversity components are given in terms of equivalent number of species, sites, regions etc.;
#' \item \strong{"normed1"} - the normed components of diversity will be returned with formula (1 – 1 / E) / (1 - 1 / Emax) [de Bello et al., 2010, Eqn 14];
#' \item \strong{"normed2"} - the normed components of diversity will be returned with formula (E – 1) / (Emax - 1) [Villéger et al., 2012].
#' }
#' For Eβ, Emax = M (the number of sites). For Eα and Eγ, Emax=S (the number of species in the data set).
#'
#' Formula type (argument 'formula'):
#' \itemize{
#' \item \strong{"QE"} - the definition of the quadratic entropy is following Rao (1982);
#' \item \strong{"EDI"} - the Euclidean Diversity Index of Champely and Chessel (2002).
#' }
#'
#' For the associated permutation test see \code{\link{randtestEqRS}}.
#' 
#' @return A data frame with each component of the selected diversity decomposition.
#' 
#' Rescaled estimates of the beta components (option = "normed1" or "normed2") reach the maximum value of 1 when sites within a region are maximally dissimilar whatever the level of diversity within sites (e.g., within any region, sites do not share species and any species from any site always is maximally dissimilar from all species in all other sites); beta components should be equal 0 when sites are identical within regions.
#' @author Sandrine Pavoine, Eric Marcon, Carlo Ricotta.
#' @references
#' Pavoine S., Marcon E., Ricotta C. (2016) ‘Equivalent numbers’ for species, phylogenetic or functional diversity in a nested hierarchy of multiple scales. Methods Ecol Evol, 7: 1152-1163. DOI:10.1111/2041-210X.12591
#' 
#' de Bello F., Lavergne S., Meynard C.N., Lepš J., Thuiller W. (2010) The partitioning of diversity: showing Theseus a way out of the labyrinth. Journal of Vegetation Science, 21: 992-1000. DOI: 10.1111/j.1654-1103.2010.01195.x
#' 
#' Villéger S., Miranda J.R., Hernandez D.F., Mouillot, D. (2012) Low functional beta-diversity despite high taxonomic beta-diversity among tropical estuarine fish communities. PLoS One, 7, e40679. DOI: 10.1371/journal.pone.0040679
#' 
#' Champely S., Chessel D. (2002) Measuring biological diversity using Euclidean metrics. Environmental and Ecological Statistics, 9, 167-177. DOI: 10.1023/A:1015170104476
#' 
#' Rao C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology, 21, 24-43. DOI: 10.1016/0040-5809(82)90004-1
#' 
#' @seealso \code{\link{randtestEqRS}}, \code{\link{EqRSintra}}, \code{\link{EqRao}}.
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
#' EqRS(ab, , stru, option="eq")
#' EqRS(ab, dTaxo, stru, formula = "QE", option="eq")
#' EqRS(ab, dSize, stru, formula = "QE", option="eq")
#' EqRS(ab, dFeed, stru, formula = "QE", option="eq")
#' EqRS(ab, dSF, stru, formula = "QE", option="eq")
#' 
#' EqRS(ab, , stru, option="normed2")
#' EqRS(ab, dTaxo, stru, formula = "QE", option="normed2")
#' EqRS(ab, dSize, stru, formula = "QE", option="normed2")
#' EqRS(ab, dFeed, stru, formula = "QE", option="normed2")
#' EqRS(ab, dSF, stru, formula = "QE", option="normed2")
#' 
#' @export

EqRS <- function(df, dis = NULL, structures = NULL, option = c("eq", "normed1", "normed2"), 
                 formula = c("QE", "EDI"), tol = 0.00000001) {
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
        forw <- cbind.data.frame(as.vector(firstw), as.vector(listw), as.vector(rep(finalw, nrow(structures))))
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
  nlev <- nrow(a$results)
  gamma <- a$results[nlev, 1]
  gamma <- 1/(1 - gamma)
  alpha <- a$results[nlev - 1, 1]
  alpha <- 1/(1 - alpha)
  beta1 <- (1 - sum(a$results[-c(1, nlev), 1]))/(1 - a$results[nlev, 1])
  if (nlev == 3) {
    if (option[1] == "eq") {
      res <- cbind.data.frame(c(beta1, alpha, gamma))
      rownames(res) <- c("Inter-sites", "Intra-sites", "Gamma")
      colnames(res) <- "Equivalent numbers"
    } else if (option[1] == "normed2") {
      beta1 <- (beta1 - 1)/(ncomm - 1)
      res <- cbind.data.frame(c(beta1))
      rownames(res) <- c("Inter-sites")
      colnames(res) <- "Normed beta diversity"
    } else if (option[1] == "normed1") {
      beta1 <- (1 - 1/beta1)/(1 - 1/ncomm)
      res <- cbind.data.frame(c(beta1))
      rownames(res) <- c("Inter-sites")
      colnames(res) <- "Normed beta diversity"
    }
    return(res)
  } else {
    v1 <- sapply(3:(nlev - 1), function(i) sum(a$results[i:(nlev - 1), 1]))
    v2 <- sapply(2:(nlev - 2), function(i) sum(a$results[i:(nlev - 1), 1]))
    resbeta <- (1 - v1)/(1 - v2)
    if (option[1] == "eq") {
      nc <- ncol(structures)
      res <- cbind.data.frame(c(beta1, resbeta, alpha, gamma))
      inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep = ""))
      intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep = ""))
      intrainter <- paste(inter[-(nc + 1)], intra[-1])
      rownames(res) <- c(inter[nc + 1], rev(intrainter), intra[1], "Gamma")
      colnames(res) <- "Equivalent numbers"
    } else if (option[1] == "normed2") {
      nc <- ncol(structures)
      beta1N <- (beta1 - 1)/(length(levels(structures[, nc])) - 1)
      fun <- function(i) {
        T <- as.factor(structures[, i])
        poidsnum <- tapply(w, structures[, i], sum)
        if (i == 1) {
          effectifs <- as.vector(table(T))
          poidsden <- poidsnum/effectifs
        } else {
          effectifs <- unlist(lapply(split(structures[, i - 1], T), function(x) length(unique(x))))
          poidsden <- poidsnum/effectifs
        }
        dfilist <- split(df, T)
        if (i > 1) {
          strulist <- split(structures, T)
          strulist <- lapply(strulist, function(x) x[, 1:(i - 1), drop = FALSE])
        }
        fun0 <- function(j) {
          tab <- dfilist[[j]]
          if (i > 1) {
            thestru <- strulist[[j]]
            for (k in 1:ncol(thestru)) thestru[, k] <- factor(thestru[, k])
          }
          options(warn = -1)
          if (i > 1) {
            ai <- apqe(as.data.frame(t(tab)), sqrt(2 * as.dist(d)), thestru)
          } else { 
            ai <- apqe(as.data.frame(t(tab)), sqrt(2 * as.dist(d)))
          }
          options(warn = op)
          nlevai <- nrow(ai$results)
          alphai <- ai$results[nlevai, 1] - ai$results[1, 1]
          return(alphai)
        }
        resi <- unlist(lapply(1:length(dfilist), fun0))
        res <- sum((1 - resi) * poidsnum)/sum((1 - resi) * poidsden)
        return(res)
      }
      M <- sapply(nc:1, fun)
      resbetaN <- (resbeta - 1)/(M - 1)
      res <- cbind.data.frame(c(beta1N, resbetaN))
      inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep = ""))
      intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep = ""))
      intrainter <- paste(inter[-(nc + 1)], intra[-1])
      rownames(res) <- c(inter[nc + 1], rev(intrainter))
      colnames(res) <- "Normed contributions to beta diversity"
    } else if (option[1] == "normed1") {
      nc <- ncol(structures)
      beta1N <- (beta1 - 1)/(length(levels(structures[, nc])) - 1)
      fun <- function(i) {
        T <- as.factor(structures[, i])
        poidsnum <- tapply(w, structures[, i], sum)
        if (i == 1) {
          effectifs <- as.vector(table(T))
          poidsden <- poidsnum/effectifs
        } else {
          effectifs <- unlist(lapply(split(structures[, i - 1], T), function(x) length(unique(x))))
          poidsden <- poidsnum/effectifs
        }
        dfilist <- split(df, T)
        if (i > 1) {
          strulist <- split(structures, T)
          strulist <- lapply(strulist, function(x) x[, 1:(i - 1), drop = FALSE])
        }
        fun0 <- function(j) {
          tab <- dfilist[[j]]
          if (i > 1) {
            thestru <- strulist[[j]]
            for (k in 1:ncol(thestru)) thestru[, k] <- factor(thestru[, k])
          }
          options(warn = -1)
          if (i > 1) {
            ai <- apqe(as.data.frame(t(tab)), sqrt(2 * as.dist(d)), thestru)
          } else {
            ai <- apqe(as.data.frame(t(tab)), sqrt(2 * as.dist(d)))
          }
          options(warn = op)
          nlevai <- nrow(ai$results)
          alphai <- ai$results[nlevai, 1] - ai$results[1, 1]
          return(alphai)
        }
        resi <- unlist(lapply(1:length(dfilist), fun0))
        res <- sum((1 - resi) * poidsnum)/sum((1 - resi) * poidsden)
        return(res)
      }
      M <- sapply(nc:1, fun)
      resbetaN <- (1 - 1/resbeta)/(1 - 1/M)
      res <- cbind.data.frame(c(beta1N, resbetaN))
      inter <- c("Inter-sites", paste("Inter-", colnames(structures), sep = ""))
      intra <- c("Intra-sites", paste("Intra-", colnames(structures), sep = ""))
      intrainter <- paste(inter[-(nc + 1)], intra[-1])
      rownames(res) <- c(inter[nc + 1], rev(intrainter))
      colnames(res) <- "Normed contributions to beta diversity"
    }
    return(res)
  }
}