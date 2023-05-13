#' Partitioning biodiversity effects into density and size components
#'
#' Function \code{densize} (DENSIty and SIZE) partitions the net, complementarity, and selection effects of biodiversity into additive components that reflect diversity-induced changes in plant density and size.
#'
#' @param init.dens A matrix or data frame consisting of initial plant density (the number of plants or seeds sown per area).
#' Rows are plots and columns are species.
#' Elements should be filled with 0 for plots without focal species.
#' The arguments \code{final.dens}, \code{final.yield}, and \code{germ.dens} follow the same format, with the plots and species arranged in the same order as in \code{init.dens}.
#' @param final.dens A matrix or data frame consisting of final plant density (the number of plants survived).
#' @param final.yield A matrix or data frame consisting of final yield (yield per area).
#' Yield can be in any unit; e.g., biomass (g), volume (m3).
#' @param germ.dens Optional. A matrix or data frame consisting of germination density (the number of seedlings emerged per area).
#' When provided, germination components are calculated in addition to density and size components.
#' @author Shinichi Tatsumi
#' @return The \code{densize} function returns a matrix with rows as plots and columns as types of components. Rows representing monoculture plots, for which biodiversity cannot be defined, are removed.
#' \itemize{
#'  \item{When the argument \code{germ.dens} was NOT provided, the function returns (1) density-mediated complementarity effects (\code{dens.compl}), (2) density-mediated selection effects (\code{dens.selec}), (3) size-mediated complementarity effects (\code{size.compl})}, and (4) size-mediated selection effects (\code{size.selec}).
#'  \item{When the argument \code{germ.dens} was provided, the function returns (1) germination-mediated complementarity effects (\code{germ.compl}), (2) germination-mediated selection effects (\code{germ.selec}), (3) post-germination density-mediated complementarity effects (\code{pg.dens.compl}), (4) post-germination density-mediated selection effects (\code{pg.dens.selec}), (5) size-mediated complementarity effects (\code{size.compl})}, and (6) size-mediated selection effects (\code{size.selec}).
#' }
#' @examples
#' dat  <- BioDivExpt.1
#' res1 <- densize(dat$InitDens, dat$FinalDens, dat$FinalYield)
#' res2 <- densize(dat$InitDens, dat$FinalDens, dat$FinalYield, dat$GermDens)
#' # res1: Density and size components
#' # res2: Germination, post-germination density, and size components
#' #  Note that the input dataset (BioDivExpt_1) has 45 rows, while
#' # each output object has 30 rows. This is because the densize
#' # function has removed rows representing monocultures (1 to 15),
#' # for which biodiversity effects cannot be defined.
#' @export

densize <- function(init.dens, final.dens, final.yield, germ.dens) {

  Nsp          <- ncol(init.dens)
  Nplot        <- nrow(init.dens)

  WhichMono    <- unlist(ifelse(rowSums(init.dens!=0)==1, apply(init.dens!=0, 1, which), rep(0, Nplot)))
  M            <- matrix(rep(sapply(1:Nsp, function(s) mean(final.yield[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
  Mmean        <- mean(M[1,])

  Dprime       <- matrix(rep(sapply(1:Nsp, function(s) mean(final.dens[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
  Dobserved    <- final.dens
  Dexpected    <- sweep(init.dens, 1, rowSums(init.dens), "/") * Dprime
  DeltaD       <- Dobserved - Dexpected

  Wprime       <- matrix(rep(sapply(1:Nsp, function(s) sum(final.yield[WhichMono==s, s]) / sum(final.dens[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
  Wobserved    <- ifelse(final.dens!=0, as.matrix(final.yield / final.dens), 0)
  Wexpected    <- Wprime
  DeltaW       <- Wobserved - Wexpected

  if(missing(germ.dens)) {

    dens.compl    <- rowSums(       Mmean/2/M  * (Wobserved + Wexpected) * DeltaD)
    dens.selec    <- rowSums((1/2 - Mmean/2/M) * (Wobserved + Wexpected) * DeltaD)
    size.compl    <- rowSums(       Mmean/2/M  * (Dobserved + Dexpected) * DeltaW)
    size.selec    <- rowSums((1/2 - Mmean/2/M) * (Dobserved + Dexpected) * DeltaW)

    Res    <- cbind(dens.compl, dens.selec, size.compl, size.selec)

  } else {

    Iprime       <- matrix(rep(sapply(1:Nsp, function(s) mean(init.dens[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
    Iobserved    <- init.dens
    Gprime       <- matrix(rep(sapply(1:Nsp, function(s) mean(germ.dens[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
    Gobserved    <- germ.dens

    gObserved    <- ifelse(init.dens!=0, as.matrix(germ.dens / init.dens), 0)
    gExpected    <- Gprime / Iprime
    Deltag       <- gObserved - gExpected

    pObserved    <- ifelse(germ.dens!=0, as.matrix(final.dens / germ.dens), 0)
    pExpected    <- Dprime / Gprime
    Deltap       <- pObserved - pExpected

    germ.compl       <- rowSums(       Mmean/2/M  * (Wobserved + Wexpected) * Iobserved * Dprime / Gprime * Deltag)
    germ.selec       <- rowSums((1/2 - Mmean/2/M) * (Wobserved + Wexpected) * Iobserved * Dprime / Gprime * Deltag)
    pg.dens.compl    <- rowSums(       Mmean/2/M  * (Wobserved + Wexpected) * Gobserved                   * Deltap)
    pg.dens.selec    <- rowSums((1/2 - Mmean/2/M) * (Wobserved + Wexpected) * Gobserved                   * Deltap)
    size.compl       <- rowSums(       Mmean/2/M  * (Dobserved + Dexpected)                               * DeltaW)
    size.selec       <- rowSums((1/2 - Mmean/2/M) * (Dobserved + Dexpected)                               * DeltaW)

    Res    <- cbind(germ.compl,    germ.selec,
                    pg.dens.compl, pg.dens.selec,
                    size.compl,    size.selec)
  }

  rownames(Res) <- 1:Nplot
  Res           <- Res[WhichMono==0,]

  return(Res)
}
