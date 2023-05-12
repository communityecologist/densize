#' Partitioning biodiversity effects into density and size components
#'
#' Function \code{densize} (DENSIty and SIZE) partitions the net, complementarity, and selection effects of biodiversity into additive components that reflect diversity-induced changes in plant density and size.
#'
#' @param init_dens Matrix or data frame consisting of initial plant (seed) density. Rows are plots and columns are species. Elements should be filled with 0 for plots without focal species. The same applies to \code{final_dens}, \code{final_yield}, and \code{germ_dens}.
#' @param final_dens Matrix or data frame consisting of final plant density. Note that the plots and species must be in the exact same order as \code{init_dens}. The same applies to \code{final_yield} and \code{germ_dens}.
#' @param final_yield Matrix or data frame consisting of final yield.
#' @param germ_dens Optional. Matrix or data frame consisting of germination density (e.g., seedling counts). When provided, germination components are calculated in addition to density and size components.
#' @author Shinichi Tatsumi
#' @examples
#' aaa
#' @export

densize <- function(init_dens, final_dens, final_yield, germ_dens) {

  Nsp          <- ncol(init_dens)
  Nplot        <- nrow(init_dens)

  WhichMono    <- unlist(ifelse(rowSums(init_dens!=0)==1, apply(init_dens!=0, 1, which), rep(0, Nplot)))
  M            <- matrix(rep(sapply(1:Nsp, function(s) mean(final_yield[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
  Mmean        <- mean(M[1,])

  Dprime       <- matrix(rep(sapply(1:Nsp, function(s) mean(final_dens[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
  Dobserved    <- final_dens
  Dexpected    <- sweep(init_dens, 1, rowSums(init_dens), "/") * Dprime
  DeltaD       <- Dobserved - Dexpected

  Wprime       <- matrix(rep(sapply(1:Nsp, function(s) sum(final_yield[WhichMono==s, s]) / sum(final_dens[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
  Wobserved    <- ifelse(final_dens!=0, as.matrix(final_yield / final_dens), 0)
  Wexpected    <- Wprime
  DeltaW       <- Wobserved - Wexpected

  if(missing(germ_dens)) {

    dens.compl    <- rowSums(       Mmean/2/M  * (Wobserved + Wexpected) * DeltaD)
    dens.selec    <- rowSums((1/2 - Mmean/2/M) * (Wobserved + Wexpected) * DeltaD)
    size.compl    <- rowSums(       Mmean/2/M  * (Dobserved + Dexpected) * DeltaW)
    size.selec    <- rowSums((1/2 - Mmean/2/M) * (Dobserved + Dexpected) * DeltaW)

    Res    <- cbind(dens.compl, dens.selec, size.compl, size.selec)

  } else {

    Iprime       <- matrix(rep(sapply(1:Nsp, function(s) mean(init_dens[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
    Iobserved    <- init_dens
    Gprime       <- matrix(rep(sapply(1:Nsp, function(s) mean(germ_dens[WhichMono==s, s])), Nplot), nrow=Nplot, byrow=T)
    Gobserved    <- germ_dens

    gObserved    <- ifelse(init_dens!=0, as.matrix(germ_dens / init_dens), 0)
    gExpected    <- Gprime / Iprime
    Deltag       <- gObserved - gExpected

    pObserved    <- ifelse(germ_dens!=0, as.matrix(final_dens / germ_dens), 0)
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
