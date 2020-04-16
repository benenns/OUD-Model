#############################################################################
## CA OAT CEA Model
## Deterministic - Parameter seetings from input file
## Feb 16, 2017
#############################################################################
## Legacy: myCAModel-0-TestParameterSetting-v10
## Updated: weighted HIV costs by Pr(on ART); new param 'probOnART'
## Updated: HIV mix corrected (same across heroin/PO)
#############################################################################
# PARAMETERS FROM FILE
#############################################################################
library(XLConnect)
#demoExcelFile <- system.file("C:/Users/Benjamin/Documents/GitHub/OUD-Model/Data/CombinedInputs-V1-Baseline.xlsx", package = "XLConnect")

WB <- loadWorkbook("C:/Users/Benjamin/Documents/GitHub/OUD-Model/Data/CombinedInputs-V1-Baseline.xlsx")

unpack.list <- function(object) {
  for(.x in names(object)){
    assign(value = object[[.x]], x=.x, envir = parent.frame())
  }
} 

unpack.list(as.list(readWorksheet(WB, sheet = "Parameters", rownames=1)))
#############################################################################
# INPUTS FROM FILE
#############################################################################
To.logic <- readWorksheet(WB, sheet = "BooleanStates", rownames=1)
Cycle.logic <- readWorksheet(WB, sheet = "BooleanCycle", rownames=1)
#Weibull <- as.matrix(readWorksheet(WB, sheet = "WeibullFor", header = FALSE) [,-1])
Weibull.Shape <- as.matrix(readWorksheet(WB, sheet = "WeibullShape", header = FALSE) [,-1])
Weibull.Scale <- as.matrix(readWorksheet(WB, sheet = "WeibullScale", header = FALSE) [,-1])
Empirical.Distribution <- readWorksheet(WB, sheet = "EmpiricalStates", rownames=1)
Death.table <- as.matrix(readWorksheet(WB, sheet = "DeathTable", header = TRUE))
Death.mult <- readWorksheet(WB, sheet = "DeathMult", rownames=1)
HIV.sero <- readWorksheet(WB, sheet = "HIVsero", rownames=1)
Cycle.Frailty <- as.matrix(readWorksheet(WB, sheet = "Frailty", rownames=1))
CycleFrailty <- Cycle.Frailty[c("FRAILTY.1", "FRAILTY.2", "FRAILTY.3"),]
Init.dist <- readWorksheet(WB, sheet = "InitDist", rownames=1)
if(switch.HIV==1){
  Init.dist.HIV.H <- cbind((Init.dist[,grep("H", names(Init.dist))] * (1 - (HIVposMix))),
                           (Init.dist[,grep("H", names(Init.dist))] * (HIVposMix)))
  check.probs(Init.dist.HIV.H / sum(Init.dist[,grep("H", names(Init.dist))]))
  
  Init.dist.HIV.P <- cbind((Init.dist[,grep("P", names(Init.dist))] * (1 - (HIVposMix))),
                           (Init.dist[,grep("P", names(Init.dist))] * (HIVposMix)))
  check.probs(Init.dist.HIV.P / sum(Init.dist[,grep("P", names(Init.dist))]))
  
  Init.dist.HIV <- cbind((Init.dist[,grep("H", names(Init.dist))] * (1 - (HIVposMix))),
                         (Init.dist[,grep("P", names(Init.dist))] * (1 - (HIVposMix))),
                         (Init.dist[,grep("H", names(Init.dist))] * (HIVposMix)),
                         (Init.dist[,grep("P", names(Init.dist))] * (HIVposMix)))
  AllHSNames.pos <- paste0(colnames(Init.dist), "+")
  AllHSNames.neg <- paste0(colnames(Init.dist), "-")
  AllHSNames.HIV <- c(AllHSNames.neg, AllHSNames.pos)
  colnames(Init.dist.HIV) <- AllHSNames.HIV
  check.probs(Init.dist.HIV)
  Init.dist <- Init.dist.HIV
}
#############################################################################
# INPUTS FROM FILE
#############################################################################
# Gender weighted death probabilities
if(switch.gender==1){
  Death.mix <- Death.table
  colnames(Death.mix) <- c("Blocks", "GenderMix.H", "GenderMix.P")
  Death.mix[,"GenderMix.H"] <- ((Death.table[, "Male"] * (1 - GenderMix.H)) + (Death.table[,"Female"] * (GenderMix.H)))
  Death.mix[,"GenderMix.P"] <- ((Death.table[, "Male"] * (1 - GenderMix.P)) + (Death.table[,"Female"] * (GenderMix.P)))
} else {Death.mix <- Death.table[,c("Blocks", "Male")]}

#############################################################################
# COSTS & QALYs
#############################################################################
# QALY inputs
StateQALYs <- readWorksheet(WB, sheet = "StateQALYs", rownames=1)
# Cost inputs
StateCosts <- readWorksheet(WB, sheet = "StateCosts", rownames=1)
CrimeCosts <- readWorksheet(WB, sheet = "CrimeCosts")

if(switch.HIV==0){
  # TX costs
  CostVec.tx <- c(as.vector(rep(StateCosts["TX", ], each=NumberOfEpisodes), mode = "numeric"), 0, 0)
  CostVec.tx.h <- c(as.vector(rep(StateCosts["TX", grep(".H", names(StateCosts), fixed = FALSE)], 
                                  each=NumberOfEpisodes), mode = "numeric"), 0)
  CostVec.tx.p <- c(as.vector(rep(StateCosts["TX", grep(".P", names(StateCosts), fixed = FALSE)], 
                                  each=NumberOfEpisodes), mode = "numeric"), 0)
  # HRU costs
  CostVec.hru <- c(as.vector(rep(StateCosts["HRU", ], each=NumberOfEpisodes), mode = "numeric"), 0, 0)
  CostVec.hru.h <- c(as.vector(rep(StateCosts["HRU", grep(".H", names(StateCosts), fixed = FALSE)], 
                                   each=NumberOfEpisodes), mode = "numeric"), 0)
  CostVec.hru.p <- c(as.vector(rep(StateCosts["HRU", grep(".P", names(StateCosts), fixed = FALSE)], 
                                   each=NumberOfEpisodes), mode = "numeric"), 0)
  # Crime costs
  CostMat.crime <- array(0, c(nrow(CrimeCosts), ((ncol(CrimeCosts) - 1) * NumberOfEpisodes) + 2))
  CostMat.crime.h <- array(0, c(nrow(CrimeCosts), ((ncol(CrimeCosts[, grep(".H", names(CrimeCosts), fixed = FALSE)])) * NumberOfEpisodes) + 1))
  CostMat.crime.p <- array(0, c(nrow(CrimeCosts), ((ncol(CrimeCosts[, grep(".P", names(CrimeCosts), fixed = FALSE)])) * NumberOfEpisodes) + 1))
  for(i in 1:nrow(CrimeCosts)){
    CostMat.crime[i,] <- c(as.vector(rep(CrimeCosts[i,-1], each=NumberOfEpisodes), mode = "numeric"), 0, 0)
    CostMat.crime.h[i,] <- c(as.vector(rep(CrimeCosts[i, grep(".H", names(CrimeCosts), fixed = FALSE)], each=NumberOfEpisodes), mode = "numeric"), 0)
    CostMat.crime.p[i,] <- c(as.vector(rep(CrimeCosts[i, grep(".P", names(CrimeCosts), fixed = FALSE)], each=NumberOfEpisodes), mode = "numeric"), 0)
  }
  
  QALYVec <- c(as.vector(rep(StateQALYs["QALYs.NEG", ], each=NumberOfEpisodes), mode = "numeric"), 0, 0)/12
  QALYVec.h <- c(as.vector(rep(StateQALYs["QALYs.NEG", grep(".H", names(StateQALYs), fixed = FALSE)], 
                               each=NumberOfEpisodes), mode = "numeric"), 0)/12
  QALYVec.p <- c(as.vector(rep(StateQALYs["QALYs.NEG", grep(".P", names(StateQALYs), fixed = FALSE)], 
                               each=NumberOfEpisodes), mode = "numeric"), 0)/12
  
  
} else if(switch.HIV==1){
  # TX costs
  CostVec.tx <- c(as.vector(rep(StateCosts["TX",], each=NumberOfEpisodes, times = 2), mode = "numeric"), 0, 0)
  CostVec.tx.h <- c(as.vector(rep(StateCosts["TX", grep(".H", names(StateCosts), fixed = FALSE)], 
                                  each=NumberOfEpisodes, times = 2), mode = "numeric"), 0)
  CostVec.tx.p <- c(as.vector(rep(StateCosts["TX", grep(".P", names(StateCosts), fixed = FALSE)], 
                                  each=NumberOfEpisodes, times = 2), mode = "numeric"), 0)
  # HRU costs
  CostVec.hru <- c(as.vector(rep(StateCosts["HRU",], each=NumberOfEpisodes, times = 2), mode = "numeric"), 0, 0)
  CostVec.hru.h <- c(as.vector(rep(StateCosts["HRU", grep(".H", names(StateCosts), fixed = FALSE)], 
                                   each=NumberOfEpisodes, times = 2), mode = "numeric"), 0)
  CostVec.hru.p <- c(as.vector(rep(StateCosts["HRU", grep(".P", names(StateCosts), fixed = FALSE)], 
                                   each=NumberOfEpisodes, times = 2), mode = "numeric"), 0)
  # HIV-specific costs
  CostVec.hivhru <- c(as.vector(rep(0, ncol(StateCosts) * NumberOfEpisodes)), 
                   as.vector(rep(StateCosts["HIV",], each=NumberOfEpisodes), mode = "numeric"), 0, 0)
  CostVec.hivhru.h <- c(as.vector(rep(0, ncol(StateCosts["HIV", grep(".H", names(StateCosts), fixed = FALSE)]) * NumberOfEpisodes)), 
                     as.vector(rep(StateCosts["HIV", grep(".H", names(StateCosts), fixed = FALSE)], 
                                   each=NumberOfEpisodes), mode = "numeric"), 0)
  CostVec.hivhru.p <- c(as.vector(rep(0, ncol(StateCosts["HIV", grep(".P", names(StateCosts), fixed = FALSE)]) * NumberOfEpisodes)), 
                     as.vector(rep(StateCosts["HIV", grep(".P", names(StateCosts), fixed = FALSE)], 
                                   each=NumberOfEpisodes), mode = "numeric"), 0)
  # ART
  CostVec.art <- c(as.vector(rep(0, ncol(StateCosts) * NumberOfEpisodes)), 
                   as.vector(rep(StateCosts["ART",] * probOnART, each=NumberOfEpisodes), mode = "numeric"), 0, 0)
  CostVec.art.h <- c(as.vector(rep(0, ncol(StateCosts["ART", grep(".H", names(StateCosts), fixed = FALSE)]) * NumberOfEpisodes)), 
                     as.vector(rep(StateCosts["ART", grep(".H", names(StateCosts), fixed = FALSE)]  * probOnART, 
                                   each=NumberOfEpisodes), mode = "numeric"), 0)
  CostVec.art.p <- c(as.vector(rep(0, ncol(StateCosts["ART", grep(".P", names(StateCosts), fixed = FALSE)]) * NumberOfEpisodes)), 
                     as.vector(rep(StateCosts["ART", grep(".P", names(StateCosts), fixed = FALSE)]  * probOnART, 
                                   each=NumberOfEpisodes), mode = "numeric"), 0)
  # Combine HIV costs
  CostVec.hiv <- CostVec.hivhru + CostVec.art
  CostVec.hiv.h <- CostVec.hivhru.h + CostVec.art.h
  CostVec.hiv.p <- CostVec.hivhru.p + CostVec.art.p
  
  # Crime costs
  CostMat.crime <- array(0, c(nrow(CrimeCosts), ((ncol(CrimeCosts) - 1) * NumberOfEpisodes) * 2 + 2))
  CostMat.crime.h <- array(0, c(nrow(CrimeCosts), ((ncol(CrimeCosts[, grep(".H", names(CrimeCosts), fixed = FALSE)])) * NumberOfEpisodes) * 2 + 1))
  CostMat.crime.p <- array(0, c(nrow(CrimeCosts), ((ncol(CrimeCosts[, grep(".P", names(CrimeCosts), fixed = FALSE)])) * NumberOfEpisodes) * 2 + 1))
  for(i in 1:nrow(CrimeCosts)){
    CostMat.crime[i,] <- c(as.vector(rep(CrimeCosts[i,-1], each=NumberOfEpisodes, times = 2), mode = "numeric"), 0, 0)
    CostMat.crime.h[i,] <- c(as.vector(rep(CrimeCosts[i, grep(".H", names(CrimeCosts), fixed = FALSE)], each=NumberOfEpisodes, times = 2), mode = "numeric"), 0)
    CostMat.crime.p[i,] <- c(as.vector(rep(CrimeCosts[i, grep(".P", names(CrimeCosts), fixed = FALSE)], each=NumberOfEpisodes, times = 2), mode = "numeric"), 0)
  }
  
  QALYVec <- c(as.vector(rep(StateQALYs["QALYs.NEG", ], each=NumberOfEpisodes), mode = "numeric"), 
               as.vector(rep(StateQALYs["QALYs.POS", ], each=NumberOfEpisodes), mode = "numeric"), 0, 0)/12
  QALYVec.h <- c(as.vector(rep(StateQALYs["QALYs.NEG", grep(".H", names(StateQALYs), fixed = FALSE)], each=NumberOfEpisodes), mode = "numeric"), 
                 as.vector(rep(StateQALYs["QALYs.POS", grep(".H", names(StateQALYs), fixed = FALSE)], each=NumberOfEpisodes), mode = "numeric"), 0)/12
  QALYVec.p <- c(as.vector(rep(StateQALYs["QALYs.NEG", grep(".P", names(StateQALYs), fixed = FALSE)], each=NumberOfEpisodes), mode = "numeric"), 
                 as.vector(rep(StateQALYs["QALYs.POS", grep(".P", names(StateQALYs), fixed = FALSE)], each=NumberOfEpisodes), mode = "numeric"), 0)/12
}



# From & To health states
From.baseline <- From.state.names(To.logic,Cycle.logic,NumberOfEpisodes)
To.baseline <- To.state.names(To.logic,Cycle.logic,NumberOfEpisodes)
