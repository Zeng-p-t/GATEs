# GATEs
Matching-Based Nonparametric Estimation of Group Average Treatment Effects

- 0-bandwidth.R: the bandwidth h is chosen by using the dpill() function in the R package KernSmooth.
- 1-MainProposedSup.R: simulation code of matching-based nonparametric estimator(MATCH).
- 2-MainProposedBCSup.R: simulation code of bias-corrected matching-based nonparametric estimator(MATCH.bc).
- 3-MainCompare.R: simulation code for comparing the proposed methods with several potentially competing methods, including the IPW, OR, AIPW and PSR methods.
- 4-MainCompareTrimming.R: simulation code for comparing the proposed methods with several potentially competing methods after trimming extreme propensity scores / weights.
- EstPSR.R: function for estimating GATE based on the PSR method.
- EstTau.R: function for estimating GATE.
- GenData.R: data-generating mechanisms for cases (C1)-(C12).
- MatchY1Y0.R: function for matching with "manhattan" distance of the proposed methods.


