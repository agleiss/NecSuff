# Degrees of necessity and of sufficiency

We suggest measures to quantify the degrees of necessity and of sufficiency (DN and DS) of prognostic factors for dichotomous, survival, ordinal and nominal outcomes. A cause, represented by certain values of prognostic factors, is considered necessary for an event if, without the cause, the event cannot develop. It is considered sufficient for an event if the event is unavoidable in the presence of the cause. Necessity and sufficiency can be seen as the two faces of causation, and this symmetry and equal relevance are reflected by the suggested measures. The measures provide an approximate, in some cases an exact, multiplicative decomposition of explained variation as defined by Schemper for dichotomous outcomes. The SAS macro 'NecSuff' and the R function 'NecSuff' provided below implement the marginal DN and DS measures for the case of a dichotomous outcome. The SAS macro 'NecSuff_Surv' implements marginal and partial DN and DS for a survival outcome. The SAS macros 'NecSuff_Ord' and 'NecSuff_Nom0' calculate marginal DN and DS for an ordinal and a nominal outcome (with reference category), respectively.
Measures for competing risks survival data are provided by the SAS macro 'NecSuff_CR'

 

## References:

Gleiss, A, Schemper, M (2019):
"Quantifying degrees of necessity and of sufficiency in cause-effect relationships with dichotomous and survival outcomes.", Statistics in Medicine. 2019;38:4733–4748.
doi:10.1002/sim.8331

Gleiss, A, Henderson, R, Schemper, M (2021):
"Degrees of necessity and of sufficiency: further results and extensions, with an application to covid-19 mortality in Austria", Statistics in Medicine. 2021;40:3352-3366.
doi:10.1002/sim.8961 

Gleiss, A, Gnant, M, Schemper, M:
"Explained variation and degrees of necessity and of sufficiency for competing risks survival data" (submitted)