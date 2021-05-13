This folder contains data and R script needed to reproduce the resutls in the paper Davies et al (2021) ' Contemporary evolution of the innate immune receptor gene TLR3 in an isolated vertebrate population'

All variable names in these data should be self-explanatory, but please don't hesitate to send an email to the corresponding author if anything is unclear

Files included:
# File_for_genepop_analysis
This file contains the data used for the genepop analysis - the first tab explains the data contained and the analysis conducted, subsequent tabs are in genepop ready format


# cleaned_TLR3_script
This R script contains all the code needed to reproduce the results from the paper (apart from the genepop analysis - see above) - using the data from the below four files


# Genotype_frequencies
This file contains information on number sampled, genotype, and allele frequency changes on Cousin between 1992-2018, for adults population and individuals hatched that year 

Columns correspond to:
year = year of datapoint (1992-2018)
Adult_population = Frequency of adult popualtion sampled (0-1)
Genotype = TLR3 genptype (AA, AC, CC)
Adult_frequency = % of TLR3 genotype for sampled adult population (0-100)
YoB_frequency = % of TLR3 genotype from birds hatched that year (0-100)
Allele = TLR3 allele (A or C)
Adult_allele_frequency = Frequency of TLR3 allele for sampled adult population (0-1)
YoB_allele_frequency = Frequency of TLR3 alleles  from birds hatched that year (0-1)


# Genotype_frequencies_islands
This file contains information on allele frequency accross the five Islands, at two different time points - the first data point for Aride and Cousine is from Cousin

Columns correspond to:
Island = Island and timepoint sampled 
Year = year of datapoint (1993-2018)
A_allele_frequency = Frequency of TLR3 A allele for sampled population (0-1)
C_allele_frequency = Frequency of TLR3 C allele for sampled population  (0-1)


# TLR3_CLEANED
This file contains all the information needed for the survival and reproductive success analysis

Columns correspond to:
BirdID = Unique identifier for each bird
SexEstimate = Sex of bird (0 = Female, 1 = Male)
HatchYear = Year of hatching (1997-2010)
CohortYear = Cohorts (97-99 or 05-10)
SeasonBorn = Season bird hatched from (SB = Major, or WB = Minor)
BirthIsland = Island bird orginiated from - should all be CN
IndependentAtFirstCapture =  Age class of bird at first capture (CH+FL (chick+fledgling), or OLDER)
Translocated = Was bird translocated to another Island at any point (Yes or No)
MinorBreedingSeason = Was their a census in the minor breeding season closest to the last seen date (MINOR) or not (NOMINOR)
AgeAtDeathCousinBiAnnualSeason = Bi-annual age of bird during the last season observed in (years = 0 - 18)
SurviveToFirstYear_BiAnnual = Did bird survive to at least one year of age (Y (yes) or N (No))
Survival_Status = Was bird still alive in during the last census (0 (Yes) or 1 (No))
Survival_Status_WS2018 = Was bird still alive in the Minor 2018 breeding season (0 (Yes) or 1 (No))
CountOfSurvivingOffspringTo1Year = Number of offspring produced over birds lifetime which surived to at least one year of age (0 - 10)
CountOfSurvivingOffspringToOFL = Number of offspring produced over birds lifetime which surived to at least three months of age (0 - 12)
TLR3 = TLR3 genotype (A:A, A:C, C:C)
A_allele_presence = Does the bird have at least one copy of the TLR3 A allele (Yes or No)?
C_allele_presence = Does the bird have at least one copy of the TLR3 C allele (Yes or No)?
Hs_obs = Individual heterozygosity (0.6-1.6)
Maternal_Hs_obs = Maternal heterozygosity (0.4 - 1.5)
MHCDiversity = Number of different MHC class I alleles present (2-8)
Aseua4 = Does bird have a copy of the MHC class I Ase-ua4 allele (0 (No) or 1 (Yes))


# HWE_graph_data
This file contains information needed to plot figure S1 - data generated from genepop analysis

Columns correspond to:
Genotype: TLR3 genotype(AA, AC, CC)
Obs_Exps_values = Is the value for Observed or expected numbers of individulas (Obs or Exp)
Obs_Exp = Is the value for observed or expected numbers of individulas (a (Obs) or b (Exp))
Category = Is analysis for all individuals (All individuals), or only those individuals which survive to at least one year of age (Survive to adulthood)
value = Number of individulas (41-273)- note for expected values number can include decimal places 
