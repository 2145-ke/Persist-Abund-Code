
## Load Libraries
library(labdsv)
library(vegan)
library(data.table)
library(lme4)
library(lmerTest)
library(ggplot2)
library(cowplot)

## Read in data
setwd("~/Dropbox/NutNet data")
original_coverdat <- read.csv("full-cover-09-April-2018.csv")
coverdat <- original_coverdat # backup
site_covars_orig <- read.csv("comb-by-plot-clim-soil-diversity-09-Apr-2018.csv")

dim(coverdat)

# Using April 9 2018 cover data downloaded from Dropbox: 186,781  rows * 18 col = 3,221,370 cells (!)

## Choose sites  ------------------------------------------------------------
setwd("~/Dropbox/NutNet Ashley")
# Tally number of years of data for each site
siteyear_tally <- data.table(coverdat)
siteyear_tally <- siteyear_tally[,list( years = length(unique(year_trt)),
                                        first_year = min(year_trt),
                                        last_year = max(year_trt)),
                                 by=list(site_code) ]

siteyear_tally$timespan <- (siteyear_tally$last_year - siteyear_tally$first_year)+1

# cumulative number of sites by number of years treatment: 
cuml_nsites <- siteyear_tally[,list(n_sites = length(unique(site_code)),
                                    sites = toString(site_code)), 
                              by = list(years, first_year, last_year,timespan)]
cuml_nsites <- cuml_nsites[order(-years),]
cuml_nsites$cuml_nsites <- cumsum(cuml_nsites$n_sites)
cuml_nsites[1:10,c("years", "first_year", "last_year", "n_sites", "cuml_nsites", "sites")]

# There are 28 sites with 8 or more years of data. All 28 sites have year 0 data. Select these sites for the analysis, and subset the cover data to these sites.
chosen_sites <- siteyear_tally$site_code[siteyear_tally$first_year == 0 &
                                           siteyear_tally$years >= 8] # select sites with 7 years of post-treatment data
chosen_sites <- droplevels(chosen_sites)
chosen_sites

# Site experiment deviations (tallied by Jon Bakker, kept for reference)
# keep this site: chosen_sites <- chosen_sites[chosen_sites$site_code != "bnch.us",] #missing two plots in Yr2
# keep this site: chosen_sites <- chosen_sites[chosen_sites$site_code != "hopl.us",] #plots 6, 16, 26 missing Yr0 data
# keep this site: chosen_sites <- chosen_sites[chosen_sites$site_code != "mcla.us",] #missing one plot in Yr2
chosen_sites <- chosen_sites[chosen_sites != "sgs.us"]  #multiple control in Yr0 but other trts missing that yr
# unique(as.character(coverdat$trt[coverdat$site_code == 'sgs.us' & coverdat$year_trt == 0]))

working_coverdat <- coverdat[coverdat$site_code %in% chosen_sites,]
working_coverdat <- droplevels(working_coverdat)
working_coverdat <- data.table(working_coverdat)

## Choose treatment(s) and years to focus on  ------------------------------------------------------------
#treatment codes
trts <- data.frame(trt = c("Control", "N", "P", "K", "PK", "NK", "NP", "NPK", "Fence", "NPK+Fence"),
                   N = c("No", "Yes", "No", "No", "No", "Yes", "Yes", "Yes", "No", "Yes"),
                   P = c("No", "No", "Yes", "No", "Yes", "No", "Yes", "Yes", "No", "Yes"),
                   K = c("No", "No", "No", "Yes", "Yes", "Yes", "No", "Yes", "No", "Yes"),
                   Fence = c("No", "No", "No", "No", "No", "No", "No", "No", "Yes", "Yes"),
                   Num.Fert = c(0, 1, 1, 1, 2, 2, 2, 3, 0, 3))
working_coverdat <- merge(x = working_coverdat, y = trts, all.x = TRUE)

# choose only control, npk, fence, and npk+fence: 
working_coverdat <- working_coverdat[working_coverdat$trt %in% c("Control", "NPK", "Fence", "NPK+Fence"), ]
# choose only the first 8 years of data (7 years of treatment):
working_coverdat <- working_coverdat[working_coverdat$year_trt <8, ] #5923 x 23
working_coverdat <- droplevels(working_coverdat)


## Taxonomic QC by site from Jon Bakker  ------------------------------------------------------------
# This resolves "problematic" species from each site. First, get rid of obvious "non-species."
working_coverdat$Taxon<-as.factor(toupper(as.character(working_coverdat$Taxon))) # convert to upper-case
nlevels(working_coverdat$Taxon)
# 1004 taxa to begin with. Get rid of non-living, non-plant "taxa":
 
working_coverdat <- working_coverdat[working_coverdat$live == 1 & working_coverdat$Family != "NULL",]
working_coverdat <- droplevels(working_coverdat)
nlevels(working_coverdat$Taxon)

# 983 taxa remaining (got rid of 19 non-live and NULL family). Get rid of non-vascular mosses:
working_coverdat <- working_coverdat[! working_coverdat$Family %in% c("Polytrichaceae", "Brachytheciaceae", "Hylocomiaceae", "Amblystegiaceae"), ]
working_coverdat <- droplevels(working_coverdat)
nlevels(working_coverdat$Taxon)
# now that we are done with tallying up the number of speices, can set this to be a character vector (probably more appropriate for analysis and QC).
working_coverdat$Taxon<-as.character(working_coverdat$Taxon) 
# 980 taxa remaining (got rid of 3 taxa of mosses).
# Now, go site-by-site to correct taxonomic errors and problematic species. 

# this chunk spits out each site's coverdat into a unique .csv file: 
# for(i in 1:length(chosen_sites)) {
#   write.csv(spplist(datafile = working_coverdat, site_code = chosen_sites[i]), file = paste(as.character(chosen_sites[i]), ".csv", sep = ""))
# }

# 1 bogong.au
Site = "bogong.au"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "ERIGERON ", x = temp_coverdat$Taxon)] <- "ERIGERON SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 2 bnch.us - need to check

# 3 burrawan.au - species list confirmed with J. Firn on 170324.

# 4 cdcr.us
Site = "cdcr.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "CHENOPODIUM ", x = temp_coverdat$Taxon)] <- "CHENOPODIUM SP."
temp_coverdat$Taxon[grep(pattern = "ASCLEPIAS ", x = temp_coverdat$Taxon)] <- "ASCLEPIAS SP."
temp_coverdat$Taxon[grep(pattern = "ERIGERON CANADENSIS", x = temp_coverdat$Taxon)] <- "CONYZA CANADENSIS"
temp_coverdat$Taxon[grep(pattern = "ERIGERON ", x = temp_coverdat$Taxon)] <- "ERIGERON SP."
temp_coverdat$Taxon[grep(pattern = "TRAGOPOGON ", x = temp_coverdat$Taxon)] <- "TRAGOPOGON DUBIUS"
temp_coverdat$Taxon[grep(pattern = "RUDBECKIA ", x = temp_coverdat$Taxon)] <- "RUDBECKIA HIRTA"
temp_coverdat$Taxon[grep(pattern = "CAREX ", x = temp_coverdat$Taxon)] <- "CAREX SP."
temp_coverdat$Taxon[grep(pattern = "CYPERUS ", x = temp_coverdat$Taxon)] <- "CYPERUS SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 5 cdpt.us - no changes

# 6 cbgb.us
Site = "cbgb.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "CIRSIUM ", x = temp_coverdat$Taxon)] <- "CIRSIUM SP."
temp_coverdat$Taxon[grep(pattern = "HELIANTHUS ", x = temp_coverdat$Taxon)] <- "HELIANTHUS SP."
temp_coverdat$Taxon[grep(pattern = "SOLANUM ", x = temp_coverdat$Taxon)] <- "SOLANUM CAROLINENSE"
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 7 cowi.ca - no changes

# 8 elliot.us
Site = "elliot.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "JUNCUS ", x = temp_coverdat$Taxon)] <- "JUNCUS SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 9 frue.ch - no changes

# 10 hall.us
Site = "hall.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "RUBUS ", x = temp_coverdat$Taxon)] <- "RUBUS SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 11 hopl.us - need to add

# 12 kiny.au
Site = "kiny.au"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "ARTHROPODIUM ", x = temp_coverdat$Taxon)] <- "ARTHROPODIUM SP."
temp_coverdat$Taxon[grep(pattern = "ATRIPLEX ", x = temp_coverdat$Taxon)] <- "ATRIPLEX SP."
temp_coverdat$Taxon[grep(pattern = "RHODANTHE ", x = temp_coverdat$Taxon)] <- "RHODANTHE SP."
temp_coverdat$Taxon[grep(pattern = "SONCHUS ", x = temp_coverdat$Taxon)] <- "SONCHUS SP."
temp_coverdat$Taxon[grep(pattern = "CRASSULA ", x = temp_coverdat$Taxon)] <- "CRASSULA SP."
temp_coverdat$Taxon[grep(pattern = "TRIFOLIUM ", x = temp_coverdat$Taxon)] <- "TRIFOLIUM SP."
temp_coverdat$Taxon[grep(pattern = "UNKNOWN ", x = temp_coverdat$Taxon)] <- "UNKNOWN "
temp_coverdat$Taxon[grep(pattern = "AVENA ", x = temp_coverdat$Taxon)] <- "AVENA SP."
temp_coverdat$Taxon[grep(pattern = "BROMUS ", x = temp_coverdat$Taxon)] <- "BROMUS SP."
temp_coverdat$Taxon[grep(pattern = "LOLIUM ", x = temp_coverdat$Taxon)] <- "LOLIUM SP."
temp_coverdat$Taxon[grep(pattern = "RYTIDOSPERMA ", x = temp_coverdat$Taxon)] <- "RYTIDOSPERMA SP."
temp_coverdat$Taxon[grep(pattern = "VULPIA ", x = temp_coverdat$Taxon)] <- "VULPIA SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# koffler.ca - email with M. Cadotte in March 2017
Site = "koffler.ca"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "CAREX ", x = temp_coverdat$Taxon)] <- "CAREX SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 13 konz.us
Site = "konz.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "SOLIDAGO ", x = temp_coverdat$Taxon)] <- "SOLIDAGO SP."
temp_coverdat$Taxon[grep(pattern = "MUHLENBERGIA ", x = temp_coverdat$Taxon)] <- "MUHLENBERGIA CUSPIDATA"
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 14 lancaster.uk - no changes

# 15 look.us
Site = "look.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "GALIUM ", x = temp_coverdat$Taxon)] <- "GALIUM SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 16 mcla.us -  need to check

# 17 mtca.au
Site = "mtca.au"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "STIPA ", x = temp_coverdat$Taxon)] <- "STIPA SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 18 sedg.us
Site = "sedg.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "TRIFOLIUM ", x = temp_coverdat$Taxon)] <- "TRIFOLIUM SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 19 sevi.us - need to check

# 20 shps.us
Site = "shps.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "ANTENNARIA ", x = temp_coverdat$Taxon)] <- "ANTENNARIA SP."
temp_coverdat$Taxon[temp_coverdat$Taxon == "PSEUDOSCLEROCHLOA RUPESTRIS"] <- "POA SECUNDA"
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 21 sgs.us - need to check

# 22 sier.us
Site = "sier.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "TORILIS ", x = temp_coverdat$Taxon)] <- "TORILIS SP."
temp_coverdat$Taxon[grep(pattern = "TRITELEIA ", x = temp_coverdat$Taxon)] <- "TRITELEIA SP."
temp_coverdat$Taxon[grep(pattern = "UNKNOWN ASTERACEAE", x = temp_coverdat$Taxon)] <- "UNKNOWN ASTERACEAE SP."
temp_coverdat$Taxon[grep(pattern = "PLAGIOBOTHRYS ", x = temp_coverdat$Taxon)] <- "PLAGIOBOTHRYS SP."
temp_coverdat$Taxon[grep(pattern = "CARDAMINE ", x = temp_coverdat$Taxon)] <- "CARDAMINE SP."
temp_coverdat$Taxon[grep(pattern = "CONVOLVULUS ", x = temp_coverdat$Taxon)] <- "CONVOLVULUS SP."
temp_coverdat$Taxon[grep(pattern = "TRIFOLIUM ", x = temp_coverdat$Taxon)] <- "TRIFOLIUM SP."
temp_coverdat$Taxon[grep(pattern = "VICIA ", x = temp_coverdat$Taxon)] <- "VICIA SP."
temp_coverdat$Taxon[grep(pattern = "ERODIUM ", x = temp_coverdat$Taxon)] <- "ERODIUM SP."
temp_coverdat$Taxon[grep(pattern = "CLARKIA ", x = temp_coverdat$Taxon)] <- "CLARKIA SP."
temp_coverdat$Taxon[grep(pattern = "LINANTHUS ", x = temp_coverdat$Taxon)] <- "LINANTHUS SP."
temp_coverdat$Taxon[grep(pattern = "GALIUM ", x = temp_coverdat$Taxon)] <- "GALIUM SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 23 smith.us
Site = "smith.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "SYMPHORICARPOS ", x = temp_coverdat$Taxon)] <- "SYMPHORICARPOS ALBUS"
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 24 spin.us - no changes

# 25 temple.us
Site = "temple.us"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "POLYTAENIA ", x = temp_coverdat$Taxon)] <- "POLYTAENIA NUTTALLII" # conversation with P. Fay
temp_coverdat$Taxon[grep(pattern = "ASCLEPIAS ", x = temp_coverdat$Taxon)] <- "ASCLEPIAS SP."
temp_coverdat$Taxon[grep(pattern = "GAILLARDIA ", x = temp_coverdat$Taxon)] <- "GAILLARDIA SP."
temp_coverdat$Taxon[grep(pattern = "LACTUCA ", x = temp_coverdat$Taxon)] <- "LACTUCA SP."
temp_coverdat$Taxon[grep(pattern = "LIATRIS ", x = temp_coverdat$Taxon)] <- "LIATRIS SP."
temp_coverdat$Taxon[grep(pattern = "SPOROBOLUS ", x = temp_coverdat$Taxon)] <- "SPOROBOLUS SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# ukul.za
Site = "ukul.za"
temp_coverdat <- working_coverdat[working_coverdat$site_code == Site,]
temp_coverdat$Taxon[grep(pattern = "PACHYCARPUS ", x = temp_coverdat$Taxon)] <- "PACHYCARPUS SP."
temp_coverdat$Taxon[grep(pattern = "ACACIA ", x = temp_coverdat$Taxon)] <- "ACACIA SP."
temp_coverdat$Taxon[grep(pattern = "CHAMAECRISTA ", x = temp_coverdat$Taxon)] <- "CHAMAECRISTA SP."
temp_coverdat$Taxon[grep(pattern = "ALBUCA", x = temp_coverdat$Taxon)] <- "ALBUCA SETOSA"
temp_coverdat$Taxon[grep(pattern = "HYPERICUM ", x = temp_coverdat$Taxon)] <- "HYPERICUM SP."
temp_coverdat$Taxon[grep(pattern = "SOLANUM ", x = temp_coverdat$Taxon)] <- "SOLANUM SP."
working_coverdat <- rbind(working_coverdat[working_coverdat$site_code != Site,], temp_coverdat)

# 26 valm.ch - no changes

length(unique(working_coverdat$Taxon))

# 948 taxa remaining: consolidated or dropped 27 taxa with QC.

## Site-specific adjustments for unique block & plot configurations: 
# cbgb.us
# six blocks. Verified with L. Biederman on 170324.
# Either delete blocks 4-6 ...or treat as two 'sites'.
levels(working_coverdat$site_code) <- c(levels(working_coverdat$site_code), "cbgb.us.A", "cbgb.us.B")
working_coverdat$site_code[(working_coverdat$site_code == "cbgb.us" & working_coverdat$block %in% c(1, 2, 3))] <- "cbgb.us.A"
working_coverdat$site_code[working_coverdat$site_code == "cbgb.us" & working_coverdat$block %in% c(4, 5, 6)] <- "cbgb.us.B"
levels(site_covars_orig$site_code) <- c(levels(site_covars_orig$site_code), "cbgb.us.A", "cbgb.us.B")
site_covars_orig$site_code[(site_covars_orig$site_code == "cbgb.us" & site_covars_orig$block %in% c(1, 2, 3))] <- "cbgb.us.A"
site_covars_orig$site_code[(site_covars_orig$site_code == "cbgb.us" & site_covars_orig$block %in% c(4, 5, 6))] <- "cbgb.us.B"
chosen_sites <- c(levels(chosen_sites), "cbgb.us.A", "cbgb.us.B")
chosen_sites <- chosen_sites[!chosen_sites == "cbgb.us"]

#cdcr.us
# This is fine for this analysis - keep extra blocks.

# five blocks (delete 4-5)
# working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "cdcr.us" & block %in% c(4, 5)), ])

# #cdpt.us: six blocks. J. Knops indicated 1,3,6 are most different, 2-5 are similar.
# #Decided to treat as two 'sites' of 3 blocks each.
# levels(working_coverdat$site_code) <- c(levels(working_coverdat$site_code),"cdpt.us.varied", "cdpt.us.similar")
# working_coverdat$site_code[working_coverdat$site_code == "cdpt.us" & working_coverdat$block %in% c(1, 3, 6)] <- "cdpt.us.varied"
# working_coverdat$site_code[working_coverdat$site_code == "cdpt.us" & working_coverdat$block %in% c(2, 4, 5)] <- "cdpt.us.similar"
# levels(site_covars_orig$site_code) <- c(levels(site_covars_orig$site_code), "cdpt.us.varied", "cdpt.us.similar")
# site_covars_orig$site_code[(site_covars_orig$site_code == "cdpt.us" & site_covars_orig$block %in% c(1, 3, 6))] <- "cdpt.us.varied"
# site_covars_orig$site_code[(site_covars_orig$site_code == "cdpt.us" & site_covars_orig$block %in% c(2, 4, 5))] <- "cdpt.us.similar"
# levels(chosen_sites$site_code) <- c(levels(chosen_sites$site_code), "cdpt.us.varied", "cdpt.us.similar")
# chosen_sites <- rbind(chosen_sites, chosen_sites[chosen_sites$site_code == "cdpt.us", ], chosen_sites[chosen_sites$site_code == "cdpt.us", ])
# chosen_sites$site_code[chosen_sites$site_code == "cdpt.us"] <- c("cdpt.us", "cdpt.us.varied", "cdpt.us.similar")

#koffler.ca: 3 control plots in each block; different subplot for plot 19 in 2014
# focus on first control plot in each block. Verified with M. Cadotte on 170324.
working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "koffler.ca" & plot %in% c(9, 11, 17, 21, 34, 36)), ])

#marc.ar: 3 control plots in blocks 1 and 2. Verified with J. Alberti on 170324.
# working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "marc.ar" & plot %in% c(6, 8, 11, 17)), ])

#mtca.au: four blocks (delete 4). Verified with S. Prober on 170324.
# working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "mtca.au" & block %in% c(4)), ])

#sedg.us: 2 control, 2 NPK (no fences) in each block
#  working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "sedg.us" & plot %in% c(7, 10, 17, 18, 27, 28)), ])

#shps.us: four blocks (delete 4?)
# working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "shps.us" & block %in% c(4)), ])

#sier.us: five blocks (delete 4-5?); no Yr0 data for blocks 4-5
working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "sier.us" & block %in% c(4, 5)), ])

#summ.za: 3 control (no fences) in each block; one control measured only in Yr0 and another not measured in Yr2
# working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "summ.za" & plot %in% c(1, 10, 15, 16, 21, 30)), ])

#temple.us: drop extra plots in block 2
# working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "temple.us" & plot %in% c(19, 20)), ])

#ukul.za: 2 control, 2 NPK (no fences) in each block
# working_coverdat <- with(working_coverdat, working_coverdat[!(site_code == "ukul.za" & plot %in% c(8, 10, 19, 20, 25, 30)), ])



## Final adjustments to cover data
working_coverdat$site_code <- factor(working_coverdat$site_code)
working_coverdat$block <- factor(working_coverdat$block)
working_coverdat$plot <- factor(working_coverdat$plot)
working_coverdat$unique_plot_ID <- paste(sub("[.].*", "",  as.character(working_coverdat$site_code)), # just site name, no country
                                         working_coverdat$block, working_coverdat$plot, sep = "_")
working_coverdat$unique_plot_year_ID <- paste(working_coverdat$unique_plot_ID, 
                                              working_coverdat$year_trt, sep = "_")

#Commented code here creates a summary table for each site based just on plots included in this analysis
# sites.working_coverdat <- unique(working_coverdat$site_code)
# for(i in 1:length(sites.working_coverdat)) {
#  write.csv(spplist(datafile = working_coverdat, site_code = sites.working_coverdat[i]), file = paste(as.character(sites.working_coverdat[i]), ".csv", sep = ""))
# }

# sum cover to account for taxonomic adjustments: 
working_coverdat<-working_coverdat[,list(max_cover=sum(max_cover)),
                                   by=list(unique_plot_year_ID, site_code, block, plot, year_trt, trt, Taxon, functional_group, local_lifespan, local_provenance, N_fixer, ps_path) ]

write.csv(working_coverdat, "~/Dropbox/NutNet Ashley/full-cover-subsetAA-25April2018.csv")


## Organize site-level covariates
site_covars <- site_covars_orig[site_covars_orig$site_code %in% chosen_sites, ]
site_covars<-data.table(site_covars)
site_covars<-site_covars[,list(managed = unique(managed), burned = unique(burned), grazed = unique(grazed),
                               anthropogenic = unique(anthropogenic), habitat = unique(habitat),
                               elevation = unique(elevation), live_mass = mean(live_mass, na.rm = TRUE), pct_N = mean(pct_N, na.rm = TRUE),
                               ppm_P = mean(ppm_P, na.rm = TRUE), ppm_K = mean(ppm_K, na.rm = TRUE),
                               PercentSand = mean(PercentSand, na.rm = TRUE)),
                         by=list(site_code) ]  

bioclim_vars <- site_covars_orig[site_covars_orig$site_code %in% chosen_sites, ]
bioclim_vars<-data.table(bioclim_vars)
bioclim_vars <- bioclim_vars[,list(
  MAT = unique(MAT_v2), MAT_RANGE = unique(MAT_RANGE_v2), ISO = unique(ISO_v2),
  TEMP_VAR = unique(TEMP_VAR_v2), MAX_TEMP = unique(MAX_TEMP_v2), MIN_TEMP = unique(MIN_TEMP_v2),
  ANN_TEMP_RANGE = unique(ANN_TEMP_RANGE_v2), TEMP_WET_Q = unique(TEMP_WET_Q_v2),
  TEMP_DRY_Q = unique(TEMP_DRY_Q_v2), TEMP_WARM_Q = unique(TEMP_WARM_Q_v2),
  TEMP_COLD_Q = unique(TEMP_COLD_Q_v2), MAP = unique(MAP_v2), MAP_WET_M = unique(MAP_WET_M_v2),
  MAP_DRY_M = unique(MAP_DRY_M_v2), MAP_VAR = unique(MAP_VAR_v2), MAP_WET_Q = unique(MAP_WET_Q_v2),
  MAP_DRY_Q = unique(MAP_DRY_Q_v2), MAP_WARM_Q = unique(MAP_WARM_Q_v2), MAP_COLD_Q = unique(MAP_COLD_Q_v2),
  RAIN_PET = unique(RAIN_PET), AI = unique(AI), PET = unique(PET), N_Dep = unique(N_Dep)),
  by = list(site_code)]

site_covars$management <- with(site_covars, managed + burned + grazed + anthropogenic) #combine multiple types of management
site_covars$management[site_covars$management > 1] <- 1
site_covars$management <- as.factor(site_covars$management)
site_covars$log.elevation <- log(site_covars$elevation)

#calculate site-level species richness
temp_coverdat<-data.table(working_coverdat)
temp_coverdat<-temp_coverdat[,list(max_cover = sum(max_cover)),
                             by = list(site_code, Taxon)]
S <- temp_coverdat[, list (S = length(Taxon)),
                   by = list(site_code)]
site_covars <- merge(x = site_covars, y = S, by = "site_code")
site_covars$log.S <- log(site_covars$S)

write.csv(bioclim_vars, "~/Dropbox/NutNet Ashley/comb-bioclim-covars-subsetAA.csv")
write.csv(site_covars, "~/Dropbox/NutNet Ashley/comb-site-covars-subsetAA.csv")
