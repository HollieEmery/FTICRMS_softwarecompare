library(tidyverse)

# create formatted data for each application ####

#function to read in .pks files
pks_to_ICBM <- function(x){
  dat <- read_table(x, skip=7,
                    col_names = c("m/z","I","rel_abund","Res.","freq","S/N"),guess_max = 26000)
  return(dat[,c("m/z","I","S/N","Res.")])
}

#read in and concatenate test data
path = "raw_2017/"
fns = dir(path, pattern = ".pks")
dat <- NA
for(i in fns){
  jnk <- pks_to_ICBM(paste0(path,i))
  jnk$index <- gsub(".pks","",i)
  dat <- rbind(dat,jnk)
}
dat <- dat[-1,] %>%
  mutate(peak_id=paste(1:nrow(.)))

dat$index[grep("SRNOM_PPL",dat$index)] <- "SRNOM"
dat$index[grep("NPEqIW",dat$index)] <- "NPEqIW"
dat$index[grep("Seawater",dat$index)] <- "BottomWater"
dat$index[grep("Shallow",dat$index)] <- "U1383C_Shallow"
dat$index[grep("Mid",dat$index)] <- "U1383C_Middle"
dat$index[grep("82a",dat$index)] <- "U1382A"
dat$index[grep("Deep",dat$index)] <- "U1383C_Deep"

dat$index <- factor(dat$index,levels = c("SRNOM","NPEqIW","BottomWater","U1382A","U1383C_Shallow","U1383C_Middle","U1383C_Deep"))

#ICBM-OCEAN format
## https://rhea.icbm.uni-oldenburg.de/geomol/
## following recommendations by Merder et al 2020 and at https://uol.de/f/5/inst/icbm/ag/mgc/icbm-ocean/ICBM-OCEAN_Userguides_v1.0.pdf
### data input format: mz	| I	| MDL_x	| Res.	| index
# or
###                    m/z	| I	| S/N	| Res.
smpls = as.character(na.omit(unique(dat$index)))
for(i in smpls){
  jnk <- filter(dat,index==i) %>% 
    select(-peak_id,-index) %>%
    write_csv(paste0("for_ICBM-OCEAN/",i,".csv"))
}

## UME ####
## http://dockersrv1.awi.de:3838/ume/
## following defaults and quick start video at https://www.awi.de/en/science/biosciences/ecological-chemistry/tools/ume.html
### in video he removes surfactants, blank, and C13 rejects
### preset filter profile "standard marine DOM": DBE max 20, N0-2, S0-1, P0, remove surfactants, remove non-verified by 13C isotope
### defaults: .5 ppm error, 199-799 m/z, database "02 NOM: N6, P3, S3"
#### database "02 NOM: N6, P3, S3" contains: 12C>1 & 13C<=1, H>1 & H<2+2C+N+P, O<=26, 14N<=6 & 15N<=1, P<=3, 32S <=3 & 34S<=1, DBE>=0 & integer, OC<=1.2, N/C<1.3, P/C<.3, S/C<.8
#### can change NPS range in post
#### autocalculates hclust, NMDS, many indices, can download both wide and long format, and can download sample-wise summaries

### data input format: file_id | peak_id | mz | i_magnitude | s_n
dat %>%
  select(file_id=index, peak_id, mz=`m/z`, i_magnitude=I, s_n = `S/N`) %>%
  write_csv("forUME.csv")

## Formularity ####
## https://pnnl-comp-mass-spec.github.io/Formularity/
## following recommendations by Leefmann et al 2018
## file format is same as MFAssignR
for(i in smpls){
  jnk <- filter(dat,index==i) %>% 
    select(1:3) %>%
    write_csv(paste0("for_Formularity/",i,".csv"))
}

## MFAssignR ####
## https://github.com/skschum/MFAssignR (forked at https://github.com/HollieEmery/MFAssignR)
## following recommendations by Schum al 2020 and in the readme on github (see above)
## data format: m/z | I | ret. time (or whatever you want)
for(i in smpls){
  jnk <- filter(dat,index==i) %>% 
    select(1:3) %>%
    write_csv(paste0("for_MFAssignR/",i,".csv"))
}

## PetroOrg ####

# import results and standardize formatting ####
# icb = ICBM-OCEAN. ume = UltraMass Explorer. pet = PetroOrg. for = Formularity.
path="output/"
dat_icb <- read_csv(paste0(path,"ICBMOCEAN/Formula_attribution_likeliestxmatchx_x_anonymous_9536108.csv"),guess_max = 10000)
dat_ume <- read_delim(paste0(path,"2022-10-03_19-45-15_UTC-UME_export-dataset_forUME/normalized_filtered_data.csv"),delim=";")
dat_pet <- read_csv(paste0(path,"PetroOrg/PetroOrg.csv"))
dat_for_s <- read_csv(paste0(path,"Formularity/s_ipdbNPEqIW.csv"))
dat_for_p <- read_csv(paste0(path,"Formularity/p_ipdbNPEqIW.csv"))

#how many peaks
dat_icb %>% filter(!is.na(Sample7_NPEqIW.csv)) %>% nrow()
dat_ume %>% filter(file_id=="NPEqIW") %>% nrow()
dat_pet %>% filter(!is.na(NPEqIW_standard)) %>% nrow()
length(unique(dat_for_p$mz))

#weighted mean mass
dat_icb %>% filter(!is.na(Sample7_NPEqIW.csv)) %>% summarize(mu=weighted.mean(mz,Sample7_NPEqIW.csv))
dat_pet %>% filter(!is.na(NPEqIW_standard)) %>% summarize(mu=weighted.mean(Theor_mz,NPEqIW_standard))
dat_ume %>% filter(file_id=="NPEqIW") %>% summarize(mu=weighted.mean(mz,i_magnitude))
weighted.mean(dat_for_p$mz[!duplicated(dat_for_p$mz)],dat_for_p$intensity[!duplicated(dat_for_p$mz)])
