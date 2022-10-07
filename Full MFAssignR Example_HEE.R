# most of this code is from the example included with MFAssignR

##########################
#Package Install Procedure
#install.packages("devtools")

# download the package from the GitHub repo: https://github.com/skschum/MFAssignR
# then set the wd to where you put it and install so it functions as a package:
#setwd("MFAssignR/")
#devtools::install("MFAssignR")
####################
library(MFAssignR)
path="for_MFAssignR/"
fns = dir(path=path)

Data <- read.csv(paste0(path,fns[2])) # just to grab NPEqIW

#########################
#Signal-to-noise estimation and check
Noise <- KMDNoise(Data) #Using new KMDNoise() noise estimation function
plot <- Noise[["KMD"]]  #Extracting the plot from the KMDNoise function
plot                    #Printing the plot
KMDN <-Noise[["Noise"]] #Extracting the estimated noise from the KMDNoise function
KMDN                    #Printing the noise
SNplot(Data, cut = KMDN * 6, mass = 401.1, window.x = 100, window.y = 100) #Reasonable settings for SNplot

#####################
#Isotope prescreening
Isotope <- IsoFiltR(Data)  #Input for IsoFiltR

Mono <- Isotope[["Mono"]]
Iso <- Isotope[["Iso"]]
##############################
#CHO formula pre-assign

Assign <- MFAssignCHO(Mono, Iso, ionMode = "neg", lowMW =50, highMW = 1000, ppm_err = .3, H_Cmin = 0.3,
                      HetCut = "off", NMScut = "on", SN = 6*KMDN)   #Standard parameters for negative mode

Unambig1 <- Assign[["Unambig"]]
Ambig1 <- Assign[["Ambig"]]
Unassigned1 <- Assign[["None"]]

MSAssign <- Assign[["MSAssign"]]
Error <- Assign[["Error"]]
MSgroups <- Assign[["MSgroups"]]
VK <- Assign[["VK"]]
MSAssign
Error
MSgroups
VK
##################################
#Highlighting possible recalibrant series
check <- RecalList(Unambig1)

##################################
#Qualitative check of recalibrant series and mass recalibration.

Test <- Recal(df = Unambig1,peaks = Mono, isopeaks = Iso, mode = "neg", SN = 6*KMDN, mzRange = 50, series1 = "O8_H_9",
              series2 = "O6_H_3", series3 = "O4_H_2", series4 = "O13_H_13", series5 = "O15_H_16")

Plot <- Test[["Plot"]]
Plot      #This plot is slow to generate
Mono2 <- Test[["Mono"]]
Iso2 <- Test[["Iso"]]
List <- Test[["RecalList"]]


#############################################
#Final formula assignment
#Parameters for both positive and negative mode formula assignment. This uses the formula extension version

Assign <- MFAssign(Mono2, Iso2, ionMode = "neg", lowMW =50, highMW = 1000,  Nx = 3, Sx = 1,  ppm_err = .3, H_Cmin = 0.3,
                   HetCut = "off", DeNovo = 400, NMScut = "on", SN = 6*KMDN)


Unambig2 <- Assign[["Unambig"]]
Ambig2 <- Assign[["Ambig"]]
Unassigned2 <- Assign[["None"]]

MSAssign <- Assign[["MSAssign"]]
Error <- Assign[["Error"]]
MSgroups <- Assign[["MSgroups"]]
VK <- Assign[["VK"]]
MSAssign
Error
MSgroups
VK


#############################################
# no recalibration or filtering
Assign <- MFAssign(Mono, Iso, ionMode = "neg", lowMW =50, highMW = 1000,  Nx = 3, Sx = 1,  ppm_err = .3, H_Cmin = 0.3,
                   HetCut = "off", DeNovo = 400, NMScut = "on")

Unambig0 <- Assign[["Unambig"]]
Ambig0 <- Assign[["Ambig"]]
Unassigned0 <- Assign[["None"]]

MSAssign <- Assign[["MSAssign"]]
Error <- Assign[["Error"]]
MSgroups <- Assign[["MSgroups"]]
VK <- Assign[["VK"]]
MSAssign
Error
MSgroups
VK

### loop through other samples

#### this doesn't work in a loop, maybe has to do with the restarting R step in the Assign function?
library(MFAssignR)
path="for_MFAssignR/"
fns = dir(path = path)
dat_amb=dat_unamb=dat_no=NA
for(i in fns){
  Data <- read.csv(i)
  #Isotope prescreening
  Isotope <- IsoFiltR(Data)  #Input for IsoFiltR
  Mono <- Isotope[["Mono"]]
  Iso <- Isotope[["Iso"]]
  
  Assign <- MFAssign(Mono, Iso, ionMode = "neg", lowMW =50, highMW = 1000,  Nx = 3, Sx = 1,  ppm_err = .3, H_Cmin = 0.3,
                     HetCut = "off", DeNovo = 400, NMScut = "on")
  dat_unamb=rbind(dat_unamb,Assign[["Unambig"]])
  dat_amb=rbind(dat_amb,Assign[["Ambig"]])
  dat_no=rbind(dat_no,Assign[["None"]])
}
