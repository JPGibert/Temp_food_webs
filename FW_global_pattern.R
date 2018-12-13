## Set to directory where Food_web_metadata.csv and struct_data.csv are
setwd("DIRECTORY!")
library("quantreg")
library("maps")
library("maptools")
library("LSD")
library('RColorBrewer')
library("MASS")
library("raster")
library("rgdal")
library("ocedata")
library("akima")          
library("fields")  
library("zoo")
library("matrixLaplacian")


##Load and prepares dataset
web <- read.csv("Food_web_metadata.csv")
struct <- read.csv("struct_data.csv")
web <- cbind(web,struct)
head(web)
#####################################################################################################
## ADD LAYERS


## ADD temperature layers for LAND
Temp <- raster("~/Desktop/JP/Papers_in_progress/JP_Global_FW_patterns/wc2.0_2.5m_bio/wc2.0_bio_2.5m_01.tif")

# Stack + project for extraction of data
Temp_stack <- stack(Temp)
proj<-"+proj=longlat +datum=WGS84 +ellps=WGS84"
projection(Temp)<-proj

lat<-web$Lat
long<-web$Long
xydat<-cbind(long,lat)
	# Extract Tº data at the coordinates sampled in the Tº layer
point.dat<-extract(Temp_stack,xydat)

web$Temp <- point.dat
colnames(web)[34] <- "Temp"

## ADD temperature layers for OCEANS
data(levitus)
class(levitus)
summary(levitus)
levitus$longitude

web$Lat[is.na(web$Temp)]
web$Long[is.na(web$Temp)][1]

# Replace NAs with spline approximated values using na.approach from package zoo.
SST2 <- na.approx(levitus$SST)

Temps_to_sub <- rep(0,length(web$Lat[is.na(web$Temp)]))
for(i in 1:length(web$Lat[is.na(web$Temp)])){
	Temps_to_sub[i] <- levitus$SST[which.min(abs(levitus$longitude-web$Long[is.na(web$Temp)][i])),which.min(abs(levitus$latitude-web$Lat[is.na(web$Temp)][i]))]
}

which(is.na(web$Temp))
web$Temp[which(is.na(web$Temp))] <- Temps_to_sub

## Missing Temperatures from online search
web[c(4,15,26,34,50,51,52,53),]
web$Temp[4] <- -0.5 #Weddel sea
web$Temp[15] <- 16.8 #Chesapeake Bay
web$Temp[26] <- 24.2 #Huizache, Mexico
web$Temp[34] <- 12 #Monterey Bay
web$Temp[50] <- 8.5 #NZ
web$Temp[51] <- 8.5 #NZ
web$Temp[52] <- 8.5 #NZ
web$Temp[53] <- 8.5 #NZ

head(web)

######################################################################################################
## STRUCTURAL EQUATION MODELING
library("semPlot")
library("OpenMx")
library("lavaan")
library("tidyverse")
library("knitr")
library("kableExtra")
library("GGally")

## Because of wildely unequal variances, needs standardizing
web$absLat <- (abs(web$Lat)-mean(abs(web$Lat)))/sd(abs(web$Lat))
web$LTemp <- (abs(web$Temp)-mean(abs(web$Temp)))/sd(abs(web$Temp))
web$LS <- (abs(web$S)-mean(abs(web$S)))/sd(abs(web$S))
web$LL <- (abs(web$L)-mean(abs(web$L)))/sd(abs(web$L))
web$Top <- (abs(web$Top)-mean(web$Top,na.rm=TRUE))/sd(web$Top,na.rm=TRUE)
web$Basal <- (abs(web$Basal)-mean(web$Basal,na.rm=TRUE))/sd(web$Basal,na.rm=TRUE)
web$C <- (abs(web$C)-mean(abs(web$C)))/sd(abs(web$C))
web$LdivS <- (abs(web$L/web$S)-mean(abs(web$L/web$S)))/sd(abs(web$L/web$S)) 
web$Typenum <- as.numeric(web$Type)

web$LdivSAgg <- (web$LinksAgg/web$SAgg-mean(web$LinksAgg/web$SAgg,na.rm=TRUE))/sd(web$LinksAgg/web$SAgg,na.rm=TRUE)

web$CAgg <- (web$LinksAgg/(web$SAgg^2)-mean(web$LinksAgg/(web$SAgg^2),na.rm=TRUE))/sd(web$LinksAgg/(web$SAgg^2),na.rm=TRUE)
web$LSAgg <- (web$SAgg-mean(web$SAgg,na.rm=TRUE))/sd(web$SAgg,na.rm=TRUE)
web$LLAgg <- (web$LinksAgg-mean(web$LinksAgg,na.rm=TRUE))/sd(web$LinksAgg,na.rm=TRUE)
web$Typenum <- as.numeric(web$Type)
#web2 <- web[-which(web$Country=="New Zealand"),]


summary(lm(web$meanTL~web$S))
plot(web$meanTL~web$S) # Makes no sense
plot(web$meanTL~web$L)
plot(web$L~web$S)

## For non aggregated Food Webs
## Model below has the highest p-val which can be used as a measure of goodness of fit
	# Also, meanTL is strongly related to maxTL amd does not add much to analysis

# Temp + Lat
model_SEM1 <-'
# Regressions
TL ~ Basal + Top + LS + LL + LTemp + absLat 
Top ~ LL + LS + LTemp + absLat
Basal ~ LS + LTemp + absLat
LL ~ Basal + Top + LTemp + absLat + LS + Typenum
LS ~ LTemp + absLat + Typenum
LTemp ~ absLat
'

# Temp
model_SEM2 <-'
# Regressions
TL ~ Basal + Top + LS + LL + LTemp
Top ~ LL + LS + LTemp
Basal ~ LS + LTemp
LL ~ Basal + Top + LTemp + LS + Typenum
LS ~ LTemp + Typenum
'

# Lat
model_SEM3 <-'
# Regressions
TL ~ Basal + Top + LS + LL + absLat
Top ~ LL + LS + absLat
Basal ~ LS + absLat
LL ~ Basal + Top + absLat + LS + Typenum
LS ~ absLat + Typenum
'

# No enviro
model_SEM4 <-'
# Regressions
TL ~ Basal + Top + LS + LL
Top ~ LL + LS
Basal ~ LS
LL ~ Basal + Top + LS + Typenum
LS ~ Typenum
'
######################################


# Temp + Lat
model_SEM1 <-'
# Regressions
TL ~ Basal + Top + LdivS + LTemp + absLat 
C ~ Basal + Top + LdivS + LTemp + absLat
LdivS ~ Basal + Top + LTemp + absLat + Typenum
Top ~ Basal + LTemp + absLat + Typenum 
Basal ~ LTemp + absLat + Typenum
LTemp ~ absLat
'

# Temp
model_SEM2 <-'
# Regressions
TL ~ Basal + Top + LdivS + LTemp 
C ~ Basal + Top + LdivS + LTemp 
LdivS ~ Basal + Top + LTemp + Typenum
Top ~ Basal + LTemp + Typenum 
Basal ~ LTemp + Typenum
'

# Lat
model_SEM3 <-'
# Regressions
TL ~ Basal + Top + LdivS + absLat 
C ~ Basal + Top + LdivS + absLat
LdivS ~ Basal + Top + absLat + Typenum
Top ~ Basal + absLat + Typenum 
Basal ~ absLat + Typenum
'

# No enviro
model_SEM4 <-'
# Regressions
TL ~ Basal + Top + LdivS  
C ~ Basal + Top + LdivS 
LdivS ~ Basal + Top + Typenum
Top ~ Basal + Typenum 
Basal ~ Typenum
'

# Temp + Lat a bit better than Temp both much better than Lat! 

fit1 <- cfa(model_SEM1, data = web) 
fit2 <- cfa(model_SEM2, data = web) 
fit3 <- cfa(model_SEM3, data = web) 
fit4 <- cfa(model_SEM4, data = web) 

varTable(fit1)

fitMeasures(fit1, c("chisq", "df", "pvalue", "cfi", "rmsea", "SRMR", "agfi", "AIC"))
fitMeasures(fit2, c("chisq", "df", "pvalue", "cfi", "rmsea", "SRMR", "agfi", "AIC"))
fitMeasures(fit3, c("chisq", "df", "pvalue", "cfi", "rmsea", "SRMR", "agfi", "AIC"))
fitMeasures(fit4, c("chisq", "df", "pvalue", "cfi", "rmsea", "SRMR", "agfi", "AIC"))
summary(fit1, standardized=T, rsq=TRUE)
summary(fit2, standardized=T, rsq=T)
summary(fit3, standardized=T, rsq=T)
summary(fit4, standardized=T, rsq=T)
semPaths(fit2, 'std', layout = 'spring')

## For aggregated Food Webs
## Model below has the highest p-val which can be used as a measure of goodness of fit
	# Also, meanTL is strongly related to maxTL amd does not add much to analysis
plot(web$TL~web$)
# Temp + Lat
model_SEM1A <-'
# Regressions
TLAgg ~ BasalAgg + TopAgg + LdivSAgg + LTemp + absLat 
CAgg ~ BasalAgg + TopAgg + LdivSAgg + LTemp + absLat 
LdivSAgg ~ BasalAgg + TopAgg + LTemp + absLat + Typenum
TopAgg ~ BasalAgg + LTemp + absLat + Typenum
BasalAgg ~ LTemp + absLat + Typenum
LTemp ~ absLat
'

# Temp
model_SEM2A <-'
# Regressions
TLAgg ~ BasalAgg + TopAgg + LdivSAgg + LTemp 
CAgg ~ BasalAgg + TopAgg + LdivSAgg + LTemp  
LdivSAgg ~ BasalAgg + TopAgg + LTemp + Typenum
TopAgg ~ BasalAgg + LTemp + Typenum
BasalAgg ~ LTemp + Typenum
'

# Lat
model_SEM3A <-'
# Regressions
TLAgg ~ BasalAgg + TopAgg + LdivSAgg + absLat 
CAgg ~ BasalAgg + TopAgg + LdivSAgg + absLat 
LdivSAgg ~ BasalAgg + TopAgg + absLat + Typenum
TopAgg ~ BasalAgg + absLat + Typenum
BasalAgg ~ absLat + Typenum
'

# No enviro
model_SEM4A <-'
# Regressions
TLAgg ~ BasalAgg + TopAgg + LdivSAgg 
CAgg ~ BasalAgg + TopAgg + LdivSAgg 
LdivSAgg ~ BasalAgg + TopAgg + Typenum
TopAgg ~ BasalAgg + Typenum
BasalAgg ~ Typenum
'

fit1A <- cfa(model_SEM1A, data = web) 
fit2A <- cfa(model_SEM2A, data = web) 
fit3A <- cfa(model_SEM3A, data = web) 
fit4A <- cfa(model_SEM4A, data = web) 

fitMeasures(fit1A, c("chisq", "df", "pvalue", "cfi", "rmsea", "SRMR", "agfi", "AIC"))
fitMeasures(fit2A, c("chisq", "df", "pvalue", "cfi", "rmsea", "SRMR", "agfi", "AIC"))
fitMeasures(fit3A, c("chisq", "df", "pvalue", "cfi", "rmsea", "SRMR", "agfi", "AIC"))
fitMeasures(fit4A, c("chisq", "df", "pvalue", "cfi", "rmsea", "SRMR", "agfi", "AIC"))
summary(fit1A, standardized=T, rsq=T)
summary(fit2A, standardized=T, rsq=T)
summary(fit3A, standardized=T, rsq=T)
summary(fit4A, standardized=T, rsq=T)
web$typenum
semPaths(fit4, 'std', layout = 'spring')

fitMeasures(fit1)
parTable(fit1)

################################################
## Figure 1

## SEM for non aggregated food webs
#define the label that will go into the nodes
lbls<-c("Trophic Level","Connectance","L/S","Top Sp","Basal Sp","Temperature","Ecosystem type")
#define the layout
ly<-matrix(c(-0.20,-1, -0.4,-1, -0.1,-0.5, -0.1,0, -0.5,0, -0.4,0.5, -0.20,0.5),ncol=2,byrow=TRUE)
#new plot
#semPaths(fit,what="std",layout=ly,residuals=FALSE,nCharNodes=0,groups=grps,color=c("blue","green","Brown"),nodeLabels=lbls,sizeMan=8,posCol=c("green","red"),edge.label.cex=1.5,legend=FALSE)

semPaths(fit2,what="stand",layout=ly,residuals=FALSE,nCharNodes=0,nodeLabels=lbls,sizeMan=8,posCol=c("green","red"),edge.label.cex=1.5,legend=FALSE, fade=FALSE,intercepts=FALSE, edge.width=0.52/0.76) # edge size scaled to maximum edge width of model two for direct comparison of widths.

fitted(fit2)

## SEM for non aggregated food webs
#define the label that will go into the nodes
dev.new()
lbls<-c("Trophic Level","Connectance","L/S","Top Sp","Basal Sp","Temperature","Ecosystem type")
#define the layout
ly<-matrix(c(-0.20,-1, -0.4,-1, -0.1,-0.5, -0.1,0, -0.5,0, -0.4,0.5, -0.20,0.5),ncol=2,byrow=TRUE)
#new plot
#semPaths(fit,what="std",layout=ly,residuals=FALSE,nCharNodes=0,groups=grps,color=c("blue","green","Brown"),nodeLabels=lbls,sizeMan=8,posCol=c("green","red"),edge.label.cex=1.5,legend=FALSE)

semPaths(fit2A,what="stand",layout=ly,residuals=FALSE,nCharNodes=0,nodeLabels=lbls,sizeMan=8,posCol=c("green","red"),edge.label.cex=1.5,legend=FALSE, fade=FALSE,intercepts=FALSE)

fitted(fit2)


###################################
## Figure 2
## Calculate direct and indirect effects

## 1) Non-aggregated
Temp->L/S 
# Direct 
-0.30
# Indirect
-0.56*(-0.51) # 0.2856

0.47*5.09/0.56


Temp->Trophic Level 
# Direct
NA
# Indirect
-0.51*(-0.60)-0.56*(-0.51)*0.39-0.30*0.39 # 0.30038

Temp->Connectance
# Direct
NA
# Indirect
-0.51*(-0.47)-0.56*(-0.51)*0.41-0.3*0.41 # 0.233796


## 2) Aggregated
Temp->L/S 
# Direct 
NA
# Indirect
-0.49*(-0.40)	# 0.196

Temp->Trophic Level 
# Direct
-0.15
# Indirect
-0.49*(-0.72)-0.49*(-0.40)*0.4 # 0.4312

Temp->Connectance
# Direct
NA
# Indirect
-0.49*(-0.63)-0.49*(-0.40)*0.28 # 0.36358














