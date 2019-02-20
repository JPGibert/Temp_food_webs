setwd("~/Desktop/JP/Papers_in_review_submitted/JP_Global_FW_patterns")
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
#library("matrixLaplacian")
library("cheddar")

##Load and prepares dataset
web <- read.csv("Food_web_metadata.csv")
struct <- read.csv("struct_data.csv")
web <- cbind(web,struct)
head(web)
#####################################################################################################
## ADD LAYERS


## ADD temperature layers for LAND
Temp <- raster("~/Desktop/JP/Papers_in_review_submitted/JP_Global_FW_patterns/wc2.0_2.5m_bio/wc2.0_bio_2.5m_01.tif")

# Stack + project for extraction of data
Temp_stack <- stack(Temp)
proj<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
projection(Temp)<-proj

lat<-web$Lat
long<-web$Long
xydat<-cbind(long,lat)
	# Extract Tº data at the coordinates sampled in the Tº layer
point.dat<-extract(Temp_stack,xydat)
web$Temp <- point.dat
#colnames(web)[34] <- "Temp"
extract(Temp_stack,xydat)


## ADD temperature layers for OCEANS
data(levitus)
class(levitus)
summary(levitus)

# Replace NAs with spline approximated values using na.approach from package zoo.
SST2 <- na.approx(levitus$SST)

Temps_to_sub <- rep(0,length(web$Lat[is.na(web$Temp)]))
for(i in 1:length(web$Lat[is.na(web$Temp)])){
	Temps_to_sub[i] <- levitus$SST[which.min(abs(levitus$longitude-web$Long[is.na(web$Temp)][i])),which.min(abs(levitus$latitude-web$Lat[is.na(web$Temp)][i]))]
}

web$Temp[which(is.na(web$Temp))] <- Temps_to_sub

## Missing Temperatures from online search (see details in Appendix 2)
web[c(6,19,32,41,57,58,59,60),]
web$Temp[6] <- -0.5 #Weddel sea (Assumes similar water surface temperature than that for the Scotia-Weddel confluence zone in Levitus)
web$Temp[19] <- 16.8 #Chesapeake Bay
web$Temp[32] <- 30 #24.2 #Huizache, Mexico (From paper, see dataset)
web$Temp[41] <- 12 #Monterey Bay
web$Temp[57] <- 9 #NZ
web$Temp[58] <- 9 #NZ
web$Temp[59] <- 9 #NZ
web$Temp[60] <- 9 #NZ

## JUST for analyses without these 8 food webs
# web <- web[-c(6,19,32,41,57,58,59,60),]

head(web)
plot(TL~meanTL, data=web, las=TRUE, ylab="maxTL", xlab="meanTL")
abline(lm(web$TL~web$meanTL))
plot(Int~Basal, data=web, las=TRUE, ylab="Fraction Intermediate", xlab="Fraction Basal")
abline(lm(Int~Basal, data=web))
summary(lm(Int~Basal, data=web))
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
web$Omnv <- (abs(web$Omnv)-mean(abs(web$Omnv),na.rm=TRUE))/sd(abs(web$Omnv),na.rm=TRUE)

web$LdivSAgg <- (web$LinksAgg/web$SAgg-mean(web$LinksAgg/web$SAgg,na.rm=TRUE))/sd(web$LinksAgg/web$SAgg,na.rm=TRUE)

web$CAgg <- (web$LinksAgg/(web$SAgg^2)-mean(web$LinksAgg/(web$SAgg^2),na.rm=TRUE))/sd(web$LinksAgg/(web$SAgg^2),na.rm=TRUE)
web$LSAgg <- (web$SAgg-mean(web$SAgg,na.rm=TRUE))/sd(web$SAgg,na.rm=TRUE)
web$LLAgg <- (web$LinksAgg-mean(web$LinksAgg,na.rm=TRUE))/sd(web$LinksAgg,na.rm=TRUE)
web$OmnvAgg <- (abs(web$OmnvAgg)-mean(abs(web$OmnvAgg),na.rm=TRUE))/sd(abs(web$OmnvAgg),na.rm=TRUE)
#web$Typenum <- as.numeric(web$Type)
#web2 <- web[-which(web$Country=="New Zealand"),]


## For non aggregated Food Webs
## Model below has the highest p-val which can be used as a measure of goodness of fit
	# Also, meanTL is strongly related to maxTL amd does not add much to analysis

# Temp + Lat
model_SEM1 <-'
# Regressions
TL ~ Basal + Top + LS + LL + LTemp + absLat
C ~ Basal + Top + LS + LL + LTemp + absLat
Omnv ~ Basal + Top + LS + LL + LTemp + absLat
Top ~ LS + Basal + LTemp + absLat + Typenum
Basal ~ LS + LTemp + absLat + Typenum
LL ~ Basal + Top + LTemp + absLat + LS + Typenum
LS ~ LTemp + absLat + Typenum
LTemp ~ absLat
'

# Temp
model_SEM2 <-'
# Regressions
TL ~ Basal + Top + LS + LL + LTemp
C ~ Basal + Top + LS + LL + LTemp
Omnv ~ Basal + Top + LS + LL + LTemp
Top ~ LS + Basal + LTemp + Typenum
Basal ~ LS + LTemp + Typenum
LL ~ Basal + Top + LTemp + LS + Typenum
LS ~ LTemp + Typenum
'

# Lat
model_SEM3 <-'
# Regressions
TL ~ Basal + Top + LS + LL + absLat
C ~ Basal + Top + LS + LL + absLat
Omnv ~ Basal + Top + LS + LL + absLat
Top ~ LS + Basal + absLat + Typenum
Basal ~ LS + absLat + Typenum
LL ~ Basal + Top + absLat + LS + Typenum
LS ~ absLat + Typenum
'

# No enviro
model_SEM4 <-'
# Regressions
TL ~ Basal + Top + LS + LL
C ~ Basal + Top + LS + LL
Omnv ~ Basal + Top + LS + LL
Top ~ LS + Basal + Typenum
Basal ~ LS + Typenum
LL ~ Basal + Top + LS + Typenum
LS ~ Typenum
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
TLAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg + LTemp + absLat 
CAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg + LTemp + absLat
OmnvAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg + LTemp + absLat 
TopAgg ~ LSAgg + BasalAgg + LTemp + absLat + Typenum
BasalAgg ~ LSAgg + LTemp + absLat + Typenum
LLAgg ~ BasalAgg + TopAgg + LSAgg + LTemp + absLat + Typenum
LSAgg ~ LTemp + absLat + Typenum
LTemp ~ absLat
'

# Temp
model_SEM2A <-'
# Regressions
TLAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg + LTemp 
CAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg + LTemp
OmnvAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg + LTemp 
TopAgg ~ LSAgg + BasalAgg + LTemp + Typenum
BasalAgg ~ LSAgg + LTemp + Typenum
LLAgg ~ BasalAgg + TopAgg + LSAgg + LTemp + Typenum
LSAgg ~ LTemp + Typenum
'

# Lat
model_SEM3A <-'
# Regressions
TLAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg + absLat 
CAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg + absLat
OmnvAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg + absLat 
TopAgg ~ LSAgg + BasalAgg + absLat + Typenum
BasalAgg ~ LSAgg + absLat + Typenum
LLAgg ~ BasalAgg + TopAgg + LSAgg + absLat + Typenum
LSAgg ~ absLat + Typenum
'

# No enviro
model_SEM4A <-'
# Regressions
TLAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg
CAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg
OmnvAgg ~ BasalAgg + TopAgg + LSAgg + LLAgg 
TopAgg ~ LSAgg + BasalAgg + Typenum
BasalAgg ~ LSAgg + Typenum
LLAgg ~ BasalAgg + TopAgg + LSAgg + Typenum
LSAgg ~ Typenum
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
semPaths(fit2, 'std', layout = 'spring')

fitMeasures(fit1)
parTable(fit1)

################################################
## Figure 2

## 1) SEM for non aggregated food webs
#define the label that will go into the nodes
lbls<-c("Trophic Level","Connectance","Omnivory","Top Sp.","Basal Sp.","N Links","Species","Temperature","Ecosystem type")
#define the layout
ly<-matrix(c(-0.10,-0.8, -0.4,-0.8, -0.7,-0.8, 0.1,0, -0.9,0, -0.3,-0.3, -0.5,0.3, -0.7,0.8, -0.1,0.8),ncol=2,byrow=TRUE)
#new plot
#semPaths(fit,what="std",layout=ly,residuals=FALSE,nCharNodes=0,groups=grps,color=c("blue","green","Brown"),nodeLabels=lbls,sizeMan=8,posCol=c("green","red"),edge.label.cex=1.5,legend=FALSE)

semPaths(fit2,what="stand",layout=ly,residuals=FALSE,nCharNodes=0,nodeLabels=lbls,sizeMan=8,posCol=c("green","red"),edge.label.cex=1.5,legend=FALSE, fade=FALSE,intercepts=FALSE,structural=FALSE, edge.width=0.52/0.76) # edge size scaled to maximum edge width of model two for direct comparison of widths.

fitted(fit2)

## 2) SEM for aggregated food webs (Appendix)
#define the label that will go into the nodes
lbls<-c("Trophic Level","Connectance","Omnivory","Top Sp.","Basal Sp.","N Links","Species","Temperature","Ecosystem type")
#define the layout
ly<-matrix(c(-0.10,-0.8, -0.4,-0.8, -0.7,-0.8, 0.1,0, -0.9,0, -0.3,-0.3, -0.5,0.3, -0.7,0.8, -0.1,0.8),ncol=2,byrow=TRUE)
#new plot
#semPaths(fit,what="std",layout=ly,residuals=FALSE,nCharNodes=0,groups=grps,color=c("blue","green","Brown"),nodeLabels=lbls,sizeMan=8,posCol=c("green","red"),edge.label.cex=1.5,legend=FALSE)

semPaths(fit2A,what="stand",layout=ly,residuals=FALSE,nCharNodes=0,nodeLabels=lbls,sizeMan=8,posCol=c("green","red"),edge.label.cex=1.5,legend=FALSE, fade=FALSE,intercepts=FALSE, edge.width=0.52/0.76) # edge size scaled to maximum edge width of model two for direct comparison of widths.
fitted(fit2)


###################################
## Figure 3
## Calculate direct and indirect effects

## BIOTIC PROPERTIES ---------------------------
## 1) Non-aggregated
Temp->Species 
# Direct 
-0.28
# Indirect
0

Temp->Links 
# Direct
-0.21
# Indirect
-0.41*(-0.5) -0.28*(0.47)*(-0.5) -0.28*0.68  # =0.0804
# Total
0.0804-0.21  # =-0.1296

Temp->Basal
# Direct
-0.41
# Indirect
-0.28*0.47  # =-0.1316
# Total 
-0.41-0.136 # =-0.546

## STRUCTURAL PROPERTIES ---------------------------
## 1) Non-aggregated
Temp->Omnivory 
# Direct 
0.21
# Indirect
-0.41*(-0.57)-0.28*(-0.25) -0.28*(0.47)*(-0.57) -0.28*0.68*0.33 -0.21*0.33 -0.28*0.47*(-0.5)*0.33 -0.41*(-0.5)*0.33 # =0.335944
# Total
0.21+0.335944 # =0.545944

0.2028*3/0.3

Temp->Connectance 
# Direct
NA
# Indirect
	# With all food webs 
-0.41*(-0.23) -0.28*(0.47)*(-0.23) -0.28*(-0.69) -0.28*(0.68)*0.44 -0.21*0.44 -0.41*(-0.5)*0.44 -0.28*0.47*(-0.5)*0.44  # =0.2607
# Total
0.2607

Temp->Trophic Level
# Direct
NA
# Indirect
-0.41*(-0.72)-0.41*(-0.5)*0.35 -0.28*0.47*(-0.72) -0.28*0.47*(-0.5)*0.35 -0.28*0.68*0.35 -0.21*0.35 # =0.344592
# Total 
0.344592



# Without 8 food webs for which manual online search for temps. was needed
## BIOTIC PROPERTIES ---------------------------
## 1) Non-aggregated
Temp->Species 
# Direct 
-0.29
# Indirect
0

Temp->Links 
# Direct
-0.30
# Indirect
-0.46*(-0.51) -0.29*(0.40)*(-0.51) -0.29*0.64  # =0.10816
# Total
0.10816-0.30  # =-0.19184

Temp->Basal
# Direct
-0.46
# Indirect
-0.29*0.40  # =-0.116
# Total 
-0.46-0.116 # =-0.576

## STRUCTURAL PROPERTIES ---------------------------
## 1) Non-aggregated
Temp->Omnivory 
# Direct 
0.27
# Indirect
-0.46*(-0.51) -0.29*(-0.23) -0.29*(0.40)*(-0.51) -0.29*0.64*0.35 -0.30*0.35 -0.29*0.40*(-0.51)*0.35 -0.46*(-0.51)*0.35 # =0.293316
# Total
0.27 + 0.293316 # =0.563316

Temp->Connectance 
# Direct
NA
# Indirect
-0.29*(-0.70) -0.29*(0.64)*0.49 -0.30*0.49 -0.46*(-0.51)*0.49 -0.29*0.40*(-0.51)*0.49  # =0.1028
# Total
0.1028

Temp->Trophic Level
# Direct
NA
# Indirect
-0.46*(-0.67) -0.46*(-0.51)*0.36 -0.29*0.40*(-0.67) -0.29*0.40*(-0.51)*0.36 -0.29*0.64*0.36 -0.30*0.36 # =0.3168
# Total 
0.3168


# AGGREGATED
##------------------------------------------------------------------------------------------------------------------
## BIOTIC PROPERTIES ---------------------------
## 1) Aggregated
Temp->Species 
# Direct 
-0.23
# Indirect
0

Temp->Links 
# Direct
0
# Indirect
-0.31*(-0.41) -0.23*(0.50)*(-0.41) -0.23*0.86  # =-0.0235
# Total
-0.0235  # =-0.0235

Temp->Basal
# Direct
-0.31
# Indirect
-0.23*0.50  # =-0.115
# Total 
-0.31-0.115 # =-0.425

## STRUCTURAL PROPERTIES ---------------------------
## 2) Aggregated
Temp->Omnivory 
# Direct 
0.12
# Indirect
-0.31*(-0.74) -0.31*(-0.41)*0.34 -0.23*0.5*(-0.74)	-0.23*(-0.19) -0.23*0.50*(-0.41)*0.34 # =0.417445
# Total
0.12+0.417445  # =0.537445

Temp->Connectance 
# Direct
NA
# Indirect
-0.31*(-0.43) -0.31*(-0.41)*0.40 -0.23*0.5*(-0.43) -0.23*(-0.6) -0.23*0.5*(-0.41)*0.4 -0.23*0.86*0.4 # =0.31133 
# Total
0.31133

Temp->Trophic Level
# Direct
-0.13
# Indirect
-0.31*(-0.84) -0.31*(-0.41)*0.39 -0.23*0.5*(-0.41)*0.39 -0.23^0.5*(-0.84) -0.23*0.86*0.39 # 0.6540653
# Total
-0.13+0.6540 # =0.524


################################################
## Figure 1
## MAP v1
size <- 1.1                               
map('world',interior=FALSE,col='gray', fill=TRUE,lty=0,xlim=c(-165,190), ylim=c(-70,130))
points(web$Lat[web$Country2!="New Zealand" & web$Country2!="Brazil" & web$Country2!="UK"]~web$Long[web$Country2!="New Zealand" & web$Country2!="Brazil" & web$Country2!="UK"], col=c("purple","orange")[web$Type2[web$Country2!="New Zealand" & web$Country2!="Brazil" & web$Country2!="UK"]], pch=16, cex=size)
points(jitter(web$Lat[web$Country2=="New Zealand"],factor=4000)~jitter(web$Long[web$Country2=="New Zealand"],factor=8000), col=c("purple","orange")[web$Type2[web$Country2=="New Zealand"]], pch=16, cex=size)
points(jitter(web$Lat[web$Country2=="Brazil"], factor=200)~web$Long[web$Country2=="Brazil"], col=c("purple","orange")[web$Type2[web$Country2=="Brazil"]], pch=16, cex=size)
points(web$Lat[web$Country2=="UK"]~jitter(web$Long[web$Country2=="UK"], factor=1500), col=c("purple","orange")[web$Type2[web$Country2=="UK"]], pch=16, cex=size)


## MAP v2
size <- 0.8                              
map('world',interior=FALSE,col='gray', fill=TRUE,lty=0,xlim=c(-165,190), ylim=c(-70,130))
points(web$Lat[web$Country2!="New Zealand" & web$Country2!="Brazil" & web$Country2!="UK"]~web$Long[web$Country2!="New Zealand" & web$Country2!="Brazil" & web$Country2!="UK"], pch=c(3,16,15,17)[web$Type[web$Country2!="New Zealand" & web$Country2!="Brazil" & web$Country2!="UK"]], cex=size)
points(jitter(web$Lat[web$Country2=="New Zealand"],factor=4000)~jitter(web$Long[web$Country2=="New Zealand"],factor=8000), pch=c(3,16,15,17)[web$Type[web$Country2=="New Zealand"]], cex=size)
points(jitter(web$Lat[web$Country2=="Brazil"], factor=200)~web$Long[web$Country2=="Brazil"], pch=c(3,16,15,17)[web$Type[web$Country2=="Brazil"]], cex=size)
points(web$Lat[web$Country2=="UK"]~jitter(web$Long[web$Country2=="UK"], factor=1500), pch=c(3,16,15,17)[web$Type[web$Country2=="UK"]], cex=size)








## THE END
