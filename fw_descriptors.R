setwd("~/Desktop/JP/Papers_in_review_submitted/JP_Global_FW_patterns/Food_webs_POSTA")

library("foodweb")

#########################################################################################################
#			Food web structure CALCULATOR
# The function below calculates all species TLs, as well as prop basal, intermediate and top
#########################################################################################################
## Normal food webs
TL_calculator <- function(MAT){
	# Find dimensions of matrix
	dims <- dim(MAT)
	# Crate array of trophic levels and basal, intermediate and top, as well as generality and vulnerability
	gen_vul <- rep(0,0)
	props <- rep(0,0,0)
	TLs <- rep(0,dims[1])
	# Find primary producers, assign TL=1 and add to basal tally
	TLs[rowSums(MAT)==0] <- 1
	props[1] <- sum(TLs)
	# Now calculate all other TLs
	for(j in 1:100){ 
		# We need to iterate several times because the first time the TLs array will only have 0s and 1s
		# and as it starts filling up the measures of TL will change.
		for(i in seq(1:dims[1])[TLs!=1]){
			TLs[i] = 1 + sum(MAT[i,]%*%t(t(TLs)) * (1/sum(MAT[i,]))) 
		} 
	}
	# Now calculate top and intermediate
	# Top predators are not eaten (so marginal column sum == 0), and can't be herbivores (thus TLs>2)
	props[3] <- sum(colSums(MAT)[TLs>2]==0) 
	props[2] <- dims[1]-props[1]-props[3]
	# Calculate average generality and vulnerability (here a standardized measure of variation in gen and vul)
	gen_vul[1] <- sd(colSums(MAT))/mean(colSums(MAT))
	gen_vul[2] <- sd(rowSums(MAT))/mean(rowSums(MAT))
	# Return TLs
	return(list(TLs, props/dims[1], gen_vul, dim(MAT)[1]))
}

###############################################################################################################
# This could be done more elegantly and faster by running the whole thing once but I wanted to make sure everything was working correctly.

## 1) (NON AGGREGATED) Calculate TLs and basal, int and top for all food webs 
setwd("~/Desktop/JP/Papers_in_review_submitted/JP_Global_FW_patterns/Food_webs_POSTA")
fws <- dir()[-1]
fws <- fws[grep("_Adj.txt",fws)]
fws <- fws[-18] # Eliminate Carpinteria,

struct_data <- matrix(rep(0,10*length(fws)), nrow=length(fws), ncol=10)
for(i in 1:length(fws)){	
	matrix <- t(as.matrix(read.table(fws[i]))) # Needs to be transposed so that 1s means i eats j.
	results <- TL_calculator(matrix) 
		# To quantify omnivory I use the number of species with non-integer trophic level. I divide that by the total number of species to get the proportion of omnivores.
	omnivores <- sum(results[[1]]!=round(results[[1]]),na.rm=TRUE)/length(results[[1]])
	struct_data[i,] <- c(max(results[[1]]),mean(results[[1]]),omnivores,results[[2]],sum(matrix),results[[3]],results[[4]])
}


struct.data <- data.frame(struct_data)
struct.data <- rbind(struct.data[1:17,],rep(NA,5),struct.data[18:dim(struct.data)[1],])
colnames(struct.data) <- c("TL","meanTL","Omnv","Basal","Int","Top","Links", "Gen", "Vul", "Taxa")


## 2) (AGGREGATED) Calculate TLs and basal, int and top for all food webs 
setwd("~/Desktop/JP/Papers_in_review_submitted/JP_Global_FW_patterns/Food_webs_POSTA")
fws <- dir()[-1]
fws <- fws[grep("_Adj.txt",fws)]
fws <- fws[-18] # Eliminate Carpinteria,

struct_data2 <- matrix(rep(0,10*length(fws)), nrow=length(fws), ncol=10)
for(i in 1:length(fws)){	
	# Load food web
	matrix <- t(as.matrix(read.table(fws[i]))) # Needs to be transposed so that 1s means i eats j.
	# Aggregate	
	# Aggregate producers by finding who eats nothing (rowSums(matrix)==0), fetching the columns of the matrix, transposing them for "duplicated" (which only works on rows), amd finding how they compare to one another
	# We find the indices of the producers that stay
	basal <- which(rowSums(matrix)==0)
	# We need to transpose the matrix for this because we are comparing basal species by who eats them 
	basal_index <- basal[!duplicated(t(matrix[,rowSums(matrix)==0]))]
	# Aggregate consumers by finding who eats something first (rowSums(matrix)>0), then removing duplicates
	non_basal <- which(rowSums(matrix)>0)
	NONbasal_index <- non_basal[!duplicated(matrix[rowSums(matrix)>0,])]
	# We concatenate new indices and make the new matrix
	matrix <- matrix[c(basal_index,NONbasal_index),c(basal_index,NONbasal_index)]
	# Measure structure
	results <- TL_calculator(matrix) 
	
		# To quantify omnivory I use the number of species with non-integer trophic level. I divide that by the total number of species to get the proportion of omnivores.
	omnivores <- sum(results[[1]]!=round(results[[1]]),na.rm=TRUE)/length(results[[1]])
	
	struct_data2[i,] <- c(max(results[[1]]),mean(results[[1]]),omnivores,results[[2]],sum(matrix),results[[3]],results[[4]])
}

struct.data2 <- data.frame(struct_data2)
struct.data2 <- rbind(struct.data2[1:17,],rep(NA,5),struct.data2[18:dim(struct.data2)[1],])
colnames(struct.data2) <- c("TLAgg","meanTLAgg","OmnvAgg","BasalAgg","IntAgg","TopAgg","LinksAgg", "GenAgg", "VulAgg", "SAgg")
struct.data <- cbind(struct.data,struct.data2)
setwd("~/Desktop/JP/Papers_in_review_submitted/JP_Global_FW_patterns/")
write.table(struct.data, "struct_data_NEW.csv",sep=",",row.names=FALSE)



##########################################################################################################
## Create all the .csv files from .txt
fws <- dir()[-1]
fws <- fws[grep("_Adj.txt",fws)]

for(i in 1:length(fws)){
	setwd("~/Desktop/JP/Papers_in_progress/JP_Global_FW_patterns/Food_webs_POSTA")
	fw_i <- read.table(fws[i])
	setwd("~/Desktop/JP/Papers_in_progress/JP_Global_FW_patterns/Fodd_webs_csv")
	write.table(fw_i,paste(substr(fws[i],1,(nchar(fws[i])-4)),".csv", sep=""), sep=",",col.names=FALSE ,row.names=FALSE)
}

# plot food web
plotweb(c("black", "blue", "red", "green", "yellow"), radii=c(20,10,10,10,10))

fws.csv <- dir()
names <- fws.csv[-c(14,41,42)]
write.table(names,"FW_names.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)


prue <- t(as.matrix(read.table("Carpinteria_Adj.txt")))



