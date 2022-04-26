#clear workspace
  rm(list=ls()) 

#set working directory - wherever the source files are save
  setwd("/Users/Evan_Craine/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/gitHub/wetChem")
  
#prepare packages
  library(agricolae)
  #install.packages("dplyr")
  library(dplyr)
  #install.packages("pastecs")
  library(pastecs)
  #install.packages(ggplot2)
  library(ggplot2)
  library(qwraps2)
  library(data.table)
  #install.packages("tables")
  library(tables) #need tables package for tabular function
  library(psych)
  
#read in data
  AdjMaster = read.csv("wetChem_adj_FINAL.csv") #AAs are in g/100g crude protein
  RawMaster = read.csv("wetChem_OGunits_FINAL.csv") #AAs are in g/100g sample
  #master = read.csv("R_WSU_quinoa_wetchem_NIR_stats.csv")
#import data set
  #wetChem_adj_FINAL$POP <- as.factor(wetChem_adj_FINAL$POP)
  #wetChem_adj_FINAL$LOCATION <- as.factor(wetChem_adj_FINAL$LOCATION)
  
  
  
#check column names
  head(AdjMaster) #0 in pop is a commercial variety/check

  
#check the structure of the data frame
  str(AdjMaster)

  
  
#pop needs to be a factor with levels
  AdjMaster$POP <- as.factor(AdjMaster$POP)
  AdjMaster$LOCATION <- as.factor(AdjMaster$LOCATION)
  str(AdjMaster) #check to make sure it is fixed

#summary stats by groups for proximates
  AdjDescrip <- describe(AdjMaster)
  
  summary(AdjMaster)
  
#df is not a data frame...?!
  is.data.frame(AdjMaster)
  as.data.frame(AdjMaster)

  
  
#$$$___descriptive statistics___$$$#
#proximates# #*********use only if you want TotalAA in g/100gsample

  #tabular operates row ~ columns; row + 1 gives you all groups at bottom; n = 1 lists the the number of samples per group in the table
  
  # remove Crude Fiver
  DescAll <- tabular( (LOCATION + POP + 1) ~ (n=1) + Format(digits=2) * (Total.AA + Crude.protein + Moisture + Crude.Fat + Ash) * (mean + sd + min + max), data=AdjMaster) #these are the columns 
  DescAll
  
  DescDets <- tabular( (1) ~ (n=1) + Format(digits=2) * (Total.AA + Crude.protein + Moisture + Crude.Fat + Ash) * (min + max + mean + sd), data=AdjMaster) #these are the columns 
  DescDets
#save as .csv  
  write.csv.tabular(DescAll, file = "ADJ_PROX_Descriptives_LOC_POP.csv", justification = "n", row.names=FALSE)
  write.csv.tabular(DescDets, file = "ADJ_PROX_Descriptives All.csv", justification = "n", row.names=FALSE)
  
#$$$___descriptive statistics___$$$#
# want descriptives for the AA profile
  DescAA <- tabular( (1) ~ Format(digits=2) * (Total.Essential.AA + Leucine + Lysine + Valine + Isoleucine + Phenylalanine + Threonine + Histidine + Methionine + Tryptophan + Total.Nonessential.AA + Glutamic.Acid + Aspartic.Acid + Arginine + Glycine + Alanine + Proline + Serine + Tyrosine + Cysteine + Taurine + Hydroxyproline + Hydroxylysine + Ornithine) * (mean + sd + min + max), data=AdjMaster) #these are the columns 
  DescAA
  
  write.csv.tabular(DescAA, file = "ADJ_AAprofile_Descriptives_FINAL.csv", justification = "n", row.names=FALSE)
  
  #want descriptives for the essential AA profile by 'groups'
  
  DescAAdeets <- tabular( (LOCATION + POP + 1) ~ (n=1) + Format(digits=2) * (Total.Essential.AA + AAA + Leucine +	Lysine + Valine + SAA + Isoleucine + Threonine + Histidine + Tryptophan) * (mean + sd), data=AdjMaster) #these are the columns 
  DescAAdeets
  
  
  
  write.csv.tabular(DescAAdeets, file = "ADJ_EssentialAA_Descriptives_FINAL.csv", justification = "n", row.names=FALSE)
  
  
  
  #$$$___ranking samples for nutritional composition___$$$#
     # copy master to new dataframe 'showdown'
     
   showdown <- copy(AdjMaster)
  str(showdown)
  # remove identifiers, Lanthionine, newEntry, locEntry, newLocation
  showdown <- showdown[,c(10:16, 18:43)] #46 variables to 33...
  head(showdown)
  str(showdown)
  
# create a vector of the colnames (i.e. nutritional attribute names)  
  attributes <- as.vector(colnames(showdown)) #this returns df column names
  attributes
  
# create new series of coloumns for attribute ranks
  
  newcolumn <- paste0(attributes, " Rank") #adds rank to each
  showdown[,newcolumn] <- NA #addes empty columns to data frame using names from newcolumn vector
  showdown[,"overall rank"] <- NA #creates the last column to sum the ranks in
  head(showdown)
  str(show) #doubles the number of columns (66)...
  
# test meat of for loop
  for (i in 34:66) { #start at first parameter rank columnn number and go to last parameter rank column number
     showdown[,34]<-rank(showdown[,(34-66)]) #replace 47 with i in for loop; this provides a ranking to the tenth, so you can easily deal with ties
     
  }
  
# for loop construction and execution
  for (i in 47:81) { #start at first parameter rank columnn and go to last parameter rank column
     showdown[,i]<-rank(showdown[,(i-34)]) #replace 47 with i in for loop
  }

# visual representation of populations in study
  hist(as.numeric(AdjMaster$POP),
       main = "Populations",
       xlab = "Population ID",
       border = "black",
       col = "blue")


# seed composition
  hist(AdjMaster$Crude.protein,
       main = "Crude Protein Content",
       xlab = "Protein (g/100g)",
       border = "black",
       col = "blue")
  
  hist(AdjMaster$Crude.Fat,
       main = "Crude Fat Content",
       xlab = "Fat (g/100g)",
       border = "black",
       col = "blue")
  
  hist(AdjMaster$Crude.Fiber, 
       main = "Crude Fiber Content",
       xlab = "Fiber (g/100g)",
       border = "black",
       col = "blue")
  
  hist(AdjMaster$Moisture, 
       main = "Moisture Content (as is)",
       xlab = "Moisture (g/100g)",
       border = "black",
       col = "blue")
  
  hist(AdjMaster$Ash,
       main = "Ash Content",
       xlab = "Ash (g/100g)",
       border = "black",
       col = "blue")

#visual representation of amino acids
  hist(AdjMaster$Total.AA,
       main = "Total Amino Acid Content",
       xlab = "Total Amino Acid (g/100g)",
       border = "black",
       col = "blue")
  
  hist(AdjMaster$Taurine,
       main = "Taurine Content",
       xlab = "Taurine (g/100g)",
       border = "black",
       col = "blue")
  
  hist(AdjMaster$Hydroxyproline,
       main = "Hydroxyproline Content",
       xlab = "Hydroxyproline (g/100g)",
       border = "black",
       col = "blue")
 
   hist(AdjMaster$Aspartic.Acid,
       main = "Aspartic Acid Content",
       xlab = "Aspartic Acid (g/100g)",
       border = "black",
       col = "blue")
  
   hist(AdjMaster$Threonine,
        main = "Threonine Content",
        xlab = "Threonine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Serine,
        main = "Serine Content",
        xlab = "Serine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Glutamic.Acid,
        main = "Glutamic Acid Content",
        xlab = "Glutamic Acid (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Proline,
        main = "Proline Content",
        xlab = "Proline (g/100g)",
        border = "black",
        col = "blue")
   
   #hist(AdjMaster$Lanthionine._,  ##### all zeroes!!!!#####
    #    main = "Lanthionine Content",
    #    xlab = "Lanthionine (g/100g)",
    #    border = "black",
    #    col = "blue")
   
   hist(AdjMaster$Glycine,
        main = "Glycine Content",
        xlab = "Glycine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Alanine,
        main = "Alanine Content",
        xlab = "Alanine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Cysteine,
        main = "Cysteine Content",
        xlab = "Cysteine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Valine,
        main = "Valine Content",
        xlab = "Valine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Methionine,
        main = "Methionine Content",
        xlab = "Methionine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Isoleucine,
        main = "Isoleucine Content",
        xlab = "Isoleucine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Leucine,
        main = "Leucine Content",
        xlab = "Leucine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Tyrosine,
        main = "Tyrosine Content",
        xlab = "Tyrosine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Phenylalanine,
        main = "Phenylalanine Content",
        xlab = "Phenylalanine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Hydroxylysine,
        main = "Hydroxylysine Content",
        xlab = "Hydroxylysine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Ornithine._,
        main = "Ornithine Content",
        xlab = "Ornithine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Lysine,
        main = "Lysine Content",
        xlab = "Lysine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Histidine,
        main = "Histidine Content",
        xlab = "Histidine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Arginine,
        main = "Arginine Content",
        xlab = "Arginine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(AdjMaster$Tryptophan,
        main = "Tryptophan Content",
        xlab = "Tryptophan (g/100g)",
        border = "black",
        col = "blue")

#********************************linear relationships********************************#
   #$$$$$test for normality$$$$$# 
#examples; HO=data are normally distributed; HA=data are not normally distributed; if p value is greater than 0.05, we accept H0; if p value is less than 0.05 we reject H0 and data are not normally distributed
   # Shapiro-Wilk normality test for mpg
   shapiro.test(my_data$mpg) # => p = 0.1229
   # Shapiro-Wilk normality test for wt
   shapiro.test(my_data$wt) # => p = 0.0
#create an empty data frame
   NormalTest <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("parameter", "statistic", "p value", "normal?"))
#shorten AdjMaster to get just the parameters we're interested in
   AdjMasterMinus <- AdjMaster[,c(10:16, 18:43)] #just the nutritional parameters
#create a vector of the parameters   
   factors <- colnames(AdjMasterMinus) #33 parameters we're interested in 
#create a for loop to test all parameters for normality and report out the results   
    for (i in 1:33) {
       Ntest <- shapiro.test(AdjMasterMinus[,i])
       NormalTest[(i), 1] <- factors[i]
       NormalTest[(i), 2] <- Ntest$statistic
       NormalTest[(i), 3] <- Ntest$p.value
       NormalTest[(i), 4] <- if (Ntest$p.value > 0.05){
         print('Yes')
       } else {
         print('No')
       }
    }  
   
   #$$$$$inspect QQ plots$$$$$# 
   library("ggpubr")
  
   ggqqplot(AdjMaster$Taurine, ylab = "Taurine") #normal
   ggqqplot(AdjMaster$Methionine, ylab = "Methionine") #nonnormal, but pretty close
   ggqqplot(AdjMaster$Hydroxylysine, ylab = "Hydroxylysine") #nonnormal
   ggqqplot(AdjMaster$Ornithine, ylab = "Ornithine") #nonnormal
   ggqqplot(AdjMaster$Arginine, ylab = "Arginine") #nonnormal
   ggqqplot(AdjMaster$Moisture, ylab = "Moisture") #nonnormal
   ggqqplot(AdjMaster$Crude.Fiber, ylab = "Crude Fiber") #nonnormal
   ggqqplot(AdjMaster$Ash, ylab = "Ash") #nonnormal
   
   #$$$$$CORRELATIONS$$$$$# 
   
   #load the required packages/functions
   library(Hmisc) #rcorr function is present in this package
   #install.packages("corrplot")
   library(corrplot)
   library(PerformanceAnalytics)
   library(RColorBrewer)

# generate correlation table with sig values = asterisk
   #significance levels
   #first run corstars function
   # x is a matrix containing the data
   # method : correlation method. "pearson"" or "spearman"" is supported
   # removeTriangle : remove upper or lower triangle
   # results :  if "html" or "latex"
   # the results will be displayed in html or latex format
   corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                       result=c("none", "html", "latex")){
      #Compute correlation matrix
      require(Hmisc)
      x <- as.matrix(x)
      correlation_matrix<-rcorr(x, type=method[1])
      R <- correlation_matrix$r # Matrix of correlation coeficients
      p <- correlation_matrix$P # Matrix of p-value 
      
      ## Define notions for significance levels; spacing is important.
      mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
      
      ## trunctuate the correlation matrix to two decimal
      R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
      
      ## build a new matrix that includes the correlations with their apropriate stars
      Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
      diag(Rnew) <- paste(diag(R), " ", sep="")
      rownames(Rnew) <- colnames(x)
      colnames(Rnew) <- paste(colnames(x), "", sep="")
      
      ## remove upper triangle of correlation matrix
      if(removeTriangle[1]=="upper"){
         Rnew <- as.matrix(Rnew)
         Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
         Rnew <- as.data.frame(Rnew)
      }
      
      ## remove lower triangle of correlation matrix
      else if(removeTriangle[1]=="lower"){
         Rnew <- as.matrix(Rnew)
         Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
         Rnew <- as.data.frame(Rnew)
      }
      
      ## remove last column and return the correlation matrix
      Rnew <- cbind(Rnew[1:length(Rnew)-1])
      if (result[1]=="none") return(Rnew)
      else{
         if(result[1]=="html") print(xtable(Rnew), type="html")
         else print(xtable(Rnew), type="latex") 
      }
   }

   # mat : is a matrix of data
   # ... : further arguments to pass to the native R cor.test function
   cor.mtest <- function(mat, ...) {
      mat <- as.matrix(mat)
      n <- ncol(mat)
      p.mat<- matrix(NA, n, n)
      diag(p.mat) <- 0
      for (i in 1:(n - 1)) {
         for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
         }
      }
      colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
      p.mat
   }
   
# RAW (not adjusted)
   RAWdataCorr = RawMaster[,c(10:16, 18:43)] #all response variables except lanthionine
   
   head(RAWdataCorr) #we don't want periods in the column names
   RAWdataCorr = RAWdataCorr[,-31] # remove crude fiber
   
   names(RAWdataCorr) <- gsub(x = names(RAWdataCorr), pattern = "\\.", replacement = " ")  #change column names strategically
   head(RAWdataCorr)
   colnames(RAWdataCorr)[colnames(RAWdataCorr)=="Crude protein"] <- "Crude Protein" #change single column name
   head(RAWdataCorr)
   
   
   
   sum(is.na(RAWdataCorr)) #no NA values
   str(RAWdataCorr) #sampleID present
   

   #In the table above correlations coefficients between the possible pairs of variables are shown.   
   #corstars(x, method = c("pearson", "spearman"), removeTriangle = c("upper",
   #                                                                 "lower"), result = c("none", "html", "latex"), labels_rows,
   #       labels_cols = 1:length(labels_rows) - 1, sig.level = 0.05,
   #      caption = c("Correlation"), filename = "")
   
   
   dataCorr <- RAWdataCorr
   head(dataCorr) #we don't want periods in the column names
   
   names(dataCorr) <- gsub(x = names(dataCorr), pattern = "\\.", replacement = " ")  #change column names strategically
   head(dataCorr)
   colnames(dataCorr)[colnames(dataCorr)=="Crude protein"] <- "Crude Protein" #change single column name
   head(dataCorr)
   
   sum(is.na(dataCorr)) #no NA values
   str(dataCorr) #sampleID present
   
   #table with significance levels
   
   CorrTable <- corstars(dataCorr, method = "spearman")
   
# save as a csv
   write.csv(CorrTable, file = "CorrTableSigRAW.csv")
   dir()
      
   
# correlation matrix for all response variables
   res <- cor(RAWdataCorr, method = "spearman")
   round(res, 2)
   
library(RColorBrewer)
library(corrplot)
# visualize the correlation matrix
   corrplot(res)
   corrplot(res, method = "circle")
   corrplot(res, method = "circle", type = "upper")
   
# changing the order of the corr matrix   
   #"AOE" is for the angular order of the eigenvectors. It is calculated from the order of the angles 
   #"FPC" for the first principal component order  
   #"hclust" for hierarchical clustering order
   #"hclust.method" for the agglomeration method to be used. "hclust.method" should be one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
   #"alphabet" for alphabetical order.   
   corrplot(res, method = "circle", type = "upper", order = "hclust")   # using res2
   
# add p values
# matrix of the p-value of the correlation
   p.mat <- cor.mtest(RAWdataCorr, method = "spearman")
   
# finalize visualizing the correlation table with RAW values
   # open a file !!! do this to save for publication quality image
   #png(height=1200, width=1200, pointsize=12, file="RAWCORR.png")
   
   # create the plot
   finalRAWcorrPLOT <-
      corrplot(res, method = "circle", type = "upper", order = "hclust",
               addrect = 4,
               col = brewer.pal(n = 8, name = "RdBu"),
               diag = FALSE, outline = TRUE,
               p.mat = p.mat, sig.level = 0.5, insig = "blank",
               tl.pos = "td", tl.col = "black", tl.srt = 45, tl.cex = 0.65, 
               cl.pos = "r", cl.lim = c(-1,1), 
               cl.ratio = 0.1, cl.align.text = "l"
               )
   tiff(filename = 'nonAdjCorrPlot180mm.tiff', height = 180, width = 180, units = "mm", bg="white", res = 300)
   #tiff("nonAdjCorrPlot85mm.tiff", height = 85, width = 85, units = "mm", pointsize = 12, compression = "lzw", res = 300)
   # close the file 7.086 = 180mm
   
   dev.off()

   dir()
   

   

#PCA to visualize correlations
   
   install.packages("vegan") #don't compile
   
   adjPCA <-  rda(ADJdataCorr, scale=TRUE)   # run a PCA
   biplot(adjPCA, display = "species") 
   
   rawPCA <-  rda(RAWdataCorr, scale=TRUE)   # run a PCA
   biplot(rawPCA, display = "species") 
   

   #$$$$$ANOVA$$$$$# 
   levels(AdjMaster$LOCATION) #show levels "Chim" "MV"   "Quil" "Sq"  
   
   library(dplyr) #summary statistics by group
   group_by(AdjMaster, LOCATION) %>% #crude protein
     summarise(
       count = n(),
       mean = mean(Crude.protein, na.rm = TRUE),
       sd = sd(Crude.protein, na.rm = TRUE)
     )
   
   group_by(AdjMaster, LOCATION) %>% #total essential amino acid content
     summarise(
       count = n(),
       mean = mean(Total.Essential.AA, na.rm = TRUE),
       sd = sd(Total.Essential.AA, na.rm = TRUE)
     )
   
   group_by(AdjMaster, LOCATION) %>% #lysine
     summarise(
       count = n(),
       mean = mean(Lysine, na.rm = TRUE),
       sd = sd(Lysine, na.rm = TRUE)
     )
   
   group_by(AdjMaster, LOCATION) %>% #Ash #sig diff likely
     summarise(
       count = n(),
       mean = mean(Ash, na.rm = TRUE),
       sd = sd(Ash, na.rm = TRUE)
     )
   
   
#visualize data
   if(!require(devtools)) install.packages("devtools")
   devtools::install_github("kassambara/ggpubr")
   
   library("ggpubr")
   ggboxplot(AdjMaster, x = "LOCATION", y = "Crude.protein", 
             color = "LOCATION", palette = "Blues",
             order = c("Chim", "Quil", "Sq", "MV"),
             ylab = "Crude Protein Content", xlab = "Location")

#complete ANOVA for crude protein
   res.aov <- aov(Crude.protein ~ LOCATION, data = AdjMaster)
   summary(res.aov)
   
   
   
#********************************PDCAAS_Infants********************************#

# clear workspace
   rm(list=ls()) 
   
# start with ADJ dataframe
   #set working directory
   #setwd("~/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #work comp
   #or#
   #setwd("~/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #home comp
   
   setwd("/Users/Evan_Craine/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/gitHub/wetChem")
   dir()
   

# read in data   
   AdjMaster = read.csv("wetChem_adj_FINAL.csv") #AAs are in g/100g crude protein
   
# check column names
   head(AdjMaster) #0 in pop is a commercial variety/check
   
   
# check the structure of the data frame
   str(AdjMaster)
   
# pop needs to be a factor with levels
   AdjMaster$POP <- as.factor(AdjMaster$POP)
   str(AdjMaster) #check to make sure it is fixed
   
# remove everything except AAs with scoring patterns
   df.PDCAAS <- AdjMaster[,c(9, 13, 21:24, 29, 30, 32:34)]
   
# reorder to match scoring pattern table
   df.PDCAAS <- df.PDCAAS[, c(1, 8, 5, 6, 7, 10, 11, 2, 9, 3)]

# calculate the Amino Acid score for each AA by age group
   # AA Score = AA content / scoring pattern ***must use same units
   
   #scoring patterns for 0.5 age group
#     Histidine	Isoleucine	Leucine	Lysine	SAA	 AAA	Threonine	Tryptophan	Valine
#0.5	   2	      3.2	      6.6	   5.7	   2.8	 5.2	3.1	      0.85	      4.3
#1-2     1.8	   3.1      	6.3	   5.2	   2.6	   4.6	2.7	0.74	4.2


# Histidine AA score at 0.5 years old group
   df.PDCAAS$His_Score <- sapply(df.PDCAAS$Histidine, function(X) (X/2))


   # Isoleucine AA score at 0.5 years old group
   df.PDCAAS$Ile_Score <- sapply(df.PDCAAS$Isoleucine, function(X) (X/3.2))

# Leucine AA score at 0.5 years old group
   df.PDCAAS$Leu_Score <- sapply(df.PDCAAS$Leucine, function(X) (X/6.6))
   
# Lysine AA score at 0.5 years old group
   df.PDCAAS$Lys_Score <- sapply(df.PDCAAS$Lysine, function(X) (X/5.7))
   
# SAA AA score at 0.5 years old group
   df.PDCAAS$SAA_Score <- sapply(df.PDCAAS$SAA, function(X) (X/2.8))

# AAA AA score at 0.5 years old group
   df.PDCAAS$AAA_Score <- sapply(df.PDCAAS$AAA, function(X) (X/5.2))
   
# Threonine AA score at 0.5 years old group
   df.PDCAAS$Thr_Score <- sapply(df.PDCAAS$Threonine, function(X) (X/3.1))
   
# Tryptophan AA score at 0.5 years old group
   df.PDCAAS$Trp_Score <- sapply(df.PDCAAS$Tryptophan, function(X) (X/0.85))
   
# Valine AA score at 0.5 years old group
   df.PDCAAS$Val_Score <- sapply(df.PDCAAS$Valine, function(X) (X/4.3))
   
# Create df with AA scores
   df.AA.Scores <- df.PDCAAS[,c(1,11:19)]

# Count the number of samples that have limiting AA content for the 0.5 age group, by AA
   df.AA.Scores$His_Lim <- sapply(df.AA.Scores$His_Score, function(X) (
      if(X>1) {
         print(".")
      } else {
         print("LIMITING")
      }
      ))
   
   df.AA.Scores$Ile_Lim <- sapply(df.AA.Scores$Ile_Score, function(X) (
      if(X>1) {
         print(".")
      } else {
         print("LIMITING")
      }
   ))
   
   df.AA.Scores$Leu_Lim <- sapply(df.AA.Scores$Leu_Score, function(X) (
      if(X>1) {
         print(".")
      } else {
         print("LIMITING")
      }
   ))
   
   df.AA.Scores$Lys_Lim <- sapply(df.AA.Scores$Lys_Score, function(X) (
      if(X>1) {
         print(".")
      } else {
         print("LIMITING")
      }
   ))
   
   df.AA.Scores$SAA_Lim <- sapply(df.AA.Scores$SAA_Score, function(X) (
      if(X>1) {
         print(".")
      } else {
         print("LIMITING")
      }
   ))
   
   df.AA.Scores$AAA_Lim <- sapply(df.AA.Scores$AAA_Score, function(X) (
      if(X>1) {
         print(".")
      } else {
         print("LIMITING")
      }
   ))
   
   df.AA.Scores$Thr_Lim <- sapply(df.AA.Scores$Thr_Score, function(X) (
      if(X>1) {
         print(".")
      } else {
         print("LIMITING")
      }
   ))
   
   df.AA.Scores$Trp_Lim <- sapply(df.AA.Scores$Trp_Score, function(X) (
      if(X>1) {
         print(".")
      } else {
         print("LIMITING")
      }
   ))
   
   df.AA.Scores$Val_Lim <- sapply(df.AA.Scores$Val_Score, function(X) (
      if(X>1) {
         print(".")
      } else {
         print("LIMITING")
      }
   ))
   
   
# Identify limiting AA by sample (i.e. row)
   df.AA.Scores$Limiting_AA <- colnames(df.AA.Scores)[apply(df.AA.Scores,1,which.min)]
   df.AA.Scores$Limiting_AA_Score <- apply(df.AA.Scores,1,which.min)

   #MARGIN = 1 is rows, MARGIN = 2
   
# calculating PDCAAS
   for (i in 1:100) {
      Lim_AA <- df.AA.Scores[i,21]
      Lim_AA_Score <- df.AA.Scores[[i, Lim_AA]] #double bracket necessary; otherwise a list
      PDCAAS <- (Lim_AA_Score*.83) #lowest AA score * digestibility value (subject to change depending on source!!!!)
      df.AA.Scores$PDCAAS[i] <- PDCAAS
      
   }
   
   is.data.frame(df.AA.Scores)
   is.list(df.AA.Scores)
   
   
# save the results
   write.csv(df.AA.Scores, "PDCAAS_infant.csv", row.names = F) 
   
   
# visualize the results
   library(ggplot2)
   #setwd("~/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #workcomputer
   #setwd("~/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #home comp
   
   setwd("/Users/Evan_Craine/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/gitHub/wetChem")
   
   
   Infant_PDCAAS <- read.csv("PDCAAS_infant_ggplot.csv")
   str(Infant_PDCAAS)
   
 
   PDCAAS_infant <- ggplot(
      data = Infant_PDCAAS, aes(
         x = newLocation,
         y = PDCAAS,
         fill = newLocation)) +
      geom_point(
         aes(color=newLocation,
             shape=newLocation),
         size = 2
      ) + 
     ylim(.65,.9) +
      facet_grid(cols = vars(POP))
  
   PDCAAS_infant #yay
   
   p2 <- PDCAAS_infant + labs(title = "PDCAA Score by Population and Location",
              #tag= "Figure XX",
              x="Location",
              #y= "",
              Legend.title = "Location") +
      theme(legend.title = element_blank(),#element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 0.5),
            legend.position = "right",
            
            plot.title = element_text(face = "bold", color = "black", size = 14, angle = 0, hjust = 0.5),
            
            axis.title.y=element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
            
            axis.title.x = element_blank(), #element_text(face = "bold", color = "black", size = 14, angle = 0),
            axis.text.x = element_blank(), #element_text(color = "black", size = 12, angle = 0), #moves x axis title 
            
            axis.ticks.y = element_line(), 
            axis.ticks.x = element_blank())
   
   p2
   ggsave("PDCAAS_infant.png", plot = p2)
   dir()
   
   !!!!!!!!!!!!!!!!
   #switch it up so locations are the column grids
   p <-ggplot(
         Infant_PDCAAS,
         aes(POP,PDCAAS)
      )
   p
   p <- p + geom_point(aes(col=POP, shape=factor(POP))) + facet_grid(cols = vars(newLocation)) + labs(title = "PDCAA Score by Location and Population",
             #tag= "Figure XX",
             x="Population",
             #y= "",
             Legend.title = "Location") +
      theme(legend.title = element_blank(),#element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 0.5),
            legend.position = "right",
            
            plot.title = element_text(face = "bold", color = "black", size = 14, angle = 0, hjust = 0.5),
            
            axis.title.y=element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
            
            axis.title.x = element_blank(), #element_text(face = "bold", color = "black", size = 14, angle = 0),
            axis.text.x = element_blank(), #element_text(color = "black", size = 12, angle = 0), #moves x axis title 
            
            axis.ticks.y = element_line(), 
            axis.ticks.x = element_blank())
   p
                  
   ggsave("PDCAAS_infant_switch.png", plot = p)
   
      
   
   #********************************PDCAAS_PreSchool (1-2)********************************#
   
   # clear workspace
   rm(list=ls()) 

# start with ADJ dataframe
#set working directory
#setwd("~/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #work comp
#or#
#setwd("~/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #home comp
setwd("/Users/Evan_Craine/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/gitHub/wetChem")

dir()

# read in data   
AdjMaster = read.csv("wetChem_adj_FINAL.csv") #AAs are in g/100g crude protein

# check column names
head(AdjMaster) #0 in pop is a commercial variety/check


# check the structure of the data frame
str(AdjMaster)

# pop needs to be a factor with levels
AdjMaster$POP <- as.factor(AdjMaster$POP)
str(AdjMaster) #check to make sure it is fixed

# remove everything except AAs with scoring patterns
df.PDCAAS <- AdjMaster[,c(9, 13, 21:24, 29, 30, 32:34)]

# reorder to match scoring pattern table
df.PDCAAS <- df.PDCAAS[, c(1, 8, 5, 6, 7, 10, 11, 2, 9, 3)]

# calculate the Amino Acid score for each AA by age group
# AA Score = AA content / scoring pattern ***must use same units

#scoring patterns for 1-2 age group
#     Histidine	Isoleucine	Leucine	Lysine	SAA	 AAA	Threonine	Tryptophan	Valine
#1-2	  1.8	        3.1	      6.3	   5.2	   2.6	 4.6	2.7	      0.74	      4.2




# Histidine AA score at 1-2 years old group
df.PDCAAS$His_Score <- sapply(df.PDCAAS$Histidine, function(X) (X/1.8))

# Isoleucine AA score at 1-2 years old group
df.PDCAAS$Ile_Score <- sapply(df.PDCAAS$Isoleucine, function(X) (X/3.1))

# Leucine AA score at 1-2 years old group
df.PDCAAS$Leu_Score <- sapply(df.PDCAAS$Leucine, function(X) (X/6.3))

# Lysine AA score at 1-2 years old group
df.PDCAAS$Lys_Score <- sapply(df.PDCAAS$Lysine, function(X) (X/5.2))

# SAA AA score at 1-2 years old group
df.PDCAAS$SAA_Score <- sapply(df.PDCAAS$SAA, function(X) (X/2.6))

# AAA AA score at 1-2 years old group
df.PDCAAS$AAA_Score <- sapply(df.PDCAAS$AAA, function(X) (X/4.6))

# Threonine AA score at 1-2 years old group
df.PDCAAS$Thr_Score <- sapply(df.PDCAAS$Threonine, function(X) (X/2.7))

# Tryptophan AA score at 1-2 years old group
df.PDCAAS$Trp_Score <- sapply(df.PDCAAS$Tryptophan, function(X) (X/0.74))

# Valine AA score at 1-2 years old group
df.PDCAAS$Val_Score <- sapply(df.PDCAAS$Valine, function(X) (X/4.2))

# Create df with AA scores
df.AA.Scores <- df.PDCAAS[,c(1,11:19)]

# Count the number of samples that have limiting AA content for the 1-2 age group, by AA
df.AA.Scores$His_Lim <- sapply(df.AA.Scores$His_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Ile_Lim <- sapply(df.AA.Scores$Ile_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Leu_Lim <- sapply(df.AA.Scores$Leu_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Lys_Lim <- sapply(df.AA.Scores$Lys_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$SAA_Lim <- sapply(df.AA.Scores$SAA_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$AAA_Lim <- sapply(df.AA.Scores$AAA_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Thr_Lim <- sapply(df.AA.Scores$Thr_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Trp_Lim <- sapply(df.AA.Scores$Trp_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Val_Lim <- sapply(df.AA.Scores$Val_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))


# Identify limiting AA by sample (i.e. row)
df.AA.Scores$Limiting_AA <- colnames(df.AA.Scores)[apply(df.AA.Scores,1,which.min)]
df.AA.Scores$Limiting_AA_Score <- apply(df.AA.Scores,1,which.min)

#MARGIN = 1 is rows, MARGIN = 2

# calculating PDCAAS
for (i in 1:100) {
   Lim_AA <- df.AA.Scores[i,21]
   Lim_AA_Score <- df.AA.Scores[[i, Lim_AA]] #double bracket necessary; otherwise a list
   PDCAAS <- (Lim_AA_Score*.83) #lowest AA score * digestibility value (subject to change depending on source!!!!)
   df.AA.Scores$PDCAAS[i] <- PDCAAS
   
}

is.data.frame(df.AA.Scores)
is.list(df.AA.Scores)


# save the results
write.csv(df.AA.Scores, "PDCAAS_preschool_1_2.csv", row.names = F) 


# visualize the results
library(ggplot2)
setwd("~/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #workcomputer
setwd("~/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #home comp

preschool_PDCAAS <- read.csv("PDCAAS_preschool_1_2_ggplot.csv")
str(preschool_PDCAAS)


PDCAAS_preschool <- ggplot(
   data = preschool_PDCAAS, aes(
      x = newLocation,
      y = PDCAAS,
      fill = newLocation)) +
   geom_point(
      aes(color=newLocation,
          shape=newLocation),
      size = 2
   ) + 
   geom_text(aes(label=entry),hjust=0, vjust=0) + #this adds entry as label next to point
   ylim(.704,.905) +
   facet_grid(cols = vars(POP))

PDCAAS_preschool #yay

p2 <- PDCAAS_preschool + labs(title = "PDCAA Score by Population and Location",
                           #tag= "Figure XX",
                           x="Location",
                           #y= "",
                           Legend.title = "Location") +
   theme(legend.title = element_blank(),#element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 0.5),
         legend.position = "right",
         
         plot.title = element_text(face = "bold", color = "black", size = 14, angle = 0, hjust = 0.5),
         
         axis.title.y=element_text(size = 12, color = "black"),
         axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
         
         axis.title.x = element_blank(), #element_text(face = "bold", color = "black", size = 14, angle = 0),
         axis.text.x = element_blank(), #element_text(color = "black", size = 12, angle = 0), #moves x axis title 
         
         axis.ticks.y = element_line(), 
         axis.ticks.x = element_blank())

p2
ggsave("PDCAAS_preschool_popgrids.png", plot = p2)
dir()

!!!!!!!!!!!!!!!!
   #switch it up so locations are the column grids

# copying above
   P <- ggplot(data = preschool_PDCAAS, aes(x = newLocation,y = PDCAAS,fill = newLocation)) +
   geom_point(aes(color=newLocation,shape=newLocation),size = 2) + 
   ylim(.704,.905) + facet_grid(cols = vars(newLocation))

P #yay

p2 <- P + labs(title = "PDCAA Score by Population and Location",
                              #tag= "Figure XX",
                              x="Location",
                              #y= "",
                              Legend.title = "Location") +
   theme(legend.title = element_blank(),#element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 0.5),
         legend.position = "right",
         
         plot.title = element_text(face = "bold", color = "black", size = 14, angle = 0, hjust = 0.5),
         
         axis.title.y=element_text(size = 12, color = "black"),
         axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
         
         axis.title.x = element_blank(), #element_text(face = "bold", color = "black", size = 14, angle = 0),
         axis.text.x = element_blank(), #element_text(color = "black", size = 12, angle = 0), #moves x axis title 
         
         axis.ticks.y = element_line(), 
         axis.ticks.x = element_blank())

p2
ggsave("PDCAAS_preschool_locgrids.png", plot = p2)
dir()
   
   
   
 #first attempt  
   p <-ggplot(preschool_PDCAAS,
      aes(POP,PDCAAS))
p
p <- p + geom_point(aes(col=POP, shape=factor(POP))) + facet_grid(cols = vars(newLocation)) + labs(title = "PDCAA Score by Location and Population",
                                                                                                   #tag= "Figure XX",
                                                                                                   x="Population",
                                                                                                   #y= "",
                                                                                                   Legend.title = "Location") +
   theme(legend.title = element_blank(),#element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 0.5),
         legend.position = "right",
         
         plot.title = element_text(face = "bold", color = "black", size = 14, angle = 0, hjust = 0.5),
         
         axis.title.y=element_text(size = 12, color = "black"),
         axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
         
         axis.title.x = element_blank(), #element_text(face = "bold", color = "black", size = 14, angle = 0),
         axis.text.x = element_blank(), #element_text(color = "black", size = 12, angle = 0), #moves x axis title 
         
         axis.ticks.y = element_line(), 
         axis.ticks.x = element_blank())
p

ggsave("PDCAAS_preschool_switch.png", plot = p)

# plot PDCAAS with labels for high and lows
library(extrafont)
font_import() # import all your fonts
fonts() #get a list of fonts
fonttable()

PrePD <- ggplot(data = preschool_PDCAAS, aes(x = entry, y = PDCAAS, fill=newLocation))+
   geom_point(aes(color=newLocation,shape=newLocation),size = 2) + 
   geom_text(aes(label=ifelse(PDCAAS>0.862,as.character(entry),'')), hjust=0.5,vjust=0)+
   geom_text(aes(label=ifelse(PDCAAS<0.78,as.character(entry),'')), hjust=0.5,vjust=0)+
theme(legend.title = element_blank(),#element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 0.5),
      legend.position = "right",
      
      plot.title = element_text(face = "bold", color = "black", size = 14, angle = 0, hjust = 0.5),
      
      axis.title.y=element_text(size = 12, color = "black", family = ),
      axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
      
      axis.title.x = element_blank(), #element_text(face = "bold", color = "black", size = 14, angle = 0),
      axis.text.x = element_blank(), #element_text(color = "black", size = 12, angle = 0), #moves x axis title 
      
      axis.ticks.y = element_line(), 
      axis.ticks.x = element_blank())
PrePD
ggplot(nba, aes(x= MIN, y= PTS, colour="green", label=Name))+
   geom_point() +
   geom_text(aes(label=ifelse(PTS>24,as.character(Name),'')),hjust=.5,vjust=0)
!!!!!!!!!!!!!!!!
   
   #********************************PDCAAS_Children (3-10)********************************#
   
   # clear workspace
   rm(list=ls()) 

# start with ADJ dataframe
#set working directory
#setwd("~/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #work comp
#or#
#setwd("~/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #home comp

setwd("/Users/Evan_Craine/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/gitHub/wetChem")
dir()


# read in data   
AdjMaster = read.csv("wetChem_adj_FINAL.csv") #AAs are in g/100g crude protein

# check column names
head(AdjMaster) #0 in pop is a commercial variety/check


# check the structure of the data frame
str(AdjMaster)

# pop needs to be a factor with levels
AdjMaster$POP <- as.factor(AdjMaster$POP)
str(AdjMaster) #check to make sure it is fixed

# remove everything except AAs with scoring patterns
df.PDCAAS <- AdjMaster[,c(9, 13, 21:24, 29, 30, 32:34)]

# reorder to match scoring pattern table
df.PDCAAS <- df.PDCAAS[, c(1, 8, 5, 6, 7, 10, 11, 2, 9, 3)]

# calculate the Amino Acid score for each AA by age group
# AA Score = AA content / scoring pattern ***must use same units

#scoring patterns for 1-2 age group
#     Histidine	Isoleucine	Leucine	Lysine	SAA	 AAA	Threonine	Tryptophan	Valine
#3-10   1.6	      3.1	      6.1	   4.8	   2.4	 4.1	2.5	      0.66	      4



# Histidine AA score at 3-10 years old group
df.PDCAAS$His_Score <- sapply(df.PDCAAS$Histidine, function(X) (X/1.6))

# Isoleucine AA score at 3-10 years old group
df.PDCAAS$Ile_Score <- sapply(df.PDCAAS$Isoleucine, function(X) (X/3.1))

# Leucine AA score at 3-10 years old group
df.PDCAAS$Leu_Score <- sapply(df.PDCAAS$Leucine, function(X) (X/6.1))

# Lysine AA score at 3-10 years old group
df.PDCAAS$Lys_Score <- sapply(df.PDCAAS$Lysine, function(X) (X/4.8))

# SAA AA score at 3-10 years old group
df.PDCAAS$SAA_Score <- sapply(df.PDCAAS$SAA, function(X) (X/2.4))

# AAA AA score at 3-10 years old group
df.PDCAAS$AAA_Score <- sapply(df.PDCAAS$AAA, function(X) (X/4.1))

# Threonine AA score at 3-10 years old group
df.PDCAAS$Thr_Score <- sapply(df.PDCAAS$Threonine, function(X) (X/2.5))

# Tryptophan AA score at 3-10 years old group
df.PDCAAS$Trp_Score <- sapply(df.PDCAAS$Tryptophan, function(X) (X/0.66))

# Valine AA score at 3-10 years old group
df.PDCAAS$Val_Score <- sapply(df.PDCAAS$Valine, function(X) (X/4))

# Create df with AA scores
df.AA.Scores <- df.PDCAAS[,c(1,11:19)]

# Count the number of samples that have limiting AA content for the 3-10 age group, by AA
df.AA.Scores$His_Lim <- sapply(df.AA.Scores$His_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Ile_Lim <- sapply(df.AA.Scores$Ile_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Leu_Lim <- sapply(df.AA.Scores$Leu_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Lys_Lim <- sapply(df.AA.Scores$Lys_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$SAA_Lim <- sapply(df.AA.Scores$SAA_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$AAA_Lim <- sapply(df.AA.Scores$AAA_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Thr_Lim <- sapply(df.AA.Scores$Thr_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Trp_Lim <- sapply(df.AA.Scores$Trp_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))

df.AA.Scores$Val_Lim <- sapply(df.AA.Scores$Val_Score, function(X) (
   if(X>1) {
      print(".")
   } else {
      print("LIMITING")
   }
))


# Identify limiting AA by sample (i.e. row)
df.AA.Scores$Limiting_AA <- colnames(df.AA.Scores)[apply(df.AA.Scores,1,which.min)]
df.AA.Scores$Limiting_AA_Score <- apply(df.AA.Scores,1,which.min)

#MARGIN = 1 is rows, MARGIN = 2

# calculating PDCAAS
for (i in 1:100) {
   Lim_AA <- df.AA.Scores[i,21]
   Lim_AA_Score <- df.AA.Scores[[i, Lim_AA]] #double bracket necessary; otherwise a list
   PDCAAS <- (Lim_AA_Score*.843) #lowest AA score * digestibility value (subject to change depending on source!!!!) (using Ranhotra et al. 1993)
   df.AA.Scores$PDCAAS[i] <- PDCAAS
   
}

is.data.frame(df.AA.Scores)
is.list(df.AA.Scores)


# save the results
write.csv(df.AA.Scores, "PDCAAS_children_3_10.csv", row.names = F) 



!!!!!!!!!!!!!!!!
      
      
# added unadjusted_wetchem_analysis.R
      # this script has code for the Principal components and windroses
      # clear workspace
      rm(list=ls()) 
   
   # set working directory
   #setwd("~/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #work comp
   #or#
   #setwd("~/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/QUINOA/wet_chem/Manuscript/ANALYSIS/R") #home comp

   
   setwd("/Users/Evan_Craine/evancraine@gmail_DROPBOX/Dropbox/evan.craine@wsu.edu/gitHub/wetChem")
   
   dir()
   
   # prepare packages
   library(agricolae)
   #install.packages("dplyr")
   library(dplyr)
   #install.packages("pastecs")
   library(pastecs)
   #install.packages(ggplot2)
   library(ggplot2)
   library(qwraps2)
   library(data.table)
   #install.packages("tables")
   library(tables) #need tables package for tabular function
   library(viridis)
   
   # color palettes for the color challenged viewer
   # The palette with grey:
   cbPaletteGray <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
   
   # The palette with black:
   cbbPalettBlack <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
   
   # To use for fills, add
   scale_fill_manual(values=cbPaletteBlack)
   
   # To use for line and point colors, add
   scale_colour_manual(values=cbPaletteBlack)
   
   # read in data
   master = read.csv("wetChem_OGunits_FINAL.csv")
   #master = read.csv("R_WSU_quinoa_wetchem_NIR_stats.csv")
   #import data set
   #wetChem_adj_FINAL$POP <- as.factor(wetChem_adj_FINAL$POP)
   #wetChem_adj_FINAL$LOCATION <- as.factor(wetChem_adj_FINAL$LOCATION)
   
   
   # check column names
   head(master) #0 in pop is a commercial variety/check
   
   #check the structure of the data frame
   str(master)
   
   # pop needs to be a factor with levels
   master$POP <- as.factor(master$POP)
   master$LOCATION <- as.factor(master$LOCATION)
   
   str(master) #check to make sure it is fixed
   
   
   
   #$$$___descriptive statistics___$$$#
   # summary stats by groups for proximates 
   
   # proximates
   
   #tabular operates row ~ columns; row + 1 gives you all groups at bottom; n = 1 lists the the number of samples per group in the table
   DescAll <- tabular( (LOCATION + POP + 1) ~ (n=1) + Format(digits=2) * (Total.AA + Crude.Protein. + Moisture + Crude.Fat + Crude.Fiber + Ash + Crude.Carbohydrates) * (mean + sd), data=master) #these are the columns 
   DescAll
   
   DescDets <- tabular( (1) ~ (n=1) + Format(digits=2) * (Total.AA + Crude.Protein. + Moisture + Crude.Fat + Crude.Fiber + Ash + Crude.Carbohydrates) * (min + max + mean + sd), data=master) #these are the columns 
   DescDets
   #save as .csv  
   write.csv.tabular(DescAll, file = "nonadJ_Descriptives.csv", justification = "n", row.names=FALSE)
   write.csv.tabular(DescDets, file = "nonadj_DetailedDescriptives_All.csv", justification = "n", row.names=FALSE)
   
   #$$$___descriptive statistics___$$$#
   #want descriptives for the AA profile
   DescAA <- tabular( (1) ~ Format(digits=2) * (Total.AA + Total.Essential.AA + Leucine + Lysine + Valine + Isoleucine + Phenylalanine + Threonine + Histidine + Methionine + Tryptophan + Total.Nonessential.AA + Glutamic.Acid + Aspartic.Acid + Arginine + Glycine + Alanine + Proline + Serine + Tyrosine + Cysteine + Taurine + Hydroxyproline + Hydroxylysine + Ornithine) * (mean + sd + min + max), data=master) #these are the columns 
   DescAA
   
   write.csv.tabular(DescAA, file = "nonadj_AAprofile_Descriptives.csv", justification = "n", row.names=FALSE)
   
   #want descriptives for the essential AA profile by 'groups'
   
   DescAAdeets <- tabular( (LOCATION + POP + 1) ~ (n=1) + Format(digits=2) * (Total.Essential.AA + Leucine +	Lysine + Valine + Isoleucine + Phenylalanine + Threonine + Histidine + Methionine + Tryptophan) * (mean + sd), data=master) #these are the columns 
   DescAAdeets
   
   
   
   write.csv.tabular(DescAAdeets, file = "nonadj_EssentialAA_Descriptives.csv", justification = "n", row.names=FALSE)
   
   
   
   #$$$___ranking samples for nutritional composition___$$$#
   #read csv for showdownOG
   
   showdown <- read.csv("showdown.csv")
   
   str(showdown)
   # keep pop, newEntry, newLocation, ID through the last parameter, remove lanthionine by skipping
   showdown <- showdown[,c(8, 44, 46, 9:11, 13:43)] #46 variables to 33...
   head(showdown)
   str(showdown)
   
   showdown$POP <- as.factor(showdown$POP)
   
   str(showdown)
   # create a vector of the colnames (i.e. nutritional attribute names)  
   attributes <- as.vector(colnames(showdown[,5:37])) #this returns df column names
   attributes
   
   # create new series of coloumns for attribute ranks
   
   newcolumn <- paste0(attributes, "_Rank") #adds rank to each
   showdown[,newcolumn] <- NA #addes empty columns to data frame using names from newcolumn vector
   showdown[,"Rank"] <- NA #creates the last column to sum the ranks in
   head(showdown)
   str(showdown) #doubles the number of columns (66)...
   
   #test meat of for loop
   for (i in 38:71) { #start at first parameter rank columnn number and go to last parameter rank column number
      showdown[,38]<-rank(showdown[,(38-33)], ties.method = "average")} #rank the first attribute column
   #put into the first attribute rank column the rank of that attribute's raw data
   #using subtraction to start at 1 and end at 33
   #pass loop meat into first attribute column number to last attribute column number - 1!!!
   #loop will end by passing rank for last attribute into last attribute rank
   
   #higher rank higher value...
   #samples with identical values share rank values for that attribute
   
   
   #for loop construction and execution
   for (i in 38:71) { #start at first parameter rank columnn and go to last parameter rank column
      showdown[,i]<-rank(showdown[,(i-33)], ties.method = "average") 
      
   }
   
   #now each rank column has the sample with the lowest content as 1, and the sample with the highest content as 100 (or the number that represents the tie)
   
   #now we need to sum the ranks and put this in the total rank column 
   showdown[,38:70]
   showdown$Rank <- rowSums( showdown[,38:70] ) #start at first rank column and end at last rank column 
   
   head(showdown)
   
   #large number good for total rank (lots of high ranks and higher content) and small number not as good...#
   
   # save this dataframe as a csv
   write.csv(showdown, file = "nonADJ_Nutritional_Composition_Rankings_FINAL.csv")
   dir()
   
   # take a look at entries ranked
   str(showdown)
   simpleRank = showdown[,c(1:4, 71)] #pull out the identifier columns from showdown & total rank column
   head(simpleRank)
   
   
   # sort by Totl Rank
   simpleRank = simpleRank[order(simpleRank$Rank, decreasing = TRUE),]
   str(simpleRank)
   
   # save ordered df as csv
   write.csv(simpleRank, file = "nonADJ_SimpleRank.csv")
   # number from highest rank to lowest 100 to 1
   # higher rank is better!
   
   
   #!!!!!!!!!!!!!!!!go and add a rank column where higher total rank = lower rank number - who's #1!!!?!?!??!
   
   # top 5
   top5 = simpleRank[1:5,]
   top5
   
   write.csv(top5, file = "nonADJ_Top5.csv")
   
   # bottom 5
   bottom5 = simpleRank[96:100,]
   bottom5
   
   write.csv(top5, file = "nonADJ_Bottom5.csv")
   
   
   # visualize ranks
   # info on prepping data for ggplot here: https://stackoverflow.com/questions/15706281/controlling-order-of-points-in-ggplot2-in-r
   # good scatter plot resource here: https://www.statmethods.net/graphs/scatterplot.html
   
   # read in *~*~new*~*~* showdown csv file as a df
   dir() #want to read in the file where total rank has been convereted to 1 to 100 - "nonADJ_SimpleRank_ordered.csv" 
   ranks <- read.csv("nonADJ_SimpleRank_ordered.csv")
   str(ranks)
   
   # pop to factor
   ranks$POP <- as.factor(ranks$POP)
   ranks$Rank <- as.numeric(ranks$Rank)
   str(ranks)
   
   ggplot() +
      geom_point(
         data = ranks,
         aes(x=ID,y=Rank),
         color = "black",
         size = 2,
         shape = 23) #first go at it
   
   rank <- ggplot(ranks, aes(x=newLocation, y=Rank, col = POP))
   rank + geom_jitter() #scree/jitter pop colored by pops (only location on x axis)
   
   rank + geom_point(position = "jitter")
   
   rank + geom_text(aes(label=ID), size=4) #practicing with overlaying categorical data
   
   
   library(ggplot2)
   
   #main grouping category is population, second is location
   #example
   #1
   ggplot(data=df, aes(x=second.cat, y=value)) + geom_point() + facet_grid(~ main.cat, scales = 'free')
   
   
   g <- ggplot(df, aes(x = x, y = y, group = group))
   g <- g + geom_line(aes(colour = group))
   g <- g + geom_point(aes(colour = group), alpha = 0.8)
   
   
   NR <- ggplot(
      data = ranks,
      aes(
         x = newLocation,
         y = Rank,
         group = newLocation,
         color = newLocation,
         shape = newLocation,
         size = 2
      )) +
      geom_point() +
      scale_y_reverse(lim = c(100,1))+ #reverse the scale so the #1 ranked sample is on top....
      facet_grid(~ POP, scales = 'free') 
   NR #BOOM WINNER WINNER CHICKEN DINNER
   
   NR2 <- ggplot(
      data = ranks,
      aes(
         x = newLocation,
         y = Rank,
         fill = newLocation)) +
      geom_point(
         aes(color=newLocation,
             shape=newLocation),
         size = 2
      ) + 
      scale_y_reverse(lim = c(100,1)) +
      facet_grid(~ POP, scales = 'free') 
   NR2 #BOOM WINNER WINNER CHICKEN DINNER
   
   NR2 + labs(title = "Rank by Population and Location",
              #tag= "Figure XX",
              x="Location",
              #y= "",
              Legend.title = "Location") +
      theme(legend.title = element_blank(),#element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 0.5),
            legend.position = "right",
            
            plot.title = element_text(face = "bold", color = "black", size = 14, angle = 0, hjust = 0.5),
            
            axis.title.y=element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
            
            axis.title.x = element_blank(), #element_text(face = "bold", color = "black", size = 14, angle = 0),
            axis.text.x = element_blank(), #element_text(color = "black", size = 12, angle = 0), #moves x axis title 
            
            axis.ticks.y = element_line(), 
            axis.ticks.x = element_blank())
   
   ggsave("Rank.png", plot = NR2)
   
   
   # look at amino acid cutoffs using adjusted data set
   #I KNOW WE"RE WORKING WITH ADJUSTED DATA IN THE UNADJ SCRIPT
   #IT"S WILD AND WONDERFUL AND WE"RE DOING IT!!!
   ADJ_master <- read.csv("wetChem_adj_FINAL.csv")   
   
   str(ADJ_master) # 'note' captures only irrigated or nonirrigated values
   ADJ_master$POP <- as.factor(ADJ_master$POP) # POP as factor
   ADJ_master <- ADJ_master[,-17] #remove Lanthionine
   ADJ_master
   
   Leucine <- ggplot(
      data = ADJ_master,
      aes(
         x = newLocation,
         y = Leucine,
         fill = newLocation,
         color = newLocation,
         shape = newLocation,
         size = 
      )) +
      geom_point() +
      facet_grid(~ POP, scales = 'free') 
   
   Leucine #add in labels
   
   ggsave("plot3.png", plot=pl, width=14, height=7, units="cm", dpi=1200 )
   
   
   Leucine + 
      labs(title = "Leucine Content by Population and Location",
           tag= "Figure XX",
           x="Location",
           y= "(g/100g crude protein)",
           Legend = "Location") +
      theme(legend.title = element_blank(),#element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 0.5),
            legend.position = "right",
            plot.title = element_text(face = "bold", color = "black", size = 14, angle = 0, hjust = 0.5),
            axis.title.y=element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
            axis.title.x = element_blank(), #element_text(face = "bold", color = "black", size = 14, angle = 0),
            axis.text.x = element_blank(), #element_text(color = "black", size = 12, angle = 0), #moves x axis title 
            axis.ticks.y = element_line(), 
            axis.ticks.x = element_blank())
   
   CrudeProtein <- ggplot(master, aes(x = newLocation, y = Crude.Protein.))
   # rewrite for trait of interest

   TEAA <- ggplot(
      data = ADJ_master,
      aes(
         x = newLocation,
         y = Total.Essential.AA,
         group = newLocation,
         color = newLocation,
         shape = newLocation
         #size = 1
      )) +
      geom_point() +
      scale_color_grey() + #remove this to get back to colors
      facet_grid(~ POP, scales = 'free') + 
      geom_text(label = ADJ_master$newEntry) #remove this to get rid of text
   
   TEAA 

   
   TEAA + labs(title = "Essential Amino Acid Content by Population and Location",
              #tag= "Figure XX",
              x="Location",
              y= "Content (g/100g protein)",
              Legend.title = "Location") +
      theme(legend.title = element_blank(),#element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 0.5),
            legend.position = "right",
            
            plot.title = element_text(face = "bold", color = "black", size = 14, angle = 0, hjust = 0.5),
            
            axis.title.y=element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
            
            axis.title.x = element_blank(), #element_text(face = "bold", color = "black", size = 14, angle = 0),
            axis.text.x = element_blank(), #element_text(color = "black", size = 12, angle = 0), #moves x axis title 
            
            axis.ticks.y = element_line(), 
            axis.ticks.x = element_blank())
   
   
   
   ggsave("Total.EAA.png", plot = TEAA)
   
# create scatter with labels as genotypes
   TEAA.scatter <- ggplot(data = ADJ_master, aes(y=Total.Essential.AA, x= newLocation)) 
   TEAA.scatter +
      geom_point() +
      geom_text(label=ADJ_master$newEntry)
   
   # visualizing variation within locations and populations
   
   
   #visual representation of populations in study
   hist(master$POP,
        main = "Populations",
        xlab = "Population ID",
        border = "black",
        col = "blue")
   
   
   #seed composition
   hist(master$Crude.Protein.,
        main = "Crude Protein Content",
        xlab = "Protein (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Crude.Fat,
        main = "Crude Fat Content",
        xlab = "Fat (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Crude.Fiber, 
        main = "Crude Fiber Content",
        xlab = "Fiber (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Moisture, 
        main = "Moisture Content (as is)",
        xlab = "Moisture (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Ash,
        main = "Ash Content",
        xlab = "Ash (g/100g)",
        border = "black",
        col = "blue")
   
   #visual representation of amino acids
   hist(master$Total.AA,
        main = "Total Amino Acid Content",
        xlab = "Total Amino Acid (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Taurine,
        main = "Taurine Content",
        xlab = "Taurine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Hydroxyproline,
        main = "Hydroxyproline Content",
        xlab = "Hydroxyproline (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Aspartic.Acid,
        main = "Aspartic Acid Content",
        xlab = "Aspartic Acid (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Threonine,
        main = "Threonine Content",
        xlab = "Threonine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Serine,
        main = "Serine Content",
        xlab = "Serine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Glutamic.Acid,
        main = "Glutamic Acid Content",
        xlab = "Glutamic Acid (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Proline,
        main = "Proline Content",
        xlab = "Proline (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Lanthionine._,  ##### all zeroes!!!!#####
        main = "Lanthionine Content",
        xlab = "Lanthionine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Glycine,
        main = "Glycine Content",
        xlab = "Glycine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Alanine,
        main = "Alanine Content",
        xlab = "Alanine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Cysteine,
        main = "Cysteine Content",
        xlab = "Cysteine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Valine,
        main = "Valine Content",
        xlab = "Valine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Methionine,
        main = "Methionine Content",
        xlab = "Methionine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Isoleucine,
        main = "Isoleucine Content",
        xlab = "Isoleucine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Leucine,
        main = "Leucine Content",
        xlab = "Leucine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Tyrosine,
        main = "Tyrosine Content",
        xlab = "Tyrosine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Phenylalanine,
        main = "Phenylalanine Content",
        xlab = "Phenylalanine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Hydroxylysine,
        main = "Hydroxylysine Content",
        xlab = "Hydroxylysine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Ornithine._,
        main = "Ornithine Content",
        xlab = "Ornithine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Lysine,
        main = "Lysine Content",
        xlab = "Lysine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Histidine,
        main = "Histidine Content",
        xlab = "Histidine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Arginine,
        main = "Arginine Content",
        xlab = "Arginine (g/100g)",
        border = "black",
        col = "blue")
   
   hist(master$Tryptophan,
        main = "Tryptophan Content",
        xlab = "Tryptophan (g/100g)",
        border = "black",
        col = "blue")
   
   #********************************linear relationships********************************#
   #$$$$$test for normality$$$$$# 
   #examples; HO=data are normally distributed; HA=data are not normally distributed; if p value is greater than 0.05, we accept H0; if p value is less than 0.05 we reject H0 and data are not normally distributed
   # Shapiro-Wilk normality test for mpg
   shapiro.test(my_data$mpg) # => p = 0.1229
   # Shapiro-Wilk normality test for wt
   shapiro.test(my_data$wt) # => p = 0.0
   #create an empty data frame
   NormalTest <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("parameter", "statistic", "p value", "normal?"))
   #shorten master to get just the attributes we're interested in
   MasterMinus <- master[,c(10:16, 18:43)] #just the nutritional attributes
   #create a vector of the attributes   
   factors <- colnames(MasterMinus) #33 attributes we're interested in 
   #create a for loop to test all attributes for normality and report out the results   
   for (i in 1:33) {
      Ntest <- shapiro.test(MasterMinus[,i])
      NormalTest[(i), 1] <- factors[i]
      NormalTest[(i), 2] <- Ntest$statistic
      NormalTest[(i), 3] <- Ntest$p.value
      NormalTest[(i), 4] <- if (Ntest$p.value > 0.05){
         print('Yes')
      } else {
         print('No')
      }
   }  
   
   #$$$$$inspect QQ plots$$$$$# 
   library("ggpubr")
   
   ggqqplot(master$Taurine, ylab = "Taurine") #normal
   ggqqplot(master$Methionine, ylab = "Methionine") #nonnormal, but pretty close
   ggqqplot(master$Hydroxylysine, ylab = "Hydroxylysine") #nonnormal
   ggqqplot(master$Ornithine, ylab = "Ornithine") #nonnormal
   ggqqplot(master$Arginine, ylab = "Arginine") #nonnormal
   ggqqplot(master$Moisture, ylab = "Moisture") #nonnormal
   ggqqplot(master$Crude.Fiber, ylab = "Crude Fiber") #nonnormal
   ggqqplot(master$Ash, ylab = "Ash") #nonnormal
   
   #$$$$$CORRELATIONS$$$$$# 
   
   #load the required packages/functions
   library(Hmisc) #rcorr function is present in this package
   install.packages("corrplot")
   library(corrplot)
   library(PerformanceAnalytics)
   
   #This section provides a simple function for formatting a correlation matrix into a table with 4 columns containing :
   
   #Column 1 : row names (variable 1 for the correlation test)
   #Column 2 : column names (variable 2 for the correlation test)
   #Column 3 : the correlation coefficients
   #Column 4 : the p-values of the correlations
   #The custom function below can be used :
   # ++++++++++++++++++++++++++++
   # flattenCorrMatrix
   # ++++++++++++++++++++++++++++
   # cormat : matrix of the correlation coefficients
   # pmat : matrix of the correlation p-values
   flattenCorrMatrix <- function(cormat, pmat) {
      ut <- upper.tri(cormat)
      data.frame(
         row = rownames(cormat)[row(cormat)[ut]],
         column = rownames(cormat)[col(cormat)[ut]],
         cor  =(cormat)[ut],
         p = pmat[ut]
      )
   }
   
   
   #prepare a simplified data frame with ID and the response variables
   #ALL nutritional components  
   dataCorr = master[,c(10:16, 18:43)] #all response variables except lanthionine
   head(dataCorr) #we don't want periods in the column names
   
   names(dataCorr) <- gsub(x = names(dataCorr), pattern = "\\.", replacement = " ")  #change column names strategically
   head(dataCorr)
   colnames(dataCorr)[colnames(dataCorr)=="Crude protein"] <- "Crude Protein" #change single column name
   head(dataCorr)
   
   
   
   sum(is.na(dataCorr)) #no NA values
   str(dataCorr) #sampleID present
   
   #correlation matrix for all response variables
   res <- cor(dataCorr)
   round(res, 2)
   #In the table above correlations coefficients between the possible pairs of variables are shown.
   
   
   #correlation matrix with significance levels (p-value)
   res2 <- rcorr(as.matrix(dataCorr)) #rcorr function provides a more sophisticated output
   res2
   #The output of the function rcorr() is a list containing the following elements : 
   #- r : the correlation matrix 
   #- n : the matrix of the number of observations used in analyzing each pair of variables 
   #- P : the p-values corresponding to the significance levels of correlations.
   # Extract the correlation coefficients
   res2$r
   # Extract p-values
   res2$P
   
   corrAll <- flattenCorrMatrix(res2$r, res2$P)
   
   write.csv(corrAll, "All Correlations.csv")
   
   #visualize correlation matrix
   corrplot(res2$r, type="upper", order="hclust", 
            p.mat = res2$P, sig.level = 0.05, insig = "blank") #insiginificant correlations are left blank
   
   
   col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
   corrplot(res2$r, method = "color", col = col(200),  
            type = "upper", order = "hclust", 
            addCoef.col = "black", # Add coefficient of correlation
            tl.col = "darkblue", tl.srt = 45, #Text label color and rotation
            # Combine with significance level
            p.mat = res2$P, sig.level = 0.05,  
            # hide correlation coefficient on the principal diagonal
            diag = FALSE 
   )
   
   corrplot(res2$r, p.mat = res2$p, method = "color",
            insig = "label_sig", pch.col = "white")
   
   
   #using ggcorrplot to create correlation matrix
   library(ggcorrplot)
   
   #ggcorrplot(corr, 
   #          method = c("square", "circle"), 
   #         type = c("full", "lower","upper"), 
   #        ggtheme = ggplot2::theme_minimal, title = "",
   #       show.legend = TRUE, legend.title = "Corr", show.diag = FALSE,
   #      colors = c("blue", "white", "red"), outline.color = "gray",
   #     hc.order = FALSE, hc.method = "complete", 
   #    lab = FALSE, lab_col = "black", lab_size = 4, p.mat = NULL, sig.level = 0.05,
   #   insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 5,
   #  tl.cex = 12, tl.col = "black", tl.srt = 45, digits = 2)
   
   #$$$$$$$$ first approach $$$$$$$$#
   
   corr <- round(cor(dataCorr), 1)
   head(corr)
   
   #visualize
   ggcorrplot(corr)
   
   # Compute a matrix of correlation p-values
   p.mat <- cor_pmat(dataCorr)
   head(p.mat[, 1:4])
   
   
   #basic
   ggcorrplot(corr, # Leave blank on no significant coefficient
              p.mat = p.mat, hc.order = TRUE,
              type = "lower", insig = "blank")
   
   
   #customized correlation matrix
   library(colorspace)
   
   colors <- rainbow_hcl(3)
   colors <- c("red", "yellow", "blue")
   
   ggcorrplot(corr, 
              method = "square", 
              type = "lower", 
              ggtheme = ggplot2::theme_minimal, title = " Correlations Between Nutritional Components",
              show.legend = TRUE, legend.title = "Pearson Correlation Coefficient", show.diag = TRUE,
              colors = "colors", outline.color = "gray",
              hc.order = TRUE, hc.method = "complete", 
              lab = TRUE, lab_col = "black", lab_size = 1, p.mat = p.mat, sig.level = 0.05,
              insig = "pch", pch.col = "black", pch.cex = 2,
              tl.cex = 4, tl.col = "black", tl.srt = 45, digits = 2)
   
   #final$$$$$$$$$$$$$$
   ggcorrplot(corr,
              method = "square",
              type = "lower",
              legend.title = "Correlation Coefficients",
              colors = c("red", "white", "blue"),
              p.mat = p.mat, hc.order = TRUE, hc.method = "complete",
              lab=TRUE, lab_col = "black", lab_size = 3, digits = 2, #size and color for the correlation coefficient labels
              insig = "pch", pch.col = "black", pch.cex = 5, #glyphs for insiginificant
              tl.cex =10, tl.col = "black",tl.srt=70) #the size, the color and the string rotation of text label (variable names).
   
   
   
   #$$$$$$$$ second approach $$$$$$$$#
   
   cormat <- round(cor(dataCorr),2)
   
   
   head(cormat)
   
   #The package reshape is required to melt the correlation matrix
   library(reshape2)
   melted_cormat <- melt(cormat)
   head(melted_cormat)
   
   #visualize
   library(ggplot2)
   ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile()
   
   # Get lower triangle of the correlation matrix
   get_lower_tri<-function(cormat){
      cormat[upper.tri(cormat)] <- NA
      return(cormat)
   }
   # Get upper triangle of the correlation matrix
   get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat)]<- NA
      return(cormat)
   }
   #usage
   upper_tri <- get_upper_tri(cormat)
   upper_tri
   
   # Melt the correlation matrix
   library(reshape2)
   melted_cormat <- melt(upper_tri, na.rm = TRUE)
   # Heatmap
   library(ggplot2)
   ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1,1), space = "Lab", 
                           name="Pearson\nCorrelation") +
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))+
      coord_fixed()   
   
   #reorder correlation matrix to reveal hidden pattern in the data
   reorder_cormat <- function(cormat){
      # Use correlation between variables as distance
      dd <- as.dist((1-cormat)/2)
      hc <- hclust(dd)
      cormat <-cormat[hc$order, hc$order]
   }
   
   
   #using a different package!
   library(psych)
   if(!require(PerformanceAnalytics)) install.packages("PerformanceAnalytics")
   library(PerformanceAnalytics)
   
   pairs.panels(dataCorr)
   chart.Correlation(dataCorr)
   
   
   #create basic correlation table
   #correlation matrix analysis
   mcor<-round(cor(dataCorr),2)
   mcor #return correlation coefficients for all comparison
   
   #hide upper or lower part of the matrix
   lower<-mcor
   lower[lower.tri(mcor, diag=TRUE)]<-""
   lower<-as.data.frame(lower)
   lower
   
   #nice correlation table in HTML format
   library(xtable)
   library(htmlTable)
   print(xtable(lower), type="html")
   
   #significance levels
   #first run corstars function
   # x is a matrix containing the data
   # method : correlation method. "pearson"" or "spearman"" is supported
   # removeTriangle : remove upper or lower triangle
   # results :  if "html" or "latex"
   # the results will be displayed in html or latex format
   corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                       result=c("none", "html", "latex")){
      #Compute correlation matrix
      require(Hmisc)
      x <- as.matrix(x)
      correlation_matrix<-rcorr(x, type=method[1])
      R <- correlation_matrix$r # Matrix of correlation coeficients
      p <- correlation_matrix$P # Matrix of p-value 
      
      ## Define notions for significance levels; spacing is important.
      mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
      
      ## trunctuate the correlation matrix to two decimal
      R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
      
      ## build a new matrix that includes the correlations with their apropriate stars
      Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
      diag(Rnew) <- paste(diag(R), " ", sep="")
      rownames(Rnew) <- colnames(x)
      colnames(Rnew) <- paste(colnames(x), "", sep="")
      
      ## remove upper triangle of correlation matrix
      if(removeTriangle[1]=="upper"){
         Rnew <- as.matrix(Rnew)
         Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
         Rnew <- as.data.frame(Rnew)
      }
      
      ## remove lower triangle of correlation matrix
      else if(removeTriangle[1]=="lower"){
         Rnew <- as.matrix(Rnew)
         Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
         Rnew <- as.data.frame(Rnew)
      }
      
      ## remove last column and return the correlation matrix
      Rnew <- cbind(Rnew[1:length(Rnew)-1])
      if (result[1]=="none") return(Rnew)
      else{
         if(result[1]=="html") print(xtable(Rnew), type="html")
         else print(xtable(Rnew), type="latex") 
      }
   }
   
   #table with significance levels
   #corstars(x, method = c("pearson", "spearman"), removeTriangle = c("upper",
   #                                                                 "lower"), result = c("none", "html", "latex"), labels_rows,
   #       labels_cols = 1:length(labels_rows) - 1, sig.level = 0.05,
   #      caption = c("Correlation"), filename = "")
   CorrTable <- corstars(dataCorr, method = "pearson")
   
   #save as a csv
   write.csv(CorrTable, file = "CorrTableSig.csv")
   dir()
   
   
   
#$$$$$$$$$$$$$$$$$_PCA_$$$$$$$$$$$$$$$$$$$$$$$$# #moved everything from here down over from unadjusted_wetchem_analysis.R on 9.12.19
   
# PCA to visualize correlations
   # read in data
   master = read.csv("wetChem_OGunits_FINAL.csv") # changed 0 to controls and 110 to Mount Vernon 9.11.19
   dataCorr = master[,c(10:16, 18:43)] #all response variables except lanthionine #master is OG chem units
   head(dataCorr) #we don't want periods in the column names
   
   names(dataCorr) <- gsub(x = names(dataCorr), pattern = "\\.", replacement = " ")  #change column names strategically
   head(dataCorr)
   colnames(dataCorr)[colnames(dataCorr)=="Crude protein"] <- "Crude Protein" #change single column name
   head(dataCorr)
   
# first clean up column names
   names(dataCorr) <- gsub(x = names(dataCorr), pattern = "\\.", replacement = " ")  #change column names strategically
   head(dataCorr)
   str(dataCorr)

# prepare DF to break up attributes 
   dataCorr #all attributes   
   allAA <- dataCorr[,1:27]
   noTotalsAA <- dataCorr[,1:22]
   essentialAA <- dataCorr[,c(4,11,13,14,19,20,22,23,24)]
   prox <- dataCorr[,c(25,28:33)]
   
# calculate PC for the different sets of attributes   
   para.pca <- prcomp(dataCorr, center = TRUE, scale. = TRUE) #all attributes
   summary(para.pca) #this returns all PC, with their SD, prop. of var., and cumulative proportion
   
   AllAA.pca <- prcomp(allAA, center = TRUE, scale. = TRUE)
   summary(AllAA.pca)
   
   EAA.pca <- prcomp(essentialAA, center = TRUE, scale. = TRUE)
   summary(EAA.pca)
   
   Prox.pca <- prcomp(prox, center = TRUE, scale. = TRUE)
   summary(Prox.pca)
   
# visualizpe PC using ggbiplot package
   #library(devtools)
# install_github("vqv/ggbiplot")
   library(plotly)
   #packageVersion('plotly')
   library(devtools)
   #install_github("vqv/ggbiplot")
   library(ggbiplot)
   library(ggrepel)
# install.packages(ggbiplot)
# !!!!!!!!1is GGfortify is 'on'? if so, the following will not work...!!!!!!!!

   #library(ggfortify)   
# basic plot with data points in space   
   ggbiplot(para.pca)
   ggbiplot(AllAA.pca)
   ggbiplot(EAA.pca)
   ggbiplot(Prox.pca)
   
   #adding labels to the data point
   #playing around... 
   ggbiplot(para.pca, labels=master$newEntry) #labels data points with entry
   ggbiplot(para.pca, labels=master$locEntry) #labels data points with entry and location
   ggbiplot(para.pca, labels=master$POP) #plot treatments in space
   ggbiplot(para.pca, ellipse = TRUE, groups=master$POP)
   #labeling data points with entry and coloring by location
   ggbiplot(para.pca, ellipse = TRUE, groups=master$newLocation, labels = master$newEntry)
   ggbiplot(AllAA.pca, ellipse = TRUE, groups=master$newLocation, labels = master$newEntry)
   ggbiplot(EAA.pca, ellipse = TRUE, groups=master$newLocation, labels = master$newEntry, repel = TRUE) #nice one!
   ggbiplot(Prox.pca, ellipse = TRUE, groups=master$newLocation, labels = master$newEntry)
   
##################these figures above are what has been saved and sent to KM on 6/4/19 for feedback###############
   #expermenting with aesthetics; lots of options
   EAA.pca <- prcomp(essentialAA, center = TRUE, scale. = TRUE)
   
   summary(EAA.pca)
   EAA.pca$sdev #standard devaiation of principal component
   EAA.pca$sdev^2 #eigenvalues
   EAA.pca$rotation
   plot(EAA.pca, type = "l", main = "Essential Amino Acid PC Scree Plot")
   
   EAA.pca <- ggbiplot(EAA.pca, #obs.scale = 1, var.scale = 1,
            groups=master$newLocation, ellipse = TRUE, circle = TRUE, repel = TRUE,
            labels = master$newEntry, labels.size = 3,
            var.axes = TRUE, varname.size = 2, varname.adjust = .5, varname.abbrev = FALSE)


   Prox.pca <- (prcomp(prox, center = TRUE, scale. = TRUE))
   summary(Prox.pca)
   Prox.pca$sdev
   Prox.pca$sdev^2
   Prox.pca$rotation
   plot(Prox.pca, type = "l")
   
   ggbiplot(Prox.pca, obs.scale = 1, var.scale = 1,
            groups=master$newLocation, ellipse = TRUE, circle = TRUE, repel = TRUE,
            labels = master$newEntry, labels.size = 3,
            var.axes = TRUE, varname.size = 3, varname.adjust = 1,varname.abbrev = FALSE)
   
# replotting to try an h
   
   #saving figures in Frontiers in Nutrition ready format
   #2.4.5. Format
   #   The following formats are accepted:
   #      TIFF (.tif) TIFF files should be saved using LZW compression or any other non-lossy compression method.
   #   JPEG (.jpg)
   #   EPS (.eps) EPS files can be uploaded upon acceptance
   #   Color Image Mode
   #   Images must be submitted in the color mode RGB.
   #   Resolution Requirements
   #   All images must be uploaded separately in the submission procedure and have a resolution of 300 dpi at final size. Check the resolution of your figure by enlarging it to 150%. If the resolution is too low, the image will appear blurry, jagged or have a stair-stepped effect.
   #   Please note saving a figure directly as an image file (JPEG, TIF) can greatly affect the resolution of your image. To avoid this, one option is to export the file as PDF, then convert into TIFF or EPS using a graphics software. EPS files can be uploaded upon acceptance.
   
   
   tiff("ProxPCA.tiff", height = 180, width = 180, units = "mm", pointsize = 12, compression = "lzw", res = 300)
   #copy from above when good to go; or use below
   #ggbiplot(Prox.pca, ellipse = TRUE, groups=master$newLocation, labels = master$newEntry, labels.size = 2.5, varname.size = 3, varname.adjust = 1)
   dev.off()
   
   
   tiff("EAAPCA.tiff", height = 180, width = 180, units = "mm", pointsize = 12, compression = "lzw", res = 300)
   ggbiplot(EAA.pca, ellipse = TRUE, groups=master$newLocation, labels = master$newEntry, labels.size = 2.5, varname.size = 3, varname.adjust = 1.5)
   dev.off()
   
   #install or load required packages
   install.packages("ggfortify")
   library("ggfortify")
   
   #clean up column names in 'master'
   names(master) <- gsub(x = names(master), pattern = "\\.", replacement = " ")  #change column names strategically
   head(master)
   
   #remove Lanthionine column
   master <- master[,c(1:16,18:46)]
   
   #PCA using autoplot (ggforitfy)
   ggplot2::autoplot(prcomp(dataCorr)) #the most basic way to visualize
   ggplot2::autoplot(prcomp(dataCorr), data=master, colour = 'newLocation') #color by location
   ggplot2::autoplot(prcomp(dataCorr), data = master, colour = 'newLocation', label = TRUE, label.size = 3) #labels using row names
   ggplot2::autoplot(prcomp(dataCorr), data = master, colour = 'newLocation', label = TRUE, label.size = 3, shape = FALSE) #removes dots
   
   ggplot2::autoplot(prcomp(dataCorr), data = master, colour = 'newLocation', label = FALSE, 
            loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE, loadings.label.size = 3) #all attributes
   
   autoplot(prcomp(allAA), data = master, colour = 'newLocation', label = FALSE, 
            loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE, loadings.label.size = 3) #all AA, including Totals
   
   autoplot(prcomp(noTotalsAA), data = master, colour = 'newLocation', label = FALSE, 
            loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE, loadings.label.size = 3) #all AA, no Totals (obvi)
   
   autoplot(prcomp(essentialAA), data = master, colour = 'newLocation', label = TRUE, 
            loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE, loadings.label.size = 3) #essential AA
   
   autoplot(prcomp(prox), data = master, colour = 'newLocation', label = n, 
            loadings = TRUE, loadings.colour = 'black', loadings.label = TRUE, loadings.label.size = 3)

#$$$$$$$$$$$$$$$$$$plotting PC with ggplot$$$$$$$$$$$$$# 
   set.seed(608)
   
   # PCA to visualize correlations
   # read in data
   master = read.csv("wetChem_OGunits_FINAL.csv") # changed 0 to controls and 110 to Mount Vernon 9.11.19
   dataCorr = master[,c(10:16, 18:43)] #all response variables except lanthionine #master is OG chem units
   head(dataCorr) #we don't want periods in the column names
   
   names(dataCorr) <- gsub(x = names(dataCorr), pattern = "\\.", replacement = " ")  #change column names strategically
   head(dataCorr)
   colnames(dataCorr)[colnames(dataCorr)=="Crude protein"] <- "Crude Protein" #change single column name
   head(dataCorr)
   
   # first clean up column names
   names(dataCorr) <- gsub(x = names(dataCorr), pattern = "\\.", replacement = " ")  #change column names strategically
   head(dataCorr)
   str(dataCorr)
   
   # prepare DF to break up attributes 
   dataCorr #all attributes   
   allAA <- dataCorr[,1:27]
   noTotalsAA <- dataCorr[,1:22]
   essentialAA <- dataCorr[,c(4,11,13,14,19,20,22,23,24)]
   prox <- dataCorr[,c(25,28:33)]
   
   # calculate PC for the different sets of attributes   
   para.pca <- prcomp(dataCorr, center = TRUE, scale. = TRUE) #all attributes
   summary(para.pca) #this returns all PC, with their SD, prop. of var., and cumulative proportion
   
   AllAA.pca <- prcomp(allAA, center = TRUE, scale. = TRUE)
   summary(AllAA.pca)
   
   EAA.pca <- prcomp(essentialAA, center = TRUE, scale. = TRUE)
   summary(EAA.pca)
   
   Prox.pca <- prcomp(prox, center = TRUE, scale. = TRUE)
   summary(Prox.pca)
   
#$$$$$$$$$$$$$$$$$$ANOVA$$$$$$$$$$$$$# 
levels(master$location) #show levels "Chim" "MV"   "Quil" "Sq"  

library(dplyr) #summary statistics by group
group_by(master, location) %>% #crude protein
summarise(
count = n(),
mean = mean(Crudeprotein, na.rm = TRUE),
sd = sd(Crudeprotein, na.rm = TRUE)
)

group_by(master, location) %>% #total essential amino acid content
summarise(
count = n(),
mean = mean(TotalEssentialAA, na.rm = TRUE),
sd = sd(TotalEssentialAA, na.rm = TRUE)
)

group_by(master, location) %>% #lysine
summarise(
count = n(),
mean = mean(Lysine, na.rm = TRUE),
sd = sd(Lysine, na.rm = TRUE)
)

group_by(master, location) %>% #Ash #sig diff likely
summarise(
count = n(),
mean = mean(Ash, na.rm = TRUE),
sd = sd(Ash, na.rm = TRUE)
)


#visualize data
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

library("ggpubr")
ggboxplot(master, x = "location", y = "Crudeprotein", 
color = "location", palette = "Blues",
order = c("Chim", "Quil", "Sq", "MV"),
ylab = "Crude Protein Content", xlab = "Location")

#complete ANOVA for crude protein
res.aov <- aov(Crudeprotein ~ location, data = master)
summary(res.aov)


#$$$$$$$$ data visualization for manuscript figures $$$$$$$$#

# circumplex/polar chart of proximates and amino acid profiles 
library(tidyverse)
library(ggplot2)

theme_update(plot.title = element_text(hjust = 0.5)) #ggplot default title is left-justified

library(RColorBrewer)

se <- function(x) sqrt(var(x)/length(x)) 

# example from online
mtcars$car = row.names(mtcars)
p = ggplot(mtcars, aes(x=car, y=mpg, fill=mpg)) +
geom_bar(width=1, stat='identity') +theme_light() +
scale_fill_gradient(low='red', high='white', limits=c(5,40)) +
theme(axis.title.y=element_text(angle=0))
p + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))

# calculate mean values for essential amino acids and proximates
EssentialAAmean <-round(apply(essentialAA, 2, mean),3)
AAs <- colnames(essentialAA)
seAA <- round(apply(essentialAA,2,se),3)
essential <- data.frame(EssentialAAmean,AAs,se)

ProxMean <- round(apply(prox, 2, mean), 3)
Proxs <- colnames(prox)
seProx <- round(apply(prox,2,se),3)
proximates <- data.frame(ProxMean,Proxs,seProx)


# determine color palette (Color)
library(RColorBrewer)

display.brewer.all(n=9, colorblindFriendly = FALSE)
brewer.pal(n = 9, name = "Set3")

# generate bar chart for essential amino acid profile
plotAA <- ggplot(data = essential, aes(x = factor(AAs, level= c("AAA", "Leucine", "Lysine",
                                              "Valine","SAA","Isoleucine","Threonine","Histidine","Tryptophan")), y = EssentialAAmean, fill = AAs)) +
geom_bar(stat = "identity", width = 1,  color = "white", position = position_dodge()) + 
geom_text(aes(label=EssentialAAmean), vjust=1.6, color="black", size=3.5) +
theme_classic() +
scale_fill_brewer(palette = "Set3") +
geom_errorbar(aes(ymin = EssentialAAmean - seAA,
   ymax = EssentialAAmean + seAA, width = .2)) +
labs(title = "Essential Amino Acid Profile",
tag= "Figure 1",
x="Essential Amino Acids",
y= "Content (g/100g sample)",
caption = "SAA = sulfur amino acids; AAA = aromatic amino acids; error bars represent one standard error of the mean") +
theme(legend.title = element_blank(),
legend.position = "none",
plot.title = element_text(face = "bold", color = "black", size = 16, angle = 0, hjust = 0.5),
axis.title.y=element_text(),
axis.text.y = element_blank(),#element_text(size = 12, color = "black"), #moves y axis title 
axis.title.x = element_text(face = "bold", color = "black", size = 14, angle = 0),
axis.text.x = element_text(color = "black", size = 12, angle = 0), #moves x axis title 
axis.ticks.y = element_blank(),
axis.ticks.x = element_line()) 
plotAA
# generate bar chart for proximate profile
plotPROX <- ggplot(data = proximates, aes(x = factor(Proxs, level= c("Crude Carbohydrates", "Crude Protein ", "Total AA",
                                                   "Moisture","Crude Fat","Ash","Crude Fiber")), y = ProxMean, fill = Proxs)) +
geom_bar(stat = "identity", width = 1,  color = "white") + 
geom_text(aes(label=ProxMean), vjust=1.6, color="black", size=3.5) +
theme_classic() +
scale_fill_brewer(palette = "Set3") +
geom_errorbar(aes(ymin = ProxMean - seProx,
   ymax = ProxMean + seProx), width = .2) +
labs(title = "Seed Composition",
tag= "Figure 2",
x="Seed Composition Components",
y= "Content (g/100g sample)",
caption = "AA = amino acids") +
theme(legend.title = element_blank(),
legend.position = "none",
plot.title = element_text(face = "bold", color = "black", size = 16, angle = 0, hjust = 0.5),
axis.title.y=element_text(),
axis.text.y = element_text(size = 12, color = "black"), #moves y axis title 
axis.title.x = element_text(face = "bold", color = "black", size = 14, angle = 0),
axis.text.x = element_text(color = "black", size = 12, angle = 0), #moves x axis title 
axis.ticks.y = element_line()) 

plotPROX

# ADJust barcharts for windrose
# generate adj bar chart for essential amino acid profile
plotAAwindrose <- ggplot(data = essential, aes(x = factor(AAs, level= c("AAA", "Leucine", "Lysine",
                                                      "Valine","SAA","Isoleucine","Threonine",
                                                      "Histidine","Tryptophan")), 
                             y = EssentialAAmean, 
                             fill = EssentialAAmean)) + #change fill to value for mean
geom_bar(stat = "identity", width = 1,  color = "white", position = position_dodge()) + 
#geom_text(aes(label=EssentialAAmean), vjust=5, color="black", size=3.5) +
theme_light() + #change theme to get gridlines
scale_fill_gradient(low="#BEBADA", high = "#FB8072", limits=c(.4,.8), name="Content (g/100g sample)") +
#scale_fill_brewer(palette = "Set3") +
#geom_errorbar(aes(ymin = EssentialAAmean - seAA,
#ymax = EssentialAAmean + seAA, width = .2)) +
labs(title = "Essential Amino Acid Profile",
tag= "Figure 1",
x="Essential Amino Acids",
y= "Content (g/100g sample)",
caption = "SAA = sulfur amino acids; AAA = aromatic amino acids; error bars represent one standard error of the mean") +
theme(legend.title = element_text(),
legend.position = "right", #use keywords such as "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center" to control the position of the legend on a given R graph.
plot.title = element_text(face = "bold", color = "black", size = 16, angle = 0, hjust = 0.5),
axis.title.y=element_blank(),#element_text(),
axis.text.y = element_blank(),#element_text(size = 12, color = "black"), #moves y axis title 
axis.title.x = element_blank(),#element_text(face = "bold", color = "black", size = 14, angle = 0),
axis.text.x = element_text(face= "bold", color = "black", size = 12, angle = 0), #moves x axis title 
axis.ticks.y = element_blank(),
axis.ticks.x = element_line()) 
plotAAwindrose

CircumAA<- plotAAwindrose + coord_polar()


# generate adj bar chart for seed composition profile
plotPROXwindrose <- ggplot(data = proximates, aes(x = factor(Proxs, level= c("Crude Carbohydrates", "Crude Protein ", "Total AA",
                                                           "Moisture","Crude Fat","Ash","Crude Fiber")), y = ProxMean, fill = ProxMean)) +
geom_bar(stat = "identity", width = 1,  color = "white") +

# geom_text(aes(label=Proxs, vjust="inward",hjust="inward"))+ adds prox label to bar
theme_light() + #change theme to get gridlines
scale_fill_gradient(low="#BEBADA", high = "#FB8072", limits=c(2,70), name="Content (g/100g sample)") +

#scale_fill_brewer(palette = "Set3") +
#geom_errorbar(aes(ymin = EssentialAAmean - seAA,
#ymax = EssentialAAmean + seAA, width = .2)) +
labs(title = "Seed Composition",
tag= "Figure 2",
x="Seed Composition Components",
y= "Content (g/100g sample)",
caption = "AA = amino acids") +
theme(legend.title = element_text(),
legend.position = "right", #use keywords such as "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center" to control the position of the legend on a given R graph.
plot.title = element_text(face = "bold", color = "black", size = 16, angle = 0, hjust = 0.5),
axis.title.y=element_blank(),#element_text(),
axis.text.y = element_blank(),#element_text(size = 12, color = "black"), #moves y axis title 
axis.title.x = element_blank(),#element_text(face = "bold", color = "black", size = 14, angle = 0),
axis.text.x = element_text(vjust=0,face= "bold", color = "black", size = 12, angle = 0), #moves x axis title 
axis.ticks.y = element_blank(),
axis.ticks.x = element_line())

plotPROXwindrose

CircumPROX<- plotPROXwindrose + coord_polar()
CircumPROX


mtcars$car = row.names(mtcars)
p = ggplot(mtcars, aes(x=car, y=mpg, fill=mpg)) +
geom_bar(binwidth=1, stat=’identity’) +theme_light() +
scale_fill_gradient(low=’red’, high=’white’, limits=c(5,40)) +
theme(axis.title.y=element_text(angle=0))
p + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))


plotAA <- ggplot(essential, aes(AAs, EssentialAAmean, fill=EssentialAAmean)) + 
geom_bar(stat = "identity", width = 1, color = "white", aes(fill=AAs, fill_palette(palette = cbbPalettBlack))) + 
geom_errorbar(aes(ymin = EssentialAAmean - se(essential$EssentialAAmean),
   ymax = EssentialAAmean + se(essential$EssentialAAmean)), width = .2) +
theme_classic() +
theme(axis.ticks = element_blank(),
axis.text.x = element_text(face = "plain", color = "black", size = 14, angle = 45),
axis.text.y = element_blank(),
axis.title = element_blank(),
axis.line = element_blank()
) +
labs(title = "Essential Amino Acid Profile",
tag = "Figure 1",
x = "Essential Amino Acids",
y = "Content (g/100g sample)") +
theme(legend.title = element_text("Content (g/100g sample)"))

# generate bar chart round II

plotAA2 = ggplot(essential, aes(x=AAs, y=EssentialAAmean)) +
geom_bar(width=1, stat="identity", color = "white") + theme_light() +
scale_fill_gradient(low="red", high="white", limits=NA) +
theme(axis.title.y=element_text(angle=0))
plotAA2
plotAA2 + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))

p = ggplot(mtcars, aes(x=car, y=mpg, fill=mpg)) +
geom_bar(binwidth=1, stat=’identity’) +theme_light() +
scale_fill_gradient(low=’red’, high=’white’, limits=c(5,40)) +
theme(axis.title.y=element_text(angle=0))
p + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))


# Revisions
AdjMaster = read.csv("wetChem_adj_FINAL.csv") #AAs are in g/100g crude protein
str(AdjMaster)
AdjMaster$POP <- as.factor(AdjMaster$POP)

# Create column for "Total Carbohydrates"
colnames(AdjMaster)

AdjMaster <- AdjMaster %>%
   mutate(Total.Carb = (100- (Crude.protein + Crude.Fat + Ash + Moisture)))

# remove all proximates and Amino ACids
DescAll <- tabular( (LOCATION + POP + 1) ~ (n=1) + Format(digits=2) * (Total.Carb +
                                                                          Crude.protein +
                                                                          Moisture +
                                                                          Crude.Fat +
                                                                          Ash +
                                                                          Total.AA +
                                                                          Total.Essential.AA +
                                                                          Histidine +
                                                                          Isoleucine +
                                                                          Leucine +
                                                                          Lysine +
                                                                          SAA +
                                                                          Methionine + 
                                                                          AAA +
                                                                          Phenylalanine + 
                                                                          Threonine +
                                                                          Tryptophan +
                                                                          Valine +
                                                                          Total.Nonessential.AA +
                                                                          Glutamic.Acid +
                                                                          Aspartic.Acid +
                                                                          Arginine +
                                                                          Glycine +
                                                                          Alanine +
                                                                          Proline +
                                                                          Serine +
                                                                          Tyrosine + 
                                                                          Cysteine + 
                                                                          Taurine +
                                                                          Hydroxyproline +
                                                                          Hydroxylysine +
                                                                          Ornithine + 
                                                                          Lanthionine) * (mean + sd), data=AdjMaster) #these are the columns 
DescAll

DescAll.m = as.matrix(DescAll)
DescAll.m.t = t(DescAll.m)

write.csv.tabular(DescAll.m.t, file = "Descriptives_All_LOC_POP.csv", justification = "n", row.names=FALSE)

# Revisions
#Figure 2
master = read.csv("wetChem_OGunits_FINAL2.csv") # changed 0 to controls and 110 to Mount Vernon 9.11.19
colnames(master)
dataCorr = master[,c(10:16, 18:40,42,43)] #all response variables except lanthionine #master is OG chem units
head(dataCorr) #we don't want periods in the column names

library(Hmisc)
res.corr <- rcorr(as.matrix(dataCorr))
head(res.corr)
summary(res.corr)

#visualize
library(ggcorrplot)
ggcorrplot(res.corr$r,
           method = "circle",
           outline.color = "white",
           type = "upper",
           hc.order = TRUE,
           hc.method= "ward",
           lab = TRUE,
           lab_col = "black",
           lab_size = 0.1,
           p.mat=res.corr$P,
           sig.level = 0.05,
           insig = "pch",
           show.diag = FALSE)
#basic
ggcorrplot(corr, # Leave blank on no significant coefficient
           p.mat = p.mat, hc.order = TRUE,
           type = "lower", insig = "blank")


#customized correlation matrix
library(colorspace)

colors <- rainbow_hcl(3)
colors <- c("red", "yellow", "blue")

ggcorrplot(corr, 
           method = "square", 
           type = "lower", 
           ggtheme = ggplot2::theme_minimal, title = " Correlations Between Nutritional Components",
           show.legend = TRUE, legend.title = "Pearson Correlation Coefficient", show.diag = TRUE,
           colors = colors, outline.color = "gray",
           hc.order = TRUE, hc.method = "complete", 
           lab = TRUE, lab_col = "black", lab_size = 1, p.mat = p.mat, sig.level = 0.05,
           insig = "pch", pch.col = "black", pch.cex = 2,
           tl.cex = 4, tl.col = "black", tl.srt = 45, digits = 2)

#final$$$$$$$$$$$$$$
ggcorrplot(corr,
           method = "square",
           type = "upper",
           legend.title = "Correlation Coefficients",
           colors = c("red", "white", "blue"),
           p.mat = p.mat, hc.order = TRUE, hc.method = "complete",
           lab=TRUE, lab_col = "black", lab_size = 3, digits = 2, #size and color for the correlation coefficient labels
           insig = "pch", pch.col = "black", pch.cex = 5, #glyphs for insiginificant
           tl.cex =10, tl.col = "black",tl.srt=70) #the size, the color and the string rotation of text label (variable names).

# revisions to Figure 2 (using)
RawMaster = read.csv("wetChem_OGunits_FINAL2.csv") #AAs are in g/100g sample
colnames(RawMaster)
# !!! make sure to run corstars above
# RAW (not adjusted)
RAWdataCorr = RawMaster[,c(10:16, 18:40,42,43)] #all response variables except lanthionine and crude fiber

head(RAWdataCorr) #we don't want periods in the column names

sum(is.na(dataCorr)) #no NA values
str(dataCorr) #sampleID present

res3 <- rcorr(as.matrix(RAWdataCorr, type = "spearman"))

library(corrplot)
library(RColorBrewer)

Figure2_V2 <-
   corrplot(res3$r, method = "circle", 
            outline = "white",
            type = "upper", 
            order = "hclust", hclust.method = "ward",
            col = brewer.pal(n = 4, name = "RdBu"),
            diag = FALSE, 
            addCoef.col = "black",
            number.cex = 0.5,
            number.digits = 2,
            p.mat = res3$P, sig.level = 0.5, insig = "pch",
            tl.pos = "td", tl.col = "black", tl.srt = 40, tl.cex = 0.65, 
            cl.pos = "r", cl.lim = c(-1,1), 
            cl.ratio = 0.1, cl.align.text = "l"
   )
tiff(filename = 'Figure2_V2.tiff', height = 180, width = 180, units = "mm", bg="white", res = 300)
Figure2_V2 <-
   corrplot(res3$r, method = "circle", 
            outline = "white",
            type = "upper", 
            order = "hclust", hclust.method = "ward",
            col = brewer.pal(n = 4, name = "RdBu"),
            diag = FALSE, 
            addCoef.col = "black",
            number.cex = 0.5,
            number.digits = 2,
            p.mat = res3$P, sig.level = 0.5, insig = "pch",
            tl.pos = "td", tl.col = "black", tl.srt = 40, tl.cex = 0.65, 
            cl.pos = "r", cl.lim = c(-1,1), 
            cl.ratio = 0.1, cl.align.text = "l"
   )
dev.off()

dir()


## Figure 3
library("FactoMineR")
library("factoextra")
# PCA to visualize correlations
# read in data
master = read.csv("wetChem_OGunits_FINAL2.csv") # changed 0 to controls and 110 to Mount Vernon 9.11.19
colnames(master)
dataCorr = master[,c(10:16, 18:43,50)] #all response variables except lanthionine #master is OG chem units
head(dataCorr) #we don't want periods in the column names

# prepare DF to break up attributes 
dataCorr #all attributes   
colnames(dataCorr)
allAA <- dataCorr[,1:27]
noTotalsAA <- dataCorr[,1:22]
essentialAA <- dataCorr[,c(4,11,13,14,19,20,22,23,24)]
prox <- dataCorr[,c(25,28:33)]


# Figure 3

# PCA for EAA
essentialAA <- dataCorr[,c(4,11,13,14,19,20,22,23,24)]
rownames(essentialAA) = master$labeltwo
labeys = master$labelthree
labeys2 = rownames(labeys[labeys!= ""])


res.pca <- PCA(essentialAA,  graph = FALSE)

fviz_pca_ind(res.pca, geom.ind = c("text","point"),col.ind = master$newLocation, repel = TRUE)
fviz_pca_var(res.pca, col.var = "black")


FigThreeLabelsworking <- fviz_pca_biplot(res.pca, repel = TRUE,
                                  geom = c("text","point"),
                                  select.row = labeys2,
                                  label = "all",
                                  col.ind = master$newLocation,
                                  col.var = "black",
                                  mean.point = FALSE
)
FigThreeLabelsworking

FigThreeLabels <- fviz_pca_biplot(res.pca, repel = TRUE,
                            geom = "point",
                            label = "sup.ind",
                            col.ind = master$newLocation,
                            col.var = "black",
                            mean.point = FALSE
)
FigThreeLabels # no labels, use preview to label

#addEllipses = TRUE, ellipse.level = 0.95)

FigThree <- fviz_pca_biplot(res.pca, 
                geom.ind = c("point"),
                col.ind = master$newLocation,
                pointsize = 2,
                labelsize = 4,
                mean.point = FALSE,
                label = "all",
                geom.var = "arrow", 
                col.var = "black",
                repel = TRUE,
                legend.title = "Location",
                legend.position = "top")
                #addEllipses = TRUE, ellipse.level = 0.95)
FigThree


#pdf("Figure3_V2.pdf")
#print(FigThree)
#dev.off()

tiff(filename = 'Figure3_V2.tiff', height = 180, width = 180, units = "mm", bg="white", res = 300)
dev.off()

# Supplementary Figure 1 - Version 2
prox <- dataCorr[,c(25,28:33)]
prox <- prox[,-5] # remove crude fiber
colnames(prox)

colnames(prox) = 



rownames(prox) = master$labeltwo
labeys = master$labelthree
labeys2 = rownames(labeys[labeys!= ""])


res.pca.prox <- PCA(prox,  graph = FALSE)

fviz_pca_ind(res.pca.prox, geom.ind = c("text","point"),col.ind = master$newLocation, repel = TRUE)
fviz_pca_var(res.pca.prox, col.var = "black")


SuppFigOneLabelsworking <- fviz_pca_biplot(res.pca.prox, repel = TRUE,
                                         geom = c("text","point"),
                                         select.row = labeys2,
                                         label = "all",
                                         col.ind = master$newLocation,
                                         col.var = "black",
                                         mean.point = FALSE
)
SuppFigOneLabelsworking

SuppFigOneLabels <- fviz_pca_biplot(res.pca.prox, repel = TRUE,
                                  geom = "point",
                                  label = "sup.ind",
                                  col.ind = master$newLocation,
                                  col.var = "black",
                                  mean.point = FALSE
)
SuppFigOneLabels # no labels, use preview to label

#addEllipses = TRUE, ellipse.level = 0.95)

SuppFigOne <- fviz_pca_biplot(res.pca.prox, 
                            geom.ind = c("point"),
                            col.ind = master$newLocation,
                            pointsize = 2,
                            labelsize = 4,
                            mean.point = FALSE,
                            label = "all",
                            geom.var = "arrow", 
                            col.var = "black",
                            repel = TRUE,
                            legend.title = "Location",
                            legend.position = "top")
#addEllipses = TRUE, ellipse.level = 0.95)
SuppFigOne

tiff(filename = 'SuppFigure1_V2.tiff', height = 180, width = 180, units = "mm", bg="white", res = 300)
dev.off()




# color by groups
set.seed(123)
res.km <- kmeans(res.pca$ind$coord, centers=3, nstart = 25)
grp <- as.factor(res.km$cluster)

# color individuals by groups
fviz_pca_biplot(res.pca, col.ind = grp, palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")

# Figure 5

# PCA for Prox
res.pca2 <- PCA(prox, graph = FALSE)

# Create Leucine content by location
library(ggplot2)
LeuLoc <- read.csv("wetChem_adj_FINAL.csv")
str(LeuLoc)


p <- ggplot(LeuLoc, aes(x = LeuRank, y=Leucine, shape = newLocation, color = newLocation)) +
   geom_point(size = 5) +
   geom_hline(yintercept=6.25, color = "black", size = 1) +
   geom_hline(yintercept=5.90, color = "black", linetype = "dashed", size = 1) +
   geom_hline(yintercept=6.6, color = "black", linetype = "dotted", size = 1) +
   labs(legend.title = "Location", x =element_blank(), y =element_blank()) +
   theme_minimal() + 
   theme(legend.position = c(0.9, 0.9)) +
   labs(color ="Location", shape = "Location") +
   theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank()) +
   
   theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
   ) + 
   theme(axis.ticks.y = element_line()) +
   theme(
      legend.position = c(1.0, 1.0),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
   )

p
ggsave(p, filename = "LeuLoc_6.46_11.79_tr.tiff",  height = 6.46, width = 11.79, units = "in", bg = "transparent")



##########################
p <- ggplot(LeuLoc, aes(x = LeuLocRank, y=Leucine, shape = newLocation, color = newLocation)) +
               geom_point() +
   geom_hline(yintercept=6.25, color = "black", size = 0.5) +
   geom_hline(yintercept=5.90, color = "black", linetype = "dashed", size = 0.5) +
   geom_hline(yintercept=6.6, color = "black", linetype = "dotted", size = 0.5) 
                   
p
a