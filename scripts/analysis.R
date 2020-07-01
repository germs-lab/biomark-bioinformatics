#devtools::install_github("sjackman/uniqtag")
library(plyr)
library(uniqtag)
library(ggplot2)

setwd("~/Downloads/biomark-private/")
raw_data_file = "~/Downloads/biomark-private/combined.cleaned.csv"
raw_data <- read.csv(raw_data_file, sep=",", header = FALSE)
header = c('ID','Name','Type','rConc','Name.1','Type','Value','Calibrated rConc','Quality','Call','Threshold','In Range','Out Range','Peak Ratio', 'Plate')
colnames(raw_data) <- header
dim(raw_data)
head(raw_data)
raw_data$full_assay <- paste(raw_data$Name,"-", raw_data$Name.1, sep="")
raw_data$unique_id <- make_unique(raw_data$full_assay)


# Metadata fie has to have a header names "sample_type" where you specify standards or env sample
# You should have metadata for both samples and standards
# It would be better if the "standard sample" had the appropriate dilution
meta_data_file = "~/Downloads/biomark-private/meta/Worle_biomark_meta.csv"
meta <- read.csv(meta_data_file, header=TRUE)


# Check the NTCs
ct_threshold_ntc = 20
ntc <- subset(raw_data, Name == "NTC")
summary(ntc$Value)
check_ntc <- subset(ntc, Value > 20 & Value < 999)
check_ntc$Value
data_nontc <- subset(raw_data, Name != "NTC")

# This is the main data table with all the metadata merged in
# It should be the same size as the data_nontc
dim(data_nontc)
merged <- merge(data_nontc, meta, by = "Name")
dim(merged)

#Standards Workflow

standards <- subset(merged, sample_type == "standard")
standards <- subset(standards, Value < 500)
gene_list <- unique(standards$Name)
standards_matrix <- matrix(, nrow = length(gene_list), ncol=4)

for (x in 1:length(gene_list)) {
   each_gene <- subset(standards, Name == gene_list[x])
   standard_only = paste(gene_list[x], '-', gene_list[x], sep='')
   each_gene_filtered = subset(each_gene, full_assay == standard_only)
   each_gene_filtered = subset(each_gene_filtered, Value < 500)
   each_gene_badfiltered = subset(each_gene, full_assay != standard_only)
   
   f <- ddply(each_gene_filtered, .(rConc, full_assay, Plate), summarise, MEAN = mean(Value), SE=sd(Value)/sqrt(length(Value)))
   limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
   f$log <- log10(f$rConc)
   p = ggplot(f, aes_string(x="log", y="MEAN", color="Plate")) + geom_point() + geom_errorbar(limits)+theme(axis.text.x = element_text(angle = 90))
   p + ggtitle(gene_list[x])
   ggsave(paste0(x,'standard.pdf'), device="pdf")
   lm_fit <- lm(MEAN ~ log, data=f)
   summary(lm_fit)
   standards_matrix[x,1] = gene_list[x]
   standards_matrix[x,2] = coefficients(lm_fit)[1]
   standards_matrix[x,3] = coefficients(lm_fit)[2]
   standards_matrix[x,4] = anova(lm_fit)$'Pr(>F)'[1]
  }

standards_matrix <- as.data.frame(standards_matrix)
colnames(standards_matrix) <- c("Name", "Intercept", "slope", "anova")
standards_matrix

#plots for our proposal
samples <- subset(merged, sample_type != "standard") #all samples
water <- subset(samples, sample_type == "water")
water <- subset(water, Value < 100)
f <- ddply(water, .(full_assay, plot, treatment, Name, Name.1, sample_day), summarise, MEAN = mean(Value), SE=sd(Value)/sqrt(length(Value)))
#genes that are detected #stats are done prior to removing bad measurements
f_detects <- subset(f, MEAN < 100)
print("Number of Genes with detects:  ")
length(unique(f_detects$Name.1))
unique(f_detects$Name.1)
f_plot <-ddply(f, .(treatment, Name.1, sample_day), summarise, MEAN2 = mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN)))
f2 <- subset(f_plot, MEAN2 < 30) #Can set as appropriate - filter by MEAN technical rep Ct
temp = subset(f2, treatment == "WCM")
f2$Name.1 <- factor(f2$Name.1, levels=temp$Name.1[order(-temp$MEAN)])
limits2<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
p = ggplot(f2, aes_string(x="Name.1", y="MEAN2", color="Name.1"))
p + geom_point() + geom_errorbar(width=0.25, limits2)+theme_classic()+ theme(text=element_text(size=15, family="Helvetica"),axis.text.x = element_text(angle = 90), legend.position="none") + facet_grid(sample_day~treatment)

samples <- subset(merged, sample_type != "standard") #all samples
water <- subset(samples, sample_type == "soil")
water <- subset(water, Value < 100)
f <- ddply(water, .(full_assay, plot, treatment, Name, Name.1, sample_day), summarise, MEAN = mean(Value), SE=sd(Value)/sqrt(length(Value)))
#genes that are detected #stats are done prior to removing bad measurements
f_detects <- subset(f, MEAN < 100)
print("Number of Genes with detects:  ")
length(unique(f_detects$Name.1))
unique(f_detects$Name.1)
f_plot <-ddply(f, .(treatment, Name.1, sample_day), summarise, MEAN2 = mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN)))
f2 <- subset(f_plot, MEAN2 < 30) #Can set as appropriate - filter by MEAN technical rep Ct
temp = subset(f2, treatment == "WCM")
f2$Name.1 <- factor(f2$Name.1, levels=temp$Name.1[order(-temp$MEAN)])
limits2<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
p = ggplot(f2, aes_string(x="Name.1", y="MEAN2", color="Name.1"))
p + geom_point() + geom_errorbar(width=0.25, limits2)+theme_classic()+ theme(text=element_text(size=15, family="Helvetica"),axis.text.x = element_text(angle = 90), legend.position="none") + facet_grid(sample_day~treatment)




samples <- subset(merged, sample_type != "standard") #all samples
water <- subset(samples, sample_type == "soil")
water <- subset(water, treatment == "WCM")
water <- subset(water, Value < 100)
water$Value_inv <- 1/water$Value
f <- ddply(water, .(treatment, Name.1, sample_day), summarise, MEAN = mean(Value_inv), SE=sd(Value_inv)/sqrt(length(Value_inv)))
temp <- subset(f, sample_day == "TB")
f$sample_day = factor(f$sample_day, levels=c("TB","T00","T02"))
f$Name.1 <- factor(f$Name.1, levels=f$Name.1[order(-f$MEAN)])
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f, aes_string(x="Name.1", y="MEAN", color="Name.1"))
p + geom_point() + geom_errorbar(width=0.25, limits)+ theme(text=element_text(size=20, family="Helvetica"),axis.text.x = element_text(angle = 90), legend.position="none") + facet_grid(sample_day~treatment)


samples <- subset(merged, sample_type != "standard") #all samples
manure <- subset(samples, sample_type == "manure")
crop <- subset(samples, treatment = "WCS")
crop <- subset(crop, sample_type == "soil")
crop <- subset(crop, soil_type == "crop")
total <- rbind(crop, manure)
total <- subset(total, Value < 100)
f <- ddply(total, .(sample_type, Name.1), summarise, MEAN = mean(Value), SE=sd(Value)/sqrt(length(Value)))
f#genes that are detected #stats are done prior to removing bad measurements

temp = subset(f2, treatment == "WCM")
f2$Name.1 <- factor(f2$Name.1, levels=temp$Name.1[order(-temp$MEAN)])
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f, aes_string(x="Name.1", y="MEAN", color = "sample_type"))
p + geom_point() + geom_errorbar(width=0.25, limits)+theme_classic()+ theme(text=element_text(size=15, family="Helvetica"),axis.text.x = element_text(angle = 90)) 




#Need to divide values by 16S
#Some code we'd run on the whole plate data to get 16S per sample for division - Abundance here is estimated
smple_v<- unique(f2$Plate)
k<- 1
f2$housekeeping_abund<- NA
housekeeper = "16S_Eub_338F_515R"
for(k in 1:length(smple_v)) {
  smple_r<-which(f2$Plate==smple_v[k]) #find rows in data frame that belong to a particular sample
  smple_abun<-mean(f2[smple_r[which(f2$Gene[smple_r]==housekeeper)],"Abundance"]) #average abundance of housekeeper gene for a given sample
  #can update function here if a different means of normalizing is desired (e.g., max)
  f2$housekeeping_abund[smple_r]<-smple_abun #assign housekeeper abundance for sample to all rows from sample
}


