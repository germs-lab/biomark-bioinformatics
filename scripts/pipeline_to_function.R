setwd("~/Desktop/biomark-private")
df <- readRDS("data/metawithbiomark.RDS") 
PC_Standards_no_NA <- readRDS("data/ampliconlengths.RDS")

# subset biomark chip to interesting things
soil <- df %>%
  filter(sample_type %in% "soil")
water <- df %>%
  filter(sample_type %in% "water") 
manure <- df %>%
  filter(sample_type %in% "manure")
control <- df %>%
  filter(sample_type %in% "control")
standard <- df %>%
  filter(sample_type %in% "standard")

m <- soil %>%
  filter(Value < 28, Call == "Pass", treatment == "WCS") %>%
  group_by(Assay) %>%
  dplyr::mutate(mean = mean(as.numeric(Value)))

options(scipen = 999)

# plot standard curve from biomark standards using expected copy numbers based on amplicon length for each Assay
plot_standard_curve_counts <- function(standards, gene){

  std <- standards %>%
    filter(Value < 990, grepl(gene, as.character(Assay)), grepl(gene, as.character(Sample_Name))) %>%
    left_join(PC_Standards_no_NA) %>%
    select(Assay, rConc, Value, amp_len) %>%
    mutate(counts = as.numeric(rConc) * (1/(10^9)) * (1/660) * ((6.023*10^23)/1) * (1/as.numeric(amp_len)) * (.0067)) %>%
    ggplot(aes(x = log10(counts + 1), y = as.numeric(Value))) +
    geom_point() + geom_smooth(method = 'lm') +
    labs(x = "Count (log10)", y = "Ct") +
    ggtitle(gene)
  
  return(std)
}

#Function to calculate mean Ct value if you wish to use it
meanCT <- function(sample_data, gene){
  soil16s <- sample_data %>%
    filter(grepl(gene, Assay), Value < 28) %>%
    select(ID, Sample_Name, Assay, Value, treatment, sample_day, sample_number, soil_type, plot) %>%
    group_by(Sample_Name) %>%
    dplyr::mutate(mean = mean(Value)) #%>% # if you load `plyr` after dplyr it will not work correctly
  #select(-Value) #%>%
  #unique()
  return(soil16s)
}

plot_counts_function <- function(stddf, countdf, gene){
  sixS <- plot_standard_curve_counts(stddf, gene)
  countdf$Value <- as.numeric(countdf$Value)
  msoil16S <- meanCT(countdf, gene)
  dt <- sixS$data
  inverse.lm <- lm(data = dt, formula = log10(counts+1) ~ as.numeric(Value))

  dt2 <- msoil16S
  val4 <- dt2$Value
  dt2$sixS <- 10 ^ predict(inverse.lm, data.frame(Value = val4), interval = "predict")[,1]
  d1 <- dt2 %>%
    filter(!sample_day %in% "T14") # drop random T14?
  return(d1)
}

list <- unique(m$Assay)
list
class(standard$Assay)
list2 <- unique(standard$Assay)
# What is in both
both <- intersect(list, list2)
both == list

intersect(list, PC_Standards_no_NA$Assay)



for (i in 1:length(list)) {
  std <- standard %>%
    filter(Value < 990, grepl("i", as.character(Assay)), grepl("i", as.character(Sample_Name)))
}
intersect(list, standard$Assay)
intersect(list, standard$Sample_Name)

test <- plot_standard_curve_counts(standard, "16S")

list

my_comparisons <- list( c("WCS", "WCSM"), c("WCS", "WCM"), c("WCM", "WCSM"))
ggboxplot(d1, x = "treatment", y = "sixS", fill = "treatment") +
  facet_grid(sample_number ~ sample_day) +
  stat_compare_means(comparisons = my_comparisons) +
  scale_y_continuous(trans = 'log10') +
  labs(x = "Depth 1 soil samples", y = "Log 10 16s gene copies") + 
  scale_fill_viridis_d() +
  theme(legend.position = "none") 
