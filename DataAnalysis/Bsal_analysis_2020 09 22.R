#####################################
##### The following R script was written by G. V. DiRenzo
### Please send questions to: grace.direnzo@gmail.com
#####################################




# Objective: 
  # 1. To recreate all analyses and figures from the Bsal experiment with 8 salamanders collected from the eastern USA




#####################################
###### Table of contents ############
#####################################

# 1. Set working directory & download libraries
# 2. Read in the data
# 3. Format data
# 4. Bd  summary
# 5. Figure S1. Bd load over time 
# 6. Bsal exposure summary 
# 7. Day 3 and 4 averages
# 8. Last swabbing event - day 81 - summary
# 9. Co-infection summary 
# 10. Body condition calculation 
# 11. Survival analysis & Fig. S1 & S2
# 12. Determine the total number of Plethodon that die & their average infection intensity
# 13. First phase: data format
# 14. First phase: Bsal infection intensity analysis, Supporting Information 1: Parts 1-4
# 15. Diagnotic plots, Supporting Information 1: Part 5 
# 16. Fig 3: 1st phase - Bsal infection intensity vs. time
# 17. Fig 5. Influence of Bd at t - 1 on Bsal at time t
# 18. First phase: Body condition analysis, Supporting Information 2: Parts 1 - 4 
# 19. Fig S3: Body condition vs. time 
# 20. Diagnotic plots, Supporting Information 2: Part 5
# 21. Second phase: data format
# 22. Second phase: Bsal infection intensity analysis, Supporting Information 3: Parts 1-4 
# 23. Diagnotic plots, Supporting Information 3: Part 5 
# 24. Fig 4: Bsal infection intensity vs. time
# 25. Second phase: Body condition analysis, Supporting Information 4: Parts 1-4 
# 26. Fig S5: Body condition vs. time 
# 27. Diagnotic plots, Supporting Information 4: Part 5
# 28. Table 1 & 2 Main text
# 29. Phylogeny: read in data
# 30. Define species to keep
# 31. Order the tips
# 32. Ancestoral state reconstruction
# 33. Fig. 6 main text
# 34. Test for phylogenetic signal


#####################################
#####################################
#####################################



# 1. Set working directory & download libraries------------------------------------------------




# Set working directory
setwd("~/Dropbox/Bsal_Summer_2015")


# download libraries
library(ggplot2)
library(plyr)
library(lme4)
library(emmeans)
library(lattice)
library(survival)
library(RVAideMemoire)
library(nlme)
library(ape)
library(picante)
library(phytools)
library(geiger)
library(AICcmodavg)




# 2. Read in the data ------------------------------------------------


# Read in data
  # This file has infection intensity of chytrids over time
  # The phylogeny and histology data will be read in later
Bsal <- read.table("./Data/GVDBsal_Experiment_MASTER_25FEB2016.txt", 
                   header = T, sep = "\t")






# 3. Format data ------------------------------------------------



# Look at data structure
str(Bsal)
# 'data.frame':	2135 obs. of  34 variables:
# $ Day              : Day of observation
# $ Month            : Month of observation
# $ Year             : Year of observation
# $ Exp_day          : Experimental day
# $ Phase            : Either phase 1 (low dose) or 2 (high dose)
# $ Site             : Collection site of individual
# $ Genus            : Genus
# $ Species          : Species
# $ ID_number        : ID number
# $ infection_tube   : The tube label during inoculation
# $ SVL              : Snout-to-vent length (mm)
# $ Tail             : Length of tail (mm)
# $ Mass             : Mass (g)
# $ Regen_Tail       : Yes or No tail regeneration?
# $ Skin_Scrapes     : Any skin scrapes present?
# $ Treatment        : Treatment: Control, Inoculated (1x), Double inoculation(2x)
# $ Incubator_number : Number of incubator where tank is housed
# $ Shelf_number     : The shelf in the incubator where tank is housed
# $ Swab_ID          : ID of swab collected
# $ Census           : Is the animal alive or dead? 0 = alive; 1 = dead
# $ Bd_hist          : 
# $ Bsal_hist        : 
# $ Specimen_ID      : Tag number assigned after euthanasia
# $ Bd_load          : Bd load from qPCR = 1st swab, 1st PCR run
# $ Bsal_load        : Bsal load from qPCR = 1st swab, 1st PCR run
# $ Repeat1_Bd       : Bd load from qPCR = 1st swab, 2nd PCR run
# $ Repeat1_Bsal     : Bsal load from qPCR = 1st swab, 2nd PCR run
# $ DoubleBd         : Bd load from qPCR = 2nd swab, 1st PCR run
# $ DoubleBsal       : Bsal load from qPCR = 2nd swab, 1st PCR run
# $ DoubleBd_repeat  : Bd load from qPCR = 2nd swab, 2nd PCR run
# $ DoubleBsal_repeat: Bsal load from qPCR = 2nd swab, 2nd PCR run
# $ Bsal             : Same as Bd_load column
# $ Bd               : Same as Bsal_load column
# $ Bd_load_t_min_1  : Bd load at the previous swabbing event


# Remove space from Month column name
colnames(Bsal)[2] <- "Month"


# Fix names
  # Change the name of D. organi to D. wrighti
  # Change the name of E. sp to E. wilderae
Bsal$Species <- as.character(Bsal$Species)
Bsal$Species <- ifelse(Bsal$Species == "organi", "wrighti", Bsal$Species)
Bsal$Species <- ifelse(Bsal$Species == "sp", "wilderae", Bsal$Species)
Bsal$Species <- as.factor(Bsal$Species)
Bsal <- droplevels(Bsal)


# Subset data from when individuals were collected in the field
field <- Bsal[Bsal$Day == 7 & Bsal$Month == 6,]

# Determine the number of animals that came in from the field Bd+ the first day
ddply(.data = field,
      .variable = c("Species"),
      .fun = summarize,
      Bd_pos = length(which(Bd_load > 0)))

# Add a dummy variable "Num" to be able to tally up values
Bsal$Num <- 1


# Determine the total number of indviduals collected in the field
collected <- ddply(.data = Bsal,
                    .variable = c("Species"),
                    .fun = summarize,
                    samp = length(unique(ID_number)))
sum(collected[,2])

# Remove 3 individuals that were exposed to Bsal but never had a positive Bsal 
# Remove the data from a few other individuals that died following the inoculations, and others died because of unforseeable circumstances
# Data from these individuals should NOT be used
  # D. fuscus number 9
  # N. viridescens number 10
  # P. cinereus number 4
  # P. cinereus number 7
  # P. cinereus number 11
  # P. cinereus number 12
  # P. cylindraceus number 10
  # P. glut number 6
  # P. montanus number 3
  # P. montanus number 7
Bsal <- Bsal[-which(Bsal$Species == "fuscus" & Bsal$ID_number == 9),]
Bsal <- Bsal[-which(Bsal$Species == "viridescens" & Bsal$ID_number == 10),]
Bsal <- Bsal[-which(Bsal$Species == "viridescens" & Bsal$ID_number == 12),]
Bsal <- Bsal[-which(Bsal$Species == "cinereus" & Bsal$ID_number == 4),]
Bsal <- Bsal[-which(Bsal$Species == "cinereus" & Bsal$ID_number == 7),]
Bsal <- Bsal[-which(Bsal$Species == "cinereus" & Bsal$ID_number == 11),]
Bsal <- Bsal[-which(Bsal$Species == "cinereus" & Bsal$ID_number == 12),]
Bsal <- Bsal[-which(Bsal$Species == "cylindraceus" & Bsal$ID_number == 6),]
Bsal <- Bsal[-which(Bsal$Species == "glutinosus" & Bsal$ID_number == 6),]
Bsal <- Bsal[-which(Bsal$Species == "montanus" & Bsal$ID_number == 3),]
Bsal <- Bsal[-which(Bsal$Species == "montanus" & Bsal$ID_number == 7),]


# Drop the levels for the salamanders that never got infected
Bsal <- droplevels(Bsal)



# Make a colum with the genus and species names
Bsal$Genus_species <- paste(Bsal$Genus, Bsal$Species, sep = " ")


# Remove entries with NA as Genus or Experimental day
Bsal <- Bsal[is.na(Bsal$Genus) == FALSE,]
Bsal <- Bsal[is.na(Bsal$Exp_day) == FALSE,]


# Add a column with log10(Bd load at t - 1) 
Bsal$log10_Bd_t_minus_1 <- log10(Bsal$Bd_load_t_min_1 + 1)


# Assign unique ID numbers that do not match across species

Bsal$ID_new <- paste(Bsal$Species, Bsal$ID_number, sep = "_")
Bsal$ID_new <- as.numeric(as.factor(Bsal$ID_new))

# Number of unique individuals
length(unique(Bsal$ID_new))


# Add a unique species ID
Bsal$Species_ID <- with(Bsal, paste(Species, ID_new, sep = "_"))




# 4. Bd  summary ------------------------------------------------



# Determine the number of individuls that tested positive for Bd at least once
# To do this, first summarize the data by genus, species, ID_number, and tally the number of Bd positives, Bd negatives, and total Bd samples
Bd_tally <- ddply(.data = Bsal, .variable = c("Genus", "Species", "ID_number"), .fun = summarize, 
      Bd_pos = length(which(Bd_load > 0)), 
      Bd_neg = length(which(Bd_load == 0)),
      sum = sum(Num))

# Convert Bd_pos to 0 and 1
Bd_tally$Bd_pos[Bd_tally$Bd_pos > 0] <- 1

# Sum up the total number of individuals with at least 1 positive
sum(Bd_tally$Bd_pos)

# Total number of individuals
nrow(Bd_tally)

# Fraction of individuals that tested Bd positive at least once:
sum(Bd_tally$Bd_pos)/ nrow(Bd_tally)




# 5. Figure S1. Bd load over time ------------------------------------------------



# Make a color color-blind friendly palette (15 colors)
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","goldenrod3")

# Bd load
ggplot(data = Bsal, aes(y = log10(Bd_load+1), x = Exp_day, col = as.factor(ID_number)))+
  facet_wrap(~ Genus_species)+ 
  geom_line(size = 1.2)+
  geom_point(size = 2)+
  scale_color_manual(values = pal)+
  geom_vline(xintercept = 42, col = "red", lwd = 1, lty = 2)+
  xlab("Experimental day")+
  ylab(expression(paste(log[10], "(", italic(Bd), " infection intensity + 1)", sep = "")))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 20, color = "black"), 
        axis.text.y = element_text(size = 20, color = "black"), 
        legend.text = element_text(size=20),
        axis.title.y = element_text(size = 24, color = "black"), 
        axis.title.x = element_text(size = 24, color = "black"),
        strip.text = element_text(size=24, face = "italic"),
        legend.title = element_text(size=24),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white")) +
  labs(col = "ID\nNumber")

# Figure S1.
# Save the plot
ggsave("./Figs_2020/FigS1.pdf", height = 10, width = 17)





# 6. Bsal exposure summary ------------------------------------------------



# Calculate the total number of indiviudals that were inoculated with Bsal and tested positive following exposure
Bsal_tally <- ddply(.data = Bsal, 
                .variable = c("Genus", "Species", "Treatment", "ID_number"), 
                .fun = summarize, 
                Bsal_pos = length(which(Bsal_load > 0)), 
                Bsal_neg = length(which(Bsal_load == 0)), 
                sum = sum(Num))

# Note that Bsal_tally has > 98 rows because some individuals are double counted that were used in the I treatment in the 1st phase, and then switched to D treatment for the second phase
  # To remedy this - we will remove the C treatment
  # And then consolidate

# Remove the control individuals
Bsal_tally <- Bsal_tally[-which(Bsal_tally$Treatment == "C"),]

# Consolidate across ID_number
Bsal_tally2 <- ddply(.data = Bsal_tally, 
              .variable = c("Genus", "Species", "ID_number"), 
              .fun = summarize, 
              Bsal_pos = sum(Bsal_pos), 
              Bsal_neg = sum(Bsal_neg))

# Convert Bsal_pos to 0 and 1
Bsal_tally2$Bsal_pos[Bsal_tally2$Bsal_pos > 0] <- 1

# Sum up the total number of individuals with at least 1 positive
sum(Bsal_tally2$Bsal_pos)
  # And 3 individuals were removed because they were always negative 




# 7. Day 3 and 4 averages ------------------------------------------------




# Average Bsal load on day 3 and 4
mean(c(Bsal[Bsal$Exp_day == 3 & Bsal$Treatment == "I",]$Bsal_load,
       Bsal[Bsal$Exp_day == 4 & Bsal$Treatment == "I",]$Bsal_load))

# Standard error of Bsal load on day 3 and 4
sd(c(Bsal[Bsal$Exp_day == 3 & Bsal$Treatment == "I",]$Bsal_load,
     Bsal[Bsal$Exp_day == 4& Bsal$Treatment == "I",]$Bsal_load))/ sqrt(length(c(Bsal[Bsal$Exp_day == 3 & Bsal$Treatment == "I",]$Bsal_load,
                                                                                Bsal[Bsal$Exp_day == 4 & Bsal$Treatment == "I",]$Bsal_load)))

# Number of salamanders
length(c(Bsal[Bsal$Exp_day == 3 & Bsal$Treatment == "I",]$Bsal_load,
         Bsal[Bsal$Exp_day == 4& Bsal$Treatment == "I",]$Bsal_load))




# Average Bd load on day 3 and 4
mean(c(Bsal[Bsal$Exp_day == 3 ,]$Bd_load,
       Bsal[Bsal$Exp_day == 4,]$Bd_load), na.rm = TRUE)

# Standard error of Bd load on day 3 and 4
sd(c(Bsal[Bsal$Exp_day == 3,]$Bd_load,
     Bsal[Bsal$Exp_day == 4,]$Bd_load), na.rm = TRUE)/ sqrt(length(c(Bsal[Bsal$Exp_day == 3,],
                                                                    Bsal[Bsal$Exp_day == 4,])))




# 8. Last swabbing event - day 81 - summary ------------------------------------------------



# Average Bd load on day 81
mean(c(Bsal[Bsal$Exp_day == 81,]$Bd_load), na.rm = TRUE)

# Standard error of Bd load day 81
sd(c(Bsal[Bsal$Exp_day == 81,]$Bd_load), na.rm = TRUE)/ sqrt(length(c(Bsal[Bsal$Exp_day == 81,]$Bd_load)))




# Subset the Bsal exposed treatments
Bsal_last <- Bsal[Bsal$Exp_day == 81 & Bsal$Treatment != "C",]
Bsal_last <- Bsal_last[is.na(Bsal_last$Bsal_load) == FALSE,]

# Average Bsal load on day 81 (last swabbing day of the experiment)
mean(Bsal_last$Bsal_load, na.rm = TRUE)

# Standard error of Bsal load on day 81 (last swabbing day of the experiment)
sd(Bsal_last$Bsal_load)/ sqrt(length(Bsal_last$Bsal_load))






# 9. Co-infection summary ------------------------------------------------



# Calculate the number of co-infected individuals

# Was Bsal ever present on the individual? Determine if Bsal load is > 1
Bsal$Bsal_PA <- ifelse(Bsal$Bsal_load > 0, 1, 0)

# Was Bd ever present on the individual? Determine if Bd load is > 1
Bsal$Bd_PA <- ifelse(Bsal$Bd_load > 0, 1, 0)

# Add up the 2 columns
  # If value = 1 -> single infection
  # If value = 2 -> co-infected
Bsal$coInfection <- Bsal$Bsal_PA + Bsal$Bd_PA

# Determine which individuals were positive for Bd and Bsal at the SAME time
coinf <- ddply(.data = Bsal, .variable = c("Genus", "Species", "ID_new", "Treatment"), 
               .fun = summarize, 
               coinfected = max(coInfection, na.rm = TRUE))
  # Again= the number of rows is > the number of individuals, because of the double counting with the 2nd phase of the experiment

# We will remove the controls - and then consolidate

# Remove the control individuals
coinf <- coinf[-which(coinf$Treatment == "C"),]

# Consolidate across ID_number
coinf2 <- ddply(.data = coinf, 
                     .variable = c("Genus", "Species", "ID_new"), 
                     .fun = summarize, 
                       coinfected = max(coinfected))

# Total number of co-infected individuals
length(which(coinf2$coinfected == 2))

# Total number of individuals
nrow(coinf2)

# Percent of individuals
length(which(coinf2$coinfected == 2))/nrow(coinf2)



# 10. Body condition calculation------------------------------------------------



# Fix typo in Mass column - entered as 0.29 and suppose to be 1.29
Bsal[Bsal$Species == "cinereus" & Bsal$Exp_day == 28 & Bsal$ID_number == 3,]$Mass <- 1.29



# Not used in the paper - just used to visualize data
## Use the data from the day animals were collected in the field
#ggplot(data = field, aes(y = Mass, x = SVL, col = Species))+
#  # facet_wrap(~ Species)+
#  geom_point()+
#  geom_smooth()


# Log transform Mass and SVL
field$ln_mass <- log(field$Mass)
field$ln_SVL <- log(field$SVL)


# Fit a linear model
# We use a SVL*Species interactions- so that all species have their own intercept & slope
SMI <- lm(ln_mass ~ ln_SVL*(Species-1), data = field)


# Log transform Mass and SVL
Bsal$ln_mass <- log(Bsal$Mass)
Bsal$ln_SVL <- log(Bsal$SVL)

# Estimate the prediction from the linear model for each individual
Bsal$BC_pred <- predict(SMI, Bsal, type = "response")

# Residuals = observed mass - predicted mass
Bsal$BC <- Bsal$BC_pred - Bsal$ln_mass


# Look at the output
summary(SMI)






# 11. Survival analysis ------------------------------------------------



# Broken down into 2 parts
  # First & second inoculations



# Make a list of all the individuals in the dataset
# First phase
first_part <- Bsal[Bsal$Exp_day < 40, ]

# Create a data frame to hold information
census_p1 <- data.frame(Species_ID = unique(first_part$Species_ID),
                        Census = NA,
                        Exp_day = NA,
                        Treatment = NA)

# Split up Species ID column
census_p1$Species <- gsub("_[0-9][0-9]", "", census_p1$Species_ID)
census_p1$Species <- gsub("_[0-9]", "", census_p1$Species)


# Now go through the data file and pull out the last record for each individual
for(i in 1:nrow(census_p1)){
  
  # Subset the data for each individual
  dat.sub <- first_part[first_part$Species_ID == census_p1$Species_ID[i],]

  # Determine last entry for that indiviudal
  last.entry <- which(dat.sub$Exp_day == max(dat.sub$Exp_day))
  
  census_p1$Exp_day[i] <- dat.sub[last.entry,]$Exp_day
  census_p1$Census[i] <- dat.sub[last.entry,]$Census
  census_p1$Treatment[i] <- as.character(dat.sub[last.entry,]$Treatment)
  
}

# Make treatment a factor
census_p1$Treatment <- as.factor(census_p1$Treatment)

# Remove data from control individuals
census_p1 <- census_p1[census_p1$Treatment != "C",]

# Average days until death for N. vir (and min and max ranges) - note that this table only goes up to experimental day 39
ddply(.data = census_p1,
      .variables = "Species",
      .fun = summarize,
      days_to_death = mean(Exp_day),
      se = sd(Exp_day)/sqrt(length(Exp_day)),
      min = min(Exp_day),
      max = max(Exp_day))




# Second phase
second_part <- Bsal[Bsal$Exp_day > 40 & Bsal$Exp_day < 82, ]

# Create a data frame to hold information
census_p2 <- data.frame(Species_ID = unique(second_part$Species_ID),
                        Census = NA,
                        Exp_day = NA,
                        Treatment = NA)

# Split up Species ID column
census_p2$Species <- gsub("_[0-9][0-9]", "", census_p2$Species_ID)
census_p2$Species <- gsub("_[0-9]", "", census_p2$Species)


# Now go through the data file and pull out the last record for each individual
for(i in 1:nrow(census_p2)){
  
  # Subset the data for each individual
  dat.sub <- second_part[second_part$Species_ID == census_p2$Species_ID[i],]
  
  # Determine last entry for that indiviudal
  last.entry <- which(dat.sub$Exp_day == max(dat.sub$Exp_day))
  
  census_p2$Exp_day[i] <- dat.sub[last.entry,]$Exp_day
  census_p2$Census[i] <- dat.sub[last.entry,]$Census
  census_p2$Treatment[i] <- as.character(dat.sub[last.entry,]$Treatment)
  
}

# Make treatment a factor
census_p2$Treatment <- as.factor(census_p2$Treatment)

# Remove data from control individuals
census_p2 <- census_p2[census_p2$Treatment != "C",]


# First phase
# Kaplan Meier Analysis
# All species together
sal.surv.p1 <- survfit(Surv(Exp_day, Census) ~ Species, data = census_p1) 


# Second phase
# Kaplan Meier Analysis
# Here- each species is pulled out and the figure is made seperately to compare between I and D treatments
census_p2_cinereus <- census_p2[census_p2$Species == "cinereus",]
census_p2_cylindraceus <- census_p2[census_p2$Species == "cylindraceus",]
census_p2_montanus <- census_p2[census_p2$Species == "montanus",]
census_p2_fuscus <- census_p2[census_p2$Species == "fuscus",]
census_p2_wrighti <- census_p2[census_p2$Species == "wrighti",]
census_p2_wilderae <- census_p2[census_p2$Species == "wilderae",]
census_p2_glutinosus <- census_p2[census_p2$Species == "glutinosus",]
census_p2_viridescens <- census_p2[census_p2$Species == "viridescens",]


sal.surv.p2.cin <- survfit(Surv(Exp_day, Census) ~ Treatment, data = census_p2_cinereus) 
sal.surv.p2.cyl <- survfit(Surv(Exp_day, Census) ~ Treatment, data = census_p2_cylindraceus) 
sal.surv.p2.mon <- survfit(Surv(Exp_day, Census) ~ Treatment, data = census_p2_montanus) 
sal.surv.p2.fus <- survfit(Surv(Exp_day, Census) ~ Treatment, data = census_p2_fuscus) 
sal.surv.p2.wri <- survfit(Surv(Exp_day, Census) ~ Treatment, data = census_p2_wrighti) 
sal.surv.p2.wil <- survfit(Surv(Exp_day, Census) ~ Treatment, data = census_p2_wilderae) 
sal.surv.p2.glu <- survfit(Surv(Exp_day, Census) ~ Treatment, data = census_p2_glutinosus) 
sal.surv.p2.vir <- survfit(Surv(Exp_day, Census) ~ Treatment, data = census_p2_viridescens) 


# First phase
pdf("./Figs_2020/FigS2.pdf",
    height = 6,
    width =7)
plot(sal.surv.p1, 
     lwd = c(2, 2),
     las = 1,
     xlab = "Experimental days",
     ylab = "Proportion alive",
     col = pal[1:8])
legend("bottomleft",
       c(expression(italic("P. cinereus")),
         expression(italic("P. cylindraceus")),
         expression(italic("D. fuscus")),
         expression(italic("P. glutinosus")),
         expression(italic("P. montanus")),
         expression(italic("N. viridescens")),
         expression(italic("E. wilderae")),
         expression(italic("D. wrighti"))),
       col = pal[1:8],
       lwd = c(2, 2)
)
dev.off()

# Second phase
pdf("./Figs_2020/FigS3.pdf",
    height = 10,
    width =5)
par(mfrow = c(4, 2))
plot(sal.surv.p2.cin, 
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(42, 84),
     lwd = c(2, 2),
     las = 1,
     xlab = "Experimental days",
     ylab = "Proportion alive",
     col = c("black", "deepskyblue3"))
legend("bottomright",
       c("D", "I"),
  col = c("black", "deepskyblue3"),
  lwd = c(2, 2)
)
text(x = 65, y = 0.8, cex = 1.5, expression(italic("P. cinereus")))


plot(sal.surv.p2.cyl, 
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(42, 84),
     lwd = c(2, 2),
     las = 1,
     xlab = "Experimental days",
     ylab = "Proportion alive",
     col = c("black", "deepskyblue3"))
legend("bottomright",
       c("D", "I"),
       col = c("black", "deepskyblue3"),
       lwd = c(2, 2)
)
text(x = 65, y = 0.8, cex = 1.5, expression(italic("P. cylindraceus")))




plot(sal.surv.p2.glu, 
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(42, 84),
     lwd = c(2, 2),
     las = 1,
     xlab = "Experimental days",
     ylab = "Proportion alive",
     col = c("black", "deepskyblue3"))
legend("bottomright",
       c("D", "I"),
       col = c("black", "deepskyblue3"),
       lwd = c(2, 2)
)
text(x = 65, y = 0.7, cex = 1.5, expression(italic("P. glutinosus")))






plot(sal.surv.p2.mon, 
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(42, 84),
     lwd = c(2, 2),
     las = 1,
     xlab = "Experimental days",
     ylab = "Proportion alive",
     col = c("black", "deepskyblue3"))
legend("bottomright",
       c("D", "I"),
       col = c("black", "deepskyblue3"),
       lwd = c(2, 2)
)
text(x = 65, y = 0.7, cex = 1.5, expression(italic("P. montanus")))




plot(sal.surv.p2.fus, 
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(42, 84),
     lwd = c(2, 2),
     las = 1,
     xlab = "Experimental days",
     ylab = "Proportion alive",
     col = c("black", "deepskyblue3"))
legend("bottomright",
       c("D", "I"),
       col = c("black", "deepskyblue3"),
       lwd = c(2, 2)
)
text(x = 50, y = 0.8, cex = 1.5, expression(italic("D. fuscus")))



plot(sal.surv.p2.wri, 
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(42, 84),
     lwd = c(2, 2),
     las = 1,
     xlab = "Experimental days",
     ylab = "Proportion alive",
     col = c("black", "deepskyblue3"))
legend("bottomright",
       c("D", "I"),
       col = c("black", "deepskyblue3"),
       lwd = c(2, 2)
)
text(x = 50, y = 0.8, cex = 1.5, expression(italic("D. wrighti")))



plot(sal.surv.p2.wil, 
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlim = c(42, 84),
     lwd = c(2, 2),
     las = 1,
     xlab = "Experimental days",
     ylab = "Proportion alive",
     col = c("black", "deepskyblue3"))
legend("bottomright",
       c("D", "I"),
       col = c("black", "deepskyblue3"),
       lwd = c(2, 2)
)
text(x = 65, y = 0.8, cex = 1.5, expression(italic("E. wilderae")))



plot(sal.surv.p2.vir, 
     cex.lab = 1.5,
     cex.axis = 1.5,
     lwd = c(2, 2),
     xlim = c(42, 84),
     las = 1,
     xlab = "Experimental days",
     ylab = "Proportion alive",
     col = c( "deepskyblue3"))
legend("bottomright",
       c("I"),
       col = c("deepskyblue3"),
       lwd = c(2, 2)
)
text(x = 65, y = 0.8, cex = 1.5, expression(italic("N. viridescens")))

dev.off()






# 12. Determine the total number of Plethodon that died & their average infection intensity ------------------------



# Remove newt data
census_p1_noNewts <- census_p1[census_p1$Species != "viridescens",]
census_p2_noNewts <- census_p2[census_p2$Species != "viridescens",]


# sum all deaths
sum(census_p1_noNewts[,2]) + sum(census_p2_noNewts[,2])

# Determine total number of individuals infected
nrow(census_p1_noNewts)



# Determine average load of Plethodon individuals
# Remove the newts
Bsal_noNewts <- Bsal[Bsal$Species != "viridescens",]

# Average Bsal load during entirety of experiment
mean(Bsal_noNewts$Bsal_load, na.rm  = TRUE)





# 13. First phase: data format ------------------------------------------------




# Then standardize experimental day
# Mean of experimental day
mE_first <- mean(first_part$Exp_day, na.rm = TRUE)

# SD of experimental day
sD_first <- sd(first_part$Exp_day, na.rm = TRUE)

# Standadization of experiemntal day
first_part$standard_day <- (first_part$Exp_day - mE_first)/sD_first


# Save the data into a new object before removing the controls
first_part_all <- first_part

# Remove all Control individuals
first_part <- first_part[first_part$Treatment != "C",]

# Remove control individuals
first_part <- droplevels(first_part)





# 14. First phase: Bsal infection intensity analysis, Supporting Information 1: Parts 1-4  ------------------------------------------------




# Need to remove data with missing observations for analysis
first_part.v1 <- first_part[is.na(first_part$Bsal_load) == FALSE,]
first_part.v1 <- first_part.v1[is.na(first_part.v1$log10_Bd_t_minus_1) == FALSE,]
first_part.v1 <- first_part.v1[first_part.v1$Exp_day > 0 ,]

# Fit a linear mixed effect model with ID number as a random effect
mod <- lme(log10(Bsal_load+1) ~ 
              standard_day*Species + log10_Bd_t_minus_1,
           random =~ 1|ID_new, 
            data = first_part.v1)

# Supporting Information S1
# Part 1 of Supporting Information 
# Look at the summary
summary(mod)

# Part 2 of Supporting Information 
intervals(mod)

# Part 3
# Determine if the interaction term is significant
anova(mod)

# Part 4 of Supporting Information 
# Compare the slopes of Bsal over time among species
lst = emmeans::emtrends(mod, ~Species, var= "standard_day")
slopes <- multcomp::cld(lst)
slopes




# 15. Diagnotic plots, Supporting Information 1: Part 5 ------------------------------------------------



# Create a dataframe with fitted values and residuals
df <- data.frame(fitted = predict(mod), 
                 res = residuals(mod))

pdf(file = "./Figs_2020/S1_diagnositic_1.pdf",
    height = 5, 
    width = 6)
# Fitted vs residual plot
plot(df$fitted, df$res, 
     ylab = "Residuals", xlab = "Fitted", 
     pch = 21, col = "black", bg = "deepskyblue3", 
     cex = 2, las = 1, cex.axis = 1.5, cex.lab = 1.5)
abline(h = 0, lty = 2)
dev.off()


pdf(file = "./Figs_2020/S1_diagnositic_2.pdf",
    height = 5, 
    width = 6)
# QQ plot
qqnorm(mod, pch = 21, 
       col = "black", fill = "deepskyblue3", 
       cex = 2, las = 1, 
       cex.axis = 1.5, cex.lab = 1.5, 
       main = "", abline = c(0, 1), lty = 2)
dev.off()





# 16. Fig 3: 1st phase - Bsal infection intensity vs. time ------------------------------------------------




# Create new dataset
plot_data <- data.frame(
            standard_day = first_part.v1$standard_day, 
            Exp_day = first_part.v1$Exp_day, 
            Species = first_part.v1$Species, 
            log10_Bd_t_minus_1 = mean(first_part.v1$log10_Bd_t_minus_1, na.rm = TRUE),
            ID_new = first_part.v1$ID_new,
            Genus_species = first_part.v1$Genus_species
)


# Create predictions for each column
plot_data$predictions <- predict(mod, plot_data, re.form = ~ (1|ID_new), allow.new.levels = TRUE)


# Figure 3 main text
# Make the plot
ggplot(data = first_part.v1, aes(y = log10(Bsal_load+1), x = Exp_day)) + 
  geom_point(col= "gray") + 
  geom_line(aes(group = as.factor(ID_new)), col = "gray")+
  geom_line(data = plot_data, aes(y = predictions, 
                                  x = Exp_day, col = as.factor(ID_new)), size = 0.5) + 
  facet_wrap(~Genus_species)+
  ylab(expression(paste(log[10], "(", italic(Bsal), " infection intensity + 1)", sep = "")))+
  xlab("Experimental day") +
  scale_color_manual(values = rep("black", times = length(unique(first_part$ID_new))))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.text.y = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.title.x = element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=17, face = "italic"),
        strip.background = element_rect(fill = "white"))

ggsave("./Figs_2020/Fig3.pdf", height = 8, width = 10)




# 17. Fig 5. Influence of Bd at t - 1 on Bsal at time t ------------------------------------------------




# Create new dataset
plot_data <- with(first_part.v1, expand.grid(
  standard_day = mean(standard_day, na.rm = TRUE), 
  Species = levels(Species), 
  log10_Bd_t_minus_1 = seq(min(log10_Bd_t_minus_1, na.rm = TRUE), 
                           max(log10_Bd_t_minus_1, na.rm = TRUE), length = 10),
 ID_new = unique(ID_new)
))

# Mean prediction
plot_data$Bsal_load <- predict(mod, plot_data)


#create design matrix
Designmat <- model.matrix(eval(eval(mod$call$fixed)[-2]), plot_data[-ncol(plot_data)])

#compute standard error for predictions
predvar <- diag(Designmat %*% mod$varFix %*% t(Designmat))
plot_data$SE <- sqrt(predvar) 
plot_data$SE2 <- sqrt(predvar+mod$sigma^2)

# Calculate lower and upper
plot_data$lower <- plot_data$Bsal_load - (2*plot_data$SE)
plot_data$upper <- plot_data$Bsal_load + (2*plot_data$SE)

# Subset the data to 1 species
plot_data.1 <- plot_data[plot_data$Species == levels(plot_data$Species)[1],]
plot_data.1 <- droplevels(plot_data.1)
plot_data.1 <- plot_data.1[plot_data.1$ID_new == unique(plot_data.1$ID_new)[1],]


# Fig 5. in main text
# Make plot
ggplot() + 
  geom_jitter(data = first_part.v1, aes(y = log10(Bsal_load+1), x = log10_Bd_t_minus_1), width = 0.05, height = 0.5) + 
  geom_line(data = plot_data.1, aes(y = Bsal_load, 
                                 x = log10_Bd_t_minus_1), size = 0.5) + 
  geom_ribbon(data = plot_data.1, aes(ymin = lower, 
                                    ymax = upper,
                                    x = log10_Bd_t_minus_1), alpha = 0.3)+ 
  ylab(expression(paste(italic(Bsal), " infection intensity at time ", italic(t), sep = "")))+
  xlab(expression(paste(italic(Bd), " infection intensity at time ", italic(t - 1), sep = ""))) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                     labels = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000))+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                     labels = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.text.y = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.title.x = element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text.align = 0)

ggsave("./Figs_2020/Fig5.pdf", height = 6, width = 8)







# 18. First phase: Body condition analysis, Supporting Information 2: Parts 1 - 4 ------------------------------------------------



# Remove missing observations
first_part_all.v1 <- first_part_all[is.na(first_part_all$BC) == FALSE,]


# Fit the model
mod.BC.1 <- lme(BC ~ 
                 Treatment*standard_day*Species,
                   random = ~1|ID_new, 
               data = first_part_all.v1,
               control=glmerControl(optimizer="bobyqa",
                                   optCtrl=list(maxfun=1e9)))

# Supporting Information S2
# Part 1 of Supporting Information 
# Look at the summary
summary(mod.BC.1)

# Part 2 of Supporting Information 
intervals(mod.BC.1)

# Part 3
# Determine if the interaction term is significant
anova(mod.BC.1)

# Part 4 of Supporting Information 
# Compare the slopes of Bsal over time among species
lst = emmeans::emtrends(mod.BC.1, ~Species, var= "standard_day")
slopes <- multcomp::cld(lst)
slopes






# 19. Fig S3: Body condition vs. time ------------------------------------------------




# Create predictions for each column
first_part_all.v1$predictions <- predict(mod.BC.1, 
                                 first_part_all.v1, 
                                 re.form = ~ (1|ID_new), 
                                 allow.new.levels = TRUE)

# Figure S4
# Make the plot
ggplot(data = first_part_all.v1, aes(y = BC, x = Exp_day)) + 
 geom_point(col= "gray") + 
 geom_line(aes(group = as.factor(ID_new)), col = "gray")+
  geom_line(data = first_part_all.v1, aes(y = predictions, 
                                  x = Exp_day, 
                                  col = as.factor(ID_new)), size = 0.5) + 
  facet_wrap(~Genus_species)+
  ylab("Body condition")+
  xlab("Experimental Day") +
  scale_color_manual(values = rep("black", times = length(unique(first_part_all.v1$ID_new))))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"), 
        axis.title.x = element_text(size = 18, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=17, face = "italic"),
        strip.background = element_rect(fill = "white"))

ggsave("./Figs_2020/FigS4.pdf", height = 8, width = 10)






# 20. Diagnotic plots, Supporting Information 2: Part 5 ------------------------------------------------




# Create a dataframe with fitted values and residuals
df <- data.frame(fitted = predict(mod.BC.1), 
                 res = residuals(mod.BC.1))

pdf(file = "./Figs_2020/S2_diagnositic_1.pdf",
    height = 5, 
    width = 6)
# Fitted vs residual plot
plot(df$fitted, df$res, 
     ylab = "Residuals", xlab = "Fitted", 
     pch = 21, col = "black", bg = "deepskyblue3", 
     cex = 2, las = 1, cex.axis = 1.5, cex.lab = 1.5)
abline(h = 0, lty = 2)
dev.off()


pdf(file = "./Figs_2020/S2_diagnositic_2.pdf",
    height = 5, 
    width = 6)
# QQ plot
qqnorm(mod.BC.1, pch = 21, 
       col = "black", fill = "deepskyblue3", 
       cex = 2, las = 1, 
       cex.axis = 1.5, cex.lab = 1.5, 
       main = "", abline = c(0, 1), lty = 2)
dev.off()








# 21. Second phase: data format ------------------------------------------------




# Then standardize experimental day
# Mean of experimental day
mE_second <- mean(second_part$Exp_day, na.rm = TRUE)

# SD of experimental day
sD_second <- sd(second_part$Exp_day, na.rm = TRUE)

# Standadization of experiemntal day
second_part$standard_day <- (second_part$Exp_day - mE_second)/sD_second


# Save the data into a new object before removing the controls
second_part_all <- second_part

# Remove all Control individuals
second_part <- second_part[second_part$Treatment != "C",]

# Remove control individuals
second_part <- droplevels(second_part)






# 22. Second phase: Bsal infection intensity analysis, Supporting Information 3: Parts 1-4  ------------------------------------------------




# Need to remove data with missing observations for analysis
second_part.v1 <- second_part[is.na(second_part$Bsal_load) == FALSE,]
second_part.v1 <- second_part.v1[is.na(second_part.v1$log10_Bd_t_minus_1) == FALSE,]

# Remove the 1 newt
second_part.v1 <- second_part.v1[second_part.v1$Species != "viridescens",]

# Remove 1 wilderae
second_part.v1 <- second_part.v1[second_part.v1$Species != "wilderae",]

# Fit a linear mixed effect model with ID number as a random effect
mod.2 <- lme(log10(Bsal_load+1) ~ 
             standard_day*Species*Treatment + log10_Bd_t_minus_1,
           random =~ 1|ID_new, 
           data = second_part.v1)

# Supporting Information S3
# Part 1 of Supporting Information 
# Look at the summary
summary(mod.2)

# Part 2 of Supporting Information 
intervals(mod.2)

# Part 3
# Determine if the interaction term is significant
anova(mod.2)

# Part 4 of Supporting Information 
# Compare the slopes of Bsal over time among species
lst = emmeans::emtrends(mod.2, ~Treatment, var= "standard_day")
slopes <- multcomp::cld(lst)
slopes






# 23. Diagnotic plots, Supporting Information 3: Part 5 ------------------------------------------------



# Create a dataframe with fitted values and residuals
df.2 <- data.frame(fitted = predict(mod.2), 
                 res = residuals(mod.2))

pdf(file = "./Figs_2020/S3_diagnositic_1.pdf",
    height = 5, 
    width = 6)
# Fitted vs residual plot
plot(df.2$fitted, df.2$res, 
     ylab = "Residuals", xlab = "Fitted", 
     pch = 21, col = "black", bg = "deepskyblue3", 
     cex = 2, las = 1, cex.axis = 1.5, cex.lab = 1.5)
abline(h = 0, lty = 2)
dev.off()


pdf(file = "./Figs_2020/S3_diagnositic_2.pdf",
    height = 5, 
    width = 6)
# QQ plot
qqnorm(mod.2, pch = 21, 
       col = "black", fill = "deepskyblue3", 
       cex = 2, las = 1, 
       cex.axis = 1.5, cex.lab = 1.5, 
       main = "", abline = c(0, 1), lty = 2)
dev.off()





# 24. Fig 4: Bsal infection intensity vs. time ------------------------------------------------




# Create new dataset
plot_data <- data.frame(
  standard_day = second_part.v1$standard_day, 
  Exp_day = second_part.v1$Exp_day, 
  Species = second_part.v1$Species, 
  log10_Bd_t_minus_1 = mean(second_part.v1$log10_Bd_t_minus_1, na.rm = TRUE),
  ID_new = second_part.v1$ID_new,
  Treatment = second_part.v1$Treatment,
  Genus_species = second_part.v1$Genus_species
)

nrow(plot_data)
nrow(second_part.v1)

# Create predictions for each column
plot_data$predictions <- predict(mod.2, plot_data)

# Rename the treatment columns
plot_data$trt <- ifelse(plot_data$Treatment == "D", "Double exposure", "Single exposure")
second_part.v1$trt <- ifelse(second_part.v1$Treatment == "D", "Double exposure", "Single exposure")

# Make the single exposure first
plot_data$trt <- factor(plot_data$trt , levels = c("Single exposure", "Double exposure"))
second_part.v1$trt <- factor(second_part.v1$trt , levels = c("Single exposure", "Double exposure"))



# Figure 4 main text
# Make the plot
ggplot(data = second_part.v1, aes(y = log10(Bsal_load+1), x = Exp_day)) + 
  geom_point( col= "gray") + 
  geom_line(aes(group = as.factor(ID_new)), col = "gray")+
  geom_line(data = plot_data, aes(y = predictions, 
                                  x = Exp_day, col = as.factor(ID_new)), size = 0.5) + 
  facet_grid(Genus_species~trt)+
  ylab(expression(paste(log[10], "(", italic(Bsal), " infection intensity + 1)", sep = "")))+
  xlab("Experimental day") +
  scale_color_manual(values = rep("black", times = length(unique(second_part.v1$ID_new))))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.text.y = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 20, color = "black"), 
        axis.title.x = element_text(size = 20, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=12, face = "italic"),
        strip.background = element_rect(fill = "white"))

ggsave("./Figs_2020/Fig4.pdf", height = 12, width = 6)






# 25. Second phase: Body condition analysis, Supporting Information 4: Parts 1-4  ------------------------------------------------



# Remove missing observations
second_part_all.v1 <- second_part_all[is.na(second_part_all$BC) == FALSE,]
# Remove the 1 newt
second_part_all.v1 <- second_part_all.v1[second_part_all.v1$Species != "viridescens",]

# Remove 1 wilderae
second_part_all.v1 <- second_part_all.v1[second_part_all.v1$Species != "wilderae",]


# Fit the model
mod.BC.2 <- lme(BC ~ 
                   Treatment*standard_day*Species,
                 random = ~1|ID_new, 
                 data = second_part_all.v1,
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=1e9)))


# Supporting Information S4
# Part 1 of Supporting Information 
# Look at the summary
summary(mod.BC.2)

# Part 2 of Supporting Information 
intervals(mod.BC.2)

# Part 3
# Determine if the interaction term is significant
anova(mod.BC.2)

# Part 4 of Supporting Information 
# Compare the slopes of Bsal over time among species
lst = emmeans::emtrends(mod.BC.2, ~Species, var= "standard_day")
slopes <- multcomp::cld(lst)
slopes






# 26. Fig S5: Body condition vs. time ------------------------------------------------


# Droplevels
second_part_all.v1 <- droplevels(second_part_all.v1)

# Create predictions for each column
second_part_all.v1$predictions <- predict(mod.BC.2, 
                                 second_part_all.v1, 
                                 re.form = ~ (1|ID_new), 
                                 allow.new.levels = TRUE)

# Figure S5
# Make the plot
ggplot(data = second_part_all.v1, aes(y = BC, x = Exp_day)) + 
  geom_point(col= "gray") + 
  geom_line(aes(group = as.factor(ID_new)), col = "gray")+
  geom_line(data = second_part_all.v1, aes(y = predictions, 
                                  x = Exp_day, 
                                  col = as.factor(ID_new)), size = 0.5) + 
  facet_grid(~Genus_species)+
  ylab("Body condition")+
  xlab("Experimental day") +
  scale_color_manual(values = rep("black", times = length(unique(second_part_all.v1$ID_new))))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 16, color = "black"), 
        axis.text.y = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"), 
        axis.title.x = element_text(size = 18, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text = element_text(size=10, face = "italic"),
        strip.background = element_rect(fill = "white"))

ggsave("./Figs_2020/FigS5.pdf", height = 5, width = 10)




# 27. Diagnotic plots, Supporting Information 4: Part 5 ------------------------------------------------




# Create a dataframe with fitted values and residuals
df <- data.frame(fitted = predict(mod.BC.2), 
                 res = residuals(mod.BC.2))

pdf(file = "./Figs_2020/S4_diagnositic_1.pdf",
    height = 5, 
    width = 6)
# Fitted vs residual plot
plot(df$fitted, df$res, 
     ylab = "Residuals", xlab = "Fitted", 
     pch = 21, col = "black", bg = "deepskyblue3", 
     cex = 2, las = 1, cex.axis = 1.5, cex.lab = 1.5)
abline(h = 0, lty = 2)
dev.off()


pdf(file = "./Figs_2020/S4_diagnositic_2.pdf",
    height = 5, 
    width = 6)
# QQ plot
qqnorm(mod.BC.2, pch = 21, 
       col = "black", fill = "deepskyblue3", 
       cex = 2, las = 1, 
       cex.axis = 1.5, cex.lab = 1.5, 
       main = "", abline = c(0, 1), lty = 2)
dev.off()







# 28. Table 1 & 2 Main text ------------------------------------------------



# Read in histopathology data
histo <- read.csv("./Data/Bsal Histopath.csv")

# Switch the order of I and D
histo$Treatment <- factor(histo$Treatment, levels = c("C", "I", "D"))

# Summary table of sample sizes
hist_tab <- ddply(.data = histo, 
      .variable = c("Species", "Treatment"), 
      .fun = summarize, 
      samp_size = length(unique(KRL..)),
      Bsal_pos = sum(Bsal.type, na.rm = TRUE),
      Bsal_neg = length(which(Bsal.type == 0)),
      Bd_pos = sum(Bd.type, na.rm = TRUE),
      Bd_neg = length(which(Bd.type == 0)),
      Bd_Bsal_pos = length(which(Bd.type == 1 & Bsal.type == 1)),
      Bd_Bsal_neg = length(which(Bd.type == 0 & Bsal.type == 0)))


# Switch the order of I and D
Bsal$Treatment <- factor(Bsal$Treatment, levels = c("C", "I", "D"))


# Table 1. Sample size table
tab_1 <- ddply(.data = Bsal, 
               .variable = c("Genus", "Species", "Treatment"), 
               .fun = summarize, 
               Exp_samp_size = length(unique(ID_number)))

tab_1 <- cbind(tab_1, hist_tab[,-c(1:2)])




# Remove control individuals
Bsal2 <- Bsal[Bsal$Treatment != "C",]


# Infection intensity when they died or on day 81
Dead <- Bsal2[Bsal2$Census == 1 | Bsal2$Exp_day == 81, ]
Dead$Num <- 1


# Switch the order of I and D
Dead$Treatment <- factor(Dead$Treatment, levels = c("C", "I", "D"))

# Table 2. Summary infections
tab_2 <- ddply(.data = Dead, 
                .variable = c("Genus", "Species", "Treatment"), 
                .fun = summarize, 
                samp_size = length(unique(ID_number)),
                Bsal = mean(Bsal_load, na.rm = T), 
                Bsal_se = sd(Bsal_load, na.rm = T)/sqrt(length(Num)),
                Bd = mean(Bd_load, na.rm = T), 
                Bd_s = sd(Bd_load, na.rm = T)/sqrt(length(Num)))


# Table 1 main text
write.csv(tab_1, file = "./Tables/Table_1.csv")

# Table 2
write.csv(tab_2, file = "./Tables/Table_2.csv")





# 29. Phylogeny: read in data ------------------------------------------------




# Load amphibian tree
amphibian_tree <- read.tree(file = "~/Dropbox/Bsal_Summer_2015/Data/Phylogeny/amphibia_species.nwk")

# Read in trait data
traits <- read.csv(file = "~/Dropbox/Bsal_Summer_2015/Data/Traits_all.csv")





# 30. Define species to keep ------------------------------------------------


# Species used in our experiment
sp1 <- which(amphibian_tree$tip.label %in% "Desmognathus_fuscus")
sp2 <- which(amphibian_tree$tip.label %in% "Notophthalmus_viridescens")
sp3 <- which(amphibian_tree$tip.label %in% "Plethodon_cinereus")
sp4 <- which(amphibian_tree$tip.label %in% "Plethodon_cylindraceus")
sp5 <- which(amphibian_tree$tip.label %in% "Plethodon_glutinosus")
sp6 <- which(amphibian_tree$tip.label %in% "Plethodon_montanus")
sp7 <- which(amphibian_tree$tip.label %in% "Desmognathus_wrighti")
sp8 <- which(amphibian_tree$tip.label %in% "Eurycea_wilderae")


# Species used in Martel et al. 2014
martel_paper <- c(
  
  which(amphibian_tree$tip.label %in% "Bombina_variegata"), # FROG
  which(amphibian_tree$tip.label %in% "Alytes_obstetricans"), # FROG
  which(amphibian_tree$tip.label %in% "Xenopus_tropicalis"), #
  which(amphibian_tree$tip.label %in% "Pelobates_fuscus"), #
  which(amphibian_tree$tip.label %in% "Pelodytes_punctatus"), #
  which(amphibian_tree$tip.label %in% "Hyla_arborea"), #
  which(amphibian_tree$tip.label %in% "Rana_catesbeiana"), #
  which(amphibian_tree$tip.label %in% "Rana_temporaria"), #
  which(amphibian_tree$tip.label %in% "Hynobius_retardatus"), #
  which(amphibian_tree$tip.label %in% "Salamandrella_keyserlingii"), #
  which(amphibian_tree$tip.label %in% "Pachyhynobius_shangchengensis"), #
  which(amphibian_tree$tip.label %in% "Siren_intermedia"), #
  which(amphibian_tree$tip.label %in% "Hydromantes_strinatii"), #
  which(amphibian_tree$tip.label %in% "Gyrinophilus_porphyriticus"), #
  which(amphibian_tree$tip.label %in% "Ambystoma_opacum"), #
  which(amphibian_tree$tip.label %in% "Salamandrina_perspicillata"), #
  which(amphibian_tree$tip.label %in% "Salamandra_salamandra"), #
  which(amphibian_tree$tip.label %in% "Pleurodeles_waltl"), #
  which(amphibian_tree$tip.label %in% "Tylototriton_wenxianensis"), #
  which(amphibian_tree$tip.label %in% "Taricha_granulosa"), #
  which(amphibian_tree$tip.label %in% "Euproctus_platycephalus"), #
  which(amphibian_tree$tip.label %in% "Lissotriton_helveticus"), #
  which(amphibian_tree$tip.label %in% "Lissotriton_italicus"), #
  which(amphibian_tree$tip.label %in% "Ichthyosaura_alpestris"), #
  which(amphibian_tree$tip.label %in% "Triturus_cristatus"), #
  which(amphibian_tree$tip.label %in% "Neurergus_crocatus"), #
  which(amphibian_tree$tip.label %in% "Cynops_cyanurus"), #
  which(amphibian_tree$tip.label %in% "Cynops_pyrrhogaster"), #
  which(amphibian_tree$tip.label %in% "Paramesotriton_deloustali"),
  which(amphibian_tree$tip.label %in% "Typhlonectes_natans"), #
  which(amphibian_tree$tip.label %in% "Discoglossus_pictus"), #
  which(amphibian_tree$tip.label %in% "Epidalea_calamita")
)
#


#---- Note these names changes

# Martel names to TREE names

# Lithobates_catesbeiana to Rana_catesbeiana
# Discoglossus_scovazzi to Discoglossus_pictus
# Xenopus_tropicalis to Silurana_tropicalis

#-------

# List all species to keep
keep <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, martel_paper)


# Total number of species in phylogeny
length(keep)

# Make a list of 1 to the number of tips
tips <- 1:length(amphibian_tree$tip.label)

# Remove all tips except the ones we want to keep
remove <- tips[-keep]

# New tree will only have our species listed
new_tree <- drop.tip(amphibian_tree, amphibian_tree$tip.label[c(remove)])

str(new_tree)




# 31. Order the tips ------------------------------------------------


#----- Order the tips

name.check(phy= new_tree, data= traits)

length(traits$tip.label)
length(new_tree$tip.label)

# Order the traits the same as the tree
traits <- traits[match(new_tree$tip.label, traits$tip.label),]




# 32. Ancestoral state reconstruction  ------------------------------------------------


#-----------------
x <- traits$Disease
names(x) <- traits$tip.label

fitER<-rerootingMethod(new_tree, x, model="ER")
fitER

fitSYM<-rerootingMethod(new_tree, x, model="SYM")
fitSYM

fitARD<-rerootingMethod(new_tree, x, model="ARD")
fitARD


# Estimate AICc values
AICcCustom(fitER$loglik, 1, return.K = FALSE, second.ord = TRUE, nobs = length(traits$Disease), c.hat = 1)
# 84.00499
AICcCustom(fitSYM$loglik, 6, return.K = FALSE, second.ord = TRUE, nobs = length(traits$Disease), c.hat = 1)
# 88.60291
AICcCustom(fitARD$loglik, 12, return.K = FALSE, second.ord = TRUE, nobs = length(traits$Disease), c.hat = 1)
# 104.7733




# 33. Fig. 6 main text ------------------------------------------------ 



# Change the tip labels
new_tree$tip.label <- as.character(traits$new.tip.label)


# Make figure 6
pdf("./Figs_2020/Fig6.pdf",
    height = 10,
    width = 12)
plot.phylo(new_tree, cex=1.2, label.offset=8, edge.width=2)
axisPhylo()
axis(side = 1, at = 105, line = 2, labels = "Evolutionary time (Mya)", cex.axis = 1.5, cex = 1.5)
rect(200, 0, 300, 40, col = "gray", border = "transparent", density = 30) # coloured
text(250, 41, expression(paste("Origin ", italic(Bsal), sep = "")), cex = 1.2)
tiplabels(pch=22, bg= as.character(traits$Color), cex=2.5)
colors <- c("red", "green", "orange", "yellow")
legend("bottomleft", c("Lethal", "Susceptible", "Tolerant", "Resistant"), 
       fill = c("red", "orange","yellow", "green"), bty = "n", title = expression(bold(Legend)), 
       cex=1.5)

#fitERpiechartsubset <- fitER$lik.anc
#fitERpiechartsubset[1:3,] <- fitERpiechartsubset[6:7,] <- fitERpiechartsubset[16,] <- 
#  fitERpiechartsubset[30:35,] <- fitERpiechartsubset[37:38,] <- c(NA,NA,NA,NA)

fitERpiechartsubset <- fitER$marginal.anc
# need to cross out nodes: 41, 42, 43, 46, 47, 56, 71, 72, 73, 75, 78, 79, 74

ROW <- match(c(41, 42, 43, 46, 47, 56, 71, 72, 73, 75, 78, 79, 74), rownames(fitER$marginal.anc))
fitERpiechartsubset[ROW,] <- NA
nodelabels(pie=fitERpiechartsubset,piecol=colors,cex=0.5)

dev.off()




# 34. Test for phylogenetic signal ------------------------------------------------ 




lambda0<-rescale(new_tree, "lambda", 0)#transforms the tree topology to one that has all internal branch lengths multiplied by 0 (i.e. lambda=0) creating one giant basal polytomy
par(mfrow=c(1,2))#remember from day1 session2 that this sets the graphical parameters so that the plotting device has 1 row and 2 colums, so we can now plot two trees next to each other.
plot(new_tree)
plot(lambda0)

#Now find the maximum likelihood estimate of lambda for susceptibility

traits.Disease <- traits$Disease
names(traits.Disease) <- traits$tip.label

bsal_lambda<-fitDiscrete(new_tree, traits.Disease, model="ER", treeTransform="lambda")


# To see if this indicates significant phylogenetic signal we can pull out the AICC of each model (better for this
# since our sample size is low), and compare the two

bsal_lambda0<-fitDiscrete(lambda0, traits.Disease, model="ER")

bsal_lambda0$opt$aicc # 105.5171
bsal_lambda$opt$aicc # 81.40654

delta.aicc <- bsal_lambda0$opt$aicc - bsal_lambda$opt$aicc
delta.aicc # is 24.11057, some phylogenetic signal exists


# End script
