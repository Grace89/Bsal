
# Load libraries
library(ape)
library(picante)
library(phytools)
library(geiger)

# Load amphibian tree
amphibian_tree <- read.tree(file = "~/Dropbox/Bsal_Summer_2015/Data/Phylogeny/amphibia_species.nwk")

str(amphibian_tree)

traits <- read.csv(file = "~/Dropbox/Bsal_Summer_2015/Data/Traits_all.csv")

is.ultrametric(amphibian_tree)
is.rooted(amphibian_tree)

sp1 <- which(amphibian_tree$tip.label %in% "Desmognathus_fuscus")
sp2 <- which(amphibian_tree$tip.label %in% "Notophthalmus_viridescens")
sp3 <- which(amphibian_tree$tip.label %in% "Plethodon_cinereus")
sp4 <- which(amphibian_tree$tip.label %in% "Plethodon_cylindraceus")
sp5 <- which(amphibian_tree$tip.label %in% "Plethodon_glutinosus")
sp6 <- which(amphibian_tree$tip.label %in% "Plethodon_montanus")
sp7 <- which(amphibian_tree$tip.label %in% "Desmognathus_wrighti")
sp8 <- which(amphibian_tree$tip.label %in% "Eurycea_wilderae")

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


#---- Names changes

# Desmognathus_organi to Desmognathus_wrighti
# Epidalea_calamita to Bufo_calamita
# Discoglossus_scovazzi to Discoglossus_pictus
# Xenopus_tropicalis to Silurana tropicalis

#-------

keep <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, martel_paper)
# Desmognathus_fuscus
# Desmognathus_fuscus
# Notophthalmus_viridescens
# Plethodon_cinereus
# Plethodon_cylindraceus
# Plethodon_glutinosus
# Plethodon_montanus

length(keep)

tips <- 1:length(amphibian_tree$tip.label)

remove <- tips[-keep]

new_tree <- drop.tip(amphibian_tree, amphibian_tree$tip.label[c(remove)])

str(new_tree)
#----- Order the tips

name.check(phy= new_tree, data= traits)

length(traits$tip.label)
length(new_tree$tip.label)

traits <- traits[match(new_tree$tip.label, traits$tip.label),]

#-----------------
x <- traits$Disease
names(x) <- traits$tip.label

fitER<-rerootingMethod(new_tree, x, model="ER")
fitER
fitSYM<-rerootingMethod(new_tree, x, model="SYM")
fitSYM
fitARD<-rerootingMethod(new_tree, x, model="ARD")
fitARD

library(AICcmodavg)

AICcCustom(fitER$loglik, 1, return.K = FALSE, second.ord = TRUE, nobs = length(traits$Disease), c.hat = 1)
# 84.00499
AICcCustom(fitSYM$loglik, 6, return.K = FALSE, second.ord = TRUE, nobs = length(traits$Disease), c.hat = 1)
# 88.60291
AICcCustom(fitARD$loglik, 12, return.K = FALSE, second.ord = TRUE, nobs = length(traits$Disease), c.hat = 1)
# 104.7733

plot.phylo(new_tree, cex=0.9, label.offset=8, edge.width=2)
axisPhylo()
axis(side = 1, at = 105, line = 2, labels = "Evolutionary time (Mya)")
rect(200, 0, 300, 40, col = "gray", border = "transparent", density = 30) # coloured
text(250, 41, expression(paste("Origin ", italic(Bsal), sep = "")))
tiplabels(pch=22, bg=as.character(traits$Color), cex=1.5)

tiplabels(pch=22, bg=as.character(traits$Color), cex=1.5)
colors <- c("red", "green", "orange", "yellow")

legend("bottomleft", c("Lethal", "Susceptible", "Tolerant", "Resistant"), 
       fill = c("red", "orange","yellow", "green"), bty = "n", title = expression(bold(Legend)), cex=0.9)

#fitERpiechartsubset <- fitER$lik.anc
#fitERpiechartsubset[1:3,] <- fitERpiechartsubset[6:7,] <- fitERpiechartsubset[16,] <- 
#  fitERpiechartsubset[30:35,] <- fitERpiechartsubset[37:38,] <- c(NA,NA,NA,NA)

fitERpiechartsubset <- fitER$marginal.anc
# need to cross out nodes: 41, 42, 43, 46, 47, 56, 71, 72, 73, 75, 78, 79, 74

ROW <- match(c(41, 42, 43, 46, 47, 56, 71, 72, 73, 75, 78, 79, 74), rownames(fitER$marginal.anc))
fitERpiechartsubset[ROW,] <- NA
nodelabels(pie=fitERpiechartsubset,piecol=colors,cex=0.4)

#dev.off()


##################### Test for phylogenetic signal ##########################

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

bsal_lambda0$opt$aicc # 107.2656
bsal_lambda$opt$aicc # 92.61347

delta.aicc <- bsal_lambda0$opt$aicc - bsal_lambda$opt$aicc
delta.aicc # is 24.11057, some phylogenetic signal exists


