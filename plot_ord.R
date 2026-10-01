library(QsRutils)
library(ade4)
library(ggordiplots)
library(GUniFrac)
load("expt.rda")
otu <- veganotu(expt)
d <- vegdist(decostand(otu, method = "log"), method = "bray")
add <- !(is.euclid(d))
ord <- cmdscale(d, k = nrow(otu)-1, eig = TRUE, add = add)
grps <- vegansam(expt)$Crop
gg_ordiplot(ord, groups = grps)


ord <- metaMDS(otu, distance = "bray")
gg_ordiplot(ord, groups = grps)

gen.uni <- GUniFrac(otu, phy_tree(expt), alpha = c(0, 0.5, 1))
str(gen.uni)
dw <- as.dist(gen.uni$unifracs[, , "d_1"])
add <- !(is.euclid(dw))
ord <- cmdscale(dw, k = nrow(otu)-1, eig = TRUE, add = add)
gg_ordiplot(ord, groups = grps)

d5 <- as.dist(gen.uni$unifracs[, , "d_0.5"])
add <- !(is.euclid(d5))
ord <- cmdscale(d5, k = nrow(otu)-1, eig = TRUE, add = add)
gg_ordiplot(ord, groups = grps)

expt.sub <- subset_samples(expt, Sampling_time == "Blooming")
unique(sample_data(expt.sub)[ , "Sampling_time"])
expt.sub <- subset_samples(expt.sub, Soil == "Bulk")
sum(taxa_sums(expt.sub)==0)
expt.sub <- prune_taxa(taxa_sums(expt.sub)>0, expt.sub)

otu <- veganotu(expt.sub)
sam <- vegansam(expt.sub)
d <- vegdist(decostand(otu, method = "log"), method = "bray")
add <- !(is.euclid(d))
ord <- cmdscale(d, k = nrow(otu)-1, eig = TRUE, add = add)
grps <- vegansam(expt.sub)$Crop
gg_ordiplot(ord, groups = grps)

rslt <- betadisper(d, group = grps, add=TRUE)
anova(rslt)
adonis(d ~ grps + sam$Season)
library(RVAideMemoire)
pairwise.perm.manova(d, grps)
