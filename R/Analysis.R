##### Contents #####
# 1. Preliminaries
# 2. Process trait data
# 3. Phylogeny
# 4. Rauniker life-form versus SR
# 5. Remove hydrophytes, helophytes, c4, and cam plants
# #. Plot light versus stomatal ratio
# #. Plot growth form versus stomatal ratio
#

##### 1. Preliminaries #####

  # Directories
  pathMain <- "~/Google Drive/StomatalRatio/BES"
  pathRawData <- paste0(pathMain, "/RawData")
  pathProcData <- paste0(pathMain, "/ProcData")
  pathR <- paste0(pathMain, "/R")
  pathMS <- paste0(pathMain, "/ms")
  pathFigures <- paste0(pathMS, "/Figures")
  
  setwd(pathMain)
  
  # Libraries
  library(vioplot)
  library(ape)
  library(rncl)
  library(phytools)
  library(rstanarm)
  library(phylolm)
  library(Rphylopars)
  
##### 2. Process data #####
  
  # Stomata
  stomata <- read.csv(paste0(pathRawData, "/stomata.csv"), stringsAsFactors = F, strip.white = T)
  stomata$sr_propAd <- apply(stomata, 1, function(X) as.numeric(X["ad_density"]) / (as.numeric(X["ab_density"]) + as.numeric(X["ad_density"])))
  stomata$sr_even <- apply(stomata, 1, function(X) min(as.numeric(X[c("ab_density", "ad_density")])) / max(as.numeric(X[c("ab_density", "ad_density")])))
  stomata <- stomata[which(!is.na(stomata$sr_propAd)), ]
  
  # OLD - Ellenberg light index
  # ellenberg_light <- read.csv(paste0(pathRawData, "/ellenberg_light.csv"), stringsAsFactors = F, 
  #                            strip.white = T)
  
  # OLD - Life-form
  # lifeform <- read.csv(paste0(pathRawData, "/life-form.csv"), stringsAsFactors = F, strip.white = T)
  
  # Photosynthetic pathway
  photo <- read.csv(paste0(pathRawData, "/photo.csv"), stringsAsFactors = F, strip.white = T)
  all(stomata$species %in% photo$species) # should be true
  
  # PLANTATT (Hill et al 2004) for lifeform and Ellenberg light indicator values
  plantatt <- read.csv(paste0(pathRawData, "/plantatt.csv"), stringsAsFactors = F, strip.white = T)
  # stomata$species[which(!stomata$species %in% plantatt$Taxon.name)] # Missing data for a few
  # Combine bulbous and nonbulbous geophytes
  plantatt$LF1[plantatt$LF1 == "Gb"] <- "Gn"

  # Combine hydrophytes and annual hydrophytes
  plantatt$LF1[plantatt$LF1 == "Hz"] <- "Hy"

  # Combine phanerophytes and nanophanerophytes
  plantatt$LF1[plantatt$LF1 == "Pn"] <- "Ph"
  
  # Combine datasets
  # x1 <- match(stomata$species, ellenberg_light$species)
  # x2 <- match(stomata$species, lifeform$species)
  xp <- match(stomata$species, photo$species)
  xa <- match(stomata$species, plantatt$Taxon.name)
  # stomata$ellenberg_light <- ellenberg_light$ellenberg_light[x1]
  # stomata$lifeform <- lifeform$lifeform[x2]
  stomata$photo <- photo$photo[xp]
  stomata$lifeform <- plantatt$LF1[xa]
  stomata$ellenberg_light <- plantatt$L[xa]
  stomata$height <- plantatt$Hght[xa]
  stomata$wh <- plantatt$W[xa]

  # Remove missing values
  stomata <- stomata[complete.cases(stomata[, c("sr_propAd", "lifeform", "ellenberg_light")]), ]
  
  rm("photo", "plantatt", "xp", "xa")
  
  # Export  
  # write.csv(stomata, paste0(pathProcData, "/stomata.csv"))

##### 3. Phylogeny #####
  
  # Import phylogeny of British Flora from Lim et al 2014
  tmp <- readLines(paste0(pathRawData, "/S15105.nex"))
  l1 <- grep("BEGIN CHARACTERS", tmp) - 1
  l2 <- grep("BEGIN TREES", tmp) - 1
  writeLines(tmp[c(1:l1, l2:length(tmp))], con = paste0(pathRawData, "/Lim_etal_2014.nex"))
  phy <- read_nexus_phylo(paste0(pathRawData, "/Lim_etal_2014.nex"))

  # Check for species missing from phylogeny
  stomata$species1 <- gsub(" ", "_", stomata$species)
  notInPhy <- which(!stomata$species1 %in% phy$tip.label)
  stomata[notInPhy, ]; length(notInPhy)
  
  # Modify tip labels
  
  # Using Cerastium brachypetalum as proxy for C. cerastoides based on similar divergence in Scheen et al. 2004 (Am J Bot, 91(6): 943-952)
  tmp <- gsub("Cerastium_brachypetalum", "Cerastium_cerastoides", tmp)
  
  # Using Coincya monensis as proxy for C. wrightii 
  tmp <- gsub("Coincya_monensis", "Coincya_wrightii", tmp)

  # Using Geranium pusillum as proxy for G. versicolor based on similar divergence in Palazzesi et al. 2012 (Bio J of Linn Soc, 107(1): 67-85)
  # tmp <- gsub("Geranium_pusillum", "Geranium_versicolor", tmp)

  # Using Gnaphalium luteoalbum as proxy for G. norvegicum 
  tmp <- gsub("Gnaphalium_luteoalbum", "Gnaphalium_norvegicum", tmp)
  
  # Lithospermum purpurocaeruleum is mispelled
  tmp <- gsub("Lithospermum_purpureocaeruleum", "Lithospermum_purpurocaeruleum", tmp)

  # Using Lythrum portula as proxy for L. hyssopifolia 
  tmp <- gsub("Lythrum_portula", "Lythrum_hyssopifolia", tmp)

  # Using Peucedanum ostruthium as proxy for P. officinale 
  tmp <- gsub("Peucedanum_ostruthium", "Peucedanum_officinale", tmp)

  # Changing Populus_nigra_sens.lat. to Populus nigra officinale 
  tmp <- gsub("Populus_nigra_sens.lat.", "Populus_nigra", tmp)

  # Changing Rosa_canina_agg. to Rosa_canina officinale 
  tmp <- gsub("Rosa_canina_agg.", "Rosa_canina", tmp)

  # Using Stachys alpina as proxy for S. germanica based on similar divergence in Salmaki et al. 2013 (Mol Phylo Evol, 69(3): 535-551)
  tmp <- gsub("Stachys_alpina", "Stachys_germanica", tmp)
  
  # Using Tripleurospermum_maritimum as proxy for T. inodorum 
  tmp <- gsub("Tripleurospermum_maritimum", "Tripleurospermum_inodorum", tmp)

  # Using Viola_arvensis as proxy for V. kitaibeliana 
  tmp <- gsub("Viola_arvensis", "Viola_kitaibeliana", tmp)
  
  # Write first modified phylogeny
  writeLines(tmp[c(1:l1, l2:length(tmp))], con = paste0(pathProcData, "/Lim_etal_2014_mod1.nex"))
  phy <- read_nexus_phylo(paste0(pathProcData, "/Lim_etal_2014_mod1.nex"))
  
  # Bind additional tips
  tips <- phy$tip.label
  nodes <- sapply(tips, function(x,y) which(y == x), y = phy$tip.label)
  edge.lengths <- setNames(phy$edge.length[sapply(nodes, function(x, y) {
    which(y == x) }, y = phy$edge[,2])], names(nodes))
  
  # Make Agrostemma githago sister to all (Silene sp., Lychnis sp.) [based on Fior et al. 2006, Am J Bot 104(2):399-411]
  n <- getMRCA(phy, phy$tip.label[c(grep("Silene", phy$tip.label), grep("Lychnis", phy$tip.label))])
  
  # Branch length is average of distance from (Silene sp., Lychnis sp.) MRCA to (Silene sp., Lychnis sp.) sp.
  el1 <- sapply(nodes[c(grep("Silene_", phy$tip.label), grep("Lychnis_", phy$tip.label))], 
                nodeheight, tree = phy)
  el2 <- nodeheight(tree = phy, node = n)
  
  phy <- bind.tip(phy, "Agrostemma_githago", edge.length = mean(el1 - el2), where = n)

  # Make Heracleum mantegazzianum as sister to H. sphondylium
  phy <- bind.tip(phy, "Heracleum_mantegazzianum", 
                  edge.length = phy$edge.length[which(phy$tip.label == "Heracleum_sphondylium")],
                  where = which(phy$tip.label == "Heracleum_sphondylium"))
  
  # Make Portulaca oleracea sister to all Claytonia (based on Angiosperm Phylogeny Website v 13 http://www.mobot.org/mobot/research/apweb/welcome.html]
  # n <- getMRCA(phy, phy$tip.label[grep("Claytonia", phy$tip.label)])
  
  # Branch length is average of distance from Claytonia MRCA to Claytonia sp.
  # el1 <- sapply(nodes[grep("Claytonia_", names(nodes))], nodeheight, tree = phy)
  # el2 <- nodeheight(tree = phy, node = n)
  
  # phy <- bind.tip(phy, "Portulaca_oleracea", edge.length = mean(el1 - el2), where = n)
  
  # Make Spartina alterniflora, Spartina maritima, and Spartina x townsendii as sister to S. anglica
  phy <- bind.tip(phy, "Spartina_alterniflora", 
                  edge.length = phy$edge.length[which(phy$tip.label == "Spartina_anglica")],
                  where = which(phy$tip.label == "Spartina_anglica"))
  
  phy <- bind.tip(phy, "Spartina_maritima", 
                  edge.length = phy$edge.length[which(phy$tip.label == "Spartina_anglica")],
                  where = which(phy$tip.label == "Spartina_anglica"))

  # phy <- bind.tip(phy, "Spartina_x_townsendii", 
  #                 edge.length = phy$edge.length[which(phy$tip.label == "Spartina_anglica")],
  #                 where = which(phy$tip.label == "Spartina_anglica"))

  # Make Tephroseris integrifolia sister to all Senecio sp.
  n <- getMRCA(phy, phy$tip.label[grep("Senecio", phy$tip.label)])
  
  # Branch length is average of distance from Senecio MRCA to Senecio sp.
  el1 <- sapply(nodes[grep("Senecio_", names(nodes))], nodeheight, tree = phy)
  el2 <- nodeheight(tree = phy, node = n)
  
  phy <- bind.tip(phy, "Tephroseris_integrifolia", edge.length = mean(el1 - el2), where = n)
  
  # Make Viola hirta as sister to V. reichenbachiana (both in Viola sect. Viola)
  phy <- bind.tip(phy, "Viola_hirta", 
                  edge.length = phy$edge.length[which(phy$tip.label == "Viola_reichenbachiana")],
                  where = which(phy$tip.label == "Viola_reichenbachiana"))

  # Write second modified phylogeny
  write.nexus(phy, file = paste0(pathProcData, "/Lim_etal_2014_mod2.nex"))

  # Check for species missing from modified phylogeny
  stomata$species1 <- gsub(" ", "_", stomata$species)
  notInPhy <- which(!stomata$species1 %in% phy$tip.label)
  stomata[notInPhy, ]

  rm("edge.lengths", "el1", "el2", "l1", "l2", "n", "nodes", "notInPhy", "tips", "tmp")
 
  # Run PATHd8 (Britton et al. 2002)
  cat(sprintf("Sequence length = 4692;
              
              %s
              
              mrca: Azolla_filiculoides, Solanum_nigrum, maxage=454;
              mrca: Osmunda_regalis, Asplenium_onopteris, minage=354;
              mrca: Pinus_sylvestris, Solanum_nigrum, minage=310;
              mrca: Nymphaea_alba, Solanum_nigrum, minage=140;
              mrca: Nymphaea_alba, Solanum_nigrum, maxage=217;
              mrca: Nymphaea_alba, Nuphar_lutea, minage=115;
              mrca: Chelidonium_majus, Solanum_nigrum, fixage=121;
              mrca: Quercus_robur, Bryonia_dioica, minage=84;
              
              name of mrca: Azolla_filiculoides, Solanum_nigrum, name=tracheophytes;
              name of mrca: Osmunda_regalis, Asplenium_onopteris, name=monilophytes;
              name of mrca: Pinus_sylvestris, Solanum_nigrum, name=seedplants;
              name of mrca: Nymphaea_alba, Solanum_nigrum, name=angiosperms;
              name of mrca: Nymphaea_alba, Nuphar_lutea, name=nymphaceae;
              name of mrca: Chelidonium_majus, Solanum_nigrum, name=eudicots;
              name of mrca: Quercus_robur, Bryonia_dioica, name=curcibitales_fagales;
              
              ", write.tree(phy)), file = paste0(pathMain, "/PATHd8/pathd8_in"))
  
  # Run PATHd8
  system("cd ~/Google\\ Drive/StomatalRatio/BES/PATHd8; ./PATHd8 pathd8_in pathd8_out")
  
  # Import ultrametric tree
  pathd8_out <- readLines(paste0(pathMain, "/PATHd8/pathd8_out"))
  phy <- read.tree(text = substr(pathd8_out[11], 14, nchar(pathd8_out[11])))

  rm("pathd8_out")
  
##### 4. Violin plot of Raunikaer life form versus stomatal ratio #####

  # chamaephyte = subshrub (woody)
  # geophyte = cryptophyte resting in dry ground
  # hemicryptophyte = perennial with overwintering buds at soil surface
  # hydrophyte = cryptophyte resting in marshy ground
  # phanerophyte = normally woody perennial
  # therophyte = annual
  
  lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "hydrophyte", "phanerophyte", "therophyte")
  names(lf) = c("Ch", "Gn", "hc", "Hy", "Ph", "Th")
  lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]
  
  pdf(paste0(pathFigures, "/FigureS_violin.pdf"), 8, 5.25)
  par(mai = c(1.75, 1.5, 0.25, 0.25), cex.lab = 1.5, las = 1, cex = 1.5, cex.axis = 0.5)
  
  vioplot(subset(stomata$sr_propAd, stomata$lifeform == names(lf[1]), na.rm = TRUE),
          subset(stomata$sr_propAd, stomata$lifeform == names(lf[2]), na.rm = TRUE),
          subset(stomata$sr_propAd, stomata$lifeform == names(lf[3]), na.rm = TRUE),
          subset(stomata$sr_propAd, stomata$lifeform == names(lf[4]), na.rm = TRUE),
          subset(stomata$sr_propAd, stomata$lifeform == names(lf[5]), na.rm = TRUE),
          subset(stomata$sr_propAd, stomata$lifeform == names(lf[6]), na.rm = TRUE),
          col = "grey", colMed = "black", names = rep("", 6))
  
  clip(grconvertX(0, "ndc", "user"), grconvertX(1, "ndc", "user"),
       grconvertY(0, "ndc", "user"), grconvertY(1, "ndc", "user"))
  # Sample sizes
  text(1:6 + 0.25, grconvertY(1.5, "inches", "user"), labels = lf, srt = 30, pos = 2)
  mtext(table(stomata$lifeform)[names(lf)], 3, at = 1:6, line = 0)
  title(ylab = "Stomatal Ratio",  line = 3.5)
  title(xlab = "Raunikaer Life Form", line = 4.5)
  title(ylab = expression(SR[propAd] == SD[adaxial] / SD[total]), line = 2, cex.lab = 1) 
  dev.off()
  
  rm("lf")
  
##### 5. Remove hydrophytes, helophytes, c4, cam, and nonangiosperm plants #####
  
  # x1 <- c("phanerophyte", "geophyte", "hemicryptophyte", "chamaephyte", "therophyte")
  # stomata <- stomata[which(stomata$lifeform %in% x1 & stomata$photo == "c3"), ]
  x1 <- c("Ph", "Ch", "hc", "Gn", "Th")
  stomata <- stomata[which(stomata$lifeform %in% x1 & stomata$photo == "c3"), ]
  
  angioPhy <- extract.clade(phy, node = length(phy$tip.label) + which(phy$node.label == "angiosperms"))
  # plot(angioPhy, show.tip.label = F, show.node.label = T)
  
  stomata <- stomata[stomata$species1 %in% angioPhy$tip.label, ]
  
  stomata$species <- stomata$species1
  stomata <- stomata[, colnames(stomata) != "species1"]
  
  # Export new dataset and tree
  write.csv(stomata, paste0(pathProcData, "/stomata.csv"))
  write.nexus(angioPhy, file = paste0(pathProcData, "/angioPhy.nex"))
  
  rm("x1", "phy")

##### 6. Raunikaer life form versus Ellenberg light indicator values #####
  
  # Make dichotomous and slightly adjust edge lengths to make exactly ultrametric
  modAngioPhy <- multi2di(angioPhy)
  modAngioPhy$edge.length[modAngioPhy$edge.length == 0] <- 0.02
  x <- diag(vcv.phylo(modAngioPhy))
  is.ultrametric(modAngioPhy)
  n <- length(modAngioPhy$tip.label)
  te <- sapply(1:n, function(x, y) which(y == x), y = modAngioPhy$edge[, 2]) # terminal edges
  modAngioPhy$edge.length[te] <- modAngioPhy$edge.length[te] + (max(x) - x)
  is.ultrametric(modAngioPhy)
  modAngioPhy <- drop.tip(modAngioPhy, 
                          modAngioPhy$tip.label[!modAngioPhy$tip.label %in% stomata$species])
  
  # Phylogenetic ANOVA
  # fitLFvEL_pp <- phylopars.lm(ellenberg_light ~ lifeform,
  #                             trait_data = stomata[, c("species", "ellenberg_light", "lifeform")],
  #                             tree = modAngioPhy, model = "OU", REML = F)
  # saveRDS(fitLFvEL_pp, file = paste0(pathR, "/fitLFvEL_pp.rds"))
  fitLFvEL_pp <- readRDS(file = paste0(pathR, "/fitLFvEL_pp.rds"))
  anova(fitLFvEL_pp)

  rownames(stomata) <- stomata$species # need to change row names for phylolm
  fitLFvEL_pl <- phylolm(ellenberg_light ~ -1 + lifeform, model = "OUrandomRoot", data = stomata, 
                         phy = modAngioPhy)
  summary(fitLFvEL_pl)
  confint(fitLFvEL_pl)
  
  lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "phanerophyte", "therophyte")
  names(lf) = c("Ch", "Gn", "hc", "Ph", "Th")
  lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]
  
  pdf(paste0(pathFigures, "/FigureS_lf-light.pdf"), 4, 7)
  par(mai = c(1.5, 1, 0, 0.25), cex.lab = 1, las = 1)
  plot(0, 0, type = "n", xlim = c(0, 9), ylim = c(0, 5), xlab = "", ylab = "", axes = F)  
  
  for (i in 0:4) abline(h = i + seq(0.2, 0.6, 0.2), col = "grey")
  hs <- tapply(stomata$ellenberg_light, stomata$lifeform, hist, plot = F, breaks = 0:9)  
  for (i in 1:5) {
    rect(0:8, rep(i - 1, 9), 1:9, i - 1 + hs[[names(lf)[i]]]$density, col = "black", border = "white", lwd = 2)
  }

  # label life forms
  s <- sprintf("%s (%s)", lf, table(stomata$lifeform)[names(lf)])
  text(0 + 0.5 * strwidth(s), 1:5, labels = s, pos = 1)
  
  # axes
  axis(1, at = 0.5:8.5, labels = 1:9, lwd = 0, line = -1.5)
  for (i in 0:4) axis(2, at = i + seq(0.2, 0.6, 0.2), labels = seq(0.2, 0.6, 0.2), lwd = 0, lwd.ticks = 1)

  # axis labels
  mtext(c("Shade\n", "Partial\nShade", "Sun\n"), 1, at = c(1.5, 4.5, 7.5), line = 1.5)
  title(xlab = "Ellenberg light\nindicator", cex.lab = 1.5, line = 5)
  title(ylab = "Proportion", cex.lab = 1.5)
  
  clip(grconvertX(0, from = "ndc", "user"), grconvertX(1, from = "ndc", "user"),
       grconvertY(0, from = "ndc", "user"), grconvertY(1, from = "ndc", "user"))
  rect(grconvertX(0, from = "npc", "user"), 0:4, 
       grconvertX(1, from = "npc", "user"), 0.75:4.75)
  
  # Means + 95% confidence intervals
  ci <- confint(fitLFvEL_pl)
  ci <- ci[paste0("lifeform", names(lf)), ] # reorder
  mu <- coef(fitLFvEL_pl)[paste0("lifeform", names(lf))]
  for (i in 1:5) {
    lines(ci[i, ], rep(i - 1 + 0.675, 2), lwd = 2)
    points(mu[i], i - 1 + 0.675, pch = 21, col = "black", bg = "white")
  }
  
  dev.off()
  
  rm("ci", "hs", "i", "mu", "n", "s", "te", "x")
  
##### 7. Raunikaer life form and Ellenberg light indicator values versus SR_even #####
  
  # Phylogenetic ANOVA (takes several minutes to run)
  # fitSR_pp <- phylopars.lm(sr_even ~ ellenberg_light * lifeform,
  #                          trait_data = stomata[, c("species", "sr_even", "ellenberg_light", "lifeform")],
  #                          tree = modAngioPhy, model = "OU", REML = F)
  # saveRDS(fitSR_pp, file = paste0(pathR, "/fitSR_pp.rds"))
  fitSR_pp <- readRDS(file = paste0(pathR, "/fitSR_pp.rds"))
  anova(fitSR_pp)
  
  # Parametric bootstrap CIs of full model for plotting
  # NOTE - this could be done with Rphylopars, which might be better because parameter estimates are slightly different for some reason
  # set.seed(36502)
  # fitSR_pl <- phylolm(sr_even ~ -1 + lifeform + ellenberg_light:lifeform, model = "OUrandomRoot",
  #                     data = stomata, phy = modAngioPhy, boot = 1e3)
  # saveRDS(fitSR_pl, file = paste0(pathR, "/fitSR_pl.rds"))
  fitSR_pl <- readRDS(file = paste0(pathR, "/fitSR_pl.rds"))
  
  # Figure
  
  pdf(paste0(pathFigures, "/Figure_SRmultReg.pdf"), 4, 7)
  par(mfrow = c(5, 1), mar = rep(0, 4), cex.lab = 1, oma = c(5, 7, 1, 1))

  for (i in 5:1) {

    plot(0, 0, type = "n", xlim = c(0, 9), ylim = c(-0.05, 1.05), xlab = "", ylab = "", axes = F,
         frame.plot = F)
    if (i %in% seq(1, 5, 2)) axis(2, at = seq(0, 1, 0.5), lwd = 0, lwd.ticks = 1, las = 1)
    
    # polygon of confidence intervals
    x <- seq(min(stomata$ellenberg_light[stomata$lifeform == names(lf)[i]]),
             max(stomata$ellenberg_light[stomata$lifeform == names(lf)[i]]), length.out = 1e3)
    bs <- sapply(x, function(X) {
      fitSR_pl$bootstrap[, sprintf('lifeform%s', names(lf)[i])] + 
        fitSR_pl$bootstrap[, sprintf('lifeform%s:ellenberg_light', names(lf)[i])] * X
    } )
    ci <- apply(bs, 2, quantile, probs = c(0.025, 0.5, 0.975))
    polygon(c(x, rev(x)), c(ci[1, ], rev(ci[3, ])), col = "grey", border = NA)
    points(x, ci[2, ], lwd = 2, type = "l") # median
    points(x, ci[1, ], lwd = 2, lty = 2, type = "l") # lower 95%
    points(x, ci[3, ], lwd = 2, lty = 2, type = "l") # upper 95%

    # plot points
    with(subset(stomata, stomata$lifeform == names(lf)[i]), 
         points(jitter(ellenberg_light), sr_even, pch = 21, bg = rgb(0, 0, 0, 0.1), cex = 2))
    
    # label life forms
    s <- sprintf("%s (%s)", lf[i], table(stomata$lifeform)[names(lf)[i]])
    text(0 + 0.5 * strwidth(s), 1.05, labels = s, pos = 1)
    
  }
  
  # frame plot
  clip(grconvertX(0, from = "ndc", "user"), grconvertX(1, from = "ndc", "user"),
       grconvertY(0, from = "ndc", "user"), grconvertY(1, from = "ndc", "user"))
  abline(h = grconvertY(0:5, "npc", "user"))
  lines(rep(grconvertX(0, "npc", "user"), 2), grconvertY(c(0, 5), "npc", "user"))
  lines(rep(grconvertX(1, "npc", "user"), 2), grconvertY(c(0, 5), "npc", "user"))

  # axis labels
  axis(1, at = 0.5:8.5, labels = 1:9, lwd = 0, line = -1)
  axis(1, at = c(1.5, 4.5, 7.5), lwd = 0, labels = c("Shade", "Partial Shade", "Sun"),
       line = 0)
  mtext(1, text = "Ellenberg light indicator", cex = 1.5, line = 3, outer = T)
  
  mtext(2, text = "Stomatal Ratio",  line = 4.5, outer = T, cex = 1.5)
  
  mtext(2, text = expression(
    SR[even] == plain(min)~group("(", list(SD[ab], SD[ad]), ")") /
      plain(max)~group("(", list(SD[ab], SD[ad]), ")")), outer = T, line = 2.5)

  dev.off()
  
  #
  
  
  
  
  tmp <- stomata[complete.cases(stomata), ]
  fit1 <- aov(sr2 ~ ellenberg_light * lifeform + log(height) + ellenberg_light * wh, data = tmp)
  fit2 <- aov(sr2 ~ ellenberg_light * wh, data = tmp)
  fit3 <- aov(sr2 ~ ellenberg_light * log(height), data = tmp)
  
  rownames(stomata) <- stomata$species1
  # stomata$ellenberg_light <- stomata$ellenberg_light - mean(stomata$ellenberg_light)
  
  tmp <- drop.tip(angioPhy, angioPhy$tip.label[!angioPhy$tip.label %in% stomata$species1])
  tmp$edge.length[tmp$edge.length == 0] <- 0.02

  fit1 <- phylolm(sr2 ~ -1 + lifeform + ellenberg_light:lifeform, model = "OUrandomRoot",
                  data = stomata[, c("sr2", "ellenberg_light", "lifeform")], phy = tmp)
  fit2 <- phylolm(sr2 ~ -1 + lifeform + ellenberg_light, model = "OUrandomRoot",
                  data = stomata[, c("sr2", "ellenberg_light", "lifeform")], phy = tmp)
  fit3 <- phylolm(sr2 ~ -1 + lifeform, model = "OUrandomRoot",
                  data = stomata[, c("sr2", "ellenberg_light", "lifeform")], phy = tmp)
  fit4 <- phylolm(sr2 ~ ellenberg_light, model = "OUrandomRoot",
                  data = stomata[, c("sr2", "ellenberg_light", "lifeform")], phy = tmp)
  sapply(list(fit1, fit2, fit3, fit4), AIC)

  fit <- phylolm(sr2 ~ LF1 * ellenberg_light, model = "OUrandomRoot",
                 data = stomata[, c("sr2", "ellenberg_light", "LF1")], phy = tmp)
  summary(fit)  
  
  fit1 <- phylolm(sr2 ~ -1 + lifeform + ellenberg_light:lifeform, model = "OUrandomRoot",
                  data = stomata[, c("sr2", "ellenberg_light", "lifeform")], phy = tmp)
  fit2 <- phylolm(sr2 ~ -1 + lifeform + ellenberg_light, model = "OUrandomRoot",
                  data = stomata[, c("sr2", "ellenberg_light", "lifeform")], phy = tmp)
  fit3 <- phylolm(sr2 ~ -1 + lifeform, model = "OUrandomRoot",
                  data = stomata[, c("sr2", "ellenberg_light", "lifeform")], phy = tmp)
  fit4 <- phylolm(sr2 ~ ellenberg_light, model = "OUrandomRoot",
                  data = stomata[, c("sr2", "ellenberg_light", "lifeform")], phy = tmp)
  sapply(list(fit1, fit2, fit3, fit4), AIC)
  
  fit <- phylolm(sr2 ~ lifeform * ellenberg_light, model = "OUrandomRoot",
                 data = stomata[, c("sr2", "ellenberg_light", "lifeform")], phy = tmp)
  summary(fit)  

##### X. path analysis - experimental #####
  
  library("phylopath")
  candidates <- list(A = DAG(logAb ~ ellenberg_light, logAd ~ ellenberg_light),
                     B = DAG(logAb ~ ellenberg_light, logAd ~ ellenberg_light, logAd ~ logAb, logAb ~ logAd))
  plot(candidates$A)
  
  library(Rphylopars)
  
  stomata$species <- stomata$species1
  all(stomata$species %in% phy$tip.label) # should be TRUE
  phy <- chronos(phy1)
  plot(phy, show.tip.label = F)
  axis(1)
  
  phy1 <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% stomata$species])
  plot(phy1, show.tip.label = F)
  ppBM <- phylopars(stomata[, c("species", "ab_density", "ad_density")], tree = phy1,
                    model = "BM", pheno_error = FALSE, pheno_correlated = FALSE)
  cov2cor(ppBM$par$phylocov)
  ppOU <- phylopars(stomata[, c("species", "ab_density", "ad_density")], tree = multi2di(phy1),
                    model = "mvOU", pheno_error = FALSE, pheno_correlated = FALSE, full_alpha = TRUE)
  
  library("mvMORPH")
  rownames(stomata) <- stomata$species1
  #phy$edge.length <- phy$edge.length * 100
  stomata$logAb <- log(stomata$ab_density + 1)
  stomata$logAd <- log(stomata$ad_density + 1)
  
  mmOU <- mvOU(multi2di(phy1), stomata[, c("ab_density", "ad_density")], scale.height = FALSE)
  mmBM <- mvBM(multi2di(phy1), stomata[, c("ab_density", "ad_density")], scale.height = FALSE)

  mmOU <- mvOU(multi2di(phy1), stomata[, c("logAb", "logAd")], scale.height = FALSE)
  mmBM <- mvBM(multi2di(phy1), stomata[, c("logAb", "logAd")], scale.height = FALSE)
  
##### 3. Plot light versus stomatal ratio #####
  
  vioplot(subset(stomata$sr2, stomata$ellenberg_light == 3, na.rm = TRUE),
          subset(stomata$sr2, stomata$ellenberg_light == 4, na.rm = TRUE),
          subset(stomata$sr2, stomata$ellenberg_light == 5, na.rm = TRUE),
          subset(stomata$sr2, stomata$ellenberg_light == 6, na.rm = TRUE),
          subset(stomata$sr2, stomata$ellenberg_light == 7, na.rm = TRUE),
          subset(stomata$sr2, stomata$ellenberg_light == 8, na.rm = TRUE),
          subset(stomata$sr2, stomata$ellenberg_light == 9, na.rm = TRUE))
  
##### 5. #####

i <- "phanerophyte"
vioplot(subset(stomata$sr2, stomata$ellenberg_light == 3 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 4 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 5 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 6 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 7 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 8 & stomata$lifeform == i))
i <- "geophyte"
vioplot(subset(stomata$sr2, stomata$ellenberg_light == 3 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 4 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 5 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 6 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 7 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 8 & stomata$lifeform == i))

i <- "hemicryptophyte"
vioplot(subset(stomata$sr2, stomata$ellenberg_light == 3 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 4 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 5 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 6 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 7 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 8 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 9 & stomata$lifeform == i))

i <- "chamaephyte"
vioplot(subset(stomata$sr2, stomata$ellenberg_light == 4 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 5 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 6 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 7 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 8 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 9 & stomata$lifeform == i))

i <- "therophyte"
vioplot(subset(stomata$sr2, stomata$ellenberg_light == 4 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 6 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 7 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 8 & stomata$lifeform == i),
        subset(stomata$sr2, stomata$ellenberg_light == 9 & stomata$lifeform == i))


##### betaMix with EM #####

setwd("~/Google Drive/StomatalRatio/Analysis/BoundedOUMix")

###	Source functions

source("~/Google Drive/StomatalRatio/Analysis/BoundedOUMix/functions.R")
stomata$sr1 <- ifelse(stomata$sr1 == 0, 0.001, stomata$sr1)
stomata$sr1 <- ifelse(stomata$sr1 == 1, 0.999, stomata$sr1)

#
# 1 component
#

{
  fit1 <- sbOUmix(stomata$sr1, k = 1)
  # logLik, AIC, BIC
  fit1$logLik; fit1$AIC; fit1$BIC
  c1 <- fit1$par
  
  # pdf("component1.pdf", 4, 4)
  # par(mar = c(2, 2, 1, 1), cex.lab = 1.5)
  pdf("component1.pdf", 3, 2)		
  par(mar = c(2, 2, 2, 0.25), cex.lab = 1.5)
  plot(0, 0, type = "n", axes = F, xlim = c(0, 1), ylim = c(0, 4))
  usr <- par("usr")
  h <- hist(stomata$sr1, plot = F, breaks = seq(0, 1, 0.05))
  clip(usr[1], usr[2], usr[3], h$density[1] - 5)
  hist(stomata$sr1, freq = F, col = rgb(0, 0, 0, 33, maxColorValue = 255), 
       las = 1, breaks = seq(0, 1, 0.05), border = "white", xlab = "", ylab = "",
       ylim = c(0, 5), add = T)
  axis(1, at = seq(0, 1, 0.2))
  axis(2, at = c(0:2, 4), las = 1, labels = c(0:2, 9), lwd = 0, lwd.ticks = 1)
  axis(3, at = 0.5, tick = F, labels = sprintf("BIC = %s", round(fit1$BIC)), line = -1.5)
  title(main = "One regime")
  
  # Polygons of densities
  X <- seq(1e-3, 1 - 1e-3, length = 1e4)
  Y <- with(c1, dsbOU(X, phi, th, cens = NULL))
  gap <- which(Y < 3 | Y > 8)
  X <- X[gap]
  Y <- Y[gap]
  Y <- ifelse(Y > 8, Y - 5, Y)
  
  clip(-1, usr[2], usr[3], h$density[1] - 5)
  polygon(c(X, rev(X)), c(Y, rep(0, length(gap))), col = rgb(1, 0, 0, 0.2), border = NULL)
  rect(-1, 3, 1, 3.1, col = "white", border = NA)
  lines(c(-0.05, -0.03), c(2.95, 3.05))
  lines(c(-0.05, -0.03), c(3.05, 3.15))
  lines(c(-0.04, -0.04), c(0, 3))
  lines(c(-0.04, -0.04), c(3.1, 4))
  
  dev.off()
  
}

#
# 2 component
#

{
  fit2 <- sbOUmix(stomata$sr1, k = 2)
  # logLik, AIC, BIC
  fit2$logLik; fit2$AIC; fit2$BIC
  c2 <- fit2$par
  
  # pdf("component2.pdf", 4, 4)
  # par(mar = rep(2, 4))
  
  pdf("component2.pdf", 3, 2)	
  par(mar = c(2, 2, 2, 0.25), cex.lab = 1.5)
  
  plot(0, 0, type = "n", axes = F, xlim = c(0, 1), ylim = c(0, 5))
  usr <- par("usr")
  h <- hist(stomata$sr1, plot = F, breaks = seq(0, 1, 0.05))
  clip(usr[1], usr[2], usr[3], h$density[1] - 5)
  hist(stomata$sr1, freq = F, col = rgb(0, 0, 0, 33, maxColorValue = 255), 
       las = 1, breaks = seq(0, 1, 0.05), border = "white", xlab = "", ylab = "",
       ylim = c(0, 4), add = T)
  axis(1, at = seq(0, 1, 0.2))
  axis(2, at = c(0:2, 4), las = 1, labels = c(0:2, 9), lwd = 0, lwd.ticks = 1)
  axis(3, at = 0.5, tick = F, labels = sprintf("BIC = %s", round(fit2$BIC)), line = -1.5)
  title(main = "Two regimes")
  
  # Polygons of densities
  X <- seq(1e-4, 1 - 1e-4, length = 1e4)
  Y1 <- with(c2, mixCoef[1] * dsbOU(X, phi[1], th[1], cens = NULL)) 
  Y2 <- with(c2, mixCoef[2] * dsbOU(X, phi[2], th[2], cens = NULL)) 
  
  gap <- which(Y1 < 3 | Y1 > 8)
  X1 <- X[gap]
  Y1 <- Y1[gap]
  Y1 <- ifelse(Y1 > 8, Y1 - 5, Y1)
  
  clip(-1, usr[2], usr[3], h$density[1] - 5)
  polygon(c(X1, rev(X1)), c(Y1, rep(0, length(gap))), col = rgb(1, 0, 0, 0.2), border = NULL)
  polygon(c(X, rev(X)), c(Y2, rep(0, 1e4)), col = rgb(0, 0, 1, 0.2), border = NULL)
  rect(-1, 3, 1, 3.1, col = "white", border = NA)
  lines(c(-0.05, -0.03), c(2.95, 3.05))
  lines(c(-0.05, -0.03), c(3.05, 3.15))
  lines(c(-0.04, -0.04), c(0, 3))
  lines(c(-0.04, -0.04), c(3.1, 5))
  dev.off()
}

#
# 3 components
#

{
  fit3 <- sbOUmix(stomata$sr1, k = 3)
  # logLik, AIC, BIC
  fit3$logLik; fit3$AIC; fit3$BIC
  c3 <- fit3$par
  
  pdf("component3.pdf", 3, 2)	
  par(mar = c(2, 2, 2, 0.25), cex.lab = 1.5)
  plot(0, 0, type = "n", axes = F, xlim = c(0, 1), ylim = c(0, 5))
  usr <- par("usr")
  h <- hist(stomata$sr1, plot = F, breaks = seq(0, 1, 0.05))
  clip(usr[1], usr[2], usr[3], h$density[1] - 5)
  hist(stomata$sr1, freq = F, col = rgb(0, 0, 0, 33, maxColorValue = 255), 
       las = 1, breaks = seq(0, 1, 0.05), border = "white", xlab = "", ylab = "",
       ylim = c(0, 4), add = T)
  axis(1, at = seq(0, 1, 0.2))
  axis(2, at = c(0:2, 4), las = 1, labels = c(0:2, 9), lwd = 0, lwd.ticks = 1)
  axis(3, at = 0.5, tick = F, labels = sprintf("BIC = %s", round(fit3$BIC)), line = -1.5)
  title(main = "Three regimes")
  
  # Polygons of densities
  X <- seq(1e-4, 1 - 1e-4, length = 1e4)
  Y1 <- with(c3, mixCoef[1] * dsbOU(X, phi[1], th[1], cens = NULL))
  Y2 <- with(c3, mixCoef[2] * dsbOU(X, phi[2], th[2], cens = NULL))
  Y3 <- with(c3, mixCoef[3] * dsbOU(X, phi[3], th[3], cens = NULL))
  
  gap <- which(Y1 < 3 | Y1 > 8)
  X1 <- X[gap]
  Y1 <- Y1[gap]
  Y1 <- ifelse(Y1 > 8, Y1 - 5, Y1)
  
  clip(-1, usr[2], usr[3], h$density[1] - 5)
  polygon(c(X1, rev(X1)), c(Y1, rep(0, length(gap))), col = rgb(1, 0, 0, 0.2), border = NULL)
  polygon(c(X, rev(X)), c(Y2, rep(0, 1e4)), col = rgb(0, 0, 1, 0.2), border = NULL)
  polygon(c(X, rev(X)), c(Y3, rep(0, 1e4)), col = rgb(1, 2/3, 0, 0.2), border = NULL)
  rect(-1, 3, 1, 3.1, col = "white", border = NA)
  lines(c(-0.05, -0.03), c(2.95, 3.05))
  lines(c(-0.05, -0.03), c(3.05, 3.15))
  lines(c(-0.04, -0.04), c(0, 3))
  lines(c(-0.04, -0.04), c(3.1, 5))
  dev.off()
}

#
# 4 components
#

{
  fit4 <- sbOUmix(stomata$sr1, k = 4)
  # logLik, AIC, BIC
  fit4$logLik; fit4$AIC; fit4$BIC
  c4 <- fit4$par
  
  pdf("component4.pdf", 4, 4)
  
  par(mar = rep(2, 4))
  hist(stomata$sr1, freq = F, col = rgb(0, 0, 0, 33, maxColorValue = 255), las = 1,
       breaks = seq(0, 1, 0.05), border = "white", xlab = "", ylab = "",
       main = "Four components")
  axis(3, at = 0.5, tick = F, labels = sprintf("BIC = %s", round(fit4$BIC)), line = -1.5)
  
  # Polygons of densities
  X <- seq(1e-4, 1 - 1e-4, length = 1e4)
  Y1 <- with(c4, mixCoef[1] * dsbOU(X, phi[1], th[1], cens = NULL))
  Y2 <- with(c4, mixCoef[2] * dsbOU(X, phi[2], th[2], cens = NULL))
  Y3 <- with(c4, mixCoef[3] * dsbOU(X, phi[3], th[3], cens = NULL))
  Y4 <- with(c4, mixCoef[4] * dsbOU(X, phi[4], th[4], cens = NULL))
  Y <- Y1 + Y2 + Y3 + Y4
  # polygon(c(X, rev(X)), c(Y, rep(0, 1e4)), col = rgb(0, 0, 0, 0.2), border = NULL)
  
  polygon(c(X, rev(X)), c(Y1, rep(0, 1e4)), col = rgb(1, 0, 0, 0.2), border = NULL)
  polygon(c(X, rev(X)), c(Y2, rep(0, 1e4)), col = rgb(1, 2/3, 0, 0.2), border = NULL)
  polygon(c(X, rev(X)), c(Y3, rep(0, 1e4)), col = rgb(0, 0, 1, 0.2), border = NULL)
  polygon(c(X, rev(X)), c(Y4, rep(0, 1e4)), col = rgb(0, 1, 1, 0.2), border = NULL)
  dev.off()
}

##### SCRATCH #####

la <- read.csv(paste0(pathRawData, "/leaf_angle.txt"), skip = 3, h = T, stringsAsFactors = F)

x <- unique(la$AccSpeciesName)
trgt <- x[which(x %in% stomata$species)]

la <- la[la$AccSpeciesName %in% trgt, ]
summary(as.factor(la$Dataset))

#### betaMix with EM ####
#### betaMix with Stan ####
dat <- list(y = stomata$sr1, N = nrow(stomata), K = 1)
betaMix1 <- stan(file = "Stan/betaMix.stan", chains = 1, iter = 2e3, data = dat)
summary(betaMix1)

dat$K <- 2
betaMix2 <- stan(file = "Stan/betaMix.stan", chains = 1, iter = 2e3, data = dat)
summary(betaMix2)

dat$K <- 3
betaMix3 <- stan(file = "Stan/betaMix.stan", chains = 1, iter = 2e4, data = dat)
summary(betaMix3)


stomata$sr1 <- ifelse(stomata$sr1 == 0, 0.001, stomata$sr1)
stomata$sr1 <- ifelse(stomata$sr1 == 1, 0.999, stomata$sr1)
stomata$ellenberg_light <- as.factor(stomata$ellenberg_light)
fit <- stan_betareg(sr1 ~ ellenberg_light + lifeform | ellenberg_light + lifeform, data = stomata, 
                    link = "logit", link.phi = "log", 
                    chains = 1, iter = 2e3)
print(fit, depth = 2)
plot(fit)
pp_check(fit)
prior_summary(fit)
summary(fit)

# Analysis using phylolm. Want to try using beta regression
library("phylolm")
stomata <- stomata[complete.cases(stomata), ]
rownames(stomata) <- stomata$species1
phy <- drop.tip(phy, phy$tip.label[!(phy$tip.label %in% stomata$species1)])
tips <- phy$tip.label
nodes <- sapply(tips, function(x,y) which(y == x), y = phy$tip.label)
edge.lengths <- setNames(phy$edge.length[sapply(nodes, function(x, y) {
  which(y == x) }, y = phy$edge[,2])], names(nodes))
sort(phy$edge.length)
phy$edge.length[which(phy$edge.length == 0)] <- 1e-6

fit <- phylostep(sr1 ~ ellenberg_light * lifeform, data = stomata, phy = phy,
                 model = "BM", direction = "both")
fit <- phylostep(sr1 ~ ellenberg_light * lifeform, data = stomata, phy = phy,
                 model = "OUrandomRoot", upper.bound = c(100, 100), direction = "both")

fit1 <- phylolm(sr1 ~ ellenberg_light * lifeform, data = stomata, phy = phy,
                model = "BM")

fit2 <- phylolm(sr1 ~ ellenberg_light + lifeform, data = stomata, phy = phy,
                model = "BM", lower.bound = c(63.51885), upper.bound = c(63.51885))

summary(fit1)$aic
summary(fit2)$aic

fit1 <- phylolm(sr1 ~ ellenberg_light * lifeform, data = stomata, phy = phy,
                model = "OUrandomRoot")
x1 <- summary(fit1)
fit2 <- phylolm(sr1 ~ ellenberg_light + lifeform, data = stomata, phy = phy,
                model = "OUrandomRoot", lower.bound = c(x1$optpar, 0), 
                upper.bound = c(x1$optpar, 100), starting.value = c(x1$optpar, x1$sigma2))

c(x1$aic, x2$aic)

fit <- aov(sr1 ~ lifeform * ellenberg_light, data = stomata)
Anova(fit)
