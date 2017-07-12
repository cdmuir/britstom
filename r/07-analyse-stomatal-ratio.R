source("R/header.R")

phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy_modified.nex"))
stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"))

##### Raunikaer life form and Ellenberg light indicator values versus sr_even -----

lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "hydrophyte", "phanerophyte", "therophyte")
names(lf) <- c("Ch", "Gn", "hc", "Hy", "Ph", "Th")
lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]

# Phylogenetic regression model comparison
stomata %<>% as.data.frame() %>% set_rownames(.$species) # need to change row names for phylolm
fitSRa <- phylolm(sr_even ~ lifeform * ellenberg_light, model = "OUrandomRoot",
                  data = stomata, phy = phy, upper.bound = 10)
fitSRb <- phylolm(sr_even ~ lifeform + ellenberg_light, model = "OUrandomRoot",
                  data = stomata, phy = phy, upper.bound = 10)
fitSRc <- phylolm(sr_even ~ ellenberg_light, model = "OUrandomRoot",
                  data = stomata, phy = phy, upper.bound = 10)
fitSRd <- phylolm(sr_even ~ lifeform, model = "OUrandomRoot",
                  data = stomata, phy = phy, upper.bound = 10)
fitSRe <- phylolm(sr_even ~ 1, model = "OUrandomRoot",
                  data = stomata, phy = phy, upper.bound = 10)
aicSR <- sapply(list(fitSRa, fitSRb, fitSRc, fitSRd, fitSRe), AIC)

# Parametric bootstrap CIs of full model for p-values
set.seed(36502)
fitSR_pl <- phylolm(sr_even ~ -1 + lifeform + ellenberg_light:lifeform, 
                    model = "OUrandomRoot", data = stomata, phy = phy, 
                    boot = 1e4, upper.bound = 10)
write_rds(fitSR_pl, path = str_c(path_objects, "/fitSR_pl.rds"))
fitSR_pl <- read_rds(str_c(path_objects, "/fitSR_pl.rds"))

# Figure

pdf(str_c(path_figures, "/figure_SRmultReg.pdf"), 4, 7)
par(mfrow = c(5, 1), mar = rep(0, 4), cex.lab = 1, oma = c(5, 7, 1, 1))

for (i in 5:1) {
  
  plot(0, 0, type = "n", xlim = c(1, 9), ylim = c(-0.05, 1.05), xlab = "", 
       ylab = "", axes = F, frame.plot = F)
  if (i %in% seq(1, 5, 2)) axis(2, at = seq(0, 1, 0.5), lwd = 0, lwd.ticks = 1, 
                                las = 1)
  
  # polygon of confidence intervals
  x <- seq(min(stomata$ellenberg_light[stomata$lifeform == names(lf)[i]]),
           max(stomata$ellenberg_light[stomata$lifeform == names(lf)[i]]), 
           length.out = 1e3)
  bs <- sapply(x, function(X) {
    fitSR_pl$bootstrap[, sprintf('lifeform%s', names(lf)[i])] + 
      fitSR_pl$bootstrap[, sprintf('lifeform%s:ellenberg_light', 
                                   names(lf)[i])] * X
  } )
  ci <- apply(bs, 2, quantile, probs = c(0.025, 0.5, 0.975))
  polygon(c(x, rev(x)), c(ci[1, ], rev(ci[3, ])), col = "grey", border = NA)
  points(x, ci[2, ], lwd = 2, type = "l") # median
  points(x, ci[1, ], lwd = 2, lty = 2, type = "l") # lower 95%
  points(x, ci[3, ], lwd = 2, lty = 2, type = "l") # upper 95%
  
  # plot points
  with(subset(stomata, stomata$lifeform == names(lf)[i]), 
       points(jitter(ellenberg_light), sr_even, pch = 21, bg = rgb(0, 0, 0, 0.1), 
              cex = 2))
  
  # label life forms, slope, and statistical significance
  s1 <- sprintf("%s (%s)", lf[i], table(stomata$lifeform)[names(lf)[i]])
  text(1 + 0.5 * strwidth(s1), 1.05, labels = s1, pos = 1, adj = 1)
  slopeCI <- confint(fitSR_pl)[sprintf("lifeform%s:ellenberg_light", names(lf)[i]), ]
  pval <- 1 - length(which(fitSR_pl$bootstrap[, sprintf("lifeform%s:ellenberg_light", 
                                                        names(lf)[i])] > 0)) / fitSR_pl$boot 
  s2 <- bquote(slope == .(round(coef(fitSR_pl)[sprintf("lifeform%s:ellenberg_light", 
                                                       names(lf)[i])], 3))^.(sigStar(pval)))
  text(1 + 0.5 * strwidth(s2), 1.05 - 1.5 * strheight(s1), labels = s2, pos = 1, adj = 1)
  
}

# frame plot
clip(grconvertX(0, from = "ndc", "user"), grconvertX(1, from = "ndc", "user"),
     grconvertY(0, from = "ndc", "user"), grconvertY(1, from = "ndc", "user"))
abline(h = grconvertY(0:5, "npc", "user"))
lines(rep(grconvertX(0, "npc", "user"), 2), grconvertY(c(0, 5), "npc", "user"))
lines(rep(grconvertX(1, "npc", "user"), 2), grconvertY(c(0, 5), "npc", "user"))

# axis labels
axis(1, at = 1:9, labels = 1:9, lwd = 0, line = -1)
axis(1, at = c(2, 5, 8), lwd = 0, labels = c("Shade", "Partial Shade", "Sun"),
     line = 0)
mtext(1, text = "Ellenberg light indicator", cex = 1.5, line = 3, outer = T)

mtext(2, text = "Stomatal Ratio",  line = 4.5, outer = T, cex = 1.5)

mtext(2, text = expression(
  SR[even] == plain(min)~group("(", list(SD[ab], SD[ad]), ")") /
    plain(max)~group("(", list(SD[ab], SD[ad]), ")")), outer = T, line = 2.5)

dev.off()

# Export objects to ms
export2ms(c("aicSR", "fitSRa", "fitSRb", "fitSRc", "fitSRd", "fitSRe"))

# NOTE: Matt suggested deleting this, just using AIC
# Parametric bootstrap statistical signicance of deltaAIC between model 
# with interaction (fitSRa) and model without interaction (fitSRb)
# set.seed(444294)
# bootModB <- predict(fitSRb) + rTrait(n = 1e3, modAngioPhy, "OU", 
#                                     list(sigma2 = fitSRb$sigma2, alpha = fitSRb$optpar))
# colnames(bootModB) <- paste0("b", 1:ncol(bootModB))

# fit model of model A to data simulated from model B
# fitA2B <- lapply(colnames(bootModB), function(col) {
#  stomata$Y <- bootModB[, col]
#  phylolm(Y ~ lifeform * ellenberg_light, model = "OUrandomRoot", data = stomata, 
#          phy = modAngioPhy, upper.bound = 10)
#  })

# fit model of model A to data simulated from model B
# fitB2B <- lapply(colnames(bootModB), function(col) {
#  stomata$Y <- bootModB[, col]
#  phylolm(Y ~ lifeform + ellenberg_light, model = "OUrandomRoot", data = stomata, 
#          phy = modAngioPhy, upper.bound = 10)
#  })

# refit models at parameter bounds
# for (i in which(round(sapply(fitA2B, function(X) X$optpar), 2) == 10)) {
#  stomata$Y <- bootModB[, colnames(bootModB)[i]]
#  fitA2B[[i]] <- phylolm(Y ~ lifeform * ellenberg_light, model = "OUrandomRoot", data = stomata,
#                        phy = modAngioPhy, lower.bound = 10, upper.bound = 200, 
#                         starting.value = 10)
# }

# refit models at parameter bounds
# for (i in which(round(sapply(fitB2B, function(X) X$optpar), 2) == 10)) {
#  stomata$Y <- bootModB[, colnames(bootModB)[i]]
#  fitB2B[[i]] <- phylolm(Y ~ lifeform + ellenberg_light, model = "OUrandomRoot", data = stomata,
#                        phy = modAngioPhy, lower.bound = 10, upper.bound = 200,
#                        starting.value = 10)
# }

# for some reason, this fixes those with alpha estimated at upper bound
# for (i in which(round(sapply(fitB2B, function(X) X$optpar), 2) == 200)) {
#   stomata$Y <- bootModB[, colnames(bootModB)[i]]
#   fitB2B[[i]] <- phylolm(Y ~ lifeform + ellenberg_light, model = "OUrandomRoot", data = stomata,
#                          phy = modAngioPhy, lower.bound = 0, upper.bound = 300, 
#                          starting.value = 1)
# }

# saveRDS(fitA2B, paste0(pathR, "/fitA2B.rds"))
# saveRDS(fitB2B, paste0(pathR, "/fitB2B.rds"))

# fitA2B <- readRDS(paste0(pathR, "/fitA2B.rds"))
# fitB2B <- readRDS(paste0(pathR, "/fitB2B.rds"))

# Distribution of deltaAIC
# obsDeltaAIC <- fitSRa$aic - fitSRb$aic
# distDeltaAIC <- sapply(fitA2B, function(X) X$aic) - sapply(fitB2B, function(X) X$aic)
# pInteraction <- length(which(obsDeltaAIC > distDeltaAIC)) / 1e3