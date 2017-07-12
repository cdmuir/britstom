source("R/header.R")

phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy_modified.nex"))
stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"))

##### Raunikaer life form versus Ellenberg light indicator values -----

# phylolm shows that there is no phylogenetic signal in this regression (very high alpha)
stomata %<>% as.data.frame() %>% set_rownames(.$species) # need to change row names for phylolm
fitLFvEL_pl <- phylolm(ellenberg_light ~ -1 + lifeform, model = "OUrandomRoot", 
                       data = stomata, phy = phy, upper.bound = 1e10)
summary(fitLFvEL_pl)

# therefore using standard ANOVA
fitLFvEL <- lm(ellenberg_light ~ lifeform, data = stomata)
aovLFvEL <- anova(fitLFvEL)

# for plotting confidence intervals
fitLFvEL <- lm(ellenberg_light ~ -1 + lifeform, data = stomata)

lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "phanerophyte", "therophyte")
names(lf) = c("Ch", "Gn", "hc", "Ph", "Th")
lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]

pdf(str_c(path_figures, "/figureS_lf-light.pdf"), 4, 7)
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
ci <- ci[str_c("lifeform", names(lf)), ] # reorder
mu <- coef(fitLFvEL_pl)[paste0("lifeform", names(lf))]
for (i in 1:5) {
  lines(ci[i, ], rep(i - 1 + 0.675, 2), lwd = 2)
  points(mu[i], i - 1 + 0.675, pch = 21, col = "black", bg = "white")
}

dev.off()

# Export objects to ms
export2ms(c("aovLFvEL"))