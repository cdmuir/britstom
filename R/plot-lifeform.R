source("R/header.R")

stomata <- read_csv(str_c(path_proc_data, "/stomata.csv"))
phy <- read.nexus(file = str_c(path_proc_data, "/Lim_etal_2014_final.nex"))

##### Plot lifeform versus stomatal ratio -----

# chamaephyte = subshrub (woody)
# geophyte = cryptophyte resting in dry ground
# hemicryptophyte = perennial with overwintering buds at soil surface
# hydrophyte = cryptophyte resting in marshy ground
# phanerophyte = normally woody perennial
# therophyte = annual

lf <- c("chamaephyte", "geophyte", "hemicryptophyte", "hydrophyte", "phanerophyte", "therophyte")
names(lf) <- c("Ch", "Gn", "hc", "Hy", "Ph", "Th")
lf <- lf[names(sort(tapply(stomata$sr_propAd, stomata$lifeform, mean)))]

pdf(str_c(path_figures, "/figureS_violin.pdf"), 8, 5.25)
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
