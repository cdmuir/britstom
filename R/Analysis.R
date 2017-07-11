
  
#### PROBABLY EVERYTHING BELOW HERE IS NO LONGER NEEDED AS OF JULY 2017  
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
