source("R/header.R")

phy <- read.nexus(file = str_c(path_proc_data, "/angio_phy_modified.nex"))
stomata <- read_csv(str_c(path_proc_data, "/stomata_filtered.csv"),
                    col_types = cols(
                      species = col_character(),
                      ab_density = col_double(),
                      ad_density = col_double(),
                      sr_propAd = col_double(),
                      sr_even = col_double(),
                      lifeform = col_character(),
                      ellenberg_light = col_integer(),
                      habit = col_character()
                    ))

##### Stomatal density and light -----

# density of stomatal density by side and L-value
dnsAb <- tapply(log10(stomata$ab_density + 1), stomata$ellenberg_light, function(X) {
  bp <- boxplot(X, plot = F)
  ret <- density(X[!X %in% bp$out], kernel = "gaussian", n = 1e3, 
                 adjust = 2, from = bp$stats[1, 1], to = bp$stats[5, 1])
  ret$stats <- bp$stats
  ret$out <- bp$out
  ret
} )

dnsAd <- tapply(log10(stomata$ad_density + 1), stomata$ellenberg_light, function(X) {
  bp <- boxplot(X, plot = F)
  ret <- density(X[!X %in% bp$out], kernel = "gaussian", n = 1e3, 
                 adjust = 2, from = bp$stats[1, 1], to = bp$stats[5, 1])
  ret$stats <- bp$stats
  ret$out <- bp$out
  ret
} )

pdf(str_c(path_figures, "/figure_SD-light.pdf"), 8, 10, useDingbats = FALSE)
par(mar = rep(0, 4), oma = c(5, 6, 1, 1), las = 1, mfrow = c(2, 1))
plot(0, 0, xlim = c(2.5, 9.5), ylim = c(0, 3), type = "n", axes = F, frame.plot = T)

mtext(1, text = "Ellenberg light indicator", outer = T, cex = 2, line = 3)
mtext(2, text = expression(paste("Stomatal density (no. ", mu, m^-2, ") [log scale]")), 
      outer = T, cex = 2, las = 0, line = 3)
axis(2, at = 0:3, labels = c(1, 10, 100, 1000))

# lines
abline(v = 3:9, col = "grey")

# adaxial
for (i in 3:9) {
  
  adDns <- dnsAd[[sprintf("%s", i)]]
  
  # scale density for plotting  
  adDns$d <- adDns$y * 0.4 / max(adDns$y)
  
  # polygon
  polygon(c(i + adDns$d, rev(i - adDns$d)), c(adDns$x, rev(adDns$x)), col = "grey")
  
  # outliers
  points(rep(i, length(adDns$out)), adDns$out, pch = 19)
  
  # bars (0.25 and 0.75 quantiles?)
  lines(rep(i, 2), adDns$stats[c(2, 4), 1], lwd = 5, lend = 1)
  
  # median
  points(i, adDns$stats[3, 1], pch = 19, cex = 1.5)
  
  # labels
  s <- "A. adaxial (upper)"
  text(grconvertX(0.025, "nfc", "user") + strwidth(s, cex = 2) / 2, 
       grconvertY(0.975, "nfc", "user"), labels = s, cex = 2, pos = 1)
  
}

plot(0, 0, xlim = c(2.5, 9.5), ylim = c(0, 3), type = "n", axes = F, frame.plot = T)

axis(1, at = 3:9)
axis(2, at = 0:3, labels = c(1, 10, 100, 1000))

# lines
abline(v = 3:9, col = "grey")

# abaxial
for (i in 3:9) {
  
  abDns <- dnsAb[[sprintf("%s", i)]]
  
  # scale density for plotting  
  abDns$d <- abDns$y * 0.4 / max(abDns$y)
  
  # polygons
  polygon(c(i - abDns$d, rev(i + abDns$d)), c(abDns$x, rev(abDns$x)), col = "grey")
  
  # outliers
  points(rep(i, length(adDns$out)), adDns$out, pch = 19)
  
  # bars (0.25 and 0.75 quantiles?)
  lines(rep(i, 2), abDns$stats[c(2, 4), 1], lwd = 5, lend = 1)
  
  # median
  points(i, abDns$stats[3, 1], pch = 19, cex = 1.5)
  
  # labels
  s <- "B. abaxial (lower)"
  text(grconvertX(0.025, "nfc", "user") + strwidth(s, cex = 2) / 2, 
       grconvertY(0.975, "nfc", "user"), labels = s, cex = 2, pos = 1)
  
}

dev.off()
