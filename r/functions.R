# R functions associated with Muir 2017

# Make cross tab table of growth form and life form
prop_table <- function(x1, x2) {
  
  tab <- table(x1, x2) %>%
    apply(1, function(X) X / sum(X)) %>% 
    as.matrix()
  tab
  
}

# Draw connector to show overlap between Raunkiaer life form and growth form
draw_connector <- function(x, y1, y2, col, nslice = 1e2) {
  
  # for debug
  # x <- c(0.2, 0.8)
  # y1 <- c(0.3, 0.5)
  # y2 <- c(0.2, 0.6)
  # nslice = 1e2
  # col <- c("#EDF8E9", "#006D2C")
  
  stopifnot(length(x) == 2L)
  stopifnot(length(y1) == 2L)
  stopifnot(length(y2) == 2L)
  stopifnot(length(col) == 2L)
  
  col <- colorRamp(col)(seq(0, 1, length.out = nslice))
  
  x_left <- seq(x[1], x[2], length.out = nslice + 1)
  nudge <- diff(x_left)[1] / 10
  x_right <- x_left[2:(nslice + 1)] %>% add(nudge)
  x_left <- x_left[1:nslice] %>% subtract(nudge)
  
  slope_bottom <- (y2[1] - y1[1]) / (x[2] - x[1])
  slope_top <- (y2[2] - y1[2]) / (x[2] - x[1])

  y_bottom_left <- slope_bottom * (x_left - x[1]) + y1[1] 
  y_bottom_right <- slope_bottom * (x_right - x[1]) + y1[1] 
  y_top_left <- slope_top * (x_left - x[1]) + y1[2]
  y_top_right <- slope_top * (x_right - x[1]) + y1[2]
  
  for (i in 1:nslice) {
    
    polygon(c(x_left[i], x_right[i], x_right[i], x_left[i]), 
            c(y_bottom_left[i], y_bottom_right[i], y_top_right[i], y_top_left[i]), 
            col = rgb(col[i, 1], col[i, 2], col[i, 3], maxColorValue = 255), 
            border = NA, lty = 0)

  }
  
}

# Convenience function to remove combinations of same growth form
remove_diag <- function(df) {
  
  df <- df[which(df[, 1] != df[, 2]), ]
  df
  
}

# Compare bootstrap samples to determine whether slopes are significantly different
compare_slopes <- function(s1, s2) {
  d <- s1 - s2
  ret <- mean(d)
  ret %<>% c(median(d))
  ret %<>% c(quantile(d, probs = c(0.025, 0.975)))
  ret %<>% set_names(c("mean", "median", "lower", "upper"))
  ret
}

# Draw 'histofan' around phylogeny
draw_histofan <- function(y, r1, r2) {
  
  # for debugging
  # y <- stomata$sr_even
  # r1 <- 1.1
  # r2 <- 1.2
  
  # rescale y to fit between radii
  y %<>% 
    divide_by(max(.) - min(.)) %>% 
    multiply_by(r2 - r1) %>% 
    subtract(min(.)) %>% 
    add(r1)
  
  n <- length(y)
  theta_mid <- seq(0, 2 * pi, length.out = n + 1)

  # Rings for guidance
  points(cos(theta_mid) * r1, sin(theta_mid) * r1, col = "grey", type = "l") # sr_even = 0
  points(cos(theta_mid) * (r1 + r2) / 2, sin(theta_mid) * (r1 + r2) / 2, 
         col = "grey", type = "l") # sr_even = 0
  points(cos(theta_mid) * r2, sin(theta_mid) * r2, col = "grey", type = "l") # sr_even = 0

  theta_mid %<>% magrittr::extract(1:n)
  theta_left <- theta_mid - pi / (n + 1)
  theta_right <- theta_mid + pi / (n + 1)
  
  X <- Y <- numeric()
  for (i in 1:n) {
    X %<>% c(cos(seq(theta_left[i], theta_right[i], length.out = 1e2)) * y[i])
    Y %<>% c(sin(seq(theta_left[i], theta_right[i], length.out = 1e2)) * y[i])
  }
  
  X %<>% c(.[1])
  Y %<>% c(.[1])
  points(X, Y, type = "l")

}

# Draw wedge for plotting trait data against radial phylogeny
draw_wedge <- function(theta1, theta2, r1, r2, col) {
  
  # for debug
  # theta1 = 0
  # theta2 = 1
  # r1 = 1
  # r2 = 1.1
  
  x <- c(cos(seq(theta1, theta2, length.out = 1e2)) * r1,
         cos(seq(theta2, theta1, length.out = 1e2)) * r2)
  
  y <- c(sin(seq(theta1, theta2, length.out = 1e2)) * r1,
         sin(seq(theta2, theta1, length.out = 1e2)) * r2)
  
  polygon(x, y, col = col, border = NA)
  
}

# Export objects to ms
export2ms <- function(x, path_export = getOption("path_export", "/export")) {
  
  for (i in 1:length(x)) {
    stopifnot(is.character(x[i]))
    tmp <- eval(parse(text = x[i]))
    path <- normalizePath(str_c(path_ms, path_export))
    if (!dir.exists(path)) dir.create(path)
    write_rds(tmp, str_c(path, "/", x[i]))
  }
  
}

# Import objects to ms
# import2ms <- function(path_export = getOption("path_export", "/export")) {
#   
#   path <- normalizePath(str_c(path_ms, path_export))
#   files <- list.files(path)
#   eval(parse(text = str_c(files, " <- read_rds('", path, "/", files, "')")),
#        envir = .GlobalEnv)
# 
# }

# Import objects to ms
import2ms <- function(path_export = getOption("path_export", "/export")) {
  
  path <- normalizePath(str_c(path_ms, path_export))
  objects <- list.files(path)
  files <- str_c(path, "/", objects) %>% 
    normalizePath() %>% 
    str_replace_all("\\\\", "\\\\\\\\")
  eval(parse(text = str_c(objects, " <- read_rds('", files, "')")),
       envir = .GlobalEnv)
  
}
# Calculate sr_even
calc_sr_even <- function(ab_density, ad_density) {
  apply(cbind(ab_density, ad_density), 1, min) / 
    apply(cbind(ab_density, ad_density), 1, max)
}
  
# Make statistical significance asterisks
sigStar <- function(pvalue)
{
  if (pvalue >= 0.05) return("n.s.")
  if (pvalue < 0.05 & pvalue >= 0.01) return("*")
  if (pvalue < 0.01 & pvalue >= 0.001) return("**")
  if (pvalue < 0.001) return("***")
}

# Determine a number's order of magnitude
oom <- function(x)
{
  # x should be a vector of numbers
  floor(log10(abs(x)))
}

# Determine number of significant digits to use for rounding in tables
sigDig <- function(x)
{
  # x should be a vector of numeric elements to be rounded to same significant digit
  ret <- round(x, -min(oom(x)))
  return(ret)
}

# Text string for p-value in tables
pval2latex <- function(x)
{
  stem <- x * 10 ^ -oom(x)
  ret <- ifelse(oom(x) > -3, 
                sprintf("%.*f", 3, x), 
                sprintf("%.*f $\\\\times10^{%s}$", 1, stem, oom(x)))
  return(ret)
}


