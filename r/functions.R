# R functions associated with Muir 2017

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
import2ms <- function(path_export = getOption("path_export", "/export")) {
  
  path <- normalizePath(str_c(path_ms, path_export))
  files <- list.files(path)
  eval(parse(text = str_c(files, " <- read_rds('", path, "/", files, "')")),
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


