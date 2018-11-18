# define inverse logit function
inv_logit <- function(x) {1/(1+exp(-x))}


# define model run wrapper
run_model <- function(model) {
  force(model)
  function(...) {
    m <- model(...)
    m$run()
    m$query$model_query_to_np()
  }
}

# define functional for rolling function apmultplications
rollapply <- function(x, n, f, ...) {
  offset <- trunc(n / 2)
  locs <- (offset + 1):(length(x) - n + offset)
  num <- vapply(
    locs,
    function(i) {f(x[(i-offset):(i+offset)], ...)},
    numeric(1)
  )

  c(rep(NA, offset), num)
}

# define interpolated ecdf for Kullback-Leibler divergence, per Perez-Cruz (2008)

ecdf_c <- function(x) {
  x <- sort(x)
  x_u <- unique(x)
  x_rle <- rle(x)$lengths
  p_e <- (cumsum(x_rle)-0.5)/length(x)
  f <- approxfun(x_u, p_e, method='linear', yleft=0, yright=1, rule=2)
  f
}
