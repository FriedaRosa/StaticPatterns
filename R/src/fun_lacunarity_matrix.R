library(dplyr)

# Helper: max_iterations_R ---------------------------------------------
max_iterations_R <- function(w_vec, mode, mat_width, mat_height) {
  iters <- length(w_vec)
  for (w in w_vec) {
    if (mode == 1) {
      r <- w
      N_r <- (mat_width - r + 1) * (mat_height - r + 1)
    } else {
      r <- (w - 1) %/% 2
      N_r <- (mat_width - 2 * r) * (mat_height - 2 * r)
    }
    iters <- iters + N_r
  }
  iters
}

# Helper: internal_lacunarity_box_masses -------------------------------
internal_lacunarity_box_masses <- function(box_masses, fun = 1, N_r) {
  if (length(box_masses) <= 1) return(NA_real_)

  if (fun == 1) {
    max_value <- max(box_masses, na.rm = TRUE)
    if (max_value < 0) {
      return(NA_real_)
    }
    box_masses_int <- as.integer(round(box_masses))
    n_S_r <- integer(max_value + 1L)
    for (val in box_masses_int) {
      if (val >= 0 && val <= max_value) {
        n_S_r[val + 1L] <- n_S_r[val + 1L] + 1
      }
    }
    Q_S_r <- n_S_r / length(box_masses_int)

    S_vals <- seq_along(Q_S_r) - 1
    Z_1 <- sum(S_vals * Q_S_r)
    Z_2 <- sum((S_vals^2) * Q_S_r)

    if (Z_1 == 0) return(NA_real_)
    lac <- Z_2 / (Z_1 * Z_1)

  } else {
    mval <- mean(box_masses, na.rm = TRUE)
    if (mval == 0) return(NA_real_)
    sdev <- stats::sd(box_masses, na.rm = TRUE)
    lac  <- 1 + (sdev^2 / (mval^2))
  }
  lac
}

rcpp_lacunarity_R <- function(mat, w_vec, fun, mode) {
  mat_height <- nrow(mat)
  mat_width  <- ncol(mat)
  out <- numeric(length(w_vec))

  # Identify all valid (non-NaN) positions
  non_na_mask <- !is.na(mat)

  for (j in seq_along(w_vec)) {
    r <- w_vec[j]
    message(paste0("Processing r = ", r))

    N_r <- (mat_width - r + 1) * (mat_height - r + 1)

    if (N_r <= 0) {
      out[j] <- NA_real_
      next
    }

    box_masses <- numeric(N_r)
    counter <- 1

    # Iterate only over valid starting positions where a full r x r window contains no NA
    for (x in 1:(mat_height - r + 1)) {
      for (y in 1:(mat_width - r + 1)) {
        # Check if the entire window contains only valid (non-NaN) values
        if (all(non_na_mask[x:(x + r - 1), y:(y + r - 1)])) {
          window <- mat[x:(x + r - 1), y:(y + r - 1), drop = FALSE]
          box_masses[counter] <- sum(window)
          counter <- counter + 1
        }
      }
    }

    # Compute lacunarity only on valid box masses
    if (counter > 1) {
      out[j] <- internal_lacunarity_box_masses(box_masses[1:(counter - 1)], N_r)
    } else {
      out[j] <- NA_real_  # If no valid windows were found
    }
  }

  out
}

###############################################################################
### 2) The user-facing "lacunarity_R()" function optimized for matrices
###############################################################################

lacunarity_R <- function(x, r_vec = NULL, r_max = NULL,
                         plot = FALSE, save_plot = FALSE,
                         progress = FALSE, ncores = 1L, test = FALSE) {

  sp_name <- names(x)

  print(paste0("is matrix?", is.matrix(x)))

  # **Continue with normal processing**
  out <- dplyr::tibble(name = character(), i = integer(), r = integer(),
                       "ln(r)" = numeric(), Lac = numeric(), "ln(Lac)" = numeric())

    mat <- x
    ras_name <- if (!is.null(sp_name)) sp_name else paste0("Raster")

    lac_fun <- 1L
    lac_vals <- rcpp_lacunarity_R(mat, r_vec, lac_fun, mode = 1)

    this_out <- dplyr::tibble(name = rep(ras_name, length(lac_vals)),
                              r = r_vec, "ln(r)" = log(r_vec),
                              Lac = lac_vals, "ln(Lac)" = log(lac_vals))
    out <- dplyr::bind_rows(out, this_out)

  dplyr::arrange(out, i)
}
