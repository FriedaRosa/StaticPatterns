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
internal_lacunarity_box_masses <- function(box_masses, fun, N_r) {
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

# Optimized rcpp_lacunarity_R (Corrected Matrix-Based Computation) -----
rcpp_lacunarity_R <- function(mat, w_vec, fun, mode, display_progress = FALSE) {
  mat_height <- nrow(mat)
  mat_width  <- ncol(mat)
  out <- numeric(length(w_vec))

  if (display_progress) {
    total_iters <- max_iterations_R(w_vec, mode, mat_width, mat_height)
    message("Total iterations (approx): ", total_iters)
  }

  for (j in seq_along(w_vec)) {
    w <- w_vec[j]
    r <- if (mode == 1) w else (w - 1) %/% 2

    N_r <- if (mode == 1) {
      (mat_width - r + 1) * (mat_height - r + 1)
    } else {
      (mat_width - 2 * r) * (mat_height - 2 * r)
    }

    if (N_r <= 0) {
      out[j] <- NA_real_
      next
    }

    box_masses <- numeric(N_r)
    counter <- 1

    for (x in 1:(mat_height - r + 1)) {
      for (y in 1:(mat_width - r + 1)) {
        window <- mat[x:(x + r - 1), y:(y + r - 1), drop = FALSE]
        #print(window)
        if (any(!is.na(window))) {
          box_masses[counter] <- if (fun == 1) sum(window, na.rm = TRUE) else max(window, na.rm = TRUE) - min(window, na.rm = TRUE)
        } else {
          box_masses[counter] <- NA_real_
        }
        counter <- counter + 1
      }
    }

    out[j] <- internal_lacunarity_box_masses(box_masses, fun, N_r)
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

  # **1. If `x` is already a list of matrices, use it directly**
  if (inherits(x, "list") && all(sapply(x, function(el) inherits(el, "matrix")))) {
    # No conversion needed
  }

  # **2. If `x` is a single matrix, wrap it in a list**
  else if (inherits(x, "matrix")) {
    x <- list(x)  # Ensure consistent format
  }

  # **3. Convert SpatRaster objects to matrices**
  else if (inherits(x, "list") && all(sapply(x, inherits, "SpatRaster"))) {
    x <- lapply(x, function(r) as.matrix(r, wide = TRUE))
  }

  else if (inherits(x, "SpatRaster")) {
    x <- list(as.matrix(x, wide = TRUE))  # Convert single SpatRaster to a list containing one matrix
  }

  # **4. If input is a folder path, read raster files and convert to matrices**
  else if (is.character(x)) {
    r_paths <- list.files(x, pattern = "\\.tif$", full.names = TRUE)
    if (length(r_paths) == 0) stop("No .tif files found.")
    x <- lapply(r_paths, function(f) as.matrix(terra::rast(f), wide = TRUE))
  }

  else {
    stop("x must be a list of matrices, a single matrix, a SpatRaster, a list of SpatRasters, or a folder path.")
  }

  # **Continue with normal processing**
  out <- dplyr::tibble(name = character(), i = integer(), r = integer(),
                       "ln(r)" = numeric(), Lac = numeric(), "ln(Lac)" = numeric())

  for (i in seq_along(x)) {
    mat <- x[[i]]
    ras_name <- if (!is.null(sp_name)) sp_name[i] else paste0("Raster_", i)

    lac_fun <- if (length(unique(mat)) <= 2) 1L else 2L
    lac_vals <- rcpp_lacunarity_R(mat, r_vec, lac_fun, mode = 1, display_progress = FALSE)

    this_out <- dplyr::tibble(name = rep(ras_name, length(lac_vals)),
                              i = i, r = r_vec, "ln(r)" = log(r_vec),
                              Lac = lac_vals, "ln(Lac)" = log(lac_vals))
    out <- dplyr::bind_rows(out, this_out)
  }

  dplyr::arrange(out, i)
}
