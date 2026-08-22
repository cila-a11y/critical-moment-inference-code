# Optional joint value--derivative multiplier band and local root certificate.
#
# This module is sourced by method.R after the primary engine has been defined.
# It deliberately does not modify, reuse, or combine decisions from the primary
# value band.  Its d=0 and d=1 slices share one augmented critical value, so all
# derivative-assisted statements are calibrated within this separate family.

augmented_feature_key <- function(moment_type, degree) {
  paste0("d", as.integer(degree), "::", moment_type)
}

power_derivative_matrix <- function(values, orders, moment_type, degree) {
  degree <- as.integer(degree)
  assert_true(length(degree) == 1L && !is.na(degree) && degree %in% 0:2,
              "The augmented derivative degree must be 0, 1, or 2.")
  out <- power_matrix(values, orders, moment_type)
  if (degree == 0L) return(out)
  log_magnitude <- log(abs(values))
  multiplier <- log_magnitude^degree
  multiplier[values == 0] <- 0
  out <- out * multiplier
  out[values == 0, ] <- 0
  if (any(!is.finite(out))) {
    stop_numerical("An augmented derivative feature is nonfinite.")
  }
  out
}

build_augmented_feature_object <- function(base_object) {
  degrees <- 0:1
  feature_x <- list()
  feature_y <- list()
  for (degree in degrees) {
    for (moment_type in base_object$moment_types) {
      key <- augmented_feature_key(moment_type, degree)
      if (degree == 0L) {
        feature_x[[key]] <- base_object$feature_x[[moment_type]]
        feature_y[[key]] <- base_object$feature_y[[moment_type]]
      } else {
        feature_x[[key]] <- power_derivative_matrix(
          base_object$x, base_object$orders, moment_type, degree
        )
        feature_y[[key]] <- power_derivative_matrix(
          base_object$y, base_object$orders, moment_type, degree
        )
      }
    }
  }
  combined_x <- do.call(cbind, feature_x)
  combined_y <- do.call(cbind, feature_y)
  mean_x <- colMeans(combined_x)
  mean_y <- colMeans(combined_y)
  centered_x <- sweep(combined_x, 2L, mean_x, FUN = "-")
  centered_y <- sweep(combined_y, 2L, mean_y, FUN = "-")
  variance <- (1 - base_object$lambda) * colMeans(centered_x^2) +
    base_object$lambda * colMeans(centered_y^2)
  raw_second_moment <- (1 - base_object$lambda) * colMeans(combined_x^2) +
    base_object$lambda * colMeans(combined_y^2)
  if (any(!is.finite(mean_x)) || any(!is.finite(mean_y)) ||
      any(!is.finite(centered_x)) || any(!is.finite(centered_y)) ||
      any(!is.finite(variance)) || any(!is.finite(raw_second_moment))) {
    stop_numerical(
      paste0(
        "An augmented feature mean, centered value, variance, or raw second ",
        "moment is nonfinite."
      )
    )
  }

  covariance_01 <- setNames(vector("list", length(base_object$moment_types)),
                            base_object$moment_types)
  for (moment_type in base_object$moment_types) {
    block_0 <- seq.int(
      (match(augmented_feature_key(moment_type, 0L), names(feature_x)) - 1L) *
        length(base_object$orders) + 1L,
      length.out = length(base_object$orders)
    )
    block_1 <- seq.int(
      (match(augmented_feature_key(moment_type, 1L), names(feature_x)) - 1L) *
        length(base_object$orders) + 1L,
      length.out = length(base_object$orders)
    )
    covariance_01[[moment_type]] <-
      (1 - base_object$lambda) * colMeans(
        centered_x[, block_0, drop = FALSE] * centered_x[, block_1, drop = FALSE]
      ) + base_object$lambda * colMeans(
        centered_y[, block_0, drop = FALSE] * centered_y[, block_1, drop = FALSE]
      )
  }

  list(
    x = base_object$x, y = base_object$y,
    orders = base_object$orders,
    moment_types = base_object$moment_types, degrees = degrees,
    n_x = base_object$n_x, n_y = base_object$n_y,
    n_eff = base_object$n_eff, lambda = base_object$lambda,
    feature_x = feature_x, feature_y = feature_y,
    combined_x = combined_x, combined_y = combined_y,
    centered_x = centered_x, centered_y = centered_y,
    delta = mean_y - mean_x, variance = variance,
    raw_second_moment = raw_second_moment,
    covariance_01 = covariance_01
  )
}

augmented_column_block <- function(object, moment_type, degree) {
  key <- augmented_feature_key(moment_type, degree)
  key_index <- match(key, names(object$feature_x))
  assert_true(!is.na(key_index), sprintf("Unknown augmented coordinate %s.", key))
  grid_size <- length(object$orders)
  seq.int((key_index - 1L) * grid_size + 1L, length.out = grid_size)
}

make_augmented_state_cache <- function(object) {
  cache <- new.env(parent = emptyenv(), hash = TRUE)
  for (moment_type in object$moment_types) {
    for (degree in object$degrees) {
      block <- augmented_column_block(object, moment_type, degree)
      for (index in seq_along(object$orders)) {
        key <- paste(
          augmented_feature_key(moment_type, degree),
          order_key(object$orders[index]), sep = "::"
        )
        assign(key, list(
          contrast = object$delta[block[index]],
          variance = object$variance[block[index]],
          raw_second_moment = object$raw_second_moment[block[index]],
          covariance_01 = object$covariance_01[[moment_type]][index]
        ), envir = cache)
      }
    }
  }
  cache
}

augmented_state_at <- function(object, order, moment_type, degree, cache = NULL) {
  degree <- as.integer(degree)
  key <- paste(augmented_feature_key(moment_type, degree), order_key(order), sep = "::")
  if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  value_x <- power_derivative_vector(object$x, order, moment_type, degree)
  value_y <- power_derivative_vector(object$y, order, moment_type, degree)
  mean_x <- mean(value_x)
  mean_y <- mean(value_y)
  variance <- (1 - object$lambda) * biased_variance(value_x) +
    object$lambda * biased_variance(value_y)
  raw_second <- (1 - object$lambda) * mean(value_x^2) +
    object$lambda * mean(value_y^2)
  state <- list(
    contrast = mean_y - mean_x, variance = variance,
    raw_second_moment = raw_second, covariance_01 = NA_real_
  )
  if (any(!is.finite(c(
    state$contrast, state$variance, state$raw_second_moment
  )))) {
    stop_numerical("An off-grid augmented feature aggregate is nonfinite.")
  }
  if (!is.null(cache)) assign(key, state, envir = cache)
  state
}

# Generic pooled empirical covariance.  The multiplier process uses the full
# cross-coordinate covariance automatically; this function exposes the exact
# finite-sample estimator for audits without materialising a large covariance
# matrix in every Monte Carlo replication.
augmented_empirical_covariance <- function(
    object, moment_type_1, degree_1, order_1,
    moment_type_2, degree_2, order_2) {
  x_1 <- power_derivative_vector(object$x, order_1, moment_type_1, degree_1)
  x_2 <- power_derivative_vector(object$x, order_2, moment_type_2, degree_2)
  y_1 <- power_derivative_vector(object$y, order_1, moment_type_1, degree_1)
  y_2 <- power_derivative_vector(object$y, order_2, moment_type_2, degree_2)
  covariance_x <- mean((x_1 - mean(x_1)) * (x_2 - mean(x_2)))
  covariance_y <- mean((y_1 - mean(y_1)) * (y_2 - mean(y_2)))
  (1 - object$lambda) * covariance_x + object$lambda * covariance_y
}

augmented_sample_local_envelopes <- function(values, left, right) {
  magnitude <- abs(values)
  log_magnitude <- log(magnitude)
  envelope_0 <- pmax(exp(left * log_magnitude), exp(right * log_magnitude))
  log_absolute <- abs(log_magnitude)
  envelope_0[magnitude == 0] <- 0
  log_absolute[magnitude == 0] <- 0
  envelope_1 <- envelope_0 * log_absolute
  envelope_2 <- envelope_1 * log_absolute
  envelope_1[magnitude == 0] <- 0
  envelope_2[magnitude == 0] <- 0
  out <- c(
    mean_h0 = mean(envelope_0), mean_h1 = mean(envelope_1),
    mean_h2 = mean(envelope_2),
    mean_h0_h1 = mean(envelope_0 * envelope_1),
    mean_h1_h2 = mean(envelope_1 * envelope_2),
    mean_h1_squared = mean(envelope_1^2),
    mean_h0_h2 = mean(envelope_0 * envelope_2)
  )
  if (any(!is.finite(out))) {
    stop_numerical("An augmented local envelope bound is nonfinite.")
  }
  out
}

augmented_local_cell_constants <- function(object, left, right, cache = NULL) {
  key <- paste("aug", order_key(left), order_key(right), sep = "::")
  if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  bound_x <- augmented_sample_local_envelopes(object$x, left, right)
  bound_y <- augmented_sample_local_envelopes(object$y, left, right)
  sample_variance_lipschitz <- function(bound, degree) {
    if (degree == 0L) {
      2 * (bound[["mean_h0_h1"]] +
             bound[["mean_h0"]] * bound[["mean_h1"]])
    } else {
      2 * (bound[["mean_h1_h2"]] +
             bound[["mean_h1"]] * bound[["mean_h2"]])
    }
  }
  sample_second_lipschitz <- function(bound, degree) {
    if (degree == 0L) 2 * bound[["mean_h0_h1"]] else
      2 * bound[["mean_h1_h2"]]
  }
  degree_constants <- lapply(0:1, function(degree) {
    derivative_name <- if (degree == 0L) "mean_h1" else "mean_h2"
    list(
      contrast_lipschitz = bound_x[[derivative_name]] + bound_y[[derivative_name]],
      variance_lipschitz =
        (1 - object$lambda) * sample_variance_lipschitz(bound_x, degree) +
        object$lambda * sample_variance_lipschitz(bound_y, degree),
      second_moment_lipschitz =
        (1 - object$lambda) * sample_second_lipschitz(bound_x, degree) +
        object$lambda * sample_second_lipschitz(bound_y, degree)
    )
  })
  names(degree_constants) <- c("d0", "d1")
  sample_covariance_lipschitz <- function(bound) {
    bound[["mean_h1_squared"]] + bound[["mean_h0_h2"]] +
      bound[["mean_h1"]]^2 + bound[["mean_h0"]] * bound[["mean_h2"]]
  }
  out <- list(
    degree = degree_constants,
    covariance_01_lipschitz =
      (1 - object$lambda) * sample_covariance_lipschitz(bound_x) +
      object$lambda * sample_covariance_lipschitz(bound_y)
  )
  if (!is.null(cache)) assign(key, out, envir = cache)
  out
}

validate_augmented_relative_guard <- function(
    object, interval, tolerance = 1e-12, initial_grid = NULL,
    max_levels = 6L, max_nodes = 28801L,
    state_cache = NULL, constant_cache = NULL) {
  # The guard is imposed separately for d=0 and d=1 as a variance/raw-second
  # ratio.  Both terms scale quadratically under a deterministic rescaling of
  # the corresponding coordinate, unlike an absolute variance cutoff.
  if (is.null(initial_grid)) {
    initial_grid <- object$orders[
      object$orders >= interval[1L] - 1e-12 &
        object$orders <= interval[2L] + 1e-12
    ]
  }
  initial_grid <- sort(unique(c(interval, initial_grid)))
  max_levels <- as.integer(max_levels)
  max_nodes <- as.integer(max_nodes)
  assert_true(length(initial_grid) >= 2L,
              "The augmented guard grid has fewer than two nodes.")
  assert_true(max_nodes >= length(initial_grid),
              "max_nodes is smaller than the initial augmented guard grid.")
  nodes <- initial_grid
  queue <- lapply(seq_len(length(initial_grid) - 1L), function(index) {
    list(left = initial_grid[index], right = initial_grid[index + 1L], level = 0L)
  })
  accepted <- list()
  failures <- list()
  minimum_ratio <- c(d0 = Inf, d1 = Inf)
  head_index <- 1L
  while (head_index <= length(queue)) {
    cell <- queue[[head_index]]
    head_index <- head_index + 1L
    radius <- (cell$right - cell$left) / 2
    constants <- augmented_local_cell_constants(
      object, cell$left, cell$right, constant_cache
    )
    ratios <- matrix(
      NA_real_, nrow = length(object$moment_types), ncol = 2L,
      dimnames = list(object$moment_types, c("d0", "d1"))
    )
    for (moment_type in object$moment_types) {
      for (degree in 0:1) {
        local <- constants$degree[[degree + 1L]]
        left_state <- augmented_state_at(
          object, cell$left, moment_type, degree, state_cache
        )
        right_state <- augmented_state_at(
          object, cell$right, moment_type, degree, state_cache
        )
        variance_lower <- min(left_state$variance, right_state$variance) -
          local$variance_lipschitz * radius
        second_upper <- max(
          left_state$raw_second_moment, right_state$raw_second_moment
        ) + local$second_moment_lipschitz * radius
        ratios[moment_type, degree + 1L] <-
          if (!is.finite(variance_lower) || !is.finite(second_upper) ||
              variance_lower <= 0 || second_upper <= 0) 0 else
            variance_lower / second_upper
      }
    }
    if (all(is.finite(ratios) & ratios > tolerance)) {
      minimum_ratio <- pmin(minimum_ratio, apply(ratios, 2L, min))
      accepted[[length(accepted) + 1L]] <- cell
      next
    }
    can_split <- cell$level < max_levels && length(nodes) < max_nodes
    if (can_split) {
      midpoint <- cell$left + (cell$right - cell$left) / 2
      if (!any(abs(nodes - midpoint) <=
               .Machine$double.eps * max(1, abs(midpoint)))) {
        nodes <- c(nodes, midpoint)
      }
      queue[[length(queue) + 1L]] <-
        list(left = cell$left, right = midpoint, level = cell$level + 1L)
      queue[[length(queue) + 1L]] <-
        list(left = midpoint, right = cell$right, level = cell$level + 1L)
    } else {
      minimum_ratio <- pmin(minimum_ratio, apply(ratios, 2L, min))
      failures[[length(failures) + 1L]] <- cell
    }
  }
  partition <- c(accepted, failures)
  grid <- sort(unique(unlist(lapply(partition, function(cell) {
    c(cell$left, cell$right)
  }))))
  list(
    valid = length(failures) == 0L,
    minimum_certified_ratio = min(minimum_ratio),
    minimum_certified_ratio_d0 = minimum_ratio[["d0"]],
    minimum_certified_ratio_d1 = minimum_ratio[["d1"]],
    grid = grid, nodes_used = length(nodes),
    maximum_level = if (length(partition)) {
      max(vapply(partition, `[[`, integer(1L), "level"))
    } else 0L,
    unresolved_cells = length(failures), limit_hit = length(failures) > 0L
  )
}

nested_augmented_multiplier_bands <- function(
    base_object, order_intervals, bootstrap_reps, alpha,
    multiplier_distribution = "rademacher", bootstrap_batch_size = 50L,
    variance_tolerance = 1e-12, guard_max_levels = 6L,
    guard_max_nodes = 28801L, bootstrap_seed = NULL) {
  intervals <- normalise_order_intervals(order_intervals, base_object$orders)
  object <- build_augmented_feature_object(base_object)
  state_cache <- make_augmented_state_cache(object)
  constant_cache <- new.env(parent = emptyenv(), hash = TRUE)
  interval_indices <- lapply(intervals, function(interval) {
    which(object$orders >= interval[1L] - 1e-12 &
            object$orders <= interval[2L] + 1e-12)
  })
  guard_results <- Map(function(interval, indices) {
    validate_augmented_relative_guard(
      object = object, interval = interval, tolerance = variance_tolerance,
      initial_grid = object$orders[indices], max_levels = guard_max_levels,
      max_nodes = guard_max_nodes, state_cache = state_cache,
      constant_cache = constant_cache
    )
  }, intervals, interval_indices)
  names(guard_results) <- names(intervals)

  bootstrap_reps <- as.integer(bootstrap_reps)
  bootstrap_batch_size <- min(as.integer(bootstrap_batch_size), bootstrap_reps)
  assert_true(bootstrap_reps >= 2L,
              "At least two augmented multiplier replications are required.")
  if (!is.null(bootstrap_seed)) {
    assign(".Random.seed", bootstrap_seed, envir = .GlobalEnv)
  }
  standard_deviation <- sqrt(object$variance)
  invalid_columns <- !is.finite(standard_deviation) | standard_deviation <= 0
  safe_standard_deviation <- standard_deviation
  safe_standard_deviation[invalid_columns] <- 1
  columns_by_interval <- lapply(interval_indices, function(indices) {
    unlist(lapply(names(object$feature_x), function(key) {
      key_index <- match(key, names(object$feature_x))
      (key_index - 1L) * length(object$orders) + indices
    }), use.names = FALSE)
  })
  suprema <- matrix(
    NA_real_, nrow = bootstrap_reps, ncol = length(intervals),
    dimnames = list(NULL, names(intervals))
  )
  next_index <- 1L
  while (next_index <= bootstrap_reps) {
    batch_size <- min(bootstrap_batch_size, bootstrap_reps - next_index + 1L)
    # Resetting to the primary bootstrap seed before entering this loop and
    # calling the same row-major generator gives exactly the same (X,Y)
    # multiplier pair in every replication as the primary band.
    multiplier_pair <- draw_multiplier_pairs(
      batch_size, object$n_x, object$n_y, multiplier_distribution
    )
    process_y <- sqrt(object$lambda / object$n_y) *
      (multiplier_pair$y %*% object$centered_y)
    process_x <- sqrt((1 - object$lambda) / object$n_x) *
      (multiplier_pair$x %*% object$centered_x)
    standardized <- sweep(
      abs(process_y - process_x), 2L, safe_standard_deviation, FUN = "/"
    )
    if (any(invalid_columns)) standardized[, invalid_columns] <- Inf
    rows <- next_index:(next_index + batch_size - 1L)
    for (interval_id in names(intervals)) {
      if (isTRUE(guard_results[[interval_id]]$valid)) {
        suprema[rows, interval_id] <- apply(
          standardized[, columns_by_interval[[interval_id]], drop = FALSE],
          1L, max
        )
      }
    }
    next_index <- next_index + batch_size
  }
  quantile_index <- min(
    bootstrap_reps, ceiling((bootstrap_reps + 1) * (1 - alpha))
  )
  fits <- lapply(names(intervals), function(interval_id) {
    guard <- guard_results[[interval_id]]
    values <- suprema[, interval_id]
    critical_value <- if (isTRUE(guard$valid)) {
      value <- sort.int(values, partial = quantile_index)[quantile_index]
      if (!is.finite(value)) {
        stop_numerical("The augmented multiplier critical value is nonfinite.")
      }
      value
    } else NA_real_
    list(
      object = object, interval_id = interval_id,
      interval = intervals[[interval_id]], indices = interval_indices[[interval_id]],
      guard = guard, guard_failed = !isTRUE(guard$valid),
      critical_value = critical_value,
      bootstrap_supremum = if (isTRUE(guard$valid)) values else
        rep(NA_real_, bootstrap_reps)
    )
  })
  names(fits) <- names(intervals)
  list(
    object = object, intervals = intervals, fits = fits,
    bootstrap_suprema = suprema, state_cache = state_cache,
    constant_cache = constant_cache, bootstrap_reps = bootstrap_reps,
    alpha = alpha, multiplier_distribution = multiplier_distribution,
    quantile_index = quantile_index
  )
}

augmented_band_at_order <- function(fit, order, moment_type, degree,
                                    state_cache = NULL) {
  state <- augmented_state_at(
    fit$object, order, moment_type, degree, state_cache
  )
  half_width <- if (fit$guard_failed || !is.finite(state$variance) ||
                    state$variance < 0) {
    Inf
  } else {
    fit$critical_value * sqrt(state$variance / fit$object$n_eff)
  }
  c(
    lower = state$contrast - half_width,
    upper = state$contrast + half_width,
    contrast_hat = state$contrast, variance_hat = state$variance,
    half_width = half_width, covariance_01 = state$covariance_01
  )
}

extract_augmented_band <- function(fit, moment_type, degree) {
  block <- augmented_column_block(fit$object, moment_type, degree)
  indices <- fit$indices
  columns <- block[indices]
  contrast <- fit$object$delta[columns]
  variance <- fit$object$variance[columns]
  half_width <- if (fit$guard_failed) rep(Inf, length(indices)) else
    fit$critical_value * sqrt(pmax(variance, 0) / fit$object$n_eff)
  list(
    orders = fit$object$orders[indices], contrast_hat = contrast,
    variance_hat = variance, half_width = half_width,
    lower = contrast - half_width, upper = contrast + half_width
  )
}

vacuous_augmented_cell_table <- function(grid, degree, reason) {
  data.frame(
    left = head(grid, -1L), right = tail(grid, -1L), level = 0L,
    degree = as.integer(degree), lower_bound = -Inf, upper_bound = Inf,
    variance_lower_bound = NA_real_, band_lipschitz = Inf,
    positive = FALSE, negative = FALSE, unresolved = TRUE,
    reason = reason, stringsAsFactors = FALSE
  )
}

adaptive_augmented_enclosure <- function(
    fit, moment_type, degree, initial_grid = NULL,
    numerical_tolerance = 0, safety_margin = 0,
    max_levels = 6L, max_nodes = 28801L,
    refine_unresolved_levels = 0L, maximum_cell_half_width = Inf,
    roundoff_inflation = FALSE, state_cache = NULL, constant_cache = NULL) {
  degree <- as.integer(degree)
  interval <- fit$interval
  if (is.null(initial_grid)) initial_grid <- fit$guard$grid
  initial_grid <- sort(unique(c(interval, initial_grid)))
  initial_grid <- initial_grid[
    initial_grid >= interval[1L] - 1e-12 &
      initial_grid <= interval[2L] + 1e-12
  ]
  if (fit$guard_failed) {
    cells <- vacuous_augmented_cell_table(
      initial_grid, degree, "augmented_relative_guard_failed"
    )
    return(list(
      cells = cells, grid = initial_grid, nodes_used = length(initial_grid),
      maximum_level = 0L, limit_hit = FALSE,
      variance_unresolved_cells = nrow(cells),
      statistically_unresolved_cells = nrow(cells)
    ))
  }
  max_levels <- as.integer(max_levels)
  max_nodes <- as.integer(max_nodes)
  refine_unresolved_levels <- as.integer(refine_unresolved_levels)
  maximum_cell_half_width <- as.numeric(maximum_cell_half_width)
  assert_true(length(maximum_cell_half_width) == 1L &&
                !is.na(maximum_cell_half_width) && maximum_cell_half_width > 0,
              "The augmented maximum cell half-width must be positive.")
  assert_true(max_nodes >= length(initial_grid),
              "max_nodes is smaller than the augmented enclosure grid.")
  nodes <- initial_grid
  queue <- lapply(seq_len(length(initial_grid) - 1L), function(index) {
    list(left = initial_grid[index], right = initial_grid[index + 1L], level = 0L)
  })
  output <- list()
  head_index <- 1L
  limit_hit <- FALSE
  split_cell <- function(cell) {
    midpoint <- cell$left + (cell$right - cell$left) / 2
    if (!any(abs(nodes - midpoint) <=
             .Machine$double.eps * max(1, abs(midpoint)))) {
      nodes <<- c(nodes, midpoint)
    }
    queue[[length(queue) + 1L]] <<-
      list(left = cell$left, right = midpoint, level = cell$level + 1L)
    queue[[length(queue) + 1L]] <<-
      list(left = midpoint, right = cell$right, level = cell$level + 1L)
  }
  append_vacuous <- function(cell, variance_lower, reason) {
    output[[length(output) + 1L]] <<- data.frame(
      left = cell$left, right = cell$right, level = cell$level,
      degree = degree, lower_bound = -Inf, upper_bound = Inf,
      variance_lower_bound = variance_lower, band_lipschitz = Inf,
      positive = FALSE, negative = FALSE, unresolved = TRUE,
      reason = reason, stringsAsFactors = FALSE
    )
  }
  while (head_index <= length(queue)) {
    cell <- queue[[head_index]]
    head_index <- head_index + 1L
    radius <- (cell$right - cell$left) / 2
    constants <- augmented_local_cell_constants(
      fit$object, cell$left, cell$right, constant_cache
    )$degree[[degree + 1L]]
    left_state <- augmented_state_at(
      fit$object, cell$left, moment_type, degree, state_cache
    )
    right_state <- augmented_state_at(
      fit$object, cell$right, moment_type, degree, state_cache
    )
    variance_lower <- min(left_state$variance, right_state$variance) -
      constants$variance_lipschitz * radius
    can_split <- cell$level < max_levels && length(nodes) < max_nodes
    if (!is.finite(variance_lower) || variance_lower <= 0) {
      if (can_split) {
        split_cell(cell)
      } else {
        limit_hit <- TRUE
        append_vacuous(cell, variance_lower, "variance_unresolved_at_limit")
      }
      next
    }
    left_band <- augmented_band_at_order(
      fit, cell$left, moment_type, degree, state_cache
    )
    right_band <- augmented_band_at_order(
      fit, cell$right, moment_type, degree, state_cache
    )
    width_lipschitz <- fit$critical_value / sqrt(fit$object$n_eff) *
      constants$variance_lipschitz / (2 * sqrt(variance_lower))
    band_lipschitz <- constants$contrast_lipschitz + width_lipschitz
    local_scale <- max(
      1, abs(left_band), abs(right_band), abs(band_lipschitz * radius),
      na.rm = TRUE
    )
    roundoff_margin <- if (isTRUE(roundoff_inflation)) {
      16 * .Machine$double.eps * local_scale
    } else 0
    margin <- safety_margin + roundoff_margin
    lower_bound <- min(left_band[["lower"]], right_band[["lower"]]) -
      band_lipschitz * radius - margin
    upper_bound <- max(left_band[["upper"]], right_band[["upper"]]) +
      band_lipschitz * radius + margin
    positive <- is.finite(lower_bound) && lower_bound > numerical_tolerance
    negative <- is.finite(upper_bound) && upper_bound < -numerical_tolerance
    unresolved <- !(positive || negative)
    endpoint_same_sign <-
      (left_band[["lower"]] > numerical_tolerance &&
         right_band[["lower"]] > numerical_tolerance) ||
      (left_band[["upper"]] < -numerical_tolerance &&
         right_band[["upper"]] < -numerical_tolerance)
    refine_for_resolution <- radius > maximum_cell_half_width
    refine_statistically <- unresolved &&
      (cell$level < refine_unresolved_levels || endpoint_same_sign)
    should_refine <- refine_for_resolution || refine_statistically
    if (should_refine && can_split) {
      split_cell(cell)
      next
    }
    if (should_refine && !can_split) limit_hit <- TRUE
    output[[length(output) + 1L]] <- data.frame(
      left = cell$left, right = cell$right, level = cell$level,
      degree = degree, lower_bound = lower_bound, upper_bound = upper_bound,
      variance_lower_bound = variance_lower, band_lipschitz = band_lipschitz,
      positive = positive, negative = negative, unresolved = unresolved,
      reason = if (positive) "certified_positive" else if (negative) {
        "certified_negative"
      } else if (should_refine && !can_split) {
        "statistically_unresolved_at_limit"
      } else "statistically_unresolved",
      stringsAsFactors = FALSE
    )
  }
  cells <- do.call(rbind, output)
  cells <- cells[order(cells$left, cells$right), , drop = FALSE]
  row.names(cells) <- NULL
  list(
    cells = cells,
    grid = sort(unique(c(cells$left, cells$right))),
    nodes_used = length(nodes),
    maximum_level = if (nrow(cells)) max(cells$level) else 0L,
    limit_hit = limit_hit,
    variance_unresolved_cells = sum(cells$reason == "variance_unresolved_at_limit"),
    statistically_unresolved_cells = sum(cells$unresolved)
  )
}

augmented_certified_runs <- function(cells, sign_name) {
  cells <- cells[order(cells$left, cells$right), , drop = FALSE]
  cells <- cells[!is.na(cells[[sign_name]]) & cells[[sign_name]], , drop = FALSE]
  if (!nrow(cells)) return(list())
  runs <- list()
  current <- cells[1L, , drop = FALSE]
  if (nrow(cells) >= 2L) {
    for (index in seq.int(2L, nrow(cells))) {
      next_cell <- cells[index, , drop = FALSE]
      tolerance <- 64 * .Machine$double.eps * max(
        1, abs(current$right[nrow(current)]), abs(next_cell$left)
      )
      if (abs(current$right[nrow(current)] - next_cell$left) <= tolerance) {
        current <- rbind(current, next_cell)
      } else {
        runs[[length(runs) + 1L]] <- current
        current <- next_cell
      }
    }
  }
  runs[[length(runs) + 1L]] <- current
  runs
}

best_augmented_anchor_pair <- function(
    fit, moment_type, run, direction, numerical_tolerance, state_cache) {
  nodes <- sort(unique(c(run$left, run$right)))
  value_bands <- t(vapply(nodes, function(order) {
    augmented_band_at_order(fit, order, moment_type, 0L, state_cache)[c("lower", "upper")]
  }, numeric(2L)))
  if (direction == "up") {
    left_indices <- which(value_bands[, "upper"] < -numerical_tolerance)
    right_indices <- which(value_bands[, "lower"] > numerical_tolerance)
  } else {
    left_indices <- which(value_bands[, "lower"] > numerical_tolerance)
    right_indices <- which(value_bands[, "upper"] < -numerical_tolerance)
  }
  if (!length(left_indices) || !length(right_indices)) return(NULL)
  candidates <- expand.grid(
    left_index = left_indices, right_index = right_indices,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  candidates <- candidates[
    nodes[candidates$left_index] < nodes[candidates$right_index], , drop = FALSE
  ]
  if (!nrow(candidates)) return(NULL)
  candidates$width <- nodes[candidates$right_index] - nodes[candidates$left_index]
  candidates$p_left <- nodes[candidates$left_index]
  candidates$p_right <- nodes[candidates$right_index]
  candidates <- candidates[order(candidates$width, candidates$p_left), , drop = FALSE]
  selected <- candidates[1L, , drop = FALSE]
  left_band <- value_bands[selected$left_index, ]
  right_band <- value_bands[selected$right_index, ]
  derivative_margin <- if (direction == "up") {
    min(run$lower_bound)
  } else {
    -max(run$upper_bound)
  }
  data.frame(
    p_left = selected$p_left, p_right = selected$p_right,
    direction = direction,
    value_left_lower = left_band[["lower"]],
    value_left_upper = left_band[["upper"]],
    value_right_lower = right_band[["lower"]],
    value_right_upper = right_band[["upper"]],
    derivative_margin = derivative_margin,
    stringsAsFactors = FALSE
  )
}

derivative_assisted_certificates <- function(
    fit, moment_type, derivative_cells, numerical_tolerance = 0,
    state_cache = NULL) {
  certificates <- list()
  for (run in augmented_certified_runs(derivative_cells, "positive")) {
    certificate <- best_augmented_anchor_pair(
      fit, moment_type, run, "up", numerical_tolerance, state_cache
    )
    if (!is.null(certificate)) certificates[[length(certificates) + 1L]] <- certificate
  }
  for (run in augmented_certified_runs(derivative_cells, "negative")) {
    certificate <- best_augmented_anchor_pair(
      fit, moment_type, run, "down", numerical_tolerance, state_cache
    )
    if (!is.null(certificate)) certificates[[length(certificates) + 1L]] <- certificate
  }
  if (!length(certificates)) {
    return(data.frame(
      p_left = numeric(0L), p_right = numeric(0L), direction = character(0L),
      value_left_lower = numeric(0L), value_left_upper = numeric(0L),
      value_right_lower = numeric(0L), value_right_upper = numeric(0L),
      derivative_margin = numeric(0L), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, certificates)
  out <- out[order(out$p_left, out$p_right, out$direction), , drop = FALSE]
  row.names(out) <- NULL
  out
}

derivative_certificate_truth <- function(certificates, truth) {
  if (!nrow(certificates)) return(logical(0L))
  vapply(seq_len(nrow(certificates)), function(index) {
    if (truth$identically_zero) return(FALSE)
    inside <- truth$roots > certificates$p_left[index] &
      truth$roots < certificates$p_right[index]
    if (sum(inside) != 1L) return(FALSE)
    root_index <- which(inside)
    truth$multiplicities[root_index] == 1L &&
      isTRUE(truth$directions[root_index] == certificates$direction[index])
  }, logical(1L))
}

summarise_derivative_assisted_row <- function(
    augmented, interval_id, moment_type, truth,
    numerical_tolerance, safety_margin, max_levels, max_nodes,
    refine_unresolved_levels, maximum_cell_half_width,
    roundoff_inflation) {
  fit <- augmented$fits[[interval_id]]
  value_band <- extract_augmented_band(fit, moment_type, 0L)
  derivative_band <- extract_augmented_band(fit, moment_type, 1L)
  value_enclosure <- adaptive_augmented_enclosure(
    fit = fit, moment_type = moment_type, degree = 0L,
    initial_grid = fit$guard$grid,
    numerical_tolerance = numerical_tolerance, safety_margin = safety_margin,
    max_levels = max_levels, max_nodes = max_nodes,
    refine_unresolved_levels = refine_unresolved_levels,
    maximum_cell_half_width = maximum_cell_half_width,
    roundoff_inflation = roundoff_inflation,
    state_cache = augmented$state_cache,
    constant_cache = augmented$constant_cache
  )
  derivative_enclosure <- adaptive_augmented_enclosure(
    fit = fit, moment_type = moment_type, degree = 1L,
    initial_grid = fit$guard$grid,
    numerical_tolerance = numerical_tolerance, safety_margin = safety_margin,
    max_levels = max_levels, max_nodes = max_nodes,
    refine_unresolved_levels = refine_unresolved_levels,
    maximum_cell_half_width = maximum_cell_half_width,
    roundoff_inflation = roundoff_inflation,
    state_cache = augmented$state_cache,
    constant_cache = augmented$constant_cache
  )
  certificates <- if (fit$guard_failed) {
    derivative_assisted_certificates(
      fit, moment_type, derivative_enclosure$cells, numerical_tolerance,
      augmented$state_cache
    )[FALSE, , drop = FALSE]
  } else {
    derivative_assisted_certificates(
      fit, moment_type, derivative_enclosure$cells, numerical_tolerance,
      augmented$state_cache
    )
  }
  truth_flags <- derivative_certificate_truth(certificates, truth)
  covariance <- fit$object$covariance_01[[moment_type]][fit$indices]
  variance_0 <- value_band$variance_hat
  variance_1 <- derivative_band$variance_hat
  denominator <- sqrt(pmax(variance_0, 0) * pmax(variance_1, 0))
  correlation <- covariance / denominator
  correlation[!is.finite(denominator) | denominator <= 0] <- NA_real_
  list(
    fit = fit, value_band = value_band, derivative_band = derivative_band,
    value_enclosure = value_enclosure,
    derivative_enclosure = derivative_enclosure,
    certificates = certificates, certificate_truth = truth_flags,
    covariance_min = if (length(covariance)) min(covariance) else NA_real_,
    covariance_max = if (length(covariance)) max(covariance) else NA_real_,
    correlation_min = if (any(is.finite(correlation))) min(correlation, na.rm = TRUE) else NA_real_,
    correlation_max = if (any(is.finite(correlation))) max(correlation, na.rm = TRUE) else NA_real_
  )
}
