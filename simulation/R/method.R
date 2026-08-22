# Statistical engine for simultaneous inference on critical moment orders.
#
# The bootstrap is evaluated only on the prespecified master grid.  Numerical
# continuum certification and root isolation evaluate new endpoints lazily and
# cache scalar states; they never allocate an n by R adaptive-grid matrix.

power_matrix <- function(values, orders, moment_type) {
  if (!all(is.finite(values))) {
    stop_numerical("The sample contains a nonfinite observation.")
  }
  assert_true(all(is.finite(orders) & orders > 0),
              "All evaluated orders must be positive.")
  log_magnitude <- log(abs(values))
  out <- exp(outer(log_magnitude, orders, `*`))
  zero_rows <- values == 0
  if (any(zero_rows)) out[zero_rows, ] <- 0
  if (identical(moment_type, "signed")) {
    out <- out * sign(values)
  } else if (!identical(moment_type, "absolute")) {
    stopf("Unknown moment type: %s", moment_type)
  }
  if (any(!is.finite(out))) {
    stop_numerical(
      "A powered observation overflowed. Reduce the order interval or inspect the data."
    )
  }
  out
}

power_vector <- function(values, order, moment_type) {
  as.vector(power_matrix(values, order, moment_type))
}

power_derivative_vector <- function(values, order, moment_type, degree = 1L) {
  degree <- as.integer(degree)
  assert_true(length(degree) == 1L && !is.na(degree) && degree >= 0L,
              "degree must be a nonnegative integer.")
  powered <- power_vector(values, order, moment_type)
  if (degree == 0L) return(powered)
  log_magnitude <- log(abs(values))
  out <- powered * log_magnitude^degree
  out[values == 0] <- 0
  if (any(!is.finite(out))) {
    stop_numerical("A derivative contribution is nonfinite.")
  }
  out
}

biased_variance <- function(values) {
  centered <- values - mean(values)
  mean(centered^2)
}

build_feature_object <- function(x, y, orders, moment_types) {
  orders <- sort(unique(as.numeric(orders)))
  moment_types <- unique(as.character(moment_types))
  assert_true(length(orders) >= 2L && all(diff(orders) > 0),
              "The master order grid must be strictly increasing.")
  assert_true(all(moment_types %in% c("absolute", "signed")),
              "Only absolute and signed moments are implemented.")

  feature_x <- lapply(moment_types, function(type) power_matrix(x, orders, type))
  feature_y <- lapply(moment_types, function(type) power_matrix(y, orders, type))
  names(feature_x) <- names(feature_y) <- moment_types
  combined_x <- do.call(cbind, feature_x)
  combined_y <- do.call(cbind, feature_y)
  mean_x <- colMeans(combined_x)
  mean_y <- colMeans(combined_y)
  centered_x <- sweep(combined_x, 2L, mean_x, FUN = "-")
  centered_y <- sweep(combined_y, 2L, mean_y, FUN = "-")
  n_x <- length(x)
  n_y <- length(y)
  lambda <- n_x / (n_x + n_y)
  variance_x <- colMeans(centered_x^2)
  variance_y <- colMeans(centered_y^2)
  variance <- (1 - lambda) * variance_x + lambda * variance_y
  raw_second_moment <- (1 - lambda) * colMeans(combined_x^2) +
    lambda * colMeans(combined_y^2)
  if (any(!is.finite(mean_x)) || any(!is.finite(mean_y)) ||
      any(!is.finite(centered_x)) || any(!is.finite(centered_y)) ||
      any(!is.finite(variance)) || any(!is.finite(raw_second_moment))) {
    stop_numerical(
      "A feature mean, centered value, variance, or raw second moment is nonfinite."
    )
  }

  list(
    x = x, y = y, orders = orders, moment_types = moment_types,
    n_x = n_x, n_y = n_y, n_eff = n_x * n_y / (n_x + n_y),
    lambda = lambda, feature_x = feature_x, feature_y = feature_y,
    combined_x = combined_x, combined_y = combined_y,
    centered_x = centered_x, centered_y = centered_y,
    delta = mean_y - mean_x, variance = variance,
    raw_second_moment = raw_second_moment
  )
}

column_block <- function(object, moment_type) {
  type_index <- match(moment_type, object$moment_types)
  assert_true(!is.na(type_index), sprintf("Moment type '%s' was not fitted.", moment_type))
  grid_size <- length(object$orders)
  seq.int((type_index - 1L) * grid_size + 1L, length.out = grid_size)
}

order_key <- function(order) formatC(order, digits = 17L, format = "fg", flag = "#")

normalise_order_intervals <- function(order_intervals, orders) {
  order_range <- range(orders)
  if (is.null(order_intervals)) {
    if (order_range[1L] <= 0.25 + 1e-12 && order_range[2L] >= 2.50 - 1e-12) {
      order_intervals <- list(K1 = c(0.25, 1.50), K2 = c(0.25, 2.00),
                              K3 = c(0.25, 2.50))
    } else {
      order_intervals <- list(K = order_range)
    }
  }
  if (is.matrix(order_intervals) || is.data.frame(order_intervals)) {
    matrix_intervals <- as.matrix(order_intervals)
    assert_true(ncol(matrix_intervals) == 2L,
                "An interval matrix must have exactly two columns.")
    order_intervals <- lapply(seq_len(nrow(matrix_intervals)),
                              function(i) as.numeric(matrix_intervals[i, ]))
  }
  if (is.numeric(order_intervals) && length(order_intervals) == 2L) {
    order_intervals <- list(K = order_intervals)
  }
  assert_true(is.list(order_intervals) && length(order_intervals) >= 1L,
              "order_intervals must be a nonempty list of endpoint pairs.")
  if (is.null(names(order_intervals)) || any(!nzchar(names(order_intervals)))) {
    names(order_intervals) <- paste0("K", seq_along(order_intervals))
  }
  assert_true(!anyDuplicated(names(order_intervals)), "Interval identifiers must be unique.")
  out <- lapply(order_intervals, function(interval) {
    interval <- as.numeric(interval)
    assert_true(length(interval) == 2L && all(is.finite(interval)) &&
                  interval[1L] < interval[2L],
                "Every order interval must have two finite increasing endpoints.")
    assert_true(interval[1L] >= order_range[1L] - 1e-12 &&
                  interval[2L] <= order_range[2L] + 1e-12,
                "Every order interval must lie inside the master grid range.")
    interval
  })
  names(out) <- names(order_intervals)
  out
}

regular_interval_grid <- function(interval, spacing = NULL, size = NULL) {
  if (!is.null(spacing)) {
    spacing <- as.numeric(spacing)
    assert_true(length(spacing) == 1L && is.finite(spacing) && spacing > 0,
                "Grid spacing must be a positive finite scalar.")
    count <- floor((interval[2L] - interval[1L]) / spacing + 1e-10)
    grid <- interval[1L] + spacing * seq.int(0L, count)
    grid <- grid[grid <= interval[2L] + 1e-12]
    if (!length(grid) || abs(tail(grid, 1L) - interval[2L]) > 1e-12) {
      grid <- c(grid, interval[2L])
    } else {
      grid[length(grid)] <- interval[2L]
    }
    return(sort(unique(c(interval[1L], grid, interval[2L]))))
  }
  size <- as.integer(size)
  assert_true(length(size) == 1L && !is.na(size) && size >= 2L,
              "Grid size must be at least two.")
  seq(interval[1L], interval[2L], length.out = size)
}

normalise_analysis_grid <- function(orders, order_intervals, label = "analysis") {
  orders <- as.numeric(orders)
  assert_true(
    length(orders) >= 2L && all(is.finite(orders)) && all(orders > 0),
    sprintf("The %s grid must contain at least two finite positive orders.", label)
  )
  endpoints <- as.numeric(unlist(order_intervals, use.names = FALSE))
  assert_true(
    length(endpoints) >= 2L && all(is.finite(endpoints)),
    "Order-interval endpoints are unavailable when constructing an analysis grid."
  )
  full_range <- range(endpoints)
  orders <- orders[
    orders >= full_range[1L] - 1e-12 & orders <= full_range[2L] + 1e-12
  ]
  grid <- sort(unique(c(orders, endpoints)))
  assert_true(
    length(grid) >= 2L && all(diff(grid) > 0),
    sprintf("The %s grid must contain at least two distinct increasing orders.", label)
  )
  grid
}

grid_on_interval <- function(orders, interval) {
  grid <- orders[
    orders >= interval[1L] - 1e-12 & orders <= interval[2L] + 1e-12
  ]
  grid <- sort(unique(c(interval, grid)))
  assert_true(length(grid) >= 2L && all(diff(grid) > 0),
              "An interval grid must contain two distinct increasing endpoints.")
  grid
}

maximum_grid_spacing <- function(orders) {
  orders <- sort(unique(as.numeric(orders)))
  if (length(orders) < 2L) NA_real_ else max(diff(orders))
}

draw_multipliers <- function(number, sample_size, distribution) {
  if (identical(distribution, "rademacher")) {
    matrix(2 * sample.int(2L, number * sample_size, replace = TRUE) - 3,
           nrow = number, ncol = sample_size)
  } else if (identical(distribution, "normal")) {
    matrix(stats::rnorm(number * sample_size), nrow = number, ncol = sample_size)
  } else {
    stopf("Unknown multiplier distribution: %s", distribution)
  }
}

draw_multiplier_pairs <- function(number, n_x, n_y, distribution) {
  number <- as.integer(number)
  total_size <- as.integer(n_x + n_y)
  if (identical(distribution, "rademacher")) {
    values <- 2 * sample.int(
      2L, number * total_size, replace = TRUE
    ) - 3
  } else if (identical(distribution, "normal")) {
    values <- stats::rnorm(number * total_size)
  } else {
    stopf("Unknown multiplier distribution: %s", distribution)
  }
  # Each row is one complete (X,Y) multiplier pair.  Row-major filling makes
  # the first b pairs invariant to both the requested total B and the size of
  # the final batch, which is required by the nested finite-B sensitivity run.
  paired <- matrix(values, nrow = number, ncol = total_size, byrow = TRUE)
  list(
    x = paired[, seq_len(n_x), drop = FALSE],
    y = paired[, n_x + seq_len(n_y), drop = FALSE]
  )
}

make_state_cache <- function(object) {
  cache <- new.env(parent = emptyenv(), hash = TRUE)
  for (moment_type in object$moment_types) {
    block <- column_block(object, moment_type)
    for (index in seq_along(object$orders)) {
      assign(paste(moment_type, order_key(object$orders[index]), sep = "::"),
             list(delta = object$delta[block[index]],
                  variance = object$variance[block[index]],
                  raw_second_moment = object$raw_second_moment[block[index]],
                  derivative = NA_real_), envir = cache)
    }
  }
  cache
}

feature_state_at <- function(object, order, moment_type, cache = NULL,
                             need_derivative = FALSE) {
  key <- paste(moment_type, order_key(order), sep = "::")
  if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
    state <- get(key, envir = cache, inherits = FALSE)
    if (!need_derivative || is.finite(state$derivative)) return(state)
  }
  value_x <- power_vector(object$x, order, moment_type)
  value_y <- power_vector(object$y, order, moment_type)
  mean_x <- mean(value_x)
  mean_y <- mean(value_y)
  variance <- (1 - object$lambda) * mean((value_x - mean_x)^2) +
    object$lambda * mean((value_y - mean_y)^2)
  raw_second <- (1 - object$lambda) * mean(value_x^2) +
    object$lambda * mean(value_y^2)
  derivative <- if (need_derivative) {
    mean(power_derivative_vector(object$y, order, moment_type)) -
      mean(power_derivative_vector(object$x, order, moment_type))
  } else {
    NA_real_
  }
  state <- list(delta = mean_y - mean_x, variance = variance,
                raw_second_moment = raw_second, derivative = derivative)
  finite_state <- c(state$delta, state$variance, state$raw_second_moment)
  if (need_derivative) finite_state <- c(finite_state, state$derivative)
  if (any(!is.finite(finite_state))) {
    stop_numerical("An off-grid feature aggregate is nonfinite.")
  }
  if (!is.null(cache)) assign(key, state, envir = cache)
  state
}

sample_local_envelope <- function(values, left, right) {
  magnitude <- abs(values)
  log_magnitude <- log(magnitude)
  left_power <- exp(left * log_magnitude)
  right_power <- exp(right * log_magnitude)
  envelope <- pmax(left_power, right_power)
  envelope[magnitude == 0] <- 0
  log_absolute <- abs(log_magnitude)
  log_absolute[magnitude == 0] <- 0
  out <- c(
    mean_h = mean(envelope),
    mean_h1 = mean(envelope * log_absolute),
    mean_h2 = mean(envelope * log_absolute^2),
    mean_h_h1 = mean(envelope^2 * log_absolute)
  )
  if (any(!is.finite(out))) {
    stop_numerical("A local envelope bound is nonfinite.")
  }
  out
}

local_cell_constants <- function(object, left, right, cache = NULL) {
  key <- paste(order_key(left), order_key(right), sep = "::")
  if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  bound_x <- sample_local_envelope(object$x, left, right)
  bound_y <- sample_local_envelope(object$y, left, right)
  variance_lipschitz_x <- 2 *
    (bound_x[["mean_h_h1"]] + bound_x[["mean_h"]] * bound_x[["mean_h1"]])
  variance_lipschitz_y <- 2 *
    (bound_y[["mean_h_h1"]] + bound_y[["mean_h"]] * bound_y[["mean_h1"]])
  out <- list(
    delta_lipschitz = bound_x[["mean_h1"]] + bound_y[["mean_h1"]],
    derivative_lipschitz = bound_x[["mean_h2"]] + bound_y[["mean_h2"]],
    variance_lipschitz = (1 - object$lambda) * variance_lipschitz_x +
      object$lambda * variance_lipschitz_y,
    second_moment_lipschitz = 2 * ((1 - object$lambda) *
      bound_x[["mean_h_h1"]] + object$lambda * bound_y[["mean_h_h1"]])
  )
  if (!is.null(cache)) assign(key, out, envir = cache)
  out
}

validate_relative_guard <- function(object, interval, moment_types,
                                    tolerance = 1e-12,
                                    initial_grid = NULL,
                                    max_levels = 6L,
                                    max_nodes = 28801L,
                                    state_cache = NULL,
                                    constant_cache = NULL) {
  if (is.null(initial_grid)) {
    initial_grid <- object$orders[object$orders >= interval[1L] - 1e-12 &
                                    object$orders <= interval[2L] + 1e-12]
  }
  initial_grid <- sort(unique(c(interval, initial_grid)))
  initial_grid <- initial_grid[initial_grid >= interval[1L] - 1e-12 &
                                 initial_grid <= interval[2L] + 1e-12]
  assert_true(length(initial_grid) >= 2L, "The guard grid has fewer than two nodes.")
  max_levels <- as.integer(max_levels)
  max_nodes <- as.integer(max_nodes)
  assert_true(max_nodes >= length(initial_grid),
              "max_nodes is smaller than the initial guard grid.")

  nodes <- initial_grid
  queue <- lapply(seq_len(length(initial_grid) - 1L), function(index) {
    list(left = initial_grid[index], right = initial_grid[index + 1L], level = 0L)
  })
  accepted <- list()
  failures <- list()
  head_index <- 1L
  minimum_ratio <- Inf

  while (head_index <= length(queue)) {
    cell <- queue[[head_index]]
    head_index <- head_index + 1L
    radius <- (cell$right - cell$left) / 2
    constants <- local_cell_constants(object, cell$left, cell$right, constant_cache)
    ratios <- vapply(moment_types, function(moment_type) {
      left_state <- feature_state_at(object, cell$left, moment_type, state_cache)
      right_state <- feature_state_at(object, cell$right, moment_type, state_cache)
      variance_lower <- min(left_state$variance, right_state$variance) -
        constants$variance_lipschitz * radius
      second_upper <- max(left_state$raw_second_moment,
                          right_state$raw_second_moment) +
        constants$second_moment_lipschitz * radius
      if (!is.finite(variance_lower) || !is.finite(second_upper) ||
          variance_lower <= 0 || second_upper <= 0) 0 else variance_lower / second_upper
    }, numeric(1L))
    if (all(is.finite(ratios) & ratios > tolerance)) {
      minimum_ratio <- min(minimum_ratio, ratios)
      accepted[[length(accepted) + 1L]] <- cell
      next
    }
    can_split <- cell$level < max_levels && length(nodes) < max_nodes
    if (can_split) {
      midpoint <- (cell$left + cell$right) / 2
      if (!any(abs(nodes - midpoint) <= .Machine$double.eps * max(1, abs(midpoint)))) {
        nodes <- c(nodes, midpoint)
      }
      queue[[length(queue) + 1L]] <-
        list(left = cell$left, right = midpoint, level = cell$level + 1L)
      queue[[length(queue) + 1L]] <-
        list(left = midpoint, right = cell$right, level = cell$level + 1L)
    } else {
      minimum_ratio <- min(minimum_ratio, ratios)
      failures[[length(failures) + 1L]] <- cell
    }
  }

  partition <- c(accepted, failures)
  grid <- sort(unique(unlist(lapply(partition, function(cell) c(cell$left, cell$right)))))
  list(
    valid = length(failures) == 0L,
    minimum_certified_ratio = minimum_ratio,
    grid = grid,
    nodes_used = length(nodes),
    maximum_level = if (length(partition)) max(vapply(partition, `[[`, integer(1L), "level")) else 0L,
    unresolved_cells = length(failures),
    limit_hit = length(failures) > 0L
  )
}

subset_feature_object <- function(object, order_indices) {
  order_indices <- as.integer(order_indices)
  full_size <- length(object$orders)
  combined_indices <- unlist(lapply(seq_along(object$moment_types), function(type_index) {
    (type_index - 1L) * full_size + order_indices
  }), use.names = FALSE)
  out <- object
  out$orders <- object$orders[order_indices]
  out$feature_x <- lapply(object$feature_x, function(value) value[, order_indices, drop = FALSE])
  out$feature_y <- lapply(object$feature_y, function(value) value[, order_indices, drop = FALSE])
  out$combined_x <- object$combined_x[, combined_indices, drop = FALSE]
  out$combined_y <- object$combined_y[, combined_indices, drop = FALSE]
  out$centered_x <- object$centered_x[, combined_indices, drop = FALSE]
  out$centered_y <- object$centered_y[, combined_indices, drop = FALSE]
  out$delta <- object$delta[combined_indices]
  out$variance <- object$variance[combined_indices]
  out$raw_second_moment <- object$raw_second_moment[combined_indices]
  out
}

nested_multiplier_bands <- function(
    x, y, orders, moment_types, order_intervals,
    bootstrap_reps, alpha,
    multiplier_distribution = "rademacher",
    bootstrap_batch_size = 50L,
    variance_tolerance = 1e-12,
    guard_max_levels = 6L,
    guard_max_nodes = 28801L) {
  order_intervals <- normalise_order_intervals(order_intervals, orders)
  master_orders <- sort(unique(c(as.numeric(orders), unlist(order_intervals))))
  object <- build_feature_object(x, y, master_orders, moment_types)
  state_cache <- make_state_cache(object)
  constant_cache <- new.env(parent = emptyenv(), hash = TRUE)

  interval_indices <- lapply(order_intervals, function(interval) {
    which(object$orders >= interval[1L] - 1e-12 &
            object$orders <= interval[2L] + 1e-12)
  })
  guard_results <- Map(function(interval, indices) {
    validate_relative_guard(
      object = object, interval = interval, moment_types = moment_types,
      tolerance = variance_tolerance, initial_grid = object$orders[indices],
      max_levels = guard_max_levels, max_nodes = guard_max_nodes,
      state_cache = state_cache, constant_cache = constant_cache
    )
  }, order_intervals, interval_indices)
  names(guard_results) <- names(order_intervals)

  bootstrap_reps <- as.integer(bootstrap_reps)
  bootstrap_batch_size <- min(as.integer(bootstrap_batch_size), bootstrap_reps)
  assert_true(bootstrap_reps >= 2L, "At least two multiplier replications are required.")
  standard_deviation <- sqrt(object$variance)
  safe_standard_deviation <- standard_deviation
  invalid_columns <- !is.finite(safe_standard_deviation) | safe_standard_deviation <= 0
  safe_standard_deviation[invalid_columns] <- 1
  suprema <- matrix(NA_real_, nrow = bootstrap_reps, ncol = length(order_intervals),
                    dimnames = list(NULL, names(order_intervals)))

  columns_by_interval <- lapply(interval_indices, function(indices) {
    unlist(lapply(seq_along(moment_types), function(type_index) {
      (type_index - 1L) * length(object$orders) + indices
    }), use.names = FALSE)
  })
  next_index <- 1L
  while (next_index <= bootstrap_reps) {
    batch_size <- min(bootstrap_batch_size, bootstrap_reps - next_index + 1L)
    multiplier_pair <- draw_multiplier_pairs(
      batch_size, object$n_x, object$n_y, multiplier_distribution
    )
    multiplier_x <- multiplier_pair$x
    multiplier_y <- multiplier_pair$y
    process_y <- sqrt(object$lambda / object$n_y) *
      (multiplier_y %*% object$centered_y)
    process_x <- sqrt((1 - object$lambda) / object$n_x) *
      (multiplier_x %*% object$centered_x)
    standardized <- abs(process_y - process_x)
    standardized <- sweep(standardized, 2L, safe_standard_deviation, FUN = "/")
    if (any(invalid_columns)) standardized[, invalid_columns] <- Inf
    rows <- next_index:(next_index + batch_size - 1L)
    for (interval_id in names(order_intervals)) {
      if (isTRUE(guard_results[[interval_id]]$valid)) {
        suprema[rows, interval_id] <-
          apply(standardized[, columns_by_interval[[interval_id]], drop = FALSE], 1L, max)
      }
    }
    next_index <- next_index + batch_size
  }

  quantile_index <- min(
    bootstrap_reps, ceiling((bootstrap_reps + 1) * (1 - alpha))
  )
  interval_fits <- lapply(names(order_intervals), function(interval_id) {
    indices <- interval_indices[[interval_id]]
    fit <- subset_feature_object(object, indices)
    guard <- guard_results[[interval_id]]
    fit$interval_id <- interval_id
    fit$interval <- order_intervals[[interval_id]]
    fit$guard <- guard
    fit$guard_failed <- !isTRUE(guard$valid)
    fit$relative_variance <- fit$variance / fit$raw_second_moment
    fit$relative_variance[fit$raw_second_moment == 0] <- 0
    if (fit$guard_failed) {
      fit$critical_value <- 0
      fit$half_width <- rep(Inf, length(fit$delta))
      fit$lower <- rep(-Inf, length(fit$delta))
      fit$upper <- rep(Inf, length(fit$delta))
      fit$bootstrap_supremum <- rep(NA_real_, bootstrap_reps)
    } else {
      values <- suprema[, interval_id]
      critical_value <- sort.int(values, partial = quantile_index)[quantile_index]
      if (!is.finite(critical_value)) {
        stop_numerical("The multiplier critical value is nonfinite.")
      }
      fit$critical_value <- critical_value
      fit$half_width <- critical_value * sqrt(pmax(fit$variance, 0) / fit$n_eff)
      fit$lower <- fit$delta - fit$half_width
      fit$upper <- fit$delta + fit$half_width
      fit$bootstrap_supremum <- values
    }
    fit
  })
  names(interval_fits) <- names(order_intervals)

  list(
    master = object, intervals = order_intervals, fits = interval_fits,
    guard_results = guard_results, bootstrap_suprema = suprema,
    state_cache = state_cache, constant_cache = constant_cache,
    bootstrap_reps = bootstrap_reps, alpha = alpha,
    multiplier_distribution = multiplier_distribution
  )
}

# Backward-compatible one-interval wrapper.
multiplier_band <- function(
    x, y, orders, moment_types, bootstrap_reps, alpha,
    multiplier_distribution = "rademacher", bootstrap_batch_size = 50L,
    variance_tolerance = 1e-12, guard_max_levels = 6L,
    guard_max_nodes = 28801L) {
  nested <- nested_multiplier_bands(
    x = x, y = y, orders = orders, moment_types = moment_types,
    order_intervals = list(K = range(orders)), bootstrap_reps = bootstrap_reps,
    alpha = alpha, multiplier_distribution = multiplier_distribution,
    bootstrap_batch_size = bootstrap_batch_size,
    variance_tolerance = variance_tolerance,
    guard_max_levels = guard_max_levels, guard_max_nodes = guard_max_nodes
  )
  fit <- nested$fits[[1L]]
  fit$nested_fit <- nested
  fit
}

extract_band <- function(fit, moment_type) {
  columns <- column_block(fit, moment_type)
  list(
    orders = fit$orders, delta_hat = fit$delta[columns],
    variance_hat = fit$variance[columns], half_width = fit$half_width[columns],
    lower = fit$lower[columns], upper = fit$upper[columns]
  )
}

empirical_delta_at <- function(x, y, order, moment_type) {
  mean(power_vector(y, order, moment_type)) -
    mean(power_vector(x, order, moment_type))
}

empirical_derivative_at <- function(x, y, order, moment_type) {
  mean(power_derivative_vector(y, order, moment_type)) -
    mean(power_derivative_vector(x, order, moment_type))
}

empirical_variance_at <- function(x, y, order, moment_type) {
  n_x <- length(x)
  n_y <- length(y)
  lambda <- n_x / (n_x + n_y)
  value_x <- power_vector(x, order, moment_type)
  value_y <- power_vector(y, order, moment_type)
  (1 - lambda) * biased_variance(value_x) + lambda * biased_variance(value_y)
}

band_at_order <- function(fit, order, moment_type, state_cache = NULL) {
  state <- feature_state_at(fit, order, moment_type, state_cache)
  half_width <- if (fit$guard_failed || !is.finite(state$variance) || state$variance < 0) {
    Inf
  } else {
    fit$critical_value * sqrt(state$variance / fit$n_eff)
  }
  c(lower = state$delta - half_width, upper = state$delta + half_width,
    delta_hat = state$delta, variance_hat = state$variance,
    half_width = half_width)
}

evaluate_fitted_band <- function(fit, orders, moment_type, chunk_size = 256L) {
  orders <- as.numeric(orders)
  chunk_size <- max(1L, as.integer(chunk_size))
  delta_hat <- variance <- numeric(length(orders))
  starts <- seq.int(1L, length(orders), by = chunk_size)
  for (start in starts) {
    indices <- start:min(length(orders), start + chunk_size - 1L)
    feature_x <- power_matrix(fit$x, orders[indices], moment_type)
    feature_y <- power_matrix(fit$y, orders[indices], moment_type)
    mean_x <- colMeans(feature_x)
    mean_y <- colMeans(feature_y)
    variance_x <- colMeans(sweep(feature_x, 2L, mean_x, FUN = "-")^2)
    variance_y <- colMeans(sweep(feature_y, 2L, mean_y, FUN = "-")^2)
    local_delta <- mean_y - mean_x
    local_variance <- (1 - fit$lambda) * variance_x + fit$lambda * variance_y
    if (any(!is.finite(local_delta)) || any(!is.finite(local_variance))) {
      stop_numerical("A fitted-band aggregate is nonfinite.")
    }
    delta_hat[indices] <- local_delta
    variance[indices] <- local_variance
  }
  half_width <- if (fit$guard_failed) rep(Inf, length(orders)) else
    fit$critical_value * sqrt(pmax(variance, 0) / fit$n_eff)
  list(orders = orders, delta_hat = delta_hat, variance_hat = variance,
       half_width = half_width, lower = delta_hat - half_width,
       upper = delta_hat + half_width)
}

vacuous_cell_table <- function(grid, reason = "vacuous") {
  data.frame(
    left = head(grid, -1L), right = tail(grid, -1L), level = 0L,
    lower_bound = -Inf, upper_bound = Inf,
    variance_lower_bound = NA_real_, band_lipschitz = Inf,
    positive = FALSE, negative = FALSE, outer = TRUE,
    limit_hit = FALSE, depth_limit_hit = FALSE,
    node_limit_hit = FALSE, variance_unresolved_at_limit = FALSE,
    statistical_unresolved_at_limit = FALSE,
    reason = reason, stringsAsFactors = FALSE
  )
}

enclosure_limit_severity <- function(cells, interval, selector) {
  selector <- as.logical(selector)
  assert_true(length(selector) == nrow(cells) && !anyNA(selector),
              "An enclosure-limit selector is malformed.")
  widths <- if (nrow(cells)) cells$right - cells$left else numeric(0L)
  selected_widths <- widths[selector]
  interval_width <- diff(as.numeric(interval))
  assert_true(length(interval_width) == 1L && is.finite(interval_width) &&
                interval_width > 0,
              "The enclosure interval must have positive finite width.")
  list(
    cells = sum(selector),
    total_width = sum(selected_widths),
    width_proportion = sum(selected_widths) / interval_width,
    maximum_cell_width = if (length(selected_widths)) {
      max(selected_widths)
    } else 0
  )
}

adaptive_continuum_enclosure <- function(
    fit, moment_type, initial_grid = NULL, numerical_tolerance = 0,
    safety_margin = 0, max_levels = 6L, max_nodes = 28801L,
    refine_unresolved_levels = 0L, maximum_cell_half_width = Inf,
    roundoff_inflation = FALSE, state_cache = NULL, constant_cache = NULL) {
  interval <- fit$interval %||% range(fit$orders)
  if (is.null(initial_grid)) initial_grid <- fit$guard$grid %||% fit$orders
  initial_grid <- sort(unique(c(interval, initial_grid)))
  initial_grid <- initial_grid[initial_grid >= interval[1L] - 1e-12 &
                                 initial_grid <= interval[2L] + 1e-12]
  if (fit$guard_failed) {
    cells <- vacuous_cell_table(initial_grid, "relative_guard_failed")
    return(list(cells = cells, grid = initial_grid, nodes_used = length(initial_grid),
                maximum_level = 0L, limit_hit = FALSE,
                depth_limit_hit = FALSE, node_limit_hit = FALSE,
                variance_limit_hit = FALSE,
                statistical_limit_hit = FALSE,
                limit_hit_cells = 0L, limit_hit_total_width = 0,
                limit_hit_width_proportion = 0,
                limit_hit_maximum_cell_width = 0,
                depth_limit_cells = 0L, depth_limit_total_width = 0,
                depth_limit_width_proportion = 0,
                depth_limit_maximum_cell_width = 0,
                node_limit_cells = 0L, node_limit_total_width = 0,
                node_limit_width_proportion = 0,
                node_limit_maximum_cell_width = 0,
                variance_limit_cells = 0L, variance_limit_total_width = 0,
                variance_limit_width_proportion = 0,
                variance_limit_maximum_cell_width = 0,
                statistical_limit_cells = 0L,
                statistical_limit_total_width = 0,
                statistical_limit_width_proportion = 0,
                statistical_limit_maximum_cell_width = 0,
                variance_unresolved_cells = nrow(cells),
                statistically_unresolved_cells = nrow(cells)))
  }
  max_levels <- as.integer(max_levels)
  max_nodes <- as.integer(max_nodes)
  refine_unresolved_levels <- as.integer(refine_unresolved_levels)
  maximum_cell_half_width <- as.numeric(maximum_cell_half_width)
  assert_true(length(maximum_cell_half_width) == 1L &&
                !is.na(maximum_cell_half_width) &&
                maximum_cell_half_width > 0,
              "maximum_cell_half_width must be positive.")
  assert_true(max_nodes >= length(initial_grid),
              "max_nodes is smaller than the initial enclosure grid.")
  nodes <- initial_grid
  queue <- lapply(seq_len(length(initial_grid) - 1L), function(index) {
    list(left = initial_grid[index], right = initial_grid[index + 1L], level = 0L)
  })
  output <- list()
  head_index <- 1L
  limit_hit <- FALSE

  append_outer <- function(cell, variance_lower, band_lipschitz, reason,
                           limit_hit = FALSE, depth_limit_hit = FALSE,
                           node_limit_hit = FALSE,
                           variance_unresolved_at_limit = FALSE,
                           statistical_unresolved_at_limit = FALSE) {
    output[[length(output) + 1L]] <<- data.frame(
      left = cell$left, right = cell$right, level = cell$level,
      lower_bound = -Inf, upper_bound = Inf,
      variance_lower_bound = variance_lower, band_lipschitz = band_lipschitz,
      positive = FALSE, negative = FALSE, outer = TRUE,
      limit_hit = limit_hit, depth_limit_hit = depth_limit_hit,
      node_limit_hit = node_limit_hit,
      variance_unresolved_at_limit = variance_unresolved_at_limit,
      statistical_unresolved_at_limit = statistical_unresolved_at_limit,
      reason = reason, stringsAsFactors = FALSE
    )
  }
  split_cell <- function(cell) {
    midpoint <- (cell$left + cell$right) / 2
    if (!any(abs(nodes - midpoint) <= .Machine$double.eps * max(1, abs(midpoint)))) {
      nodes <<- c(nodes, midpoint)
    }
    queue[[length(queue) + 1L]] <<-
      list(left = cell$left, right = midpoint, level = cell$level + 1L)
    queue[[length(queue) + 1L]] <<-
      list(left = midpoint, right = cell$right, level = cell$level + 1L)
  }

  while (head_index <= length(queue)) {
    cell <- queue[[head_index]]
    head_index <- head_index + 1L
    radius <- (cell$right - cell$left) / 2
    constants <- local_cell_constants(fit, cell$left, cell$right, constant_cache)
    left_state <- feature_state_at(fit, cell$left, moment_type, state_cache)
    right_state <- feature_state_at(fit, cell$right, moment_type, state_cache)
    variance_lower <- min(left_state$variance, right_state$variance) -
      constants$variance_lipschitz * radius
    can_split <- cell$level < max_levels && length(nodes) < max_nodes
    if (!is.finite(variance_lower) || variance_lower <= 0) {
      if (can_split) {
        split_cell(cell)
      } else {
        limit_hit <- TRUE
        append_outer(
          cell, variance_lower, Inf, "variance_unresolved_at_limit",
          limit_hit = TRUE,
          depth_limit_hit = cell$level >= max_levels,
          node_limit_hit = length(nodes) >= max_nodes,
          variance_unresolved_at_limit = TRUE
        )
      }
      next
    }

    left_band <- band_at_order(fit, cell$left, moment_type, state_cache)
    right_band <- band_at_order(fit, cell$right, moment_type, state_cache)
    width_lipschitz <- fit$critical_value / sqrt(fit$n_eff) *
      constants$variance_lipschitz / (2 * sqrt(variance_lower))
    band_lipschitz <- constants$delta_lipschitz + width_lipschitz
    local_scale <- max(
      1, abs(left_band), abs(right_band),
      abs(band_lipschitz * radius), na.rm = TRUE
    )
    roundoff_margin <- if (isTRUE(roundoff_inflation)) {
      16 * .Machine$double.eps * local_scale
    } else 0
    total_margin <- safety_margin + roundoff_margin
    lower_bound <- min(left_band[["lower"]], right_band[["lower"]]) -
      band_lipschitz * radius - total_margin
    upper_bound <- max(left_band[["upper"]], right_band[["upper"]]) +
      band_lipschitz * radius + total_margin
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
    depth_limit_hit <- should_refine && !can_split &&
      cell$level >= max_levels
    node_limit_hit <- should_refine && !can_split &&
      length(nodes) >= max_nodes
    statistical_unresolved_at_limit <- should_refine && !can_split &&
      unresolved
    output[[length(output) + 1L]] <- data.frame(
      left = cell$left, right = cell$right, level = cell$level,
      lower_bound = lower_bound, upper_bound = upper_bound,
      variance_lower_bound = variance_lower, band_lipschitz = band_lipschitz,
      positive = positive, negative = negative, outer = unresolved,
      limit_hit = should_refine && !can_split,
      depth_limit_hit = depth_limit_hit,
      node_limit_hit = node_limit_hit,
      variance_unresolved_at_limit = FALSE,
      statistical_unresolved_at_limit = statistical_unresolved_at_limit,
      reason = if (positive) "certified_positive" else if (negative) {
        "certified_negative"
      } else if (should_refine && !can_split) {
        "statistically_unresolved_at_limit"
      } else {
        "statistically_unresolved"
      }, stringsAsFactors = FALSE
    )
  }

  cells <- do.call(rbind, output)
  cells <- cells[order(cells$left, cells$right), , drop = FALSE]
  row.names(cells) <- NULL
  grid <- sort(unique(c(cells$left, cells$right)))
  overall_severity <- enclosure_limit_severity(
    cells, interval, cells$limit_hit
  )
  depth_severity <- enclosure_limit_severity(
    cells, interval, cells$depth_limit_hit
  )
  node_severity <- enclosure_limit_severity(
    cells, interval, cells$node_limit_hit
  )
  variance_severity <- enclosure_limit_severity(
    cells, interval, cells$variance_unresolved_at_limit
  )
  statistical_severity <- enclosure_limit_severity(
    cells, interval, cells$statistical_unresolved_at_limit
  )
  list(
    cells = cells, grid = grid, nodes_used = length(nodes),
    maximum_level = if (nrow(cells)) max(cells$level) else 0L,
    limit_hit = limit_hit,
    depth_limit_hit = any(cells$depth_limit_hit),
    node_limit_hit = any(cells$node_limit_hit),
    variance_limit_hit = any(cells$variance_unresolved_at_limit),
    statistical_limit_hit = any(cells$statistical_unresolved_at_limit),
    limit_hit_cells = overall_severity$cells,
    limit_hit_total_width = overall_severity$total_width,
    limit_hit_width_proportion = overall_severity$width_proportion,
    limit_hit_maximum_cell_width = overall_severity$maximum_cell_width,
    depth_limit_cells = depth_severity$cells,
    depth_limit_total_width = depth_severity$total_width,
    depth_limit_width_proportion = depth_severity$width_proportion,
    depth_limit_maximum_cell_width = depth_severity$maximum_cell_width,
    node_limit_cells = node_severity$cells,
    node_limit_total_width = node_severity$total_width,
    node_limit_width_proportion = node_severity$width_proportion,
    node_limit_maximum_cell_width = node_severity$maximum_cell_width,
    variance_limit_cells = variance_severity$cells,
    variance_limit_total_width = variance_severity$total_width,
    variance_limit_width_proportion = variance_severity$width_proportion,
    variance_limit_maximum_cell_width = variance_severity$maximum_cell_width,
    statistical_limit_cells = statistical_severity$cells,
    statistical_limit_total_width = statistical_severity$total_width,
    statistical_limit_width_proportion =
      statistical_severity$width_proportion,
    statistical_limit_maximum_cell_width =
      statistical_severity$maximum_cell_width,
    variance_unresolved_cells = sum(cells$reason == "variance_unresolved_at_limit"),
    statistically_unresolved_cells = sum(cells$outer)
  )
}

scenario_absolute_moment <- function(scenario, sample_id, order) {
  moment_function <- scenario[[paste0("moment_", sample_id)]]
  if (!is.function(moment_function)) return(NA_real_)
  arguments <- list(order = order)
  if ("moment_type" %in% names(formals(moment_function))) {
    arguments$moment_type <- "absolute"
  }
  value <- tryCatch(do.call(moment_function, arguments), error = function(error) NA_real_)
  if (length(value) != 1L || !is.finite(value) || value < 0) NA_real_ else value
}

population_delta_lipschitz <- function(scenario, left, right) {
  domain <- as.numeric(scenario$moment_domain %||% c(0, Inf))
  if (length(domain) != 2L || is.na(domain[1L]) || is.na(domain[2L])) {
    domain <- c(0, Inf)
  }
  lower_room <- left - domain[1L]
  upper_room <- if (is.finite(domain[2L])) domain[2L] - right else Inf
  epsilon <- min(0.25, lower_room / 2, upper_room / 2)
  if (!is.finite(epsilon) || epsilon <= 0) return(Inf)
  lower_order <- left - epsilon
  upper_order <- right + epsilon
  moments <- c(
    scenario_absolute_moment(scenario, "x", lower_order),
    scenario_absolute_moment(scenario, "x", upper_order),
    scenario_absolute_moment(scenario, "y", lower_order),
    scenario_absolute_moment(scenario, "y", upper_order)
  )
  if (any(!is.finite(moments))) return(Inf)
  # For z <= 1, z^p |log z| <= z^(left-epsilon)/epsilon;
  # for z >= 1 it is <= z^(right+epsilon)/epsilon.  Summing the
  # two population bounds for X and Y controls |Delta'(p)| on the cell.
  sum(moments) / epsilon
}

certified_gap_lower_bound <- function(endpoint_margins, margin_lipschitz,
                                      radius, safety_margin = 0) {
  assert_true(length(endpoint_margins) == 2L &&
                all(is.finite(endpoint_margins)),
              "A cell gap needs two finite endpoint margins.")
  assert_true(length(margin_lipschitz) == 1L &&
                is.finite(margin_lipschitz) && margin_lipschitz >= 0,
              "The gap Lipschitz constant must be nonnegative and finite.")
  assert_true(length(radius) == 1L && is.finite(radius) && radius >= 0,
              "The cell radius must be nonnegative and finite.")
  assert_true(length(safety_margin) == 1L && is.finite(safety_margin) &&
                safety_margin >= 0,
              "The cell safety margin must be nonnegative and finite.")
  min(endpoint_margins) - margin_lipschitz * radius - safety_margin
}

certify_continuum_truth_coverage <- function(
    fit, scenario, moment_type, initial_grid = NULL,
    numerical_tolerance = 0, safety_margin = 0,
    max_levels = 6L, max_nodes = 28801L,
    roundoff_inflation = FALSE, state_cache = NULL,
    constant_cache = NULL, stop_on_endpoint_violation = TRUE) {
  interval <- fit$interval %||% range(fit$orders)
  if (is.null(initial_grid)) initial_grid <- fit$guard$grid %||% fit$orders
  initial_grid <- sort(unique(c(interval, initial_grid)))
  initial_grid <- initial_grid[initial_grid >= interval[1L] - 1e-12 &
                                 initial_grid <= interval[2L] + 1e-12]
  if (fit$guard_failed) {
    return(list(
      covered = TRUE, nodes_used = length(initial_grid), maximum_level = 0L,
      unresolved_cells = 0L, limit_hit = FALSE,
      endpoint_violation = FALSE
    ))
  }
  max_levels <- as.integer(max_levels)
  max_nodes <- as.integer(max_nodes)
  assert_true(max_nodes >= length(initial_grid),
              "max_nodes is smaller than the initial truth-coverage grid.")
  nodes <- initial_grid
  queue <- lapply(seq_len(length(initial_grid) - 1L), function(index) {
    list(left = initial_grid[index], right = initial_grid[index + 1L], level = 0L)
  })
  unresolved <- 0L
  limit_hit <- FALSE
  maximum_level <- 0L
  any_endpoint_violation <- FALSE
  head_index <- 1L

  while (head_index <= length(queue)) {
    cell <- queue[[head_index]]
    head_index <- head_index + 1L
    maximum_level <- max(maximum_level, cell$level)
    radius <- (cell$right - cell$left) / 2
    left_band <- band_at_order(fit, cell$left, moment_type, state_cache)
    right_band <- band_at_order(fit, cell$right, moment_type, state_cache)
    true_left <- scenario$delta(cell$left, moment_type)
    true_right <- scenario$delta(cell$right, moment_type)
    endpoint_violation <-
      left_band[["lower"]] > true_left + numerical_tolerance ||
      left_band[["upper"]] < true_left - numerical_tolerance ||
      right_band[["lower"]] > true_right + numerical_tolerance ||
      right_band[["upper"]] < true_right - numerical_tolerance
    any_endpoint_violation <- any_endpoint_violation || endpoint_violation
    # This scorer is used only because the population contrast is known in a
    # simulation.  A failed endpoint comparison is already an exact witness
    # that the simultaneous band does not cover that contrast, so further
    # subdivision cannot change the scientific coverage outcome.  The
    # optional exhaustive path is retained only for a focused regression test.
    if (endpoint_violation && isTRUE(stop_on_endpoint_violation)) {
      return(list(
        covered = FALSE, nodes_used = length(nodes),
        maximum_level = maximum_level, unresolved_cells = 0L,
        limit_hit = FALSE, endpoint_violation = TRUE
      ))
    }

    constants <- local_cell_constants(fit, cell$left, cell$right, constant_cache)
    left_state <- feature_state_at(fit, cell$left, moment_type, state_cache)
    right_state <- feature_state_at(fit, cell$right, moment_type, state_cache)
    variance_lower <- min(left_state$variance, right_state$variance) -
      constants$variance_lipschitz * radius
    population_lipschitz <- population_delta_lipschitz(
      scenario, cell$left, cell$right
    )
    certified <- FALSE
    if (is.finite(variance_lower) && variance_lower > 0 &&
        is.finite(population_lipschitz)) {
      width_lipschitz <- fit$critical_value / sqrt(fit$n_eff) *
        constants$variance_lipschitz / (2 * sqrt(variance_lower))
      band_lipschitz <- constants$delta_lipschitz + width_lipschitz
      margin_lipschitz <- band_lipschitz + population_lipschitz
      local_scale <- max(
        1, abs(left_band), abs(right_band), abs(true_left), abs(true_right),
        abs(margin_lipschitz * radius), na.rm = TRUE
      )
      roundoff_margin <- if (isTRUE(roundoff_inflation)) {
        16 * .Machine$double.eps * local_scale
      } else 0
      margin <- safety_margin + roundoff_margin
      # Coverage is a pointwise statement.  Comparing separate extrema of the
      # truth and the band is invalid because those extrema need not occur at
      # the same order.  Instead certify the two gap functions
      # Delta(p)-ell(p) and u(p)-Delta(p).  Each is Lipschitz with constant at
      # most population_lipschitz + band_lipschitz, and every point in the
      # cell lies within radius of at least one endpoint.
      lower_endpoint_margins <- c(
        true_left - left_band[["lower"]],
        true_right - right_band[["lower"]]
      )
      upper_endpoint_margins <- c(
        left_band[["upper"]] - true_left,
        right_band[["upper"]] - true_right
      )
      certified_lower_margin <- certified_gap_lower_bound(
        lower_endpoint_margins, margin_lipschitz, radius, margin
      )
      certified_upper_margin <- certified_gap_lower_bound(
        upper_endpoint_margins, margin_lipschitz, radius, margin
      )
      certified <- is.finite(certified_lower_margin) &&
        is.finite(certified_upper_margin) &&
        certified_lower_margin >= 0 && certified_upper_margin >= 0
    }
    if (certified) next

    can_split <- cell$level < max_levels && length(nodes) < max_nodes
    if (can_split) {
      midpoint <- cell$left + (cell$right - cell$left) / 2
      nodes <- c(nodes, midpoint)
      queue[[length(queue) + 1L]] <-
        list(left = cell$left, right = midpoint, level = cell$level + 1L)
      queue[[length(queue) + 1L]] <-
        list(left = midpoint, right = cell$right, level = cell$level + 1L)
    } else {
      unresolved <- unresolved + 1L
      limit_hit <- TRUE
    }
  }
  list(
    covered = unresolved == 0L,
    nodes_used = length(nodes), maximum_level = maximum_level,
    unresolved_cells = unresolved, limit_hit = limit_hit,
    endpoint_violation = any_endpoint_violation
  )
}

# Compatibility wrapper.  Unlike the first pilot, failure of one local cell
# does not make every cell vacuous.
continuum_cell_classification <- function(
    fit, band = NULL, moment_type, numerical_tolerance = 0,
    safety_margin = 0, max_levels = 6L, max_nodes = 28801L,
    refine_unresolved_levels = 0L, maximum_cell_half_width = Inf,
    roundoff_inflation = FALSE) {
  enclosure <- adaptive_continuum_enclosure(
    fit = fit, moment_type = moment_type,
    initial_grid = if (is.null(band)) fit$orders else band$orders,
    numerical_tolerance = numerical_tolerance, safety_margin = safety_margin,
    max_levels = max_levels, max_nodes = max_nodes,
    refine_unresolved_levels = refine_unresolved_levels,
    maximum_cell_half_width = maximum_cell_half_width,
    roundoff_inflation = roundoff_inflation
  )
  enclosure$cells
}

points_in_cell_union <- function(points, cells) {
  if (!length(points)) return(logical(0L))
  vapply(points, function(point) {
    any(cells$left <= point & point <= cells$right)
  }, logical(1L))
}

distance_to_interval_union <- function(point, cells) {
  if (!nrow(cells)) return(Inf)
  if (any(cells$left <= point & point <= cells$right)) return(0)
  min(pmin(abs(point - cells$left), abs(point - cells$right)))
}

hausdorff_finite_to_intervals <- function(roots, cells) {
  roots <- sort(unique(roots[is.finite(roots)]))
  if (!length(roots) && !nrow(cells)) return(0)
  if (!length(roots) || !nrow(cells)) return(Inf)
  root_to_set <- max(vapply(roots, distance_to_interval_union, numeric(1L), cells = cells))
  midpoints <- if (length(roots) >= 2L) (head(roots, -1L) + tail(roots, -1L)) / 2 else numeric(0L)
  set_to_root <- max(vapply(seq_len(nrow(cells)), function(index) {
    candidates <- c(cells$left[index], cells$right[index],
                    midpoints[midpoints >= cells$left[index] & midpoints <= cells$right[index]])
    max(vapply(candidates, function(value) min(abs(value - roots)), numeric(1L)))
  }, numeric(1L)))
  max(root_to_set, set_to_root)
}

outer_set_summary <- function(classification, true_root = NA_real_, true_roots = NULL) {
  if (is.null(true_roots)) true_roots <- true_root[is.finite(true_root)]
  outer_cells <- classification[classification$outer, c("left", "right"), drop = FALSE]
  outer <- classification$outer
  starts <- which(outer & !c(FALSE, head(outer, -1L)))
  ends <- which(outer & !c(tail(outer, -1L), FALSE))
  components <- if (length(starts)) data.frame(
    left = classification$left[starts], right = classification$right[ends]
  ) else data.frame(left = numeric(0L), right = numeric(0L))
  contains <- points_in_cell_union(true_roots, outer_cells)
  list(
    component_count = nrow(components),
    total_length = sum(outer_cells$right - outer_cells$left),
    contains_root = if (length(true_roots) == 1L) contains else NA,
    contains_all_roots = if (length(true_roots)) all(contains) else NA,
    covered_root_count = sum(contains),
    components = components,
    component_left = components$left, component_right = components$right,
    hausdorff_distance = hausdorff_finite_to_intervals(true_roots, components)
  )
}

tangency_global_geometry_summary <- function(
    outer, interval, root, n_eff, applicable, continuum_covered) {
  empty <- list(
    applicable = FALSE, coverage_condition_met = FALSE,
    outer_nonempty = NA, outer_contains_root = NA,
    component_count = NA_integer_, total_length = NA_real_,
    hausdorff = NA_real_, left_radius = NA_real_, right_radius = NA_real_,
    scaled_total_length = NA_real_, scaled_hausdorff = NA_real_,
    scaled_left_radius = NA_real_, scaled_right_radius = NA_real_,
    touches_K_left = NA, touches_K_right = NA, touches_K_boundary = NA,
    root_component_present = NA, root_component_left = NA_real_,
    root_component_right = NA_real_, root_component_length = NA_real_,
    root_component_left_radius = NA_real_,
    root_component_right_radius = NA_real_,
    scaled_root_component_length = NA_real_,
    scaled_root_component_left_radius = NA_real_,
    scaled_root_component_right_radius = NA_real_,
    root_component_touches_K_left = NA,
    root_component_touches_K_right = NA,
    root_component_touches_K_boundary = NA
  )
  if (!isTRUE(applicable)) return(empty)
  interval <- as.numeric(interval)
  assert_true(length(interval) == 2L && all(is.finite(interval)) &&
                interval[2L] > interval[1L],
              "Tangency geometry needs a finite nondegenerate K interval.")
  assert_true(length(root) == 1L && is.finite(root) &&
                root >= interval[1L] && root <= interval[2L],
              "Tangency geometry needs one root inside K.")
  assert_true(length(n_eff) == 1L && is.finite(n_eff) && n_eff > 0,
              "Tangency geometry needs a positive effective sample size.")
  components <- outer$components
  if (is.null(components)) {
    components <- data.frame(
      left = outer$component_left, right = outer$component_right
    )
  }
  assert_true(
    is.data.frame(components) &&
      all(c("left", "right") %in% names(components)) &&
      all(is.finite(c(components$left, components$right))) &&
      all(components$left <= components$right),
    "The numerical outer components are malformed."
  )
  tolerance <- 128 * .Machine$double.eps * max(1, abs(interval), abs(root))
  nonempty <- nrow(components) > 0L
  contains_root <- if (nonempty) {
    any(components$left <= root + tolerance &
          components$right >= root - tolerance)
  } else FALSE
  root_indices <- if (nonempty) which(
    components$left <= root + tolerance &
      components$right >= root - tolerance
  ) else integer(0L)
  assert_true(length(root_indices) <= 1L,
              "More than one disjoint outer component contains a root.")
  root_component <- if (length(root_indices)) {
    components[root_indices, , drop = FALSE]
  } else data.frame(left = numeric(0L), right = numeric(0L))
  global_left <- if (nonempty) min(components$left) else NA_real_
  global_right <- if (nonempty) max(components$right) else NA_real_
  left_radius <- if (nonempty) max(0, root - global_left) else NA_real_
  right_radius <- if (nonempty) max(0, global_right - root) else NA_real_
  root_left_radius <- if (nrow(root_component)) {
    max(0, root - root_component$left[1L])
  } else NA_real_
  root_right_radius <- if (nrow(root_component)) {
    max(0, root_component$right[1L] - root)
  } else NA_real_
  root_length <- if (nrow(root_component)) {
    root_component$right[1L] - root_component$left[1L]
  } else NA_real_
  scale <- n_eff^(1 / 4)
  condition_met <- isTRUE(continuum_covered)
  value_if_covered <- function(value) {
    if (condition_met) value else if (is.logical(value)) NA else NA_real_
  }
  touches_left <- nonempty && global_left <= interval[1L] + tolerance
  touches_right <- nonempty && global_right >= interval[2L] - tolerance
  root_touches_left <- nrow(root_component) > 0L &&
    root_component$left[1L] <= interval[1L] + tolerance
  root_touches_right <- nrow(root_component) > 0L &&
    root_component$right[1L] >= interval[2L] - tolerance
  list(
    applicable = TRUE, coverage_condition_met = condition_met,
    outer_nonempty = value_if_covered(nonempty),
    outer_contains_root = value_if_covered(contains_root),
    component_count = value_if_covered(as.integer(nrow(components))),
    total_length = value_if_covered(outer$total_length),
    hausdorff = value_if_covered(outer$hausdorff_distance),
    left_radius = value_if_covered(left_radius),
    right_radius = value_if_covered(right_radius),
    scaled_total_length = value_if_covered(scale * outer$total_length),
    scaled_hausdorff = value_if_covered(scale * outer$hausdorff_distance),
    scaled_left_radius = value_if_covered(scale * left_radius),
    scaled_right_radius = value_if_covered(scale * right_radius),
    touches_K_left = value_if_covered(touches_left),
    touches_K_right = value_if_covered(touches_right),
    touches_K_boundary = value_if_covered(touches_left || touches_right),
    root_component_present = value_if_covered(nrow(root_component) > 0L),
    root_component_left = value_if_covered(if (nrow(root_component)) {
      root_component$left[1L]
    } else NA_real_),
    root_component_right = value_if_covered(if (nrow(root_component)) {
      root_component$right[1L]
    } else NA_real_),
    root_component_length = value_if_covered(root_length),
    root_component_left_radius = value_if_covered(root_left_radius),
    root_component_right_radius = value_if_covered(root_right_radius),
    scaled_root_component_length = value_if_covered(scale * root_length),
    scaled_root_component_left_radius = value_if_covered(
      scale * root_left_radius
    ),
    scaled_root_component_right_radius = value_if_covered(
      scale * root_right_radius
    ),
    root_component_touches_K_left = value_if_covered(root_touches_left),
    root_component_touches_K_right = value_if_covered(root_touches_right),
    root_component_touches_K_boundary = value_if_covered(
      root_touches_left || root_touches_right
    )
  )
}

intersect_outer_cells <- function(classification, interval) {
  cells <- classification[classification$outer, c("left", "right"), drop = FALSE]
  if (!nrow(cells)) return(cells)
  cells$left <- pmax(cells$left, interval[1L])
  cells$right <- pmin(cells$right, interval[2L])
  cells[cells$left <= cells$right, , drop = FALSE]
}

tangency_localization_summary <- function(
    classification, neighborhood, root, n_eff, applicable) {
  if (!isTRUE(applicable)) {
    return(list(
      applicable = FALSE, local_outer_nonempty = NA,
      local_outer_total_length = NA_real_, scaled_radius = NA_real_,
      scaled_hausdorff = NA_real_, touches_left_boundary = NA,
      touches_right_boundary = NA, censored_by_window = NA
    ))
  }
  local_cells <- intersect_outer_cells(classification, neighborhood)
  nonempty <- nrow(local_cells) > 0L
  if (!nonempty) {
    return(list(
      applicable = TRUE, local_outer_nonempty = FALSE,
      local_outer_total_length = 0, scaled_radius = NA_real_,
      scaled_hausdorff = Inf, touches_left_boundary = FALSE,
      touches_right_boundary = FALSE, censored_by_window = FALSE
    ))
  }
  scale <- n_eff^(1 / 4)
  radius <- max(abs(c(local_cells$left, local_cells$right) - root))
  hausdorff <- hausdorff_finite_to_intervals(root, local_cells)
  tolerance <- 128 * .Machine$double.eps * max(
    1, abs(neighborhood), abs(root)
  )
  touches_left <- min(local_cells$left) <= neighborhood[1L] + tolerance
  touches_right <- max(local_cells$right) >= neighborhood[2L] - tolerance
  list(
    applicable = TRUE, local_outer_nonempty = TRUE,
    local_outer_total_length = sum(local_cells$right - local_cells$left),
    scaled_radius = scale * radius,
    scaled_hausdorff = scale * hausdorff,
    touches_left_boundary = touches_left,
    touches_right_boundary = touches_right,
    censored_by_window = touches_left || touches_right
  )
}

certified_sign_sequence <- function(orders, lower, upper, numerical_tolerance = 0) {
  status <- integer(length(orders))
  status[lower > numerical_tolerance] <- 1L
  status[upper < -numerical_tolerance] <- -1L
  indices <- which(status != 0L)
  if (!length(indices)) {
    return(list(count = 0L, indices = integer(0L), orders = numeric(0L),
                signs = integer(0L), directions = character(0L)))
  }
  keep <- c(TRUE, diff(status[indices]) != 0L)
  compressed_indices <- indices[keep]
  signs <- status[compressed_indices]
  directions <- if (length(signs) >= 2L) {
    ifelse(head(signs, -1L) < tail(signs, -1L), "up", "down")
  } else character(0L)
  list(count = max(0L, length(signs) - 1L), indices = compressed_indices,
       orders = orders[compressed_indices], signs = signs, directions = directions)
}

certified_alternation_count <- function(lower, upper, numerical_tolerance = 0) {
  certified_sign_sequence(seq_along(lower), lower, upper, numerical_tolerance)$count
}

find_certified_brackets <- function(orders, lower, upper, numerical_tolerance = 0) {
  sequence <- certified_sign_sequence(orders, lower, upper, numerical_tolerance)
  if (sequence$count == 0L) {
    return(data.frame(lower = numeric(0L), upper = numeric(0L),
                      direction = character(0L), stringsAsFactors = FALSE))
  }
  data.frame(
    lower = head(sequence$orders, -1L), upper = tail(sequence$orders, -1L),
    direction = sequence$directions, stringsAsFactors = FALSE
  )
}

find_certified_bracket <- function(orders, lower, upper, numerical_tolerance = 0) {
  brackets <- find_certified_brackets(orders, lower, upper, numerical_tolerance)
  if (!nrow(brackets)) {
    return(list(available = FALSE, lower = NA_real_, upper = NA_real_,
                direction = NA_character_))
  }
  index <- which.min(brackets$upper - brackets$lower)
  list(available = TRUE, lower = brackets$lower[index], upper = brackets$upper[index],
       direction = brackets$direction[index])
}

mellin_budget_decision <- function(
    alternation_count, structural_budget,
    true_distinct_root_count = NA_real_,
    true_total_multiplicity = true_distinct_root_count) {
  applicable <- length(structural_budget) == 1L && !is.na(structural_budget) &&
    is.finite(structural_budget) && structural_budget >= 0L
  rejected <- applicable && alternation_count > structural_budget
  exact_count_certified <- applicable && structural_budget > 0L &&
    alternation_count == structural_budget
  exact_count_statement_true <- applicable &&
    is.finite(true_distinct_root_count) &&
    is.finite(true_total_multiplicity) &&
    true_distinct_root_count == structural_budget &&
    true_total_multiplicity == structural_budget
  list(
    applicable = applicable, rejected = rejected,
    exact_count_certified = exact_count_certified,
    simplicity_certified = exact_count_certified,
    exact_count_statement_true = exact_count_statement_true
  )
}

evaluate_curve_state_chunked <- function(x, y, orders, moment_type,
                                         chunk_size = 256L) {
  orders <- as.numeric(orders)
  delta <- derivative <- variance <- numeric(length(orders))
  starts <- seq.int(1L, length(orders), by = max(1L, as.integer(chunk_size)))
  n_x <- length(x)
  n_y <- length(y)
  lambda <- n_x / (n_x + n_y)
  log_x <- log(abs(x))
  log_y <- log(abs(y))
  log_x[x == 0] <- 0
  log_y[y == 0] <- 0
  for (start in starts) {
    indices <- start:min(length(orders), start + chunk_size - 1L)
    feature_x <- power_matrix(x, orders[indices], moment_type)
    feature_y <- power_matrix(y, orders[indices], moment_type)
    mean_x <- colMeans(feature_x)
    mean_y <- colMeans(feature_y)
    local_delta <- mean_y - mean_x
    local_derivative <- colMeans(feature_y * log_y) - colMeans(feature_x * log_x)
    local_variance <- (1 - lambda) *
      colMeans(sweep(feature_x, 2L, mean_x, FUN = "-")^2) +
      lambda * colMeans(sweep(feature_y, 2L, mean_y, FUN = "-")^2)
    if (any(!is.finite(local_delta)) || any(!is.finite(local_derivative)) ||
        any(!is.finite(local_variance))) {
      stop_numerical("An empirical root-curve aggregate is nonfinite.")
    }
    delta[indices] <- local_delta
    derivative[indices] <- local_derivative
    variance[indices] <- local_variance
  }
  list(delta = delta, derivative = derivative, variance = variance)
}

isolate_empirical_roots <- function(
    x, y, interval, moment_type, initial_scan_size = 2001L,
    max_levels = 24L, max_evaluations = 50000L,
    root_tolerance = 1e-10, numerical_tolerance = 1e-12,
    chunk_size = 256L) {
  interval <- as.numeric(interval)
  assert_true(length(interval) == 2L && interval[1L] < interval[2L],
              "A root-isolation interval must have increasing endpoints.")
  initial_scan_size <- as.integer(initial_scan_size)
  max_levels <- as.integer(max_levels)
  max_evaluations <- as.integer(max_evaluations)
  assert_true(initial_scan_size >= 3L && max_evaluations >= initial_scan_size,
              "The root-isolation evaluation budget is too small.")
  object <- list(
    x = x, y = y, lambda = length(x) / (length(x) + length(y)),
    n_x = length(x), n_y = length(y),
    n_eff = length(x) * length(y) / (length(x) + length(y)),
    orders = numeric(0L), moment_types = moment_type
  )
  constant_cache <- new.env(parent = emptyenv(), hash = TRUE)
  state_cache <- new.env(parent = emptyenv(), hash = TRUE)
  initial_grid <- seq(interval[1L], interval[2L], length.out = initial_scan_size)
  initial_state <- evaluate_curve_state_chunked(
    x, y, initial_grid, moment_type, chunk_size = chunk_size
  )
  for (index in seq_along(initial_grid)) {
    assign(paste(moment_type, order_key(initial_grid[index]), sep = "::"),
           list(delta = initial_state$delta[index],
                variance = initial_state$variance[index],
                raw_second_moment = NA_real_,
                derivative = initial_state$derivative[index]), envir = state_cache)
  }
  evaluated_nodes <- initial_grid
  queue <- lapply(seq_len(initial_scan_size - 1L), function(index) {
    list(left = initial_grid[index], right = initial_grid[index + 1L], level = 0L)
  })
  roots <- numeric(0L)
  brackets <- list()
  unresolved <- list()
  evaluation_limit_hit <- FALSE
  head_index <- 1L

  state_at <- function(order) {
    key <- paste(moment_type, order_key(order), sep = "::")
    if (exists(key, envir = state_cache, inherits = FALSE)) {
      return(get(key, envir = state_cache, inherits = FALSE))
    }
    state <- feature_state_at(object, order, moment_type, state_cache,
                              need_derivative = TRUE)
    evaluated_nodes <<- c(evaluated_nodes, order)
    state
  }
  objective <- function(order) state_at(order)$delta

  bisect_bracket <- function(left, right, left_value, right_value) {
    if (left_value == 0) {
      return(list(root = left, lower = left, upper = left,
                  bisections = 0L, tolerance_met = TRUE,
                  evaluation_limit_hit = FALSE))
    }
    if (right_value == 0) {
      return(list(root = right, lower = right, upper = right,
                  bisections = 0L, tolerance_met = TRUE,
                  evaluation_limit_hit = FALSE))
    }
    assert_true(
      (left_value < 0 && right_value > 0) ||
        (left_value > 0 && right_value < 0),
      "Bisection requires a strict sign-changing bracket."
    )
    lower <- left
    upper <- right
    lower_value <- left_value
    upper_value <- right_value
    bisections <- 0L
    hit_evaluation_limit <- FALSE
    for (step in seq_len(max_levels)) {
      if (upper - lower <= root_tolerance) break
      if (length(unique(order_key(evaluated_nodes))) >= max_evaluations) {
        hit_evaluation_limit <- TRUE
        break
      }
      midpoint <- lower + (upper - lower) / 2
      midpoint_value <- objective(midpoint)
      bisections <- bisections + 1L
      if (midpoint_value == 0) {
        lower <- upper <- midpoint
        lower_value <- upper_value <- 0
        break
      }
      if ((lower_value < 0 && midpoint_value > 0) ||
          (lower_value > 0 && midpoint_value < 0)) {
        upper <- midpoint
        upper_value <- midpoint_value
      } else {
        lower <- midpoint
        lower_value <- midpoint_value
      }
    }
    list(
      root = lower + (upper - lower) / 2,
      lower = lower, upper = upper, bisections = bisections,
      tolerance_met = upper - lower <= root_tolerance,
      evaluation_limit_hit = hit_evaluation_limit
    )
  }

  boundary_state <- c(state_at(interval[1L])$delta, state_at(interval[2L])$delta)
  boundary_root <- any(abs(boundary_state) <= numerical_tolerance)

  while (head_index <= length(queue)) {
    cell <- queue[[head_index]]
    head_index <- head_index + 1L
    radius <- (cell$right - cell$left) / 2
    constants <- local_cell_constants(object, cell$left, cell$right, constant_cache)
    left_state <- state_at(cell$left)
    right_state <- state_at(cell$right)
    delta_lower <- min(left_state$delta, right_state$delta) -
      constants$delta_lipschitz * radius
    delta_upper <- max(left_state$delta, right_state$delta) +
      constants$delta_lipschitz * radius
    if (delta_lower > numerical_tolerance || delta_upper < -numerical_tolerance) next

    derivative_lower <- min(left_state$derivative, right_state$derivative) -
      constants$derivative_lipschitz * radius
    derivative_upper <- max(left_state$derivative, right_state$derivative) +
      constants$derivative_lipschitz * radius
    increasing <- derivative_lower > numerical_tolerance
    decreasing <- derivative_upper < -numerical_tolerance
    monotone <- increasing || decreasing
    opposite <- (left_state$delta <= 0 && right_state$delta >= 0) ||
      (left_state$delta >= 0 && right_state$delta <= 0)

    if (monotone && opposite) {
      isolated <- bisect_bracket(
        cell$left, cell$right, left_state$delta, right_state$delta
      )
      evaluation_limit_hit <- evaluation_limit_hit ||
        isolated$evaluation_limit_hit
      if (is.finite(isolated$root) && !isolated$evaluation_limit_hit) {
        roots <- c(roots, isolated$root)
        brackets[[length(brackets) + 1L]] <- data.frame(
          lower = isolated$lower, upper = isolated$upper,
          direction = if (increasing) "up" else "down",
          bisections = isolated$bisections,
          tolerance_met = isolated$tolerance_met,
          stringsAsFactors = FALSE
        )
        next
      }
    }
    if (monotone && !opposite) next

    can_split <- cell$level < max_levels &&
      length(unique(order_key(evaluated_nodes))) < max_evaluations
    if (can_split) {
      midpoint <- (cell$left + cell$right) / 2
      state_at(midpoint)
      queue[[length(queue) + 1L]] <-
        list(left = cell$left, right = midpoint, level = cell$level + 1L)
      queue[[length(queue) + 1L]] <-
        list(left = midpoint, right = cell$right, level = cell$level + 1L)
    } else {
      unresolved[[length(unresolved) + 1L]] <- cell
    }
  }

  bracket_table <- if (length(brackets)) do.call(rbind, brackets) else
    data.frame(lower = numeric(0L), upper = numeric(0L),
               direction = character(0L), bisections = integer(0L),
               tolerance_met = logical(0L), stringsAsFactors = FALSE)
  if (length(roots)) {
    ordering <- order(roots)
    roots <- roots[ordering]
    bracket_table <- bracket_table[ordering, , drop = FALSE]
    duplicate_tolerance <- max(2 * root_tolerance, 64 * .Machine$double.eps)
    keep <- c(TRUE, diff(roots) > duplicate_tolerance)
    roots <- roots[keep]
    bracket_table <- bracket_table[keep, , drop = FALSE]
  }
  tolerance_limit_hit <- nrow(bracket_table) > 0L &&
    !all(bracket_table$tolerance_met)
  limit_hit <- length(unresolved) > 0L || evaluation_limit_hit ||
    tolerance_limit_hit
  list(
    success = !boundary_root && !limit_hit && length(roots) == 1L,
    roots = roots, root_count = length(roots), brackets = bracket_table,
    boundary_root = boundary_root, unresolved_cells = length(unresolved),
    evaluations = length(unique(order_key(evaluated_nodes))),
    maximum_level = if (length(queue)) max(vapply(queue, `[[`, integer(1L), "level")) else 0L,
    maximum_bisections_used = if (nrow(bracket_table)) {
      max(bracket_table$bisections)
    } else 0L,
    all_root_tolerances_met = if (nrow(bracket_table)) {
      all(bracket_table$tolerance_met)
    } else NA,
    limit_hit = limit_hit
  )
}

# Compatibility root finder.  It now returns only roots certified by the
# interval algorithm, rather than roots inferred solely from a point scan.
find_empirical_roots <- function(x, y, interval, moment_type, scan_size = 2001L,
                                 max_levels = 24L, max_evaluations = 50000L) {
  isolation <- isolate_empirical_roots(
    x = x, y = y, interval = interval, moment_type = moment_type,
    initial_scan_size = scan_size, max_levels = max_levels,
    max_evaluations = max_evaluations
  )
  isolation$roots
}

wald_from_isolation <- function(object, isolation, moment_type, true_root,
                                alpha = 0.05, slope_tolerance = 1e-10,
                                variance_tolerance = 1e-12,
                                state_cache = NULL) {
  root_hat <- if (isTRUE(isolation$success) && length(isolation$roots) == 1L) {
    isolation$roots[1L]
  } else NA_real_
  derivative_hat <- variance_hat <- standard_error <- lower <- upper <- NA_real_
  reported <- FALSE
  if (is.finite(root_hat)) {
    state <- feature_state_at(object, root_hat, moment_type, state_cache,
                              need_derivative = TRUE)
    derivative_hat <- state$derivative
    variance_hat <- state$variance
    if (is.finite(derivative_hat) && abs(derivative_hat) > slope_tolerance &&
        is.finite(variance_hat) && variance_hat > variance_tolerance) {
      standard_error <- sqrt(variance_hat / object$n_eff) / abs(derivative_hat)
      quantile <- stats::qnorm(1 - alpha / 2)
      lower <- root_hat - quantile * standard_error
      upper <- root_hat + quantile * standard_error
      reported <- TRUE
    }
  }
  conditional_cover <- if (reported && is.finite(true_root)) {
    lower <= true_root && true_root <= upper
  } else NA
  formal_cover <- if (is.finite(true_root)) {
    if (reported) isTRUE(conditional_cover) else TRUE
  } else NA
  report_and_cover <- if (is.finite(true_root)) reported && isTRUE(conditional_cover) else NA
  list(
    reported = reported, root_hat = root_hat,
    derivative_hat = derivative_hat, variance_hat = variance_hat,
    standard_error = standard_error, lower = lower, upper = upper,
    length = if (reported) upper - lower else NA_real_,
    conditional_cover = conditional_cover,
    formal_unconditional_cover = formal_cover,
    report_and_cover = report_and_cover
  )
}

joint_root_multiplier_refinement <- function(
    object, root_details, moment_type, bootstrap_reps, alpha,
    multiplier_distribution, bootstrap_batch_size, bootstrap_seed,
    slope_tolerance = 1e-10, variance_tolerance = 1e-12) {
  detail_ids <- names(root_details)
  if (is.null(detail_ids)) detail_ids <- paste0("Q", seq_along(root_details))
  applicable <- length(root_details) > 0L && all(vapply(root_details, function(detail) {
    !detail$truth$identically_zero && length(detail$truth$roots) == 1L &&
      detail$truth$multiplicities[1L] == 1L &&
      !is.na(detail$truth$directions[1L])
  }, logical(1L)))
  empty_result <- function(reported = FALSE) {
    truth_finite <- applicable && all(vapply(root_details, function(detail) {
      is.finite(detail$truth$roots[1L])
    }, logical(1L)))
    list(
      applicable = applicable, reported = reported, critical_value = NA_real_,
      root_ids = detail_ids, root_hats = numeric(0L), lower = numeric(0L),
      upper = numeric(0L), conditional_cover = if (reported) FALSE else NA,
      formal_unconditional_cover = if (truth_finite) TRUE else NA,
      report_and_cover = if (truth_finite) FALSE else NA
    )
  }
  if (!applicable) return(empty_result())
  ready <- all(vapply(root_details, function(detail) {
    isTRUE(detail$isolation$success) && length(detail$isolation$roots) == 1L
  }, logical(1L)))
  if (!ready || is.null(bootstrap_seed)) return(empty_result())

  root_hats <- vapply(root_details, function(detail) {
    detail$isolation$roots[1L]
  }, numeric(1L))
  feature_x <- power_matrix(object$x, root_hats, moment_type)
  feature_y <- power_matrix(object$y, root_hats, moment_type)
  mean_x <- colMeans(feature_x)
  mean_y <- colMeans(feature_y)
  centered_x <- sweep(feature_x, 2L, mean_x, FUN = "-")
  centered_y <- sweep(feature_y, 2L, mean_y, FUN = "-")
  variance <- (1 - object$lambda) * colMeans(centered_x^2) +
    object$lambda * colMeans(centered_y^2)
  derivative <- vapply(root_hats, function(order) {
    empirical_derivative_at(object$x, object$y, order, moment_type)
  }, numeric(1L))
  regular <- all(is.finite(variance) & variance > variance_tolerance) &&
    all(is.finite(derivative) & abs(derivative) > slope_tolerance)
  if (!regular) return(empty_result())

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  assign(".Random.seed", bootstrap_seed, envir = .GlobalEnv)

  bootstrap_reps <- as.integer(bootstrap_reps)
  bootstrap_batch_size <- min(as.integer(bootstrap_batch_size), bootstrap_reps)
  maximum_statistics <- numeric(bootstrap_reps)
  next_index <- 1L
  while (next_index <= bootstrap_reps) {
    batch_size <- min(bootstrap_batch_size, bootstrap_reps - next_index + 1L)
    multiplier_pair <- draw_multiplier_pairs(
      batch_size, object$n_x, object$n_y, multiplier_distribution
    )
    process <- sqrt(object$lambda / object$n_y) *
      (multiplier_pair$y %*% centered_y) -
      sqrt((1 - object$lambda) / object$n_x) *
      (multiplier_pair$x %*% centered_x)
    standardized <- sweep(abs(process), 2L, sqrt(variance), FUN = "/")
    rows <- next_index:(next_index + batch_size - 1L)
    maximum_statistics[rows] <- apply(standardized, 1L, max)
    next_index <- next_index + batch_size
  }
  quantile_index <- min(
    bootstrap_reps, ceiling((bootstrap_reps + 1) * (1 - alpha))
  )
  critical_value <- sort.int(
    maximum_statistics, partial = quantile_index
  )[quantile_index]
  standard_error <- sqrt(variance / object$n_eff) / abs(derivative)
  lower <- root_hats - critical_value * standard_error
  upper <- root_hats + critical_value * standard_error
  true_roots <- vapply(root_details, function(detail) detail$truth$roots[1L], numeric(1L))
  conditional_cover <- all(lower <= true_roots & true_roots <= upper)
  list(
    applicable = TRUE, reported = TRUE, critical_value = critical_value,
    root_ids = detail_ids, root_hats = root_hats, lower = lower, upper = upper,
    conditional_cover = conditional_cover,
    formal_unconditional_cover = conditional_cover,
    report_and_cover = conditional_cover
  )
}

pointwise_band <- function(fit, orders, moment_type, alpha = 0.05,
                           chunk_size = 256L) {
  evaluated <- evaluate_fitted_band(fit, orders, moment_type, chunk_size)
  half_width <- stats::qnorm(1 - alpha / 2) *
    sqrt(pmax(evaluated$variance_hat, 0) / fit$n_eff)
  evaluated$half_width <- half_width
  evaluated$lower <- evaluated$delta_hat - half_width
  evaluated$upper <- evaluated$delta_hat + half_width
  evaluated
}

moment_weight_diagnostic <- function(values, order, top_fraction = 0.01) {
  magnitude <- abs(values)
  log_weight <- order * log(magnitude)
  log_weight[magnitude == 0] <- -Inf
  if (all(!is.finite(log_weight))) {
    return(c(ess = 0, ess_ratio = 0, maximum_share = NA_real_,
             top_share = NA_real_))
  }
  maximum <- max(log_weight[is.finite(log_weight)])
  weight <- exp(log_weight - maximum)
  weight[!is.finite(weight)] <- 0
  shares <- weight / sum(weight)
  top_number <- max(1L, ceiling(top_fraction * length(values)))
  c(
    ess = 1 / sum(shares^2),
    ess_ratio = 1 / (length(values) * sum(shares^2)),
    maximum_share = max(shares),
    top_share = sum(sort(shares, decreasing = TRUE)[seq_len(top_number)])
  )
}

extract_sample_values <- function(generated) {
  if (is.list(generated) && !is.null(generated$values)) {
    list(values = as.numeric(generated$values), metadata = generated[setdiff(names(generated), "values")])
  } else {
    list(values = as.numeric(generated), metadata = attributes(generated) %||% list())
  }
}

metadata_number <- function(metadata, candidates) {
  search_one_level <- function(value) {
    if (!is.list(value)) return(NA_real_)
    for (candidate in candidates) {
      candidate_value <- value[[candidate]]
      if (!is.null(candidate_value) && length(candidate_value) == 1L) {
        numeric_value <- suppressWarnings(as.numeric(candidate_value))
        if (length(numeric_value) == 1L && is.finite(numeric_value)) {
          return(numeric_value)
        }
      }
    }
    NA_real_
  }
  direct <- search_one_level(metadata)
  if (is.finite(direct)) return(direct)
  if (is.list(metadata)) {
    for (component in metadata) {
      nested <- search_one_level(component)
      if (is.finite(nested)) return(nested)
    }
  }
  NA_real_
}

scenario_abs_quantiles <- function(scenario, sample_id, probabilities) {
  function_names <- if (identical(sample_id, "x")) {
    c("magnitude_quantile_x", "q_abs_x", "abs_quantile_x", "q_x")
  } else c("magnitude_quantile_y", "q_abs_y", "abs_quantile_y", "q_y")
  for (name in function_names) {
    if (is.function(scenario[[name]])) return(as.numeric(scenario[[name]](probabilities)))
  }
  if (identical(scenario$family, "lognormal") && !is.null(scenario$parameters)) {
    suffix <- if (identical(sample_id, "x")) "x" else "y"
    return(stats::qlnorm(
      probabilities,
      meanlog = scenario$parameters[[paste0("meanlog_", suffix)]],
      sdlog = scenario$parameters[[paste0("sdlog_", suffix)]]
    ))
  }
  rep(NA_real_, length(probabilities))
}

tail_diagnostics <- function(x, y, upper_order, scenario = NULL,
                             quantile_probabilities = c(0.95, 0.99),
                             top_fraction = 0.01,
                             metadata_x = list(), metadata_y = list()) {
  x_p1 <- moment_weight_diagnostic(x, 1, top_fraction)
  y_p1 <- moment_weight_diagnostic(y, 1, top_fraction)
  x_upper <- moment_weight_diagnostic(x, upper_order, top_fraction)
  y_upper <- moment_weight_diagnostic(y, upper_order, top_fraction)
  qx <- if (is.null(scenario)) rep(NA_real_, length(quantile_probabilities)) else
    scenario_abs_quantiles(scenario, "x", quantile_probabilities)
  qy <- if (is.null(scenario)) rep(NA_real_, length(quantile_probabilities)) else
    scenario_abs_quantiles(scenario, "y", quantile_probabilities)
  counts_x <- ifelse(is.finite(qx), vapply(qx, function(q) sum(abs(x) > q), integer(1L)), NA)
  counts_y <- ifelse(is.finite(qy), vapply(qy, function(q) sum(abs(y) > q), integer(1L)), NA)
  list(
    x_ess_ratio_p1 = x_p1[["ess_ratio"]], y_ess_ratio_p1 = y_p1[["ess_ratio"]],
    x_max_share_p1 = x_p1[["maximum_share"]], y_max_share_p1 = y_p1[["maximum_share"]],
    x_top1pct_share_p1 = x_p1[["top_share"]], y_top1pct_share_p1 = y_p1[["top_share"]],
    x_ess_ratio_upper = x_upper[["ess_ratio"]],
    y_ess_ratio_upper = y_upper[["ess_ratio"]],
    x_max_share_upper = x_upper[["maximum_share"]],
    y_max_share_upper = y_upper[["maximum_share"]],
    x_top1pct_share_upper = x_upper[["top_share"]],
    y_top1pct_share_upper = y_upper[["top_share"]],
    x_count_above_q95 = counts_x[1L], y_count_above_q95 = counts_y[1L],
    x_count_above_q99 = counts_x[min(2L, length(counts_x))],
    y_count_above_q99 = counts_y[min(2L, length(counts_y))],
    x_latent_high_count = metadata_number(metadata_x, c("latent_high_count", "high_count")),
    y_latent_high_count = metadata_number(metadata_y, c("latent_high_count", "high_count"))
  )
}

value_by_moment_type <- function(value, moment_type, default = NULL) {
  if (is.null(value)) return(default)
  if (is.function(value)) return(value(moment_type))
  if (is.list(value)) {
    if (!is.null(value[[moment_type]])) return(value[[moment_type]])
    if (length(value) == 1L) return(value[[1L]])
    return(default)
  }
  if (!is.null(names(value)) && moment_type %in% names(value)) return(value[[moment_type]])
  value
}

scenario_roots_for_type <- function(scenario, moment_type) {
  root_table <- value_by_moment_type(scenario$root_truth, moment_type, NULL)
  if (is.data.frame(root_table) && "order" %in% names(root_table)) {
    roots <- as.numeric(root_table$order)
    return(roots[is.finite(roots)])
  }
  candidates <- c("positive_roots", "roots", "critical_orders")
  for (candidate in candidates) {
    if (!is.null(scenario[[candidate]])) {
      roots <- value_by_moment_type(scenario[[candidate]], moment_type, numeric(0L))
      roots <- as.numeric(roots)
      return(sort(unique(roots[is.finite(roots)])))
    }
  }
  numeric(0L)
}

scenario_identically_zero_for_type <- function(scenario, moment_type) {
  value <- scenario$identically_zero_by_type %||% scenario$identically_zero
  isTRUE(value_by_moment_type(value, moment_type, FALSE))
}

scenario_budget_for_type <- function(scenario, moment_type) {
  value <- value_by_moment_type(
    scenario$structural_budgets %||% scenario$structural_budget,
    moment_type, NA_integer_
  )
  if (length(value) != 1L || is.na(value)) NA_integer_ else as.integer(value)
}

scenario_directions_for_type <- function(scenario, moment_type, number_of_roots) {
  root_table <- value_by_moment_type(scenario$root_truth, moment_type, NULL)
  value <- if (is.data.frame(root_table) && "direction" %in% names(root_table)) {
    root_table$direction
  } else {
    value_by_moment_type(
      scenario$root_directions %||% scenario$root_direction,
      moment_type, NA_character_
    )
  }
  value <- as.character(value)
  if (!length(value)) value <- NA_character_
  if (length(value) == 1L && number_of_roots > 1L) value <- rep(value, number_of_roots)
  length(value) <- number_of_roots
  value
}

scenario_multiplicities_for_type <- function(scenario, moment_type, number_of_roots) {
  root_table <- value_by_moment_type(scenario$root_truth, moment_type, NULL)
  value <- if (is.data.frame(root_table) && "multiplicity" %in% names(root_table)) {
    root_table$multiplicity
  } else {
    value_by_moment_type(
      scenario$root_multiplicities %||% scenario$multiplicities,
      moment_type, rep(1L, number_of_roots)
    )
  }
  value <- as.integer(value)
  if (!length(value)) value <- integer(0L)
  if (length(value) == 1L && number_of_roots > 1L) value <- rep(value, number_of_roots)
  length(value) <- number_of_roots
  value[is.na(value) | value < 1L] <- 1L
  value
}

scenario_truth_for_interval <- function(scenario, moment_type, interval) {
  if (!is.null(scenario$root_truth) &&
      exists("truth_on_interval", mode = "function", inherits = TRUE)) {
    truth_function <- get("truth_on_interval", mode = "function", inherits = TRUE)
    truth_arguments <- list(
      scenario = scenario, p_min = interval[1L], p_max = interval[2L]
    )
    if ("moment_type" %in% names(formals(truth_function))) {
      truth_arguments$moment_type <- moment_type
    }
    raw <- do.call(truth_function, truth_arguments)
    roots <- as.numeric(raw$roots %||% numeric(0L))
    multiplicities <- as.integer(
      raw$root_multiplicities %||% rep.int(1L, length(roots))
    )
    directions <- as.character(
      raw$root_directions %||% rep.int(NA_character_, length(roots))
    )
    identically_zero <- isTRUE(raw$identically_zero)
    interior <- roots > interval[1L] + 1e-12 & roots < interval[2L] - 1e-12
    return(list(
      roots = roots, directions = directions, multiplicities = multiplicities,
      distinct_root_count = if (identically_zero) Inf else
        as.numeric(raw$root_count %||% length(roots)),
      root_count = if (identically_zero) Inf else
        as.numeric(raw$total_multiplicity %||% sum(multiplicities)),
      has_any_root = isTRUE(raw$has_any_root),
      has_reversal = isTRUE(raw$has_reversal),
      unique_root = if (!identically_zero && length(roots) == 1L &&
                          isTRUE(interior[1L]) && multiplicities[1L] == 1L) {
        roots[1L]
      } else NA_real_,
      identically_zero = identically_zero
    ))
  }
  all_roots <- scenario_roots_for_type(scenario, moment_type)
  in_interval <- all_roots >= interval[1L] - 1e-12 & all_roots <= interval[2L] + 1e-12
  roots <- all_roots[in_interval]
  interior <- roots > interval[1L] + 1e-12 & roots < interval[2L] - 1e-12
  all_directions <- scenario_directions_for_type(scenario, moment_type, length(all_roots))
  all_multiplicities <- scenario_multiplicities_for_type(
    scenario, moment_type, length(all_roots)
  )
  directions <- all_directions[in_interval]
  multiplicities <- all_multiplicities[in_interval]
  identically_zero <- scenario_identically_zero_for_type(scenario, moment_type)
  list(
    roots = roots, directions = directions, multiplicities = multiplicities,
    distinct_root_count = if (identically_zero) Inf else length(roots),
    root_count = if (identically_zero) Inf else sum(multiplicities),
    has_any_root = identically_zero || length(roots) > 0L,
    has_reversal = any(interior & !is.na(directions)),
    unique_root = if (!identically_zero && length(roots) == 1L && interior &&
                        multiplicities[1L] == 1L) roots else NA_real_,
    identically_zero = identically_zero
  )
}

normalise_root_intervals <- function(root_intervals = NULL, root_interval = NULL) {
  if (is.null(root_intervals)) root_intervals <- root_interval
  if (is.null(root_intervals)) return(list())
  if (is.numeric(root_intervals) && length(root_intervals) == 2L) {
    root_intervals <- list(Q0 = as.numeric(root_intervals))
  }
  assert_true(is.list(root_intervals), "root_intervals must be a list of endpoint pairs.")
  # An empty list is the deliberate representation of "no prespecified
  # root-isolation interval" (for example, the absolute component of the
  # mixed-sign scenario).  In R, paste0("Q", integer(0)) returns "Q", so the
  # generic naming branch below cannot be applied to a zero-length list.
  if (!length(root_intervals)) return(list())
  if (is.null(names(root_intervals)) || any(!nzchar(names(root_intervals)))) {
    names(root_intervals) <- paste0("Q", seq_along(root_intervals))
  }
  root_intervals <- lapply(root_intervals, function(interval) {
    interval <- as.numeric(interval)
    assert_true(length(interval) == 2L && all(is.finite(interval)) &&
                  interval[1L] < interval[2L],
                "Each root-isolation interval must have increasing finite endpoints.")
    interval
  })
  root_intervals
}

normalise_root_intervals_by_type <- function(
    root_intervals = NULL, root_interval = NULL, scenario = NULL,
    moment_types) {
  source <- root_intervals
  if (is.null(source) && !is.null(scenario$root_intervals)) {
    source <- scenario$root_intervals
  }
  if (is.null(source)) source <- root_interval

  by_type <- setNames(vector("list", length(moment_types)), moment_types)
  is_type_specific <- is.list(source) && !is.null(names(source)) &&
    any(moment_types %in% names(source))
  for (moment_type in moment_types) {
    value <- if (is_type_specific) source[[moment_type]] else source
    by_type[[moment_type]] <- normalise_root_intervals(value, NULL)
  }
  by_type
}

bracket_truth_flags <- function(brackets, truth) {
  if (!nrow(brackets)) {
    return(list(covers = logical(0L), direction_correct = logical(0L)))
  }
  covers <- direction_correct <- logical(nrow(brackets))
  for (index in seq_len(nrow(brackets))) {
    inside <- which(truth$roots >= brackets$lower[index] - 1e-12 &
                      truth$roots <= brackets$upper[index] + 1e-12)
    covers[index] <- length(inside) > 0L
    direction_correct[index] <- if (!length(inside)) FALSE else
      any(!is.na(truth$directions[inside]) &
            truth$directions[inside] == brackets$direction[index])
  }
  list(covers = covers, direction_correct = direction_correct)
}

false_cell_sign_statement <- function(cells, scenario, moment_type, truth,
                                      numerical_tolerance = 0) {
  signed_cells <- cells[cells$positive | cells$negative, , drop = FALSE]
  if (!nrow(signed_cells)) return(FALSE)
  if (truth$identically_zero) return(TRUE)
  for (index in seq_len(nrow(signed_cells))) {
    cell <- signed_cells[index, ]
    if (any(truth$roots >= cell$left - 1e-12 & truth$roots <= cell$right + 1e-12)) {
      return(TRUE)
    }
    probes <- c(cell$left, (cell$left + cell$right) / 2, cell$right)
    values <- scenario$delta(probes, moment_type)
    if (cell$positive && any(values <= numerical_tolerance)) return(TRUE)
    if (cell$negative && any(values >= -numerical_tolerance)) return(TRUE)
  }
  FALSE
}

collapse_numeric <- function(values) {
  if (!length(values)) "" else paste(format(values, digits = 16L, scientific = FALSE),
                                      collapse = ";")
}

collapse_character <- function(values) {
  if (!length(values)) "" else paste(values, collapse = ";")
}

# Keep the optional augmented construction in a separate module.  It is
# loaded here, after all primary numerical helpers have been defined, and
# does not alter the primary band or any primary decision.
source_cmo("R/augmented.R", local = environment())

analyse_one_replication <- function(
    scenario, n_x, n_y, orders, moment_types,
    bootstrap_reps, alpha, multiplier_distribution,
    bootstrap_batch_size, variance_tolerance, numerical_tolerance,
    audit_grid_size, root_interval = NULL, root_scan_size = 2001L,
    bootstrap_seed = NULL, order_intervals = NULL, root_intervals = NULL,
    guard_max_levels = 6L, enclosure_max_levels = 6L,
    enclosure_max_nodes = 28801L, enclosure_refine_unresolved_levels = 0L,
    enclosure_safety_margin = 1e-12,
    root_max_levels = 24L, root_max_evaluations = 50000L,
    root_tolerance = 1e-10, slope_tolerance = 1e-10,
    root_variance_tolerance = 1e-12,
    evaluation_chunk_size = 256L,
    maximum_enclosure_bisections = NULL,
    maximum_continuum_enclosure_bisections = NULL,
    maximum_enclosure_nodes = NULL,
    maximum_enclosure_subintervals = NULL,
    maximum_enclosure_half_width = NULL,
    enclosure_absolute_tolerance = NULL,
    maximum_root_bisections = NULL,
    maximum_root_evaluations = NULL,
    root_absolute_tolerance = NULL,
    tangency_neighborhood = NULL,
    audit_grid_spacing = NULL,
    maximum_truth_bisection_levels = 20L,
    maximum_truth_subintervals = NULL,
    contrast_variance_tolerance = NULL,
    augmented_variance_tolerance = NULL,
    tail_diagnostic_probabilities = c(0.95, 0.99),
    tail_top_fraction = 0.01,
    roundoff_inflation_required = FALSE,
    run_root_multiplier_refinement = TRUE,
    run_wald_refinement = TRUE,
    run_pointwise_ablation = TRUE,
    run_grid_only_ablation = TRUE,
    run_derivative_assisted_ablation = FALSE,
    enclosure_orders = NULL) {
  replication_started <- proc.time()[["elapsed"]]
  invisible(gc(reset = TRUE))
  if (!is.null(maximum_enclosure_bisections)) {
    enclosure_max_levels <- as.integer(maximum_enclosure_bisections)
    guard_max_levels <- as.integer(maximum_enclosure_bisections)
  }
  # The legacy control above intentionally keeps its historical joint effect
  # on the relative-variance guard and the continuum enclosure.  Paired
  # depth diagnostics need to vary only the latter, so this more specific
  # control is applied afterwards and takes precedence for the enclosure.
  if (!is.null(maximum_continuum_enclosure_bisections)) {
    enclosure_max_levels <- as.integer(
      maximum_continuum_enclosure_bisections
    )
  }
  if (!is.null(maximum_enclosure_nodes)) {
    enclosure_max_nodes <- as.integer(maximum_enclosure_nodes)
    if (!is.null(maximum_enclosure_subintervals)) {
      enclosure_max_nodes <- min(
        enclosure_max_nodes,
        as.integer(maximum_enclosure_subintervals) + 1L
      )
    }
  } else if (!is.null(maximum_enclosure_subintervals)) {
    enclosure_max_nodes <- as.integer(maximum_enclosure_subintervals) + 1L
  }
  if (!is.null(maximum_root_bisections)) {
    root_max_levels <- as.integer(maximum_root_bisections)
  }
  if (!is.null(maximum_root_evaluations)) {
    root_max_evaluations <- as.integer(maximum_root_evaluations)
  }
  if (!is.null(root_absolute_tolerance)) {
    root_tolerance <- as.numeric(root_absolute_tolerance)
  }
  if (!is.null(contrast_variance_tolerance)) {
    variance_tolerance <- as.numeric(contrast_variance_tolerance)
  }
  if (!is.null(enclosure_absolute_tolerance)) {
    enclosure_safety_margin <- max(
      as.numeric(enclosure_safety_margin),
      as.numeric(enclosure_absolute_tolerance)
    )
  }
  assert_true(length(tail_diagnostic_probabilities) >= 2L &&
                all(is.finite(tail_diagnostic_probabilities)) &&
                all(tail_diagnostic_probabilities > 0 &
                      tail_diagnostic_probabilities < 1),
              "tail_diagnostic_probabilities must contain at least two probabilities in (0,1).")
  assert_true(length(tail_top_fraction) == 1L && is.finite(tail_top_fraction) &&
                tail_top_fraction > 0 && tail_top_fraction <= 1,
              "tail_top_fraction must lie in (0,1].")
  sample_x <- scenario$draw_x %||% scenario$r_x
  sample_y <- scenario$draw_y %||% scenario$r_y
  assert_true(is.function(sample_x) && is.function(sample_y),
              "The scenario must provide draw_x/draw_y or r_x/r_y generators.")
  generated_x <- extract_sample_values(sample_x(n_x))
  generated_y <- extract_sample_values(sample_y(n_y))
  x <- generated_x$values
  y <- generated_y$values
  assert_true(length(x) == n_x && length(y) == n_y,
              "A data generator returned the wrong sample size.")
  bootstrap_seed_used <- if (is.null(bootstrap_seed)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    bootstrap_seed
  }
  assign(".Random.seed", bootstrap_seed_used, envir = .GlobalEnv)

  intervals <- normalise_order_intervals(order_intervals, orders)
  bootstrap_orders <- normalise_analysis_grid(
    orders, intervals, label = "bootstrap"
  )
  enclosure_master_orders <- if (is.null(enclosure_orders)) NULL else
    normalise_analysis_grid(
      enclosure_orders,
      intervals,
      label = "enclosure"
    )
  root_intervals_by_type <- normalise_root_intervals_by_type(
    root_intervals = root_intervals, root_interval = root_interval,
    scenario = scenario, moment_types = moment_types
  )
  nested <- nested_multiplier_bands(
    x = x, y = y, orders = bootstrap_orders, moment_types = moment_types,
    order_intervals = intervals, bootstrap_reps = bootstrap_reps,
    alpha = alpha, multiplier_distribution = multiplier_distribution,
    bootstrap_batch_size = bootstrap_batch_size,
    variance_tolerance = variance_tolerance,
    guard_max_levels = guard_max_levels,
    guard_max_nodes = enclosure_max_nodes
  )

  augmented_variance_tolerance <- as.numeric(
    augmented_variance_tolerance %||% variance_tolerance
  )
  assert_true(length(augmented_variance_tolerance) == 1L &&
                is.finite(augmented_variance_tolerance) &&
                augmented_variance_tolerance >= 0,
              "augmented_variance_tolerance must be a nonnegative finite scalar.")
  augmented <- if (isTRUE(run_derivative_assisted_ablation)) {
    nested_augmented_multiplier_bands(
      base_object = nested$master, order_intervals = intervals,
      bootstrap_reps = bootstrap_reps, alpha = alpha,
      multiplier_distribution = multiplier_distribution,
      bootstrap_batch_size = bootstrap_batch_size,
      variance_tolerance = augmented_variance_tolerance,
      guard_max_levels = guard_max_levels,
      guard_max_nodes = enclosure_max_nodes,
      bootstrap_seed = bootstrap_seed_used
    )
  } else NULL

  # Root isolation is independent of K and is therefore performed once for
  # each moment-type/isolation-interval pair and reused by all nested rows.
  root_analyses <- setNames(vector("list", length(moment_types)), moment_types)
  for (moment_type in moment_types) {
    type_root_intervals <- root_intervals_by_type[[moment_type]]
    root_analyses[[moment_type]] <- setNames(
      vector("list", length(type_root_intervals)), names(type_root_intervals)
    )
    for (root_id in names(type_root_intervals)) {
      root_range <- type_root_intervals[[root_id]]
      isolation <- isolate_empirical_roots(
        x = x, y = y, interval = root_range, moment_type = moment_type,
        initial_scan_size = root_scan_size, max_levels = root_max_levels,
        max_evaluations = root_max_evaluations,
        root_tolerance = root_tolerance,
        numerical_tolerance = numerical_tolerance,
        chunk_size = evaluation_chunk_size
      )
      root_truth <- scenario_truth_for_interval(scenario, moment_type, root_range)
      true_root <- if (!root_truth$identically_zero && length(root_truth$roots) == 1L) {
        root_truth$roots[1L]
      } else NA_real_
      wald <- wald_from_isolation(
        object = nested$master, isolation = isolation,
        moment_type = moment_type, true_root = true_root, alpha = alpha,
        slope_tolerance = slope_tolerance,
        variance_tolerance = root_variance_tolerance,
        state_cache = nested$state_cache
      )
      root_analyses[[moment_type]][[root_id]] <-
        list(interval = root_range, truth = root_truth,
             isolation = isolation, wald = wald)
    }
  }

  rows <- list()
  row_counter <- 0L
  root_multiplier_cache <- new.env(parent = emptyenv(), hash = TRUE)
  truth_max_nodes <- if (!is.null(maximum_truth_subintervals)) {
    as.integer(maximum_truth_subintervals) + 1L
  } else {
    enclosure_max_nodes
  }
  for (interval_id in names(intervals)) {
    fit <- nested$fits[[interval_id]]
    interval <- intervals[[interval_id]]
    enclosure_initial_grid <- if (is.null(enclosure_master_orders)) {
      grid_on_interval(fit$guard$grid %||% fit$orders, interval)
    } else {
      grid_on_interval(enclosure_master_orders, interval)
    }
    audit_orders <- regular_interval_grid(
      interval, spacing = audit_grid_spacing, size = audit_grid_size
    )
    for (moment_type in moment_types) {
      row_counter <- row_counter + 1L
      band <- extract_band(fit, moment_type)
      truth <- scenario_truth_for_interval(scenario, moment_type, interval)
      true_delta <- scenario$delta(band$orders, moment_type)
      derivative_assisted <- if (!is.null(augmented)) {
        augmented_for_enclosure <- augmented
        if (!is.null(enclosure_master_orders)) {
          augmented_for_enclosure$fits[[interval_id]]$guard$grid <-
            enclosure_initial_grid
        }
        summarise_derivative_assisted_row(
          augmented = augmented_for_enclosure, interval_id = interval_id,
          moment_type = moment_type, truth = truth,
          numerical_tolerance = numerical_tolerance,
          safety_margin = enclosure_safety_margin,
          max_levels = enclosure_max_levels,
          max_nodes = enclosure_max_nodes,
          refine_unresolved_levels = enclosure_refine_unresolved_levels,
          maximum_cell_half_width = maximum_enclosure_half_width %||% Inf,
          roundoff_inflation = isTRUE(roundoff_inflation_required)
        )
      } else NULL
      audit_band <- evaluate_fitted_band(
        fit, audit_orders, moment_type, chunk_size = evaluation_chunk_size
      )
      audit_true_delta <- scenario$delta(audit_orders, moment_type)

      enclosure <- adaptive_continuum_enclosure(
        fit = fit, moment_type = moment_type,
        initial_grid = enclosure_initial_grid,
        numerical_tolerance = numerical_tolerance,
        safety_margin = enclosure_safety_margin,
        max_levels = enclosure_max_levels, max_nodes = enclosure_max_nodes,
        refine_unresolved_levels = enclosure_refine_unresolved_levels,
        maximum_cell_half_width = maximum_enclosure_half_width %||% Inf,
        roundoff_inflation = isTRUE(roundoff_inflation_required),
        state_cache = nested$state_cache,
        constant_cache = nested$constant_cache
      )
      continuum_truth <- certify_continuum_truth_coverage(
        fit = fit, scenario = scenario, moment_type = moment_type,
        initial_grid = enclosure_initial_grid,
        numerical_tolerance = numerical_tolerance,
        safety_margin = enclosure_safety_margin,
        max_levels = maximum_truth_bisection_levels,
        max_nodes = truth_max_nodes,
        roundoff_inflation = isTRUE(roundoff_inflation_required),
        state_cache = nested$state_cache,
        constant_cache = nested$constant_cache
      )
      cells <- enclosure$cells
      outer <- outer_set_summary(cells, true_roots = truth$roots)
      tangency_interval <- tangency_neighborhood %||% c(0.50, 1.45)
      tangency_applicable <- length(truth$roots) == 1L &&
        length(truth$multiplicities) == 1L && truth$multiplicities[1L] == 2L &&
        tangency_interval[1L] >= interval[1L] - 1e-12 &&
        tangency_interval[2L] <= interval[2L] + 1e-12
      tangency <- tangency_localization_summary(
        classification = cells, neighborhood = tangency_interval,
        root = if (length(truth$roots)) truth$roots[1L] else NA_real_,
        n_eff = fit$n_eff, applicable = tangency_applicable
      )
      tangency_global <- tangency_global_geometry_summary(
        outer = outer, interval = interval,
        root = if (length(truth$roots)) truth$roots[1L] else NA_real_,
        n_eff = fit$n_eff, applicable = tangency_applicable,
        continuum_covered = continuum_truth$covered
      )
      brackets <- find_certified_brackets(
        band$orders, band$lower, band$upper, numerical_tolerance
      )
      bracket_truth <- bracket_truth_flags(brackets, truth)
      bracket_root_matches <- if (nrow(brackets) && length(truth$roots)) {
        lapply(seq_len(nrow(brackets)), function(index) {
          which(truth$roots >= brackets$lower[index] - 1e-12 &
                  truth$roots <= brackets$upper[index] + 1e-12)
        })
      } else list()
      complete_distinct_bracket_coverage <- if (length(truth$roots) >= 1L &&
                                                  nrow(brackets) >= 1L) {
        matched <- unique(unlist(bracket_root_matches, use.names = FALSE))
        length(matched) == length(truth$roots) && all(bracket_truth$covers)
      } else NA
      alternation_count <- nrow(brackets)
      budget <- scenario_budget_for_type(scenario, moment_type)
      budget_decision <- mellin_budget_decision(
        alternation_count, budget, truth$distinct_root_count,
        truth$root_count
      )
      continuum_no_root <- all(!cells$outer)
      node_status_false <- any(
        (band$lower > numerical_tolerance & true_delta <= 0) |
          (band$upper < -numerical_tolerance & true_delta >= 0)
      )
      cell_status_false <- false_cell_sign_statement(
        cells, scenario, moment_type, truth, numerical_tolerance
      )
      false_sign <- node_status_false || cell_status_false
      reversal <- nrow(brackets) > 0L
      false_reversal <- reversal && (
        truth$identically_zero ||
          !all(bracket_truth$covers & bracket_truth$direction_correct)
      )
      false_direction <- reversal && !all(bracket_truth$direction_correct)
      false_bracket <- reversal && !all(bracket_truth$covers)
      false_no_root <- continuum_no_root && truth$has_any_root
      exact_count <- budget_decision$exact_count_certified
      false_exact <- exact_count && !budget_decision$exact_count_statement_true
      outer_truth_coverage <- if (truth$identically_zero) {
        all(cells$outer)
      } else if (length(truth$roots)) {
        isTRUE(outer$contains_all_roots)
      } else TRUE
      false_outer <- !outer_truth_coverage
      ideal_outer_contains_root <- if (length(truth$roots) == 1L) {
        at_root <- band_at_order(fit, truth$roots[1L], moment_type,
                                 nested$state_cache)
        at_root[["lower"]] <= 0 && 0 <= at_root[["upper"]]
      } else NA

      # Grid-only and deliberately unadjusted pointwise ablations.  Their
      # switches suppress only the ablation outputs; the primary procedure
      # may still use directly evaluated grid nodes as valid sign anchors.
      if (isTRUE(run_grid_only_ablation)) {
        grid_excludes_zero <- band$lower > numerical_tolerance |
          band$upper < -numerical_tolerance
        grid_no_root <- all(grid_excludes_zero)
        grid_only_reversal <- nrow(find_certified_brackets(
          band$orders, band$lower, band$upper, numerical_tolerance
        )) > 0L
      } else {
        grid_no_root <- NA
        grid_only_reversal <- NA
      }
      if (isTRUE(run_pointwise_ablation)) {
        pointwise <- pointwise_band(
          fit, band$orders, moment_type, alpha,
          chunk_size = evaluation_chunk_size
        )
        pointwise_true <- scenario$delta(pointwise$orders, moment_type)
        pointwise_false_sign <- any(
          (pointwise$lower > numerical_tolerance & pointwise_true <= 0) |
            (pointwise$upper < -numerical_tolerance & pointwise_true >= 0)
        )
        pointwise_brackets <- find_certified_brackets(
          pointwise$orders, pointwise$lower, pointwise$upper,
          numerical_tolerance
        )
        pointwise_no_root <- all(
          pointwise$lower > numerical_tolerance |
            pointwise$upper < -numerical_tolerance
        )
        pointwise_root_coverage <- if (
          !truth$identically_zero && length(truth$roots) > 0L
        ) {
          pointwise_at_roots <- pointwise_band(
            fit, truth$roots, moment_type, alpha,
            chunk_size = evaluation_chunk_size
          )
          all(
            pointwise_at_roots$lower <= 0 &
              pointwise_at_roots$upper >= 0
          )
        } else NA
        pointwise_grid_curve_coverage <- all(
          pointwise$lower <= pointwise_true &
            pointwise_true <= pointwise$upper
        )
        pointwise_reversal <- nrow(pointwise_brackets) > 0L
      } else {
        pointwise_false_sign <- NA
        pointwise_no_root <- NA
        pointwise_root_coverage <- NA
        pointwise_grid_curve_coverage <- NA
        pointwise_reversal <- NA
      }

      # Reuse root results only when the prespecified Q is wholly contained in K.
      type_root_intervals <- root_intervals_by_type[[moment_type]]
      eligible_root_ids <- names(type_root_intervals)[vapply(type_root_intervals, function(value) {
        value[1L] >= interval[1L] - 1e-12 && value[2L] <= interval[2L] + 1e-12
      }, logical(1L))]
      root_details <- root_analyses[[moment_type]][eligible_root_ids]
      truth_bearing <- root_details[vapply(root_details, function(detail) {
        !detail$truth$identically_zero && length(detail$truth$roots) == 1L &&
          detail$truth$multiplicities[1L] == 1L &&
          !is.na(detail$truth$directions[1L])
      }, logical(1L))]
      multiplier_key <- paste(
        moment_type, paste(eligible_root_ids, collapse = ","), sep = "::"
      )
      if (isTRUE(run_root_multiplier_refinement)) {
        if (!exists(multiplier_key, envir = root_multiplier_cache, inherits = FALSE)) {
          assign(
            multiplier_key,
            joint_root_multiplier_refinement(
              object = nested$master, root_details = root_details,
              moment_type = moment_type, bootstrap_reps = bootstrap_reps,
              alpha = alpha, multiplier_distribution = multiplier_distribution,
              bootstrap_batch_size = bootstrap_batch_size,
              bootstrap_seed = bootstrap_seed_used,
              slope_tolerance = slope_tolerance,
              variance_tolerance = root_variance_tolerance
            ),
            envir = root_multiplier_cache
          )
        }
        root_multiplier <- get(
          multiplier_key, envir = root_multiplier_cache, inherits = FALSE
        )
      } else {
        root_multiplier <- list(
          applicable = FALSE, reported = NA,
          critical_value = NA_real_, root_ids = eligible_root_ids,
          root_hats = numeric(0L), lower = numeric(0L), upper = numeric(0L),
          conditional_cover = NA, formal_unconditional_cover = NA,
          report_and_cover = NA
        )
      }
      wald_applicable <- isTRUE(run_wald_refinement) && length(truth_bearing) > 0L
      wald_reported <- if (wald_applicable) {
        all(vapply(truth_bearing, function(detail) detail$wald$reported, logical(1L)))
      } else NA
      wald_formal_cover <- if (wald_applicable) {
        all(vapply(truth_bearing, function(detail) {
          isTRUE(detail$wald$formal_unconditional_cover)
        }, logical(1L)))
      } else NA
      wald_report_and_cover <- if (wald_applicable) {
        all(vapply(truth_bearing, function(detail) {
          isTRUE(detail$wald$report_and_cover)
        }, logical(1L)))
      } else NA
      wald_conditional_cover <- if (isTRUE(wald_reported)) {
        all(vapply(truth_bearing, function(detail) {
          isTRUE(detail$wald$conditional_cover)
        }, logical(1L)))
      } else NA
      isolated_roots <- unlist(lapply(root_details, function(detail) detail$isolation$roots))
      isolation_bracket_ids <- unlist(Map(function(detail, detail_id) {
        rep.int(detail_id, nrow(detail$isolation$brackets))
      }, root_details, names(root_details)), use.names = FALSE)
      isolation_bracket_lower <- unlist(lapply(root_details, function(detail) {
        detail$isolation$brackets$lower
      }), use.names = FALSE)
      isolation_bracket_upper <- unlist(lapply(root_details, function(detail) {
        detail$isolation$brackets$upper
      }), use.names = FALSE)
      isolation_bracket_direction <- unlist(lapply(root_details, function(detail) {
        detail$isolation$brackets$direction
      }), use.names = FALSE)
      marginal_wald_ids <- names(truth_bearing)
      marginal_wald_hats <- vapply(truth_bearing, function(detail) {
        detail$wald$root_hat
      }, numeric(1L))
      marginal_wald_lower <- vapply(truth_bearing, function(detail) {
        detail$wald$lower
      }, numeric(1L))
      marginal_wald_upper <- vapply(truth_bearing, function(detail) {
        detail$wald$upper
      }, numeric(1L))
      marginal_true_roots <- vapply(truth_bearing, function(detail) {
        detail$truth$roots[1L]
      }, numeric(1L))
      marginal_root_errors <- marginal_wald_hats - marginal_true_roots
      isolation_success <- if (length(root_details)) {
        all(vapply(root_details, function(detail) detail$isolation$success, logical(1L)))
      } else NA
      # Root-estimation diagnostics are defined by the prespecified regular
      # root collection and its numerical isolation, independently of whether
      # a Wald interval was requested.  Requiring every eligible Q to be a
      # truth-bearing simple-root interval prevents silently reporting only a
      # convenient subset of the prespecified collection.
      root_collection_applicable <- length(root_details) > 0L &&
        length(truth_bearing) == length(root_details) &&
        is.finite(truth$distinct_root_count) &&
        length(truth_bearing) == as.integer(truth$distinct_root_count)
      root_collection_joint_isolation_success <- if (
          root_collection_applicable) {
        all(vapply(truth_bearing, function(detail) {
          isTRUE(detail$isolation$success) &&
            length(detail$isolation$roots) == 1L &&
            is.finite(detail$isolation$roots[1L])
        }, logical(1L)))
      } else NA
      root_collection_true_roots <- if (root_collection_applicable) {
        vapply(truth_bearing, function(detail) {
          detail$truth$roots[1L]
        }, numeric(1L))
      } else numeric(0L)
      root_collection_hats <- if (
          isTRUE(root_collection_joint_isolation_success)) {
        vapply(truth_bearing, function(detail) {
          detail$isolation$roots[1L]
        }, numeric(1L))
      } else numeric(0L)
      root_collection_errors <- if (length(root_collection_hats)) {
        root_collection_hats - root_collection_true_roots
      } else numeric(0L)
      single_detail <- if (wald_applicable && length(truth_bearing) == 1L) {
        truth_bearing[[1L]]
      } else NULL
      root_hat <- if (!is.null(single_detail)) single_detail$wald$root_hat else NA_real_
      true_root <- if (!is.null(single_detail)) single_detail$truth$roots[1L] else NA_real_

      tail <- tail_diagnostics(
        x, y, upper_order = interval[2L], scenario = scenario,
        quantile_probabilities = tail_diagnostic_probabilities,
        top_fraction = tail_top_fraction,
        metadata_x = generated_x$metadata, metadata_y = generated_y$metadata
      )
      if (!is.null(derivative_assisted)) {
        augmented_fit <- derivative_assisted$fit
        augmented_certificates <- derivative_assisted$certificates
        augmented_certificate_truth <- derivative_assisted$certificate_truth
        augmented_reported <- !isTRUE(augmented_fit$guard_failed)
        augmented_certified <- nrow(augmented_certificates) > 0L
        augmented_false_report <- augmented_certified &&
          !all(augmented_certificate_truth)
        augmented_true_value <- scenario$delta(
          derivative_assisted$value_band$orders, moment_type
        )
        augmented_true_derivative <- scenario$delta_derivative(
          derivative_assisted$derivative_band$orders, moment_type
        )
        augmented_value_grid_coverage <- all(
          derivative_assisted$value_band$lower <= augmented_true_value &
            augmented_true_value <= derivative_assisted$value_band$upper
        )
        augmented_derivative_grid_coverage <- all(
          derivative_assisted$derivative_band$lower <= augmented_true_derivative &
            augmented_true_derivative <= derivative_assisted$derivative_band$upper
        )
        augmented_joint_grid_coverage <- augmented_value_grid_coverage &&
          augmented_derivative_grid_coverage
        augmented_value_maximum_half_width <- max(
          derivative_assisted$value_band$half_width
        )
        augmented_derivative_maximum_half_width <- max(
          derivative_assisted$derivative_band$half_width
        )
        augmented_width_inflation <- if (
          is.finite(augmented_value_maximum_half_width) &&
            is.finite(max(band$half_width)) && max(band$half_width) > 0
        ) augmented_value_maximum_half_width / max(band$half_width) else NA_real_
        augmented_critical_value_ratio <- if (
          is.finite(augmented_fit$critical_value) &&
            is.finite(fit$critical_value) && fit$critical_value > 0
        ) augmented_fit$critical_value / fit$critical_value else NA_real_
        augmented_enclosure_limit <-
          isTRUE(derivative_assisted$value_enclosure$limit_hit) ||
          isTRUE(derivative_assisted$derivative_enclosure$limit_hit)
        augmented_failure <- isTRUE(augmented_fit$guard_failed) ||
          augmented_enclosure_limit
        augmented_failure_reason <- if (isTRUE(augmented_fit$guard_failed)) {
          "relative_variance_guard_failed"
        } else if (augmented_enclosure_limit) {
          "local_enclosure_limit_hit"
        } else ""
        augmented_derivative_positive_K <- augmented_reported && all(
          derivative_assisted$derivative_enclosure$cells$positive
        )
        augmented_derivative_negative_K <- augmented_reported && all(
          derivative_assisted$derivative_enclosure$cells$negative
        )
        augmented_value_left <- augmented_band_at_order(
          augmented_fit, interval[1L], moment_type, 0L,
          augmented$state_cache
        )
        augmented_value_right <- augmented_band_at_order(
          augmented_fit, interval[2L], moment_type, 0L,
          augmented$state_cache
        )
        augmented_exact_root_on_K <-
          (augmented_derivative_positive_K &&
             augmented_value_left[["upper"]] < -numerical_tolerance &&
             augmented_value_right[["lower"]] > numerical_tolerance) ||
          (augmented_derivative_negative_K &&
             augmented_value_left[["lower"]] > numerical_tolerance &&
             augmented_value_right[["upper"]] < -numerical_tolerance)
        augmented_exact_root_on_K_direction <- if (
          augmented_derivative_positive_K &&
            augmented_value_left[["upper"]] < -numerical_tolerance &&
            augmented_value_right[["lower"]] > numerical_tolerance
        ) {
          "up"
        } else if (
          augmented_derivative_negative_K &&
            augmented_value_left[["lower"]] > numerical_tolerance &&
            augmented_value_right[["upper"]] < -numerical_tolerance
        ) {
          "down"
        } else ""
        augmented_exact_root_on_K_statement_true <- if (
          isTRUE(augmented_exact_root_on_K)
        ) {
          !truth$identically_zero && length(truth$roots) == 1L &&
            length(truth$multiplicities) == 1L &&
            truth$multiplicities[1L] == 1L &&
            length(truth$directions) == 1L &&
            identical(
              truth$directions[1L], augmented_exact_root_on_K_direction
            )
        } else NA
        augmented_false_exact_root_on_K <- isTRUE(augmented_exact_root_on_K) &&
          !isTRUE(augmented_exact_root_on_K_statement_true)
      } else {
        augmented_fit <- NULL
        augmented_certificates <- data.frame()
        augmented_certificate_truth <- logical(0L)
        augmented_reported <- NA
        augmented_certified <- NA
        augmented_false_report <- NA
        augmented_value_grid_coverage <- NA
        augmented_derivative_grid_coverage <- NA
        augmented_joint_grid_coverage <- NA
        augmented_value_maximum_half_width <- NA_real_
        augmented_derivative_maximum_half_width <- NA_real_
        augmented_width_inflation <- NA_real_
        augmented_critical_value_ratio <- NA_real_
        augmented_failure <- NA
        augmented_failure_reason <- ""
        augmented_derivative_positive_K <- NA
        augmented_derivative_negative_K <- NA
        augmented_exact_root_on_K <- NA
        augmented_exact_root_on_K_direction <- ""
        augmented_exact_root_on_K_statement_true <- NA
        augmented_false_exact_root_on_K <- NA
      }
      any_primary_false <- false_sign || false_reversal || false_direction ||
        false_bracket || false_no_root || false_exact ||
        isTRUE(budget_decision$rejected) || false_outer

      row <- data.frame(
        interval_id = interval_id, moment_type = moment_type,
        n_x = n_x, n_y = n_y, n_eff = fit$n_eff,
        p_min = interval[1L], p_max = interval[2L],
        grid_size = length(band$orders),
        bootstrap_grid_size = length(band$orders),
        bootstrap_grid_spacing = maximum_grid_spacing(band$orders),
        enclosure_grid_size = length(enclosure_initial_grid),
        enclosure_grid_spacing = maximum_grid_spacing(enclosure_initial_grid),
        configured_guard_maximum_levels = guard_max_levels,
        configured_continuum_enclosure_maximum_levels =
          enclosure_max_levels,
        bootstrap_reps = bootstrap_reps,
        alpha = alpha, guard_failed = fit$guard_failed,
        guard_minimum_ratio = fit$guard$minimum_certified_ratio,
        guard_nodes_used = fit$guard$nodes_used,
        guard_maximum_level = fit$guard$maximum_level,
        critical_value = if (fit$guard_failed) NA_real_ else fit$critical_value,
        maximum_half_width = if (fit$guard_failed) Inf else max(band$half_width),
        integrated_half_width = if (fit$guard_failed) Inf else
          sum(diff(band$orders) * (head(band$half_width, -1L) +
                                    tail(band$half_width, -1L)) / 2),
        band_covers_truth_on_grid = all(band$lower <= true_delta & true_delta <= band$upper),
        band_covers_truth_on_audit_grid = all(
          audit_band$lower <= audit_true_delta & audit_true_delta <= audit_band$upper
        ),
        audit_grid_size = length(audit_orders),
        continuum_band_covers_truth = continuum_truth$covered,
        continuum_truth_nodes_used = continuum_truth$nodes_used,
        continuum_truth_maximum_level = continuum_truth$maximum_level,
        continuum_truth_unresolved_cells = continuum_truth$unresolved_cells,
        continuum_truth_limit_hit = continuum_truth$limit_hit,
        continuum_truth_endpoint_violation = continuum_truth$endpoint_violation,
        enclosure_nodes_used = enclosure$nodes_used,
        enclosure_maximum_level = enclosure$maximum_level,
        enclosure_limit_hit = enclosure$limit_hit,
        enclosure_depth_limit_hit = enclosure$depth_limit_hit,
        enclosure_node_limit_hit = enclosure$node_limit_hit,
        enclosure_variance_limit_hit = enclosure$variance_limit_hit,
        enclosure_statistical_limit_hit = enclosure$statistical_limit_hit,
        enclosure_limit_hit_cells = enclosure$limit_hit_cells,
        enclosure_limit_hit_total_width = enclosure$limit_hit_total_width,
        enclosure_limit_hit_width_proportion =
          enclosure$limit_hit_width_proportion,
        enclosure_limit_hit_maximum_cell_width =
          enclosure$limit_hit_maximum_cell_width,
        enclosure_depth_limit_cells = enclosure$depth_limit_cells,
        enclosure_depth_limit_total_width =
          enclosure$depth_limit_total_width,
        enclosure_depth_limit_width_proportion =
          enclosure$depth_limit_width_proportion,
        enclosure_depth_limit_maximum_cell_width =
          enclosure$depth_limit_maximum_cell_width,
        enclosure_node_limit_cells = enclosure$node_limit_cells,
        enclosure_node_limit_total_width = enclosure$node_limit_total_width,
        enclosure_node_limit_width_proportion =
          enclosure$node_limit_width_proportion,
        enclosure_node_limit_maximum_cell_width =
          enclosure$node_limit_maximum_cell_width,
        enclosure_variance_limit_cells = enclosure$variance_limit_cells,
        enclosure_variance_limit_total_width =
          enclosure$variance_limit_total_width,
        enclosure_variance_limit_width_proportion =
          enclosure$variance_limit_width_proportion,
        enclosure_variance_limit_maximum_cell_width =
          enclosure$variance_limit_maximum_cell_width,
        enclosure_statistical_limit_cells =
          enclosure$statistical_limit_cells,
        enclosure_statistical_limit_total_width =
          enclosure$statistical_limit_total_width,
        enclosure_statistical_limit_width_proportion =
          enclosure$statistical_limit_width_proportion,
        enclosure_statistical_limit_maximum_cell_width =
          enclosure$statistical_limit_maximum_cell_width,
        enclosure_variance_unresolved_cells = enclosure$variance_unresolved_cells,
        enclosure_unresolved_cells = enclosure$statistically_unresolved_cells,
        roundoff_inflation_requested = isTRUE(roundoff_inflation_required),
        roundoff_inflation_margin_applied = isTRUE(roundoff_inflation_required),
        continuum_no_root_certified = continuum_no_root,
        grid_only_no_root = grid_no_root,
        grid_continuum_no_root_disagreement = if (
          is.na(grid_no_root)
        ) NA else grid_no_root != continuum_no_root,
        reversal_certified = reversal,
        grid_only_reversal_certified = grid_only_reversal,
        certified_direction = collapse_character(brackets$direction),
        direction_correct = if (reversal) all(bracket_truth$direction_correct) else NA,
        false_direction_report = false_direction,
        false_sign_report = false_sign,
        false_reversal_report = false_reversal,
        false_no_root_report = false_no_root,
        budget_rejected = budget_decision$rejected,
        false_budget_rejection = budget_decision$rejected,
        zero_budget_diagnostic_applicable = is.finite(truth$root_count) &&
          truth$has_reversal && truth$root_count >= 1,
        zero_budget_rejected = alternation_count > 0L,
        exact_count_certified = exact_count,
        simplicity_certified = budget_decision$simplicity_certified,
        exact_count_statement_true = budget_decision$exact_count_statement_true,
        false_exact_count_report = false_exact,
        certified_alternation_count = alternation_count,
        bracket_count = nrow(brackets),
        bracket_lower = if (nrow(brackets)) min(brackets$lower) else NA_real_,
        bracket_upper = if (nrow(brackets)) max(brackets$upper) else NA_real_,
        bracket_length = if (nrow(brackets)) sum(brackets$upper - brackets$lower) else NA_real_,
        bracket_covers_root = if (nrow(brackets)) all(bracket_truth$covers) else NA,
        brackets_cover_all_distinct_roots = complete_distinct_bracket_coverage,
        false_bracket_report = false_bracket,
        true_roots = collapse_numeric(truth$roots),
        true_root_count = truth$root_count,
        true_distinct_root_count = truth$distinct_root_count,
        true_root = true_root,
        outer_band_contains_root = ideal_outer_contains_root,
        outer_set_covers_truth = outer_truth_coverage,
        numerical_outer_component_count = outer$component_count,
        numerical_outer_component_lefts = collapse_numeric(
          outer$component_left
        ),
        numerical_outer_component_rights = collapse_numeric(
          outer$component_right
        ),
        numerical_outer_total_length = outer$total_length,
        numerical_outer_retained_proportion = outer$total_length /
          (interval[2L] - interval[1L]),
        numerical_outer_contains_root = if (length(truth$roots) == 1L) {
          outer$contains_root
        } else NA,
        numerical_outer_contains_all_roots = outer$contains_all_roots,
        numerical_outer_hausdorff = outer$hausdorff_distance,
        false_outer_report = false_outer,
        tangency_localization_applicable = tangency$applicable,
        tangency_local_outer_nonempty = tangency$local_outer_nonempty,
        tangency_local_outer_total_length = tangency$local_outer_total_length,
        # Compatibility fields above retain their historical names.  The
        # aliases below make explicit that this geometry is intersected with
        # a fixed local window and can therefore be boundary-censored.
        tangency_truncated_local_outer_nonempty =
          tangency$local_outer_nonempty,
        tangency_truncated_local_outer_total_length =
          tangency$local_outer_total_length,
        tangency_truncated_local_touches_left_boundary =
          tangency$touches_left_boundary,
        tangency_truncated_local_touches_right_boundary =
          tangency$touches_right_boundary,
        tangency_truncated_local_censored_by_window =
          tangency$censored_by_window,
        tangency_scaled_outer_radius = if (isTRUE(continuum_truth$covered)) {
          tangency$scaled_radius
        } else NA_real_,
        tangency_scaled_outer_hausdorff = if (isTRUE(continuum_truth$covered)) {
          tangency$scaled_hausdorff
        } else if (isTRUE(tangency$applicable) &&
                     !isTRUE(tangency$local_outer_nonempty)) {
          Inf
        } else NA_real_,
        tangency_truncated_local_scaled_outer_radius = if (
          isTRUE(continuum_truth$covered)
        ) tangency$scaled_radius else NA_real_,
        tangency_truncated_local_scaled_outer_hausdorff = if (
          isTRUE(continuum_truth$covered)
        ) {
          tangency$scaled_hausdorff
        } else if (isTRUE(tangency$applicable) &&
                     !isTRUE(tangency$local_outer_nonempty)) {
          Inf
        } else NA_real_,
        tangency_global_geometry_applicable = tangency_global$applicable,
        tangency_global_geometry_conditioned_on_continuum_band_coverage =
          tangency_global$coverage_condition_met,
        tangency_global_outer_nonempty_conditional_on_continuum_band_coverage =
          tangency_global$outer_nonempty,
        tangency_global_outer_contains_root_conditional_on_continuum_band_coverage =
          tangency_global$outer_contains_root,
        tangency_global_outer_component_count_conditional_on_continuum_band_coverage =
          tangency_global$component_count,
        tangency_global_outer_total_length_conditional_on_continuum_band_coverage =
          tangency_global$total_length,
        tangency_global_outer_hausdorff_conditional_on_continuum_band_coverage =
          tangency_global$hausdorff,
        tangency_global_outer_left_radius_conditional_on_continuum_band_coverage =
          tangency_global$left_radius,
        tangency_global_outer_right_radius_conditional_on_continuum_band_coverage =
          tangency_global$right_radius,
        tangency_global_scaled_outer_total_length_conditional_on_continuum_band_coverage =
          tangency_global$scaled_total_length,
        tangency_global_scaled_outer_hausdorff_conditional_on_continuum_band_coverage =
          tangency_global$scaled_hausdorff,
        tangency_global_scaled_outer_left_radius_conditional_on_continuum_band_coverage =
          tangency_global$scaled_left_radius,
        tangency_global_scaled_outer_right_radius_conditional_on_continuum_band_coverage =
          tangency_global$scaled_right_radius,
        tangency_global_outer_touches_K_left_conditional_on_continuum_band_coverage =
          tangency_global$touches_K_left,
        tangency_global_outer_touches_K_right_conditional_on_continuum_band_coverage =
          tangency_global$touches_K_right,
        tangency_global_outer_touches_K_boundary_conditional_on_continuum_band_coverage =
          tangency_global$touches_K_boundary,
        tangency_global_root_component_present_conditional_on_continuum_band_coverage =
          tangency_global$root_component_present,
        tangency_global_root_component_left_conditional_on_continuum_band_coverage =
          tangency_global$root_component_left,
        tangency_global_root_component_right_conditional_on_continuum_band_coverage =
          tangency_global$root_component_right,
        tangency_global_root_component_length_conditional_on_continuum_band_coverage =
          tangency_global$root_component_length,
        tangency_global_root_component_left_radius_conditional_on_continuum_band_coverage =
          tangency_global$root_component_left_radius,
        tangency_global_root_component_right_radius_conditional_on_continuum_band_coverage =
          tangency_global$root_component_right_radius,
        tangency_global_scaled_root_component_length_conditional_on_continuum_band_coverage =
          tangency_global$scaled_root_component_length,
        tangency_global_scaled_root_component_left_radius_conditional_on_continuum_band_coverage =
          tangency_global$scaled_root_component_left_radius,
        tangency_global_scaled_root_component_right_radius_conditional_on_continuum_band_coverage =
          tangency_global$scaled_root_component_right_radius,
        tangency_global_root_component_touches_K_left_conditional_on_continuum_band_coverage =
          tangency_global$root_component_touches_K_left,
        tangency_global_root_component_touches_K_right_conditional_on_continuum_band_coverage =
          tangency_global$root_component_touches_K_right,
        tangency_global_root_component_touches_K_boundary_conditional_on_continuum_band_coverage =
          tangency_global$root_component_touches_K_boundary,
        root_isolation_applicable = length(root_details) > 0L,
        root_isolation_success = isolation_success,
        root_isolation_limit_hit = if (length(root_details)) {
          any(vapply(root_details, function(detail) detail$isolation$limit_hit, logical(1L)))
        } else NA,
        empirical_root_count = length(isolated_roots),
        isolated_roots = collapse_numeric(isolated_roots),
        root_isolation_interval_ids = collapse_character(names(root_details)),
        root_isolation_bracket_ids = collapse_character(isolation_bracket_ids),
        root_isolation_bracket_lower = collapse_numeric(isolation_bracket_lower),
        root_isolation_bracket_upper = collapse_numeric(isolation_bracket_upper),
        root_isolation_bracket_direction = collapse_character(
          isolation_bracket_direction
        ),
        root_isolation_maximum_bisections = if (length(root_details)) {
          max(vapply(root_details, function(detail) {
            detail$isolation$maximum_bisections_used
          }, integer(1L)))
        } else NA_integer_,
        root_isolation_all_tolerances_met = if (length(root_details)) {
          all(vapply(root_details, function(detail) {
            isTRUE(detail$isolation$all_root_tolerances_met)
          }, logical(1L)))
        } else NA,
        root_hat = root_hat,
        root_error = if (is.finite(root_hat) && is.finite(true_root)) root_hat - true_root else NA_real_,
        root_error_conditional_on_isolation = if (
          is.finite(root_hat) && is.finite(true_root)
        ) root_hat - true_root else NA_real_,
        derivative_hat = if (!is.null(single_detail)) single_detail$wald$derivative_hat else NA_real_,
        variance_at_root_hat = if (!is.null(single_detail)) single_detail$wald$variance_hat else NA_real_,
        root_standard_error = if (!is.null(single_detail)) single_detail$wald$standard_error else NA_real_,
        wald_reported = wald_reported,
        wald_interval_ids = collapse_character(marginal_wald_ids),
        wald_root_hats = collapse_numeric(marginal_wald_hats),
        wald_marginal_true_roots = collapse_numeric(marginal_true_roots),
        wald_marginal_lowers = collapse_numeric(marginal_wald_lower),
        wald_marginal_uppers = collapse_numeric(marginal_wald_upper),
        wald_marginal_root_errors = collapse_numeric(marginal_root_errors),
        wald_collection_bias = if (wald_applicable &&
                                      all(is.finite(marginal_root_errors))) {
          mean(marginal_root_errors)
        } else NA_real_,
        wald_collection_bias_conditional_on_isolation = if (
          wald_applicable && isTRUE(isolation_success) &&
            length(marginal_root_errors) > 0L &&
            all(is.finite(marginal_root_errors))
        ) {
          mean(marginal_root_errors)
        } else NA_real_,
        wald_collection_rmse = if (wald_applicable &&
                                      all(is.finite(marginal_root_errors))) {
          sqrt(mean(marginal_root_errors^2))
        } else NA_real_,
        wald_collection_rmse_conditional_on_isolation = if (
          wald_applicable && isTRUE(isolation_success) &&
            length(marginal_root_errors) > 0L &&
            all(is.finite(marginal_root_errors))
        ) {
          sqrt(mean(marginal_root_errors^2))
        } else NA_real_,
        root_collection_applicable = root_collection_applicable,
        root_collection_required_root_count = if (
          root_collection_applicable
        ) length(truth_bearing) else 0L,
        root_collection_joint_isolation_success =
          root_collection_joint_isolation_success,
        root_collection_root_ids = if (root_collection_applicable) {
          collapse_character(names(truth_bearing))
        } else "",
        root_collection_true_roots = collapse_numeric(
          root_collection_true_roots
        ),
        root_collection_root_hats = collapse_numeric(root_collection_hats),
        root_collection_errors = collapse_numeric(root_collection_errors),
        root_collection_mean_error_conditional_on_joint_isolation = if (
          length(root_collection_errors) > 0L
        ) mean(root_collection_errors) else NA_real_,
        root_collection_mean_squared_error_conditional_on_joint_isolation = if (
          length(root_collection_errors) > 0L
        ) mean(root_collection_errors^2) else NA_real_,
        wald_collection_total_length = if (isTRUE(wald_reported)) {
          sum(marginal_wald_upper - marginal_wald_lower)
        } else NA_real_,
        wald_collection_is_joint = FALSE,
        wald_lower = if (!is.null(single_detail)) single_detail$wald$lower else NA_real_,
        wald_upper = if (!is.null(single_detail)) single_detail$wald$upper else NA_real_,
        wald_length = if (!is.null(single_detail)) single_detail$wald$length else NA_real_,
        wald_conditional_cover = wald_conditional_cover,
        wald_covers_root = wald_conditional_cover,
        wald_formal_unconditional_cover = wald_formal_cover,
        wald_report_and_cover = wald_report_and_cover,
        root_multiplier_applicable = root_multiplier$applicable,
        root_multiplier_reported = root_multiplier$reported,
        root_multiplier_critical_value = root_multiplier$critical_value,
        root_multiplier_interval_ids = collapse_character(root_multiplier$root_ids),
        root_multiplier_root_hats = collapse_numeric(root_multiplier$root_hats),
        root_multiplier_lowers = collapse_numeric(root_multiplier$lower),
        root_multiplier_uppers = collapse_numeric(root_multiplier$upper),
        root_multiplier_total_length = if (isTRUE(root_multiplier$reported)) {
          sum(root_multiplier$upper - root_multiplier$lower)
        } else NA_real_,
        root_multiplier_conditional_cover = root_multiplier$conditional_cover,
        root_multiplier_formal_unconditional_cover =
          root_multiplier$formal_unconditional_cover,
        root_multiplier_report_and_cover = root_multiplier$report_and_cover,
        derivative_assisted_requested = isTRUE(run_derivative_assisted_ablation),
        derivative_assisted_implemented = TRUE,
        derivative_assisted_applicable = isTRUE(run_derivative_assisted_ablation),
        derivative_assisted_reported = augmented_reported,
        derivative_assisted_failure = augmented_failure,
        derivative_assisted_failure_reason = augmented_failure_reason,
        derivative_assisted_same_multiplier_seed = if (
          isTRUE(run_derivative_assisted_ablation)
        ) TRUE else NA,
        derivative_assisted_critical_value = if (!is.null(augmented_fit)) {
          augmented_fit$critical_value
        } else NA_real_,
        derivative_assisted_critical_value_ratio = augmented_critical_value_ratio,
        derivative_assisted_value_maximum_half_width =
          augmented_value_maximum_half_width,
        derivative_assisted_derivative_maximum_half_width =
          augmented_derivative_maximum_half_width,
        derivative_assisted_value_width_inflation = augmented_width_inflation,
        derivative_assisted_guard_failed = if (!is.null(augmented_fit)) {
          augmented_fit$guard_failed
        } else NA,
        derivative_assisted_guard_minimum_ratio = if (!is.null(augmented_fit)) {
          augmented_fit$guard$minimum_certified_ratio
        } else NA_real_,
        derivative_assisted_guard_minimum_ratio_d0 = if (!is.null(augmented_fit)) {
          augmented_fit$guard$minimum_certified_ratio_d0
        } else NA_real_,
        derivative_assisted_guard_minimum_ratio_d1 = if (!is.null(augmented_fit)) {
          augmented_fit$guard$minimum_certified_ratio_d1
        } else NA_real_,
        derivative_assisted_guard_nodes_used = if (!is.null(augmented_fit)) {
          augmented_fit$guard$nodes_used
        } else NA_integer_,
        derivative_assisted_guard_maximum_level = if (!is.null(augmented_fit)) {
          augmented_fit$guard$maximum_level
        } else NA_integer_,
        derivative_assisted_value_enclosure_nodes_used = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$value_enclosure$nodes_used else NA_integer_,
        derivative_assisted_derivative_enclosure_nodes_used = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$derivative_enclosure$nodes_used else NA_integer_,
        derivative_assisted_value_unresolved_cells = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$value_enclosure$statistically_unresolved_cells else
          NA_integer_,
        derivative_assisted_derivative_unresolved_cells = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$derivative_enclosure$statistically_unresolved_cells else
          NA_integer_,
        derivative_assisted_value_variance_unresolved_cells = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$value_enclosure$variance_unresolved_cells else
          NA_integer_,
        derivative_assisted_derivative_variance_unresolved_cells = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$derivative_enclosure$variance_unresolved_cells else
          NA_integer_,
        derivative_assisted_enclosure_limit_hit = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$value_enclosure$limit_hit ||
          derivative_assisted$derivative_enclosure$limit_hit else NA,
        derivative_assisted_derivative_positive_throughout_K = if (
          !is.null(derivative_assisted) && augmented_reported
        ) augmented_derivative_positive_K else NA,
        derivative_assisted_derivative_negative_throughout_K = if (
          !is.null(derivative_assisted) && augmented_reported
        ) augmented_derivative_negative_K else NA,
        derivative_assisted_certificate_count = if (
          !is.null(derivative_assisted)
        ) nrow(augmented_certificates) else NA_integer_,
        # Backward-compatible name: "exact" is within the serialized local
        # anchor interval.  The next two fields make local versus whole-K
        # scope explicit.
        derivative_assisted_exact_root_certified = augmented_certified,
        derivative_assisted_local_exact_root_certified = augmented_certified,
        derivative_assisted_exact_root_on_K_certified = augmented_exact_root_on_K,
        derivative_assisted_exact_root_on_K_direction =
          augmented_exact_root_on_K_direction,
        derivative_assisted_exact_root_on_K_statement_true =
          augmented_exact_root_on_K_statement_true,
        derivative_assisted_false_exact_root_on_K_report =
          augmented_false_exact_root_on_K,
        derivative_assisted_certificate_scope = if (
          isTRUE(augmented_certified)
        ) "local_anchor_interval" else "",
        derivative_assisted_certificate_reported = augmented_certified,
        derivative_assisted_certificate_directions = if (
          nrow(augmented_certificates)
        ) collapse_character(augmented_certificates$direction) else "",
        derivative_assisted_certificate_left_anchors = if (
          nrow(augmented_certificates)
        ) collapse_numeric(augmented_certificates$p_left) else "",
        derivative_assisted_certificate_right_anchors = if (
          nrow(augmented_certificates)
        ) collapse_numeric(augmented_certificates$p_right) else "",
        derivative_assisted_certificate_derivative_margins = if (
          nrow(augmented_certificates)
        ) collapse_numeric(augmented_certificates$derivative_margin) else "",
        derivative_assisted_certificate_statement_true = if (
          nrow(augmented_certificates)
        ) all(augmented_certificate_truth) else NA,
        derivative_assisted_false_exact_root_report = augmented_false_report,
        derivative_assisted_value_grid_coverage = augmented_value_grid_coverage,
        derivative_assisted_derivative_grid_coverage =
          augmented_derivative_grid_coverage,
        derivative_assisted_joint_grid_coverage = augmented_joint_grid_coverage,
        derivative_assisted_covariance_01_min = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$covariance_min else NA_real_,
        derivative_assisted_covariance_01_max = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$covariance_max else NA_real_,
        derivative_assisted_correlation_01_min = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$correlation_min else NA_real_,
        derivative_assisted_correlation_01_max = if (
          !is.null(derivative_assisted)
        ) derivative_assisted$correlation_max else NA_real_,
        pointwise_root_coverage = pointwise_root_coverage,
        pointwise_grid_curve_coverage = pointwise_grid_curve_coverage,
        pointwise_false_sign_report = pointwise_false_sign,
        pointwise_reversal_report = pointwise_reversal,
        pointwise_no_root_report = pointwise_no_root,
        any_primary_false_report = any_primary_false,
        stringsAsFactors = FALSE
      )
      row <- cbind(row, as.data.frame(tail, stringsAsFactors = FALSE))
      rows[[row_counter]] <- row
    }
  }
  result <- do.call(rbind, rows)
  row.names(result) <- NULL
  result$joint_band_covers_truth_on_grid <- FALSE
  result$joint_band_covers_truth_on_audit_grid <- FALSE
  result$joint_continuum_band_covers_truth <- FALSE
  result$joint_any_primary_false_report <- FALSE
  result$joint_pointwise_grid_curve_coverage <- FALSE
  result$joint_derivative_assisted_grid_coverage <- NA
  result$joint_derivative_assisted_false_exact_root_report <- NA
  # Backward-compatible joint_* fields are joint across moment types within a
  # fixed interval.  They are not simultaneous outcomes across K1/K2/K3.
  for (interval_id in unique(result$interval_id)) {
    indices <- which(result$interval_id == interval_id)
    result$joint_band_covers_truth_on_grid[indices] <-
      all(result$band_covers_truth_on_grid[indices])
    result$joint_band_covers_truth_on_audit_grid[indices] <-
      all(result$band_covers_truth_on_audit_grid[indices])
    result$joint_continuum_band_covers_truth[indices] <-
      all(result$continuum_band_covers_truth[indices])
    result$joint_any_primary_false_report[indices] <-
      any(result$any_primary_false_report[indices])
    result$joint_pointwise_grid_curve_coverage[indices] <-
      all(result$pointwise_grid_curve_coverage[indices])
    if (all(!is.na(result$derivative_assisted_joint_grid_coverage[indices]))) {
      result$joint_derivative_assisted_grid_coverage[indices] <-
        all(result$derivative_assisted_joint_grid_coverage[indices])
    }
    if (all(!is.na(
      result$derivative_assisted_false_exact_root_report[indices]
    )) && all(!is.na(
      result$derivative_assisted_false_exact_root_on_K_report[indices]
    ))) {
      result$joint_derivative_assisted_false_exact_root_report[indices] <-
        any(
          result$derivative_assisted_false_exact_root_report[indices] |
            result$derivative_assisted_false_exact_root_on_K_report[indices]
        )
    }
  }
  # Bridge.6 adds the complementary family across all configured intervals,
  # separately for each moment type.  The same scalar event is repeated on the
  # corresponding K rows so downstream group summaries remain self-contained;
  # the summarizer also emits a deduplicated ALL_K table with one unit event.
  result$configured_interval_count <- NA_integer_
  result$across_interval_joint_band_covers_truth_on_grid <- FALSE
  result$across_interval_joint_band_covers_truth_on_audit_grid <- FALSE
  result$across_interval_joint_continuum_band_covers_truth <- FALSE
  result$across_interval_joint_any_primary_false_report <- FALSE
  for (moment_type in unique(result$moment_type)) {
    indices <- which(result$moment_type == moment_type)
    expected_interval_ids <- sort(names(intervals))
    observed_interval_ids <- sort(unique(as.character(
      result$interval_id[indices]
    )))
    assert_true(
      identical(observed_interval_ids, expected_interval_ids) &&
        length(indices) == length(expected_interval_ids),
      sprintf(
        "Moment type %s is missing or duplicates a configured interval.",
        moment_type
      )
    )
    result$configured_interval_count[indices] <-
      length(expected_interval_ids)
    result$across_interval_joint_band_covers_truth_on_grid[indices] <-
      all(result$band_covers_truth_on_grid[indices])
    result$across_interval_joint_band_covers_truth_on_audit_grid[indices] <-
      all(result$band_covers_truth_on_audit_grid[indices])
    result$across_interval_joint_continuum_band_covers_truth[indices] <-
      all(result$continuum_band_covers_truth[indices])
    result$across_interval_joint_any_primary_false_report[indices] <-
      any(result$any_primary_false_report[indices])
  }
  result$elapsed_seconds <- proc.time()[["elapsed"]] - replication_started
  memory_snapshot <- gc()
  result$peak_r_heap_mb <- sum(memory_snapshot[, 6L])
  result
}
