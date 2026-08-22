# Exact data-generating mechanisms for the frozen Monte Carlo design.
#
# Every scenario supplies both a numeric sampler (r_x/r_y, retained for the
# worker API) and a diagnostic sampler (draw_x/draw_y) that preserves latent
# component counts. Population moments, derivatives, roots, multiplicities,
# structural Mellin budgets, tail-crossing scales, and magnitude quantiles are
# analytic and are used only to score the Monte Carlo experiment.

check_orders <- function(order, upper = Inf) {
  assert_true(
    is.numeric(order) && length(order) >= 1L &&
      all(is.finite(order)) && all(order > -1) && all(order < upper),
    sprintf("Moment orders must be finite and lie in (-1, %s).", upper)
  )
  as.numeric(order)
}

check_moment_type <- function(moment_type, allowed) {
  assert_true(
    length(moment_type) == 1L && moment_type %in% allowed,
    sprintf(
      "Unknown moment type '%s'; allowed types are %s.",
      paste(moment_type, collapse = ","), paste(allowed, collapse = ", ")
    )
  )
  moment_type
}

uniform_multiplier_bundle <- function(order, eta = 0.10) {
  order <- check_orders(order)
  lower <- 1 - eta
  upper <- 1 + eta
  exponent <- order + 1
  numerator <- upper^exponent - lower^exponent
  numerator_1 <- upper^exponent * log(upper) -
    lower^exponent * log(lower)
  numerator_2 <- upper^exponent * log(upper)^2 -
    lower^exponent * log(lower)^2
  denominator <- 2 * eta
  moment <- numerator / (denominator * exponent)
  derivative <- (
    numerator_1 / exponent - numerator / exponent^2
  ) / denominator
  second_derivative <- (
    numerator_2 / exponent -
      2 * numerator_1 / exponent^2 +
      2 * numerator / exponent^3
  ) / denominator
  list(
    moment = moment,
    derivative = derivative,
    second_derivative = second_derivative
  )
}

uniform_multiplier_moment <- function(order, eta = 0.10) {
  uniform_multiplier_bundle(order, eta)$moment
}

uniform_multiplier_moment_derivative <- function(order, eta = 0.10) {
  uniform_multiplier_bundle(order, eta)$derivative
}

uniform_multiplier_moment_second_derivative <- function(order, eta = 0.10) {
  uniform_multiplier_bundle(order, eta)$second_derivative
}

finite_scale_moment_bundle <- function(order, scales, probabilities, eta = 0.10) {
  order <- check_orders(order)
  assert_true(
    length(scales) == length(probabilities) &&
      all(is.finite(scales) & scales > 0) &&
      all(is.finite(probabilities) & probabilities >= 0) &&
      abs(sum(probabilities) - 1) <= 1e-12,
    "Invalid finite-scale mixture."
  )
  base <- uniform_multiplier_bundle(order, eta)
  log_scales <- log(scales)
  scale_moment <- vapply(order, function(p) {
    sum(probabilities * scales^p)
  }, numeric(1L))
  scale_derivative <- vapply(order, function(p) {
    sum(probabilities * scales^p * log_scales)
  }, numeric(1L))
  scale_second <- vapply(order, function(p) {
    sum(probabilities * scales^p * log_scales^2)
  }, numeric(1L))
  list(
    moment = base$moment * scale_moment,
    derivative = base$derivative * scale_moment +
      base$moment * scale_derivative,
    second_derivative = base$second_derivative * scale_moment +
      2 * base$derivative * scale_derivative +
      base$moment * scale_second
  )
}

finite_scale_cdf <- function(value, scales, probabilities, eta = 0.10) {
  value <- as.numeric(value)
  vapply(value, function(z) {
    if (is.na(z)) return(NA_real_)
    if (!is.finite(z)) return(if (z < 0) 0 else 1)
    if (z <= 0) return(0)
    component_cdf <- pmin(
      1,
      pmax(0, (z / scales - (1 - eta)) / (2 * eta))
    )
    sum(probabilities * component_cdf)
  }, numeric(1L))
}

finite_scale_quantile <- function(probability, scales, probabilities, eta = 0.10) {
  probability <- as.numeric(probability)
  assert_true(
    all(is.finite(probability) & probability >= 0 & probability <= 1),
    "Quantile probabilities must lie in [0,1]."
  )
  ordering <- order(scales)
  scales <- scales[ordering]
  probabilities <- probabilities[ordering]
  assert_true(
    all(head(scales * (1 + eta), -1L) < tail(scales * (1 - eta), -1L)),
    "The analytic mixture quantile requires disjoint ordered supports."
  )
  cumulative <- cumsum(probabilities)
  previous <- c(0, head(cumulative, -1L))
  vapply(probability, function(probability_now) {
    if (probability_now <= 0) return(scales[1L] * (1 - eta))
    if (probability_now >= 1) return(tail(scales, 1L) * (1 + eta))
    component <- which(cumulative >= probability_now)[1L]
    local_probability <-
      (probability_now - previous[component]) / probabilities[component]
    scales[component] * ((1 - eta) + 2 * eta * local_probability)
  }, numeric(1L))
}

draw_finite_scale <- function(n, scales, probabilities, eta = 0.10) {
  n <- as.integer(n)
  assert_true(length(n) == 1L && !is.na(n) && n >= 1L, "n must be positive.")
  component <- sample.int(
    length(scales), size = n, replace = TRUE, prob = probabilities
  )
  multiplier <- stats::runif(n, min = 1 - eta, max = 1 + eta)
  values <- multiplier * scales[component]
  component_counts <- tabulate(component, nbins = length(scales))
  list(
    values = values,
    latent = list(
      component = component - 1L,
      component_counts = component_counts,
      high_count = sum(component > 1L),
      highest_count = component_counts[length(component_counts)]
    )
  )
}

# A four-component signed mixture has a one-dimensional contrast space after
# imposing the probability constraint at p=0 and two further root constraints.
# Cofactors give that null vector without adding a package dependency.  The
# normalization is L1 because the bridge.6 probability-floor construction is
# invariant to the arbitrary scale of the null vector.
finite_mixture_null_contrast <- function(
    scales, root_orders, root_multiplicities) {
  scales <- as.numeric(scales)
  root_orders <- as.numeric(root_orders)
  root_multiplicities <- as.integer(root_multiplicities)
  assert_true(
    length(scales) == 4L && all(is.finite(scales) & scales > 0) &&
      all(diff(scales) > 0),
    "The signal-variance construction requires four increasing scales."
  )
  assert_true(
    length(root_orders) == length(root_multiplicities) &&
      all(is.finite(root_orders) & root_orders > 0) &&
      all(root_multiplicities %in% c(1L, 2L)) &&
      sum(root_multiplicities) == 2L,
    "The positive-root constraints must have total multiplicity two."
  )
  rows <- list(rep(1, length(scales)))
  for (index in seq_along(root_orders)) {
    root <- root_orders[index]
    rows[[length(rows) + 1L]] <- scales^root
    if (root_multiplicities[index] == 2L) {
      rows[[length(rows) + 1L]] <- scales^root * log(scales)
    }
  }
  constraint_matrix <- do.call(rbind, rows)
  assert_true(
    identical(dim(constraint_matrix), c(3L, 4L)),
    "The signal-variance constraint matrix must be 3 by 4."
  )
  cofactors <- vapply(seq_len(4L), function(column) {
    (-1)^(column + 1L) * det(constraint_matrix[, -column, drop = FALSE])
  }, numeric(1L))
  assert_true(
    all(is.finite(cofactors)) && sum(abs(cofactors)) > 1e-12,
    "The signal-variance constraint matrix is rank deficient."
  )
  contrast <- cofactors / sum(abs(cofactors))
  # The last coefficient fixes the otherwise arbitrary sign.  With ordered
  # supports this also orients d=p_Y-p_X as positive below the first root.
  if (tail(contrast, 1L) < 0) contrast <- -contrast
  reference_order <- root_orders[1L] / 2
  assert_true(
    tail(contrast, 1L) > 0 &&
      sum(contrast * scales^reference_order) > 0 &&
      max(abs(as.vector(constraint_matrix %*% contrast))) <= 5e-11,
    "The finite-mixture null contrast does not satisfy its root constraints."
  )
  contrast
}

finite_mixture_tail_crossings <- function(
    scales, probability_contrast, eta = 0.10, tolerance = 1e-12) {
  ordering <- order(scales)
  scales <- as.numeric(scales[ordering])
  probability_contrast <- as.numeric(probability_contrast[ordering])
  assert_true(
    length(scales) == length(probability_contrast) &&
      abs(sum(probability_contrast)) <= tolerance &&
      all(head(scales * (1 + eta), -1L) <
            tail(scales * (1 - eta), -1L)),
    "Tail crossings require a zero-sum contrast on disjoint supports."
  )
  crossings <- numeric(0L)
  for (index in seq_along(scales)) {
    high_mass <- if (index < length(scales)) {
      sum(probability_contrast[(index + 1L):length(scales)])
    } else 0
    component_mass <- probability_contrast[index]
    if (abs(component_mass) <= tolerance) next
    survival_fraction <- -high_mass / component_mass
    if (survival_fraction > tolerance &&
        survival_fraction < 1 - tolerance) {
      crossings <- c(
        crossings,
        scales[index] * (
          (1 + eta) - 2 * eta * survival_fraction
        )
      )
    }
  }
  sort(crossings)
}

build_signal_variance_power_spec <- function(
    kind, scales, root_orders, probability_floor = 0.02, eta = 0.10) {
  assert_true(
    length(kind) == 1L && kind %in% c("tangency", "two_root"),
    "Unknown signal-variance power-mixture kind."
  )
  multiplicities <- if (kind == "tangency") 2L else c(1L, 1L)
  normalized_contrast <- finite_mixture_null_contrast(
    scales, root_orders, multiplicities
  )
  assert_true(
    is.finite(probability_floor) && probability_floor > 0 &&
      probability_floor < 1 / 4,
    "The component probability floor must lie in (0,1/4)."
  )
  # With ||q||_1=1 and sum(q)=0, this is the largest contrast that leaves the
  # stated floor in both samples.  It avoids a nominally large contrast whose
  # high-order variance is generated by nearly absent components.
  probability_contrast <-
    2 * (1 - 4 * probability_floor) * normalized_contrast
  average_probabilities <-
    probability_floor + abs(probability_contrast) / 2
  probabilities_x <- average_probabilities - probability_contrast / 2
  probabilities_y <- average_probabilities + probability_contrast / 2
  assert_true(
    all(probabilities_x >= probability_floor - 1e-12) &&
      all(probabilities_y >= probability_floor - 1e-12) &&
      abs(sum(probabilities_x) - 1) <= 1e-12 &&
      abs(sum(probabilities_y) - 1) <= 1e-12,
    "The signal-variance probability-floor construction is invalid."
  )
  tail_crossings <- finite_mixture_tail_crossings(
    scales, probability_contrast, eta
  )
  assert_true(
    length(tail_crossings) == 2L,
    "The selected signal-variance DGP must have two tail sign changes."
  )
  list(
    kind = kind,
    eta = eta,
    scales = as.numeric(scales),
    root_orders = as.numeric(root_orders),
    root_multiplicities = multiplicities,
    probability_floor = probability_floor,
    normalized_contrast = normalized_contrast,
    probability_contrast = probability_contrast,
    average_probabilities = average_probabilities,
    probabilities_x = probabilities_x,
    probabilities_y = probabilities_y,
    tail_crossings = tail_crossings
  )
}

empty_root_table <- function() {
  data.frame(
    order = numeric(0L), multiplicity = integer(0L),
    direction = character(0L), stringsAsFactors = FALSE
  )
}

root_table <- function(order, multiplicity, direction) {
  data.frame(
    order = as.numeric(order), multiplicity = as.integer(multiplicity),
    direction = as.character(direction), stringsAsFactors = FALSE
  )
}

new_scenario <- function(
    scenario_id,
    label,
    family,
    moment_types,
    parameters,
    draw_x,
    draw_y,
    moment_x,
    moment_y,
    moment_derivative_x,
    moment_derivative_y,
    moment_second_derivative_x,
    moment_second_derivative_y,
    root_truth,
    identically_zero,
    structural_budgets,
    tail_crossing_scales,
    magnitude_quantile_x,
    magnitude_quantile_y,
    magnitude_cdf_x,
    magnitude_cdf_y,
    root_intervals = list(),
    latent_high_probabilities = list(x = NA_real_, y = NA_real_),
    moment_domain = c(-1, Inf)) {
  assert_true(length(unique(moment_types)) == length(moment_types),
              "Scenario moment types must be unique.")
  assert_true(all(moment_types %in% c("absolute", "signed")),
              "Unknown scenario moment type.")
  assert_true(setequal(names(root_truth), moment_types),
              "root_truth must be named by moment type.")
  assert_true(setequal(names(identically_zero), moment_types),
              "identically_zero must be named by moment type.")
  assert_true(setequal(names(structural_budgets), moment_types),
              "structural_budgets must be named by moment type.")
  assert_true(setequal(names(tail_crossing_scales), moment_types),
              "tail_crossing_scales must be named by moment type.")

  delta <- function(order, moment_type = moment_types[1L]) {
    check_moment_type(moment_type, moment_types)
    moment_y(order, moment_type) - moment_x(order, moment_type)
  }
  delta_derivative <- function(order, moment_type = moment_types[1L]) {
    check_moment_type(moment_type, moment_types)
    moment_derivative_y(order, moment_type) -
      moment_derivative_x(order, moment_type)
  }
  delta_second_derivative <- function(order, moment_type = moment_types[1L]) {
    check_moment_type(moment_type, moment_types)
    moment_second_derivative_y(order, moment_type) -
      moment_second_derivative_x(order, moment_type)
  }
  primary_type <- moment_types[1L]
  primary_roots <- root_truth[[primary_type]]
  primary_direction <- if (nrow(primary_roots) == 1L) {
    primary_roots$direction[1L]
  } else {
    NA_character_
  }
  primary_crossings <- tail_crossing_scales[[primary_type]]

  list(
    scenario_id = scenario_id,
    label = label,
    family = family,
    parameters = parameters,
    moment_types = moment_types,
    primary_moment_type = primary_type,
    moment_domain = moment_domain,
    root_truth = root_truth,
    root_intervals = root_intervals,
    identically_zero_by_type = identically_zero,
    structural_budgets = structural_budgets,
    tail_crossing_scales = tail_crossing_scales,
    latent_high_probabilities = latent_high_probabilities,
    latent_count_truth = function(n_x, n_y) {
      probability_x <- latent_high_probabilities$x
      probability_y <- latent_high_probabilities$y
      list(
        expected_high_x = if (is.finite(probability_x)) n_x * probability_x else NA_real_,
        expected_high_y = if (is.finite(probability_y)) n_y * probability_y else NA_real_,
        probability_zero_high_x = if (is.finite(probability_x)) {
          (1 - probability_x)^n_x
        } else {
          NA_real_
        },
        probability_zero_high_y = if (is.finite(probability_y)) {
          (1 - probability_y)^n_y
        } else {
          NA_real_
        }
      )
    },
    structural_budget = unname(structural_budgets[[primary_type]]),
    identically_zero = unname(identically_zero[[primary_type]]),
    positive_roots = primary_roots$order,
    root_direction = primary_direction,
    tail_crossing_scale = if (length(primary_crossings) == 1L) {
      unname(primary_crossings)
    } else {
      NA_real_
    },
    draw_x = draw_x,
    draw_y = draw_y,
    r_x = function(n) draw_x(n)$values,
    r_y = function(n) draw_y(n)$values,
    moment_x = moment_x,
    moment_y = moment_y,
    moment_derivative_x = moment_derivative_x,
    moment_derivative_y = moment_derivative_y,
    moment_second_derivative_x = moment_second_derivative_x,
    moment_second_derivative_y = moment_second_derivative_y,
    delta = delta,
    delta_derivative = delta_derivative,
    delta_second_derivative = delta_second_derivative,
    magnitude_quantile_x = magnitude_quantile_x,
    magnitude_quantile_y = magnitude_quantile_y,
    magnitude_cdf_x = magnitude_cdf_x,
    magnitude_cdf_y = magnitude_cdf_y
  )
}

lognormal_moment <- function(order, meanlog, sdlog) {
  order <- check_orders(order)
  exp(meanlog * order + 0.5 * sdlog^2 * order^2)
}

lognormal_moment_derivative <- function(order, meanlog, sdlog) {
  order <- check_orders(order)
  lognormal_moment(order, meanlog, sdlog) *
    (meanlog + sdlog^2 * order)
}

lognormal_moment_second_derivative <- function(order, meanlog, sdlog) {
  order <- check_orders(order)
  moment <- lognormal_moment(order, meanlog, sdlog)
  score <- meanlog + sdlog^2 * order
  moment * (score^2 + sdlog^2)
}

make_lognormal_scenario <- function(
    scenario_id,
    label,
    meanlog_x,
    sdlog_x,
    meanlog_y,
    sdlog_y,
    structural_budget,
    root_order = numeric(0L),
    root_direction = character(0L),
    identically_zero = FALSE) {
  default_root_interval <- c(0.75, 1.90)
  roots <- if (length(root_order) == 0L) empty_root_table() else {
    root_table(root_order, rep.int(1L, length(root_order)), root_direction)
  }
  crossing <- if (sdlog_x != sdlog_y) {
    exp((sdlog_y * meanlog_x - sdlog_x * meanlog_y) /
          (sdlog_y - sdlog_x))
  } else {
    numeric(0L)
  }
  moment_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    lognormal_moment(order, meanlog_x, sdlog_x)
  }
  moment_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    lognormal_moment(order, meanlog_y, sdlog_y)
  }
  derivative_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    lognormal_moment_derivative(order, meanlog_x, sdlog_x)
  }
  derivative_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    lognormal_moment_derivative(order, meanlog_y, sdlog_y)
  }
  second_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    lognormal_moment_second_derivative(order, meanlog_x, sdlog_x)
  }
  second_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    lognormal_moment_second_derivative(order, meanlog_y, sdlog_y)
  }
  new_scenario(
    scenario_id = scenario_id,
    label = label,
    family = "lognormal",
    moment_types = "absolute",
    parameters = list(
      meanlog_x = meanlog_x, sdlog_x = sdlog_x,
      meanlog_y = meanlog_y, sdlog_y = sdlog_y
    ),
    draw_x = function(n) list(
      values = stats::rlnorm(n, meanlog_x, sdlog_x), latent = list()
    ),
    draw_y = function(n) list(
      values = stats::rlnorm(n, meanlog_y, sdlog_y), latent = list()
    ),
    moment_x = moment_x,
    moment_y = moment_y,
    moment_derivative_x = derivative_x,
    moment_derivative_y = derivative_y,
    moment_second_derivative_x = second_x,
    moment_second_derivative_y = second_y,
    root_truth = list(absolute = roots),
    identically_zero = c(absolute = identically_zero),
    structural_budgets = c(absolute = structural_budget),
    tail_crossing_scales = list(absolute = crossing),
    magnitude_quantile_x = function(probability) {
      stats::qlnorm(probability, meanlog_x, sdlog_x)
    },
    magnitude_quantile_y = function(probability) {
      stats::qlnorm(probability, meanlog_y, sdlog_y)
    },
    magnitude_cdf_x = function(value) stats::plnorm(value, meanlog_x, sdlog_x),
    magnitude_cdf_y = function(value) stats::plnorm(value, meanlog_y, sdlog_y),
    root_intervals = if (
        length(root_order) == 1L &&
          root_order > default_root_interval[1L] &&
          root_order < default_root_interval[2L]) {
      list(absolute = list(default_root_interval))
    } else {
      list(absolute = list())
    },
    latent_high_probabilities = list(x = NA_real_, y = NA_real_),
    moment_domain = c(-1, Inf)
  )
}

make_rare_scenario <- function(kind, probability) {
  eta <- 0.10
  scale <- 4
  suffix <- if (abs(probability - 0.01) < 1e-14) "P01" else "P05"
  if (kind == "crossing") {
    probabilities_x <- c(1 - 5 * probability, 5 * probability)
    scales_x <- c(1, scale)
    probabilities_y <- c(1 - probability, probability)
    scales_y <- c(1, scale^2)
    label <- sprintf("Continuous rare-tail crossing, pi=%.2f", probability)
    roots <- root_table(1, 1L, "up")
    zero <- FALSE
    budget <- 1L
    crossing <- scale * (1 + 3 * eta / 5)
  } else if (kind == "no_root") {
    probabilities_x <- c(1 - 2 * probability, 2 * probability)
    scales_x <- c(1, scale)
    probabilities_y <- c(1 - probability, probability)
    scales_y <- c(1, scale^2)
    label <- sprintf("Continuous rare-tail no-root, pi=%.2f", probability)
    roots <- empty_root_table()
    zero <- FALSE
    budget <- 1L
    crossing <- scale
  } else if (kind == "equal") {
    probabilities_x <- c(1 - probability, probability)
    scales_x <- c(1, scale^2)
    probabilities_y <- probabilities_x
    scales_y <- scales_x
    label <- sprintf("Continuous rare-tail equality, pi=%.2f", probability)
    roots <- empty_root_table()
    zero <- TRUE
    budget <- NA_integer_
    crossing <- numeric(0L)
  } else {
    stopf("Unknown rare-tail kind: %s", kind)
  }
  scenario_id <- sprintf("RARE_%s_%s", toupper(kind), suffix)
  moment_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales_x, probabilities_x, eta)$moment
  }
  moment_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales_y, probabilities_y, eta)$moment
  }
  derivative_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales_x, probabilities_x, eta)$derivative
  }
  derivative_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales_y, probabilities_y, eta)$derivative
  }
  second_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales_x, probabilities_x, eta)$second_derivative
  }
  second_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales_y, probabilities_y, eta)$second_derivative
  }
  new_scenario(
    scenario_id = scenario_id,
    label = label,
    family = "continuous_rare_tail",
    moment_types = "absolute",
    parameters = list(
      kind = kind, probability = probability, eta = eta, scale = scale,
      scales_x = scales_x, probabilities_x = probabilities_x,
      scales_y = scales_y, probabilities_y = probabilities_y
    ),
    draw_x = function(n) draw_finite_scale(n, scales_x, probabilities_x, eta),
    draw_y = function(n) draw_finite_scale(n, scales_y, probabilities_y, eta),
    moment_x = moment_x,
    moment_y = moment_y,
    moment_derivative_x = derivative_x,
    moment_derivative_y = derivative_y,
    moment_second_derivative_x = second_x,
    moment_second_derivative_y = second_y,
    root_truth = list(absolute = roots),
    identically_zero = c(absolute = zero),
    structural_budgets = c(absolute = budget),
    tail_crossing_scales = list(absolute = crossing),
    magnitude_quantile_x = function(probability_now) {
      finite_scale_quantile(probability_now, scales_x, probabilities_x, eta)
    },
    magnitude_quantile_y = function(probability_now) {
      finite_scale_quantile(probability_now, scales_y, probabilities_y, eta)
    },
    magnitude_cdf_x = function(value) {
      finite_scale_cdf(value, scales_x, probabilities_x, eta)
    },
    magnitude_cdf_y = function(value) {
      finite_scale_cdf(value, scales_y, probabilities_y, eta)
    },
    root_intervals = if (kind == "crossing") {
      list(absolute = list(c(0.50, 1.45)))
    } else {
      list(absolute = list())
    },
    latent_high_probabilities = list(
      x = tail(probabilities_x, 1L), y = tail(probabilities_y, 1L)
    ),
    moment_domain = c(-1, Inf)
  )
}

make_mixed_sign_scenario <- function() {
  eta <- 0.10
  scale <- 4
  rho <- 0.05
  beta <- c(4, -5, 1)
  scales <- scale^(0:2)
  component_probabilities <- rep(1 / 3, 3L)
  positive_probabilities_y <- (1 + 3 * rho * beta) / 2

  draw_x <- function(n) {
    component <- sample.int(3L, n, replace = TRUE) - 1L
    multiplier <- stats::runif(n, 1 - eta, 1 + eta)
    sign_value <- ifelse(stats::rbinom(n, 1L, 0.5) == 1L, 1, -1)
    list(
      values = sign_value * multiplier * scale^component,
      latent = list(
        component = component,
        component_counts = tabulate(component + 1L, nbins = 3L),
        positive_count = sum(sign_value > 0),
        high_count = sum(component > 0L),
        highest_count = sum(component == 2L)
      )
    )
  }
  draw_y <- function(n) {
    component <- sample.int(3L, n, replace = TRUE) - 1L
    multiplier <- stats::runif(n, 1 - eta, 1 + eta)
    positive_probability <- positive_probabilities_y[component + 1L]
    sign_value <- ifelse(stats::runif(n) <= positive_probability, 1, -1)
    list(
      values = sign_value * multiplier * scale^component,
      latent = list(
        component = component,
        component_counts = tabulate(component + 1L, nbins = 3L),
        positive_count = sum(sign_value > 0),
        positive_count_by_component = tabulate(
          component[sign_value > 0] + 1L, nbins = 3L
        ),
        high_count = sum(component > 0L),
        highest_count = sum(component == 2L)
      )
    )
  }
  absolute_bundle <- function(order) {
    finite_scale_moment_bundle(order, scales, component_probabilities, eta)
  }
  signed_y_bundle <- function(order) {
    order <- check_orders(order)
    base <- uniform_multiplier_bundle(order, eta)
    log_scales <- log(scales)
    polynomial <- vapply(order, function(p) {
      rho * sum(beta * scales^p)
    }, numeric(1L))
    polynomial_1 <- vapply(order, function(p) {
      rho * sum(beta * scales^p * log_scales)
    }, numeric(1L))
    polynomial_2 <- vapply(order, function(p) {
      rho * sum(beta * scales^p * log_scales^2)
    }, numeric(1L))
    list(
      moment = base$moment * polynomial,
      derivative = base$derivative * polynomial + base$moment * polynomial_1,
      second_derivative = base$second_derivative * polynomial +
        2 * base$derivative * polynomial_1 + base$moment * polynomial_2
    )
  }
  moment_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, c("absolute", "signed"))
    if (moment_type == "absolute") absolute_bundle(order)$moment else {
      rep(0, length(order))
    }
  }
  moment_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, c("absolute", "signed"))
    if (moment_type == "absolute") absolute_bundle(order)$moment else {
      signed_y_bundle(order)$moment
    }
  }
  derivative_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, c("absolute", "signed"))
    if (moment_type == "absolute") absolute_bundle(order)$derivative else {
      rep(0, length(order))
    }
  }
  derivative_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, c("absolute", "signed"))
    if (moment_type == "absolute") absolute_bundle(order)$derivative else {
      signed_y_bundle(order)$derivative
    }
  }
  second_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, c("absolute", "signed"))
    if (moment_type == "absolute") absolute_bundle(order)$second_derivative else {
      rep(0, length(order))
    }
  }
  second_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, c("absolute", "signed"))
    if (moment_type == "absolute") absolute_bundle(order)$second_derivative else {
      signed_y_bundle(order)$second_derivative
    }
  }
  new_scenario(
    scenario_id = "MIXED_SIGN",
    label = "Mixed-sign directional-tail crossing",
    family = "mixed_sign_finite_mixture",
    moment_types = c("absolute", "signed"),
    parameters = list(
      eta = eta, scale = scale, rho = rho, beta = beta,
      positive_probabilities_y = positive_probabilities_y
    ),
    draw_x = draw_x,
    draw_y = draw_y,
    moment_x = moment_x,
    moment_y = moment_y,
    moment_derivative_x = derivative_x,
    moment_derivative_y = derivative_y,
    moment_second_derivative_x = second_x,
    moment_second_derivative_y = second_y,
    root_truth = list(
      absolute = empty_root_table(),
      signed = root_table(1, 1L, "up")
    ),
    identically_zero = c(absolute = TRUE, signed = FALSE),
    structural_budgets = c(absolute = NA_integer_, signed = 1L),
    tail_crossing_scales = list(absolute = numeric(0L), signed = 4.24),
    magnitude_quantile_x = function(probability) {
      finite_scale_quantile(probability, scales, component_probabilities, eta)
    },
    magnitude_quantile_y = function(probability) {
      finite_scale_quantile(probability, scales, component_probabilities, eta)
    },
    magnitude_cdf_x = function(value) {
      finite_scale_cdf(value, scales, component_probabilities, eta)
    },
    magnitude_cdf_y = function(value) {
      finite_scale_cdf(value, scales, component_probabilities, eta)
    },
    root_intervals = list(
      absolute = list(), signed = list(c(0.50, 1.45))
    ),
    latent_high_probabilities = list(x = 2 / 3, y = 2 / 3),
    moment_domain = c(-1, Inf)
  )
}

make_polynomial_scenario <- function(kind, strength = "baseline") {
  eta <- 0.10
  assert_true(
    length(strength) == 1L && strength %in% c("baseline", "strong"),
    "Polynomial-mixture strength must be 'baseline' or 'strong'."
  )
  strong <- identical(strength, "strong")
  if (kind == "tangency") {
    scale <- 4
    epsilon <- if (strong) 0.018 else 0.01
    beta <- c(-16, 24, -9, 1)
    scenario_id <- if (strong) "TANGENCY_STRONG" else "TANGENCY"
    label <- if (strong) {
      "Stronger bounded tangency with one double root"
    } else {
      "Bounded tangency with one double root"
    }
    roots <- root_table(1, 2L, "tangent")
    root_intervals <- list(absolute = list())
    tail_scales <- c(
      4 * ((1 - eta) + 2 * eta * (2 / 3)),
      16 * ((1 - eta) + 2 * eta * (8 / 9))
    )
  } else if (kind == "two_root") {
    scale <- 2
    epsilon <- if (strong) 0.03 else 0.02
    beta <- c(-8, 14, -7, 1)
    scenario_id <- if (strong) "TWO_ROOT_STRONG" else "TWO_ROOT"
    label <- if (strong) {
      "Stronger bounded mechanism with two simple roots"
    } else {
      "Bounded mechanism with two simple roots"
    }
    roots <- root_table(c(1, 2), c(1L, 1L), c("down", "up"))
    root_intervals <- list(
      absolute = list(c(0.50, 1.45), c(1.65, 2.35))
    )
    tail_scales <- c(
      2 * ((1 - eta) + 2 * eta * (4 / 7)),
      4 * ((1 - eta) + 2 * eta * (6 / 7))
    )
  } else {
    stopf("Unknown polynomial-mixture kind: %s", kind)
  }
  probabilities_x <- 1 / 4 - epsilon * beta / 2
  probabilities_y <- 1 / 4 + epsilon * beta / 2
  assert_true(
    all(probabilities_x > 0) && all(probabilities_y > 0) &&
      abs(sum(probabilities_x) - 1) <= 1e-14 &&
      abs(sum(probabilities_y) - 1) <= 1e-14,
    sprintf("Invalid polynomial-mixture probabilities for %s.", scenario_id)
  )
  scales <- scale^(0:3)
  moment_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_x, eta)$moment
  }
  moment_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_y, eta)$moment
  }
  derivative_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_x, eta)$derivative
  }
  derivative_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_y, eta)$derivative
  }
  second_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_x, eta)$second_derivative
  }
  second_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_y, eta)$second_derivative
  }
  new_scenario(
    scenario_id = scenario_id,
    label = label,
    family = "bounded_polynomial_mixture",
    moment_types = "absolute",
    parameters = list(
      kind = kind, eta = eta, scale = scale, epsilon = epsilon, beta = beta,
      scales = scales, probabilities_x = probabilities_x,
      probabilities_y = probabilities_y
    ),
    draw_x = function(n) draw_finite_scale(n, scales, probabilities_x, eta),
    draw_y = function(n) draw_finite_scale(n, scales, probabilities_y, eta),
    moment_x = moment_x,
    moment_y = moment_y,
    moment_derivative_x = derivative_x,
    moment_derivative_y = derivative_y,
    moment_second_derivative_x = second_x,
    moment_second_derivative_y = second_y,
    root_truth = list(absolute = roots),
    identically_zero = c(absolute = FALSE),
    structural_budgets = c(absolute = 2L),
    tail_crossing_scales = list(absolute = tail_scales),
    magnitude_quantile_x = function(probability) {
      finite_scale_quantile(probability, scales, probabilities_x, eta)
    },
    magnitude_quantile_y = function(probability) {
      finite_scale_quantile(probability, scales, probabilities_y, eta)
    },
    magnitude_cdf_x = function(value) {
      finite_scale_cdf(value, scales, probabilities_x, eta)
    },
    magnitude_cdf_y = function(value) {
      finite_scale_cdf(value, scales, probabilities_y, eta)
    },
    root_intervals = root_intervals,
    latent_high_probabilities = list(
      x = 1 - probabilities_x[1L], y = 1 - probabilities_y[1L]
    ),
    moment_domain = c(-1, Inf)
  )
}

make_signal_variance_power_scenario <- function(kind) {
  if (identical(kind, "tangency")) {
    scenario_id <- "TANGENCY_SV_POWER"
    label <- "Signal-variance-selected bounded tangency power design"
    scales <- c(1, 6, 24, 64)
    root_orders <- 0.60
    roots <- root_table(0.60, 2L, "tangent")
    root_intervals <- list(absolute = list())
    selection_objective <- "maximin_standardized_signal_at_p_root_plus_or_minus_0.25"
  } else if (identical(kind, "two_root")) {
    scenario_id <- "TWO_ROOT_SV_POWER"
    label <- "Signal-variance-selected bounded two-root power design"
    scales <- c(1, 8, 28, 64)
    root_orders <- c(0.75, 1.75)
    roots <- root_table(c(0.75, 1.75), c(1L, 1L), c("down", "up"))
    root_intervals <- list(
      absolute = list(Q1 = c(0.45, 1.05), Q2 = c(1.40, 2.10))
    )
    selection_objective <- paste(
      "maximin_of_the_maximum_standardized_signal_in_the_three_sign_regions"
    )
  } else {
    stopf("Unknown signal-variance power-mixture kind: %s", kind)
  }
  specification <- build_signal_variance_power_spec(
    kind = kind, scales = scales, root_orders = root_orders,
    probability_floor = 0.02, eta = 0.10
  )
  probabilities_x <- specification$probabilities_x
  probabilities_y <- specification$probabilities_y
  eta <- specification$eta
  moment_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_x, eta)$moment
  }
  moment_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_y, eta)$moment
  }
  derivative_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_x, eta)$derivative
  }
  derivative_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(order, scales, probabilities_y, eta)$derivative
  }
  second_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(
      order, scales, probabilities_x, eta
    )$second_derivative
  }
  second_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    finite_scale_moment_bundle(
      order, scales, probabilities_y, eta
    )$second_derivative
  }
  new_scenario(
    scenario_id = scenario_id,
    label = label,
    family = "bounded_signal_variance_power_mixture",
    moment_types = "absolute",
    parameters = c(
      specification,
      list(
        analytic_screen_id = "signal_variance_screen_v001",
        selection_objective = selection_objective
      )
    ),
    draw_x = function(n) draw_finite_scale(n, scales, probabilities_x, eta),
    draw_y = function(n) draw_finite_scale(n, scales, probabilities_y, eta),
    moment_x = moment_x,
    moment_y = moment_y,
    moment_derivative_x = derivative_x,
    moment_derivative_y = derivative_y,
    moment_second_derivative_x = second_x,
    moment_second_derivative_y = second_y,
    root_truth = list(absolute = roots),
    identically_zero = c(absolute = FALSE),
    structural_budgets = c(absolute = 2L),
    tail_crossing_scales = list(absolute = specification$tail_crossings),
    magnitude_quantile_x = function(probability) {
      finite_scale_quantile(probability, scales, probabilities_x, eta)
    },
    magnitude_quantile_y = function(probability) {
      finite_scale_quantile(probability, scales, probabilities_y, eta)
    },
    magnitude_cdf_x = function(value) {
      finite_scale_cdf(value, scales, probabilities_x, eta)
    },
    magnitude_cdf_y = function(value) {
      finite_scale_cdf(value, scales, probabilities_y, eta)
    },
    root_intervals = root_intervals,
    latent_high_probabilities = list(
      x = 1 - probabilities_x[1L], y = 1 - probabilities_y[1L]
    ),
    moment_domain = c(-1, Inf)
  )
}

pareto_moment_bundle <- function(order, lower, shape) {
  order <- check_orders(order, upper = shape)
  moment <- shape * lower^order / (shape - order)
  score <- log(lower) + 1 / (shape - order)
  list(
    moment = moment,
    derivative = moment * score,
    second_derivative = moment * (score^2 + 1 / (shape - order)^2)
  )
}

make_pareto_scenario <- function() {
  lower_x <- 1
  shape_x <- 9
  lower_y <- 81 / 88
  shape_y <- 11 / 2
  draw_pareto <- function(n, lower, shape) {
    uniform <- pmax(stats::runif(n), .Machine$double.xmin)
    list(values = lower * uniform^(-1 / shape), latent = list())
  }
  moment_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    pareto_moment_bundle(order, lower_x, shape_x)$moment
  }
  moment_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    pareto_moment_bundle(order, lower_y, shape_y)$moment
  }
  derivative_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    pareto_moment_bundle(order, lower_x, shape_x)$derivative
  }
  derivative_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    pareto_moment_bundle(order, lower_y, shape_y)$derivative
  }
  second_x <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    pareto_moment_bundle(order, lower_x, shape_x)$second_derivative
  }
  second_y <- function(order, moment_type = "absolute") {
    check_moment_type(moment_type, "absolute")
    pareto_moment_bundle(order, lower_y, shape_y)$second_derivative
  }
  new_scenario(
    scenario_id = "PARETO",
    label = "Heavy-tailed Pareto stress design",
    family = "pareto",
    moment_types = "absolute",
    parameters = list(
      lower_x = lower_x, shape_x = shape_x,
      lower_y = lower_y, shape_y = shape_y
    ),
    draw_x = function(n) draw_pareto(n, lower_x, shape_x),
    draw_y = function(n) draw_pareto(n, lower_y, shape_y),
    moment_x = moment_x,
    moment_y = moment_y,
    moment_derivative_x = derivative_x,
    moment_derivative_y = derivative_y,
    moment_second_derivative_x = second_x,
    moment_second_derivative_y = second_y,
    root_truth = list(absolute = root_table(1, 1L, "up")),
    identically_zero = c(absolute = FALSE),
    structural_budgets = c(absolute = 1L),
    tail_crossing_scales = list(absolute = (81 / 88)^(-11 / 7)),
    magnitude_quantile_x = function(probability) {
      lower_x * (1 - probability)^(-1 / shape_x)
    },
    magnitude_quantile_y = function(probability) {
      lower_y * (1 - probability)^(-1 / shape_y)
    },
    magnitude_cdf_x = function(value) {
      ifelse(value < lower_x, 0, 1 - (lower_x / value)^shape_x)
    },
    magnitude_cdf_y = function(value) {
      ifelse(value < lower_y, 0, 1 - (lower_y / value)^shape_y)
    },
    root_intervals = list(absolute = list(c(0.50, 1.45))),
    latent_high_probabilities = list(x = NA_real_, y = NA_real_),
    moment_domain = c(-1, shape_y)
  )
}

build_scenario_registry <- function() {
  list(
    LN_EQUAL = make_lognormal_scenario(
      "LN_EQUAL", "Lognormal equality", 0, 1 / 5, 0, 1 / 5,
      NA_integer_, identically_zero = TRUE
    ),
    LN_MODERATE = make_lognormal_scenario(
      "LN_MODERATE", "Lognormal moderate upcrossing",
      0, 1 / 5, -2 / 25, 2 / 5, 1L, 4 / 3, "up"
    ),
    LN_HIGH_VARIANCE = make_lognormal_scenario(
      "LN_HIGH_VARIANCE", "Lognormal high-variance upcrossing",
      0, 3 / 5, -3 / 10, 9 / 10, 1L, 4 / 3, "up"
    ),
    LN_WEAK = make_lognormal_scenario(
      "LN_WEAK", "Lognormal weak upcrossing",
      0, 7 / 10, -29 / 600, 3 / 4, 1L, 4 / 3, "up"
    ),
    LN_NO_ROOT = make_lognormal_scenario(
      "LN_NO_ROOT", "Lognormal no positive root",
      0, 3 / 5, -1 / 5, 3 / 5, 0L
    ),
    LN_ROOT_OUTSIDE = make_lognormal_scenario(
      "LN_ROOT_OUTSIDE", "Lognormal root outside K3",
      0, 3 / 5, -27 / 40, 9 / 10, 1L, 3, "up"
    ),
    RARE_CROSSING_P01 = make_rare_scenario("crossing", 0.01),
    RARE_NO_ROOT_P01 = make_rare_scenario("no_root", 0.01),
    RARE_EQUAL_P01 = make_rare_scenario("equal", 0.01),
    RARE_CROSSING_P05 = make_rare_scenario("crossing", 0.05),
    RARE_NO_ROOT_P05 = make_rare_scenario("no_root", 0.05),
    RARE_EQUAL_P05 = make_rare_scenario("equal", 0.05),
    MIXED_SIGN = make_mixed_sign_scenario(),
    TANGENCY = make_polynomial_scenario("tangency"),
    TWO_ROOT = make_polynomial_scenario("two_root"),
    PARETO = make_pareto_scenario(),
    TANGENCY_STRONG = make_polynomial_scenario("tangency", "strong"),
    TWO_ROOT_STRONG = make_polynomial_scenario("two_root", "strong"),
    TANGENCY_SV_POWER = make_signal_variance_power_scenario("tangency"),
    TWO_ROOT_SV_POWER = make_signal_variance_power_scenario("two_root")
  )
}

validate_scenario_registry <- function(registry, tolerance = 2e-7) {
  expected <- c(
    "LN_EQUAL", "LN_MODERATE", "LN_HIGH_VARIANCE", "LN_WEAK",
    "LN_NO_ROOT", "LN_ROOT_OUTSIDE",
    "RARE_CROSSING_P01", "RARE_NO_ROOT_P01", "RARE_EQUAL_P01",
    "RARE_CROSSING_P05", "RARE_NO_ROOT_P05", "RARE_EQUAL_P05",
    "MIXED_SIGN", "TANGENCY", "TWO_ROOT", "PARETO",
    "TANGENCY_STRONG", "TWO_ROOT_STRONG",
    "TANGENCY_SV_POWER", "TWO_ROOT_SV_POWER"
  )
  assert_true(identical(names(registry), expected),
              "The scenario registry is not the frozen bridge.6 20-DGP design.")
  probe_orders <- c(0.25, 0.80, 1.20, 1.70, 2.50)
  for (scenario in registry) {
    for (moment_type in scenario$moment_types) {
      direct_delta <- scenario$moment_y(probe_orders, moment_type) -
        scenario$moment_x(probe_orders, moment_type)
      assert_true(
        max(abs(direct_delta - scenario$delta(probe_orders, moment_type))) <=
          tolerance * max(1, max(abs(direct_delta))),
        sprintf("Moment contrast identity failed for %s/%s.",
                scenario$scenario_id, moment_type)
      )
      step <- 1e-5
      finite_difference <- (
        scenario$delta(probe_orders + step, moment_type) -
          scenario$delta(probe_orders - step, moment_type)
      ) / (2 * step)
      analytic <- scenario$delta_derivative(probe_orders, moment_type)
      assert_true(
        max(abs(finite_difference - analytic)) <=
          2e-5 * max(1, max(abs(analytic))),
        sprintf("First derivative validation failed for %s/%s.",
                scenario$scenario_id, moment_type)
      )
      truth <- scenario$root_truth[[moment_type]]
      if (nrow(truth) > 0L) {
        root_values <- scenario$delta(truth$order, moment_type)
        assert_true(
          max(abs(root_values)) <= tolerance * max(
            1,
            max(abs(scenario$moment_x(truth$order, moment_type))),
            max(abs(scenario$moment_y(truth$order, moment_type)))
          ),
          sprintf("Exact root validation failed for %s/%s.",
                  scenario$scenario_id, moment_type)
        )
        for (index in seq_len(nrow(truth))) {
          derivative <- scenario$delta_derivative(truth$order[index], moment_type)
          if (truth$multiplicity[index] == 1L) {
            expected_sign <- if (truth$direction[index] == "up") 1 else -1
            assert_true(
              is.finite(derivative) && expected_sign * derivative > 1e-10,
              sprintf("Simple-root direction failed for %s/%s.",
                      scenario$scenario_id, moment_type)
            )
          } else if (truth$multiplicity[index] == 2L) {
            curvature <- scenario$delta_second_derivative(
              truth$order[index], moment_type
            )
            assert_true(abs(derivative) <= 1e-8 && curvature > 0,
                        "Tangency multiplicity validation failed.")
          }
        }
      }
      if (isTRUE(scenario$identically_zero_by_type[[moment_type]])) {
        assert_true(
          max(abs(scenario$delta(probe_orders, moment_type))) <= 1e-12,
          sprintf("Equality-null validation failed for %s/%s.",
                  scenario$scenario_id, moment_type)
        )
      }
    }
    probabilities <- c(0.95, 0.99)
    for (sample_name in c("x", "y")) {
      quantile_function <- scenario[[paste0("magnitude_quantile_", sample_name)]]
      cdf_function <- scenario[[paste0("magnitude_cdf_", sample_name)]]
      quantiles <- quantile_function(probabilities)
      assert_true(all(is.finite(quantiles)) && all(diff(quantiles) > 0),
                  sprintf("Magnitude quantiles failed for %s/%s.",
                          scenario$scenario_id, sample_name))
      assert_true(
        max(abs(cdf_function(quantiles) - probabilities)) <= 2e-10,
        sprintf("Magnitude CDF/quantile identity failed for %s/%s.",
                scenario$scenario_id, sample_name)
      )
    }
  }
  lognormal_root_scenarios <- c(
    "LN_MODERATE", "LN_HIGH_VARIANCE", "LN_WEAK"
  )
  for (scenario_id in lognormal_root_scenarios) {
    scenario <- registry[[scenario_id]]
    intervals <- scenario$root_intervals$absolute
    roots <- scenario$root_truth$absolute$order
    assert_true(
      length(intervals) == 1L &&
        max(abs(intervals[[1L]] - c(0.75, 1.90))) <= 1e-14 &&
        length(roots) == 1L &&
        roots > intervals[[1L]][1L] && roots < intervals[[1L]][2L],
      sprintf(
        "Permanent lognormal root interval failed for %s.", scenario_id
      )
    )
  }
  assert_true(
    length(registry$LN_ROOT_OUTSIDE$root_intervals$absolute) == 0L,
    "The outside-K3 lognormal must not acquire a root-inference interval."
  )

  polynomial_pairs <- list(
    c("TANGENCY", "TANGENCY_STRONG"),
    c("TWO_ROOT", "TWO_ROOT_STRONG")
  )
  for (scenario_pair in polynomial_pairs) {
    baseline <- registry[[scenario_pair[1L]]]
    strong <- registry[[scenario_pair[2L]]]
    scale_factor <- strong$parameters$epsilon / baseline$parameters$epsilon
    assert_true(
      identical(strong$moment_types, baseline$moment_types) &&
        identical(strong$root_truth, baseline$root_truth) &&
        identical(strong$structural_budgets, baseline$structural_budgets) &&
        identical(strong$root_intervals, baseline$root_intervals) &&
        identical(strong$tail_crossing_scales,
                  baseline$tail_crossing_scales) &&
        identical(strong$parameters$scale, baseline$parameters$scale) &&
        identical(strong$parameters$beta, baseline$parameters$beta) &&
        max(abs(
          strong$delta(probe_orders) -
            scale_factor * baseline$delta(probe_orders)
        )) <= tolerance * max(1, max(abs(strong$delta(probe_orders)))),
      sprintf("Strong-companion identity failed for %s.", scenario_pair[2L])
    )
  }

  polynomial_weights <- list(
    TANGENCY = list(
      epsilon = 0.01,
      x = c(0.330, 0.130, 0.295, 0.245),
      y = c(0.170, 0.370, 0.205, 0.255)
    ),
    TWO_ROOT = list(
      epsilon = 0.02,
      x = c(0.330, 0.110, 0.320, 0.240),
      y = c(0.170, 0.390, 0.180, 0.260)
    ),
    TANGENCY_STRONG = list(
      epsilon = 0.018,
      x = c(0.394, 0.034, 0.331, 0.241),
      y = c(0.106, 0.466, 0.169, 0.259)
    ),
    TWO_ROOT_STRONG = list(
      epsilon = 0.03,
      x = c(0.370, 0.040, 0.355, 0.235),
      y = c(0.130, 0.460, 0.145, 0.265)
    )
  )
  for (scenario_id in names(polynomial_weights)) {
    scenario <- registry[[scenario_id]]
    weights <- polynomial_weights[[scenario_id]]
    assert_true(
      abs(scenario$parameters$epsilon - weights$epsilon) <= 1e-14 &&
        max(abs(scenario$parameters$probabilities_x - weights$x)) <= 1e-14 &&
        max(abs(scenario$parameters$probabilities_y - weights$y)) <= 1e-14,
      sprintf("Polynomial-mixture weights failed for %s.", scenario_id)
    )
    order <- probe_orders
    z <- scenario$parameters$scale^order
    polynomial <- if (scenario$parameters$kind == "tangency") {
      (z - 1) * (z - 4)^2
    } else {
      (z - 1) * (z - 2) * (z - 4)
    }
    expected_delta <- uniform_multiplier_moment(order, 0.10) *
      scenario$parameters$epsilon * polynomial
    assert_true(
      max(abs(scenario$delta(order) - expected_delta)) <=
        2e-12 * max(1, max(abs(expected_delta))),
      sprintf("Polynomial contrast formula failed for %s.", scenario_id)
    )
  }
  assert_true(
    abs(registry$RARE_CROSSING_P01$delta_derivative(1) -
          0.12 * log(4)) <= 1e-12,
    "Rare-crossing slope formula failed."
  )
  assert_true(
    abs(registry$MIXED_SIGN$delta_derivative(1, "signed") -
          0.6 * log(4)) <= 1e-12,
    "Mixed-sign slope formula failed."
  )
  for (scenario_id in c("TANGENCY", "TANGENCY_STRONG")) {
    scenario <- registry[[scenario_id]]
    expected_curvature <- 96 * scenario$parameters$epsilon * log(4)^2
    assert_true(
      abs(scenario$delta_second_derivative(1) - expected_curvature) <= 1e-10,
      sprintf("Tangency curvature formula failed for %s.", scenario_id)
    )
  }
  for (scenario_id in c("TWO_ROOT", "TWO_ROOT_STRONG")) {
    scenario <- registry[[scenario_id]]
    epsilon <- scenario$parameters$epsilon
    expected_slopes <- c(
      -4 * epsilon * log(2),
      24 * epsilon * uniform_multiplier_moment(2, 0.10) * log(2)
    )
    assert_true(
      max(abs(scenario$delta_derivative(c(1, 2)) - expected_slopes)) <=
        1e-11,
      sprintf("Two-root slope formulas failed for %s.", scenario_id)
    )
  }
  selected_power_specs <- list(
    TANGENCY_SV_POWER = list(
      kind = "tangency", scales = c(1, 6, 24, 64), roots = 0.60,
      multiplicities = 2L, directions = "tangent"
    ),
    TWO_ROOT_SV_POWER = list(
      kind = "two_root", scales = c(1, 8, 28, 64),
      roots = c(0.75, 1.75), multiplicities = c(1L, 1L),
      directions = c("down", "up")
    )
  )
  for (scenario_id in names(selected_power_specs)) {
    scenario <- registry[[scenario_id]]
    expected_spec <- selected_power_specs[[scenario_id]]
    probabilities <- c(
      scenario$parameters$probabilities_x,
      scenario$parameters$probabilities_y
    )
    assert_true(
      identical(scenario$parameters$kind, expected_spec$kind) &&
        max(abs(scenario$parameters$scales - expected_spec$scales)) <= 1e-14 &&
        max(abs(
          scenario$root_truth$absolute$order - expected_spec$roots
        )) <= 1e-12 &&
        identical(
          scenario$root_truth$absolute$multiplicity,
          as.integer(expected_spec$multiplicities)
        ) &&
        identical(
          scenario$root_truth$absolute$direction,
          expected_spec$directions
        ) &&
        abs(sum(abs(scenario$parameters$normalized_contrast)) - 1) <= 1e-12 &&
        tail(scenario$parameters$normalized_contrast, 1L) > 0 &&
        min(probabilities) >= 0.02 - 1e-12 &&
        length(scenario$tail_crossing_scales$absolute) == 2L &&
        scenario$structural_budgets[["absolute"]] == 2L,
      sprintf("Signal-variance-selected DGP contract failed for %s.", scenario_id)
    )
  }
  assert_true(
    abs(registry$PARETO$delta(1)) <= 1e-12 &&
      registry$PARETO$delta_derivative(1) > 0,
    "Pareto root formula failed."
  )
  invisible(TRUE)
}

.dgp_registry_cache <- new.env(parent = emptyenv())

scenario_registry <- function() {
  if (!exists("registry", envir = .dgp_registry_cache, inherits = FALSE)) {
    registry <- build_scenario_registry()
    validate_scenario_registry(registry)
    assign("registry", registry, envir = .dgp_registry_cache)
  }
  get("registry", envir = .dgp_registry_cache, inherits = FALSE)
}

get_scenario <- function(scenario_id) {
  aliases <- c(
    LN_UP_MODERATE = "LN_MODERATE",
    LN_UP_STRONG = "LN_HIGH_VARIANCE",
    LN_UP_WEAK = "LN_WEAK"
  )
  if (scenario_id %in% names(aliases)) scenario_id <- aliases[[scenario_id]]
  registry <- scenario_registry()
  if (!scenario_id %in% names(registry)) {
    stopf(
      "Unknown scenario '%s'. Available scenarios: %s",
      scenario_id, paste(names(registry), collapse = ", ")
    )
  }
  registry[[scenario_id]]
}

truth_on_interval <- function(
    scenario,
    p_min,
    p_max,
    moment_type = scenario$primary_moment_type) {
  check_moment_type(moment_type, scenario$moment_types)
  zero <- isTRUE(scenario$identically_zero_by_type[[moment_type]])
  if (zero) {
    return(list(
      moment_type = moment_type,
      roots = numeric(0L),
      root_multiplicities = integer(0L),
      root_directions = character(0L),
      root_count = Inf,
      total_multiplicity = Inf,
      has_any_root = TRUE,
      has_reversal = FALSE,
      unique_root = NA_real_,
      identically_zero = TRUE,
      structural_budget = NA_integer_
    ))
  }
  truth <- scenario$root_truth[[moment_type]]
  keep <- truth$order >= p_min & truth$order <= p_max
  interval_truth <- truth[keep, , drop = FALSE]
  interior <- interval_truth$order > p_min & interval_truth$order < p_max
  reversal <- interval_truth$direction %in% c("up", "down") & interior
  list(
    moment_type = moment_type,
    roots = interval_truth$order,
    root_multiplicities = interval_truth$multiplicity,
    root_directions = interval_truth$direction,
    root_count = nrow(interval_truth),
    total_multiplicity = sum(interval_truth$multiplicity),
    has_any_root = nrow(interval_truth) > 0L,
    has_reversal = any(reversal),
    unique_root = if (nrow(interval_truth) == 1L && interior) {
      interval_truth$order
    } else {
      NA_real_
    },
    identically_zero = FALSE,
    structural_budget = unname(scenario$structural_budgets[[moment_type]])
  )
}

# Validate once on source. Subsequent get_scenario() calls use the cached,
# already checked registry and incur no validation cost inside Monte Carlo.
scenario_registry()
