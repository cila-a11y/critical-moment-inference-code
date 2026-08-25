args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else "."

band <- utils::read.csv(file.path(root, "moment_band.csv"), check.names = FALSE)
critical <- band$critical_value[[1L]]
output <- file.path(root, "figures", "figure_03_hillstrom_moment_contrast_jns.pdf")

grDevices::cairo_pdf(output, width = 8.25, height = 3.35, family = "sans", onefile = TRUE)
graphics::par(
  mfrow = c(1L, 2L),
  mar = c(4.15, 4.25, 2.3, 0.8),
  mgp = c(2.4, 0.7, 0),
  tcl = -0.25,
  cex = 0.82,
  cex.axis = 0.9,
  cex.lab = 0.95,
  cex.main = 1.0,
  las = 1
)

blue <- "#0868AC"
band_blue <- "#9ECAE1"
positive_green <- "#DFF0D8"
threshold_red <- "#CB181D"

graphics::plot(
  band$order,
  band$delta_hat,
  type = "n",
  xlim = range(band$order),
  ylim = range(c(0, band$upper)) * c(1, 1.04),
  xlab = expression(paste("Moment order ", italic(p))),
  ylab = expression(hat(Delta)[a * "," * N](italic(p))),
  main = "(a) Moment contrast",
  font.main = 2,
  adj = 0
)
graphics::rect(
  min(band$order),
  par("usr")[[3L]],
  max(band$order),
  par("usr")[[4L]],
  col = positive_green,
  border = NA
)
graphics::polygon(
  c(band$order, rev(band$order)),
  c(band$lower, rev(band$upper)),
  col = band_blue,
  border = NA
)
graphics::lines(band$order, band$delta_hat, col = blue, lwd = 2)
graphics::abline(h = 0, lty = 2, lwd = 0.8)
graphics::abline(v = 1, col = "grey45", lty = 3, lwd = 0.8)
graphics::box()
graphics::legend(
  "topleft",
  legend = c(
    expression(hat(Delta)[a * "," * N](italic(p))),
    "Implemented 95% simultaneous band",
    "Whole-cell positive classification",
    expression(italic(p) == 1 ~ "reference")
  ),
  col = c(blue, band_blue, positive_green, "grey45"),
  lty = c(1, NA, NA, 3),
  lwd = c(2, NA, NA, 0.8),
  pch = c(NA, 15, 15, NA),
  pt.cex = c(NA, 1.45, 1.45, NA),
  bty = "o",
  bg = "white",
  box.col = "grey70",
  cex = 0.67,
  inset = 0.015,
  y.intersp = 0.9,
  x.intersp = 0.6
)

graphics::plot(
  band$order,
  band$studentized,
  type = "l",
  col = "#08519C",
  lwd = 2,
  xlim = range(band$order),
  ylim = c(-5.6, max(band$studentized) * 1.05),
  xlab = expression(paste("Moment order ", italic(p))),
  ylab = "Studentized contrast",
  main = "(b) Studentized contrast",
  font.main = 2,
  adj = 0
)
graphics::abline(h = 0, col = "black", lwd = 0.8)
graphics::abline(h = c(-critical, critical), col = threshold_red, lty = 2, lwd = 1)
graphics::abline(v = 1, col = "grey45", lty = 3, lwd = 0.8)
graphics::legend(
  "bottomleft",
  legend = c(
    "Studentized contrast",
    bquote("Simultaneous thresholds: " %+-% .(formatC(critical, digits = 3L, format = "f"))),
    expression(italic(p) == 1 ~ "reference")
  ),
  col = c("#08519C", threshold_red, "grey45"),
  lty = c(1, 2, 3),
  lwd = c(2, 1, 0.8),
  bty = "o",
  bg = "white",
  box.col = "grey70",
  cex = 0.72,
  inset = 0.015,
  y.intersp = 0.9,
  x.intersp = 0.7
)

grDevices::dev.off()
