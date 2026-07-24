# resazurin-family-diagnostics.R
# One diagnostic plot per family showing raw curves with annotated curve
# features (peak, initial slope, AUC, trajectory class) alongside feature
# distributions across individuals. Useful for outlier identification and
# ground-truthing feature extraction.

library(tidyverse)
library(patchwork)
library(here)

# ── paths ──────────────────────────────────────────────────────────────────────
metabolism_path    <- here("output", "metabolism.csv")
curve_feature_path <- here("output", "resazurin-curves", "curve_features.csv")
out_dir            <- here("output", "family-diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── load data ──────────────────────────────────────────────────────────────────
metabolism <- read_csv(metabolism_path, show_col_types = FALSE)
curve_features <- read_csv(curve_feature_path, show_col_types = FALSE) %>%
  filter(value_metric == "corrected_fc")   # one metric only

# ── colour palette for trajectory classes ──────────────────────────────────────
traj_cols <- c(
  "monotonic_increase"   = "#009E73",
  "rise_then_depress"    = "#0072B2",
  "depress_then_recover" = "#E69F00",
  "stable"               = "#999999",
  "complex"              = "#CC79A7"
)

# features shown in the distribution panel, in display order
DIST_FEATURES <- c(
  "auc_total", "initial_slope", "peak_value",
  "metabolic_scope", "metabolic_depression_index",
  "resilience_ratio"
)
DIST_LABELS <- c(
  auc_total                = "AUC total",
  initial_slope            = "Init. slope",
  peak_value               = "Peak value",
  metabolic_scope          = "Met. scope",
  metabolic_depression_index = "Depress. idx",
  resilience_ratio         = "Resilience"
)

# ── helper: trapezoid AUC ─────────────────────────────────────────────────────
trap_auc <- function(t, v) {
  ord <- order(t)
  t <- t[ord]; v <- v[ord]
  sum(diff(t) * (head(v, -1) + tail(v, -1)) / 2)
}

# ── plot function ──────────────────────────────────────────────────────────────
plot_family <- function(fam_id) {

  # --- data for this family --------------------------------------------------
  raw <- metabolism %>%
    filter(family == fam_id) %>%
    left_join(
      curve_features %>% select(trace_id, trajectory_class),
      by = c("unique" = "trace_id")
    )

  feats <- curve_features %>%
    filter(family_id_group == as.character(fam_id))

  if (nrow(raw) == 0 || nrow(feats) == 0) {
    message("Skipping family ", fam_id, " — no data")
    return(invisible(NULL))
  }

  phenotype_grp <- feats$phenotype[1]
  n_indiv       <- n_distinct(raw$unique)

  # mean curve
  mean_curve <- raw %>%
    group_by(timepoint) %>%
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(timepoint)

  # feature annotations on the mean curve
  peak_idx <- which.max(mean_curve$mean_val)
  peak_t   <- mean_curve$timepoint[peak_idx]
  peak_v   <- mean_curve$mean_val[peak_idx]

  v_t0 <- mean_curve$mean_val[mean_curve$timepoint == min(mean_curve$timepoint)]
  v_t1_row <- mean_curve %>% filter(timepoint == 1)
  v_t1 <- if (nrow(v_t1_row) > 0) v_t1_row$mean_val[1] else NA_real_

  auc_mean <- trap_auc(mean_curve$timepoint, mean_curve$mean_val)

  # per-individual feature points: peak and trough
  feat_pts <- feats %>%
    select(trace_id, time_to_peak, peak_value,
           time_to_trough, trough_value,
           time_to_min_slope, trajectory_class) %>%
    left_join(
      curve_features %>% select(trace_id, trajectory_class) %>% distinct(),
      by = "trace_id"
    ) %>%
    select(-trajectory_class.y) %>%
    rename(trajectory_class = trajectory_class.x)

  # summary stats for subtitle
  mean_auc   <- round(mean(feats$auc_total,    na.rm = TRUE), 1)
  mean_slope <- round(mean(feats$initial_slope, na.rm = TRUE), 2)
  mean_peak  <- round(mean(feats$peak_value,    na.rm = TRUE), 2)
  mean_scope <- round(mean(feats$metabolic_scope, na.rm = TRUE), 2)

  traj_tbl <- feats %>%
    count(trajectory_class) %>%
    mutate(pct = round(100 * n / sum(n))) %>%
    arrange(desc(n))
  traj_text <- paste(
    paste0(traj_tbl$trajectory_class, " (", traj_tbl$pct, "%)"),
    collapse = " | "
  )

  y_max <- max(raw$value, na.rm = TRUE) * 1.05

  # ── PANEL A: curves + feature annotations ──────────────────────────────────
  p_curves <- ggplot() +
    # AUC shading under the mean curve
    geom_area(
      data = mean_curve,
      aes(x = timepoint, y = mean_val),
      fill = "#0072B2", alpha = 0.12
    ) +
    # individual traces, coloured by trajectory class
    geom_line(
      data = raw,
      aes(x = timepoint, y = value,
          group = unique, colour = trajectory_class),
      alpha = 0.45, linewidth = 0.45
    ) +
    # individual peak dots
    geom_point(
      data = feat_pts,
      aes(x = time_to_peak, y = peak_value),
      colour = "#D55E00", shape = 21,
      fill = "#D55E00", size = 1.8, alpha = 0.55
    ) +
    # initial slope segment on mean (t_start → t=1)
    {if (!is.na(v_t1))
      geom_segment(
        aes(x = min(mean_curve$timepoint), xend = 1,
            y = v_t0, yend = v_t1),
        colour = "#009E73", linewidth = 2, linetype = "dashed"
      )
    } +
    # mean curve (bold black)
    geom_line(
      data = mean_curve,
      aes(x = timepoint, y = mean_val),
      colour = "black", linewidth = 1.4
    ) +
    # mean peak (large orange dot)
    geom_point(
      aes(x = peak_t, y = peak_v),
      colour = "#D55E00", size = 4.5
    ) +
    # mean peak label
    annotate(
      "text", x = peak_t, y = peak_v * 1.08,
      label = sprintf("peak = %.2f\nt = %g h", peak_v, peak_t),
      size = 3, colour = "#D55E00", hjust = 0.5, lineheight = 0.9
    ) +
    # initial slope label
    {if (!is.na(v_t1))
      annotate(
        "text", x = 0.5,
        y = (v_t0 + v_t1) / 2 - y_max * 0.06,
        label = sprintf("slope = %.2f", mean_slope),
        size = 3, colour = "#009E73", hjust = 0.5
      )
    } +
    # AUC label
    annotate(
      "text", x = max(mean_curve$timepoint) * 0.55,
      y = mean_curve$mean_val[which.max(mean_curve$mean_val)] * 0.3,
      label = sprintf("AUC = %.1f", mean_auc),
      size = 3, colour = "#0072B2", alpha = 0.8
    ) +
    scale_colour_manual(
      values = traj_cols, name = "Trajectory class",
      na.value = "grey60",
      guide = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1))
    ) +
    scale_x_continuous(breaks = 0:4) +
    coord_cartesian(ylim = c(NA, y_max)) +
    labs(
      title    = sprintf("Family %s  ·  %s  ·  n = %d",
                         fam_id, phenotype_grp, n_indiv),
      subtitle = sprintf(
        "Mean: AUC = %s  ·  Slope = %s  ·  Peak = %s  ·  Scope = %s\nTrajectories: %s",
        mean_auc, mean_slope, mean_peak, mean_scope, traj_text
      ),
      x = "Timepoint (h)",
      y = "Resazurin fold change"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 8, colour = "grey30"),
      legend.position  = "bottom",
      legend.text      = element_text(size = 8),
      legend.title     = element_text(size = 8.5),
      legend.key.width = unit(1.2, "cm")
    )

  # ── PANEL B: feature distributions per individual ──────────────────────────
  feat_long <- feats %>%
    select(trace_id, all_of(DIST_FEATURES)) %>%
    pivot_longer(-trace_id, names_to = "feature", values_to = "val") %>%
    filter(is.finite(val)) %>%
    mutate(
      feature = factor(feature, levels = DIST_FEATURES),
      feature_label = DIST_LABELS[as.character(feature)]
    )

  p_dist <- ggplot(feat_long, aes(x = val, y = 0)) +
    geom_boxplot(
      aes(x = val),
      orientation = "y",
      outlier.shape = NA,
      fill = "grey92", colour = "grey50",
      width = 0.55
    ) +
    geom_jitter(
      aes(x = val), height = 0.18,
      colour = "#0072B2", alpha = 0.65, size = 1.6
    ) +
    facet_wrap(
      ~ feature_label, ncol = 1,
      scales = "free_x", strip.position = "left"
    ) +
    labs(x = "Value", y = NULL, title = "Feature distributions") +
    theme_classic(base_size = 9) +
    theme(
      strip.text          = element_text(size = 8, hjust = 1),
      strip.placement     = "outside",
      strip.background    = element_blank(),
      strip.clip          = "off",
      axis.text.y         = element_blank(),
      axis.ticks.y        = element_blank(),
      panel.spacing.y     = unit(0.3, "lines"),
      plot.title          = element_text(size = 10, face = "bold")
    )

  # ── combine ────────────────────────────────────────────────────────────────
  p_combined <- p_curves + p_dist +
    plot_layout(widths = c(2.8, 1))

  fname <- sprintf(
    "family_%03d_%s.png",
    as.integer(fam_id),
    gsub("-", "_", phenotype_grp)
  )
  ggsave(
    file.path(out_dir, fname),
    p_combined, width = 12, height = 6, dpi = 150
  )
  message("  saved: ", fname)
}

# ── run for all families ───────────────────────────────────────────────────────
families <- sort(unique(metabolism$family))
message("Generating diagnostic plots for ", length(families), " families...")
walk(families, plot_family)
message("\nDone. Plots written to: ", out_dir)
