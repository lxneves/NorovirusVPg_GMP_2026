# Dot score
# --- Score SD error bars ---
etd + 
  geom_errorbar(
    data = score_etd %>%
      mutate(
        ymin = to_y(pmax(mean_score - sd_score, 0)),
        ymax = to_y(pmin(mean_score + sd_score, 300))
      ),
    aes(
      x = `fragmentation energy (ETD)`,
      ymin = ymin,
      ymax = ymax,
      group = 1
    ),
    inherit.aes = FALSE,
    width = 3,
    color = "grey20",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  
  # --- Score mean dot ---
  geom_point(
    data = score_etd %>% mutate(y = to_y(mean_score)),
    aes(
      x = `fragmentation energy (ETD)`,
      y = y
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = "grey20",
    size = 3,
    stroke = 0.8,
    show.legend = FALSE
  ) +
  
  # --- Score mean label above dot ---
  geom_text(
    data = score_etd %>% mutate(y = to_y(mean_score)),
    aes(
      x = `fragmentation energy (ETD)`,
      y = y,
      label = round(mean_score, 1)
    ),
    inherit.aes = FALSE,
    vjust = -1.0,
    size = 3.1,
    color = "grey20",
    show.legend = FALSE
  )


etd +
  # --- PSM score SD ribbon (mapped to left axis, no legend) ---
  geom_ribbon(
    data = score_etd %>%
      mutate(
        ymin = to_y(pmax(mean_score - sd_score, 0)),
        ymax = to_y(pmin(mean_score + sd_score, 250))
      ),
    aes(x = `fragmentation energy (ETD)`, ymin = ymin, ymax = ymax, group = 1),
    inherit.aes = FALSE, fill = "orange", alpha = 0.1, linewidth = 0,
    show.legend = FALSE
  ) +
  
  # --- PSM score mean (dashed; no legend) ---
  geom_line(
    data = score_etd %>% mutate(y = to_y(pmin(pmax(mean_score, 0), 250))),
    aes(x = `fragmentation energy (ETD)`, y = y, group = 1),
    inherit.aes = FALSE, color = "red4", linewidth = 1, linetype = "dashed",
    show.legend = FALSE
  ) +
  
  # --- Secondary axis with 0..300 labeling ---
  scale_y_continuous(
    name = "Log10-average intensity",
    sec.axis = sec_axis(~ to_score(.), name = "PSM score")
  )
