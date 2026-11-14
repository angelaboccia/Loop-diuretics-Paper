# You need to run helpers.R and 1. Cases selection.R scripts before running this script
# Disproportionality analysis -----
### Indication -----
Indi_diur <- copy(Indi)
indi <- "cardiac failure"
Indi_diur[indi_pt %in% indi_pt_selected, indi_pt := indi]

setDT(Indi_diur)

### Reactions -----
reac1 <- "Cardiac Arrhythmias"
reac2 <- "Acute renal failure"
reac3 <- "Syncope"

Reac_diur <- copy(Reac)
Reac_diur[pt %in% event_arr,  pt := reac1]  # AR
Reac_diur[pt %in% event_ren, pt := reac2]   # REN
Reac_diur[pt %in% event_syn, pt := reac3]   # SYN

### Drugs ----------------------------------------------------------------------
# DPA ------------------------------------------------------------
## Drug manipulation -----
Drug_temp <- copy(Drug_GDMT_source)
Drug_temp[substance %in% diuretics, substance := "Diuretics"]
Drug_temp[substance %in% ace_inhibitors, substance := "ACEi"]
Drug_temp[substance %in% arb,       substance := "ARB"]
Drug_temp[substance %in% arni,      substance := "ARNI"]
Drug_temp[substance %in% beta_blockers, substance := "Beta-blockers"]
Drug_temp[substance %in% mra,       substance := "MRA"]
Drug_temp[substance %in% sglt,      substance := "SGLT-2"]

Drug_temp[primaryid %in% pids_diuretic_ace  & substance %in% c("Diuretics", "ACEi"),          substance := "Diuretics + ACEi"]
Drug_temp[primaryid %in% pids_diuretic_arb  & substance %in% c("Diuretics", "ARB"),           substance := "Diuretics + ARB"]
Drug_temp[primaryid %in% pids_diuretic_arni & substance %in% c("Diuretics", "ARNI"),          substance := "Diuretics + ARNI"]
Drug_temp[primaryid %in% pids_diuretic_beta & substance %in% c("Diuretics", "Beta-blockers"), substance := "Diuretics + Beta-blockers"]
Drug_temp[primaryid %in% pids_diuretic_mra  & substance %in% c("Diuretics", "MRA"),           substance := "Diuretics + MRA"]
Drug_temp[primaryid %in% pids_diuretic_sglt & substance %in% c("Diuretics", "SGLT-2"),        substance := "Diuretics + SGLT-2"]

drug_selected_ace       <- list("Diuretics" = as.list("Diuretics"), "ACEi" = as.list("ACEi"))
drug_selected_ace_diur  <- list("Diuretics + ACEi" = as.list("Diuretics + ACEi"))
drug_selected_arb       <- list("Diuretics" = as.list("Diuretics"), "ARB" = as.list("ARB"))
drug_selected_arb_diur  <- list("Diuretics + ARB" = as.list("Diuretics + ARB"))
drug_selected_arni      <- list("Diuretics" = as.list("Diuretics"), "ARNI" = as.list("ARNI"))
drug_selected_arni_diur <- list("Diuretics + ARNI" = as.list("Diuretics + ARNI"))
drug_selected_beta      <- list("Diuretics" = as.list("Diuretics"), "Beta-blockers" = as.list("Beta-blockers"))
drug_selected_beta_diur <- list("Diuretics + Beta-blockers" = as.list("Diuretics + Beta-blockers"))
drug_selected_mra       <- list("Diuretics" = as.list("Diuretics"), "MRA" = as.list("MRA"))
drug_selected_mra_diur  <- list("Diuretics + MRA" = as.list("Diuretics + MRA"))
drug_selected_sglt      <- list("Diuretics" = as.list("Diuretics"), "SGLT-2" = as.list("SGLT-2"))
drug_selected_sglt_diur <- list("Diuretics + SGLT-2" = as.list("Diuretics + SGLT-2"))

## Reaction manipulation -----
Reac_temp <- copy(Reac_diur)
reac_selected <- list(
  "Cardiac Arrhythmias" = as.list("Cardiac Arrhythmias"),
  "Acute renal failure" = as.list("Acute renal failure"),
  "Syncope" = as.list("Syncope"))

## Restrictions -----
# Isolate single drug classes: exclude reports belonging to >1 class among
# diuretics, ACEi, ARB, ARNI, beta-blockers, MRA and SGLT-2
class_pids <- list(
  "Diuretics"     = pids_diuretics,
  "ACEi"          = pids_ace,
  "ARB"           = pids_arb,
  "ARNI"          = pids_arni,
  "Beta-blockers" = pids_beta,
  "MRA"           = pids_mra,
  "SGLT-2"        = pids_sglt)

class_membership <- rbindlist(lapply(names(class_pids), function(x)
  data.table(primaryid = class_pids[[x]], class = x)))

pids_intersections <- class_membership[, .(n_classes = uniqueN(class)), by = primaryid][
  n_classes > 1, primaryid]

pids_inter      <- setdiff(Demo$primaryid, pids_intersections)
pids_inter_indi <- intersect(pids_inter, pids_indi)

# Diuretic + GDMT combinations: exclude reports containing >1 GDMT class.
# GDMT_class_count is already created in script 1.
pids_more_than_one_class <- unique(GDMT_class_count[n_classes > 1, primaryid])
pids_mono      <- setdiff(Demo$primaryid, pids_more_than_one_class)
pids_mono_indi <- intersect(pids_mono, pids_indi)

restrictions <- list(
  "Full classes" = as.list(Demo$primaryid),
  "HF indication" = as.list(pids_indi),
  "Isolated classes" = as.list(pids_inter),
  "Isolated classes & HF indication" = as.list(pids_inter_indi)
)

restrictions_c <- list(
  "Full classes" = as.list(Demo$primaryid),
  "HF indication" = as.list(pids_indi),
  "One GDMT class" = as.list(pids_mono),
  "One GDMT class & HF indication" = as.list(pids_mono_indi)
)

## DPAs -----
add_dpa_bonferroni <- function(disproportionality_df, index = "IC", group = "nested") {
  
  df <- copy(disproportionality_df)
  lower_col <- paste0(index, "_lower")
  upper_col <- paste0(index, "_upper")
  median_col <- paste0(index, "_median")
  threshold <- ifelse(index == "IC", 0, 1)
  
  if (!all(c(lower_col, upper_col, median_col) %in% names(df))) return(df)
  
  df[, SE_temp := (get(upper_col) - get(lower_col)) / (2 * 1.96)]
  if (group %in% names(df)) {
    df[, alpha_bonf := .05 / .N, by = group]
  } else {
    df[, alpha_bonf := .05 / .N]
  }
  df[, z_bonf := stats::qnorm(1 - alpha_bonf / 2)]
  df[, Bonferroni_lower := get(median_col) - z_bonf * SE_temp]
  df[, Bonferroni_upper := get(median_col) + z_bonf * SE_temp]
  df[, Bonferroni := Bonferroni_lower > threshold]
  df[, c("SE_temp", "alpha_bonf", "z_bonf") := NULL]
  df
}
run_dpa_by_availability <- function(drug_selected, reac_selected, temp_drug, temp_reac,
                                    restrictions, availability_years, demo = Demo) {
  
  out <- data.table()
  demo_year <- copy(demo)[, year := as.numeric(substr(init_fda_dt, 1, 4))]
  
  for (drug_name in names(drug_selected)) {
    for (n in seq_along(restrictions)) {
      t <- restrictions[n]
      t_name <- names(t)
      t_pids <- unlist(t)
      start_year <- availability_years[[drug_name]]
      year_pids <- unique(demo_year[year >= start_year, primaryid])
      
      df <- disproportionality_analysis(
        drug_selected = setNames(list(drug_selected[[drug_name]]), drug_name),
        reac_selected = reac_selected,
        temp_drug = temp_drug,
        temp_reac = temp_reac,
        restriction = intersect(t_pids, year_pids)
      )[, `:=`(nested = t_name, start_year = start_year)]
      
      out <- rbindlist(list(out, df), fill = TRUE)
    }
  }
  
  add_dpa_bonferroni(out, index = "IC", group = "nested")
}
render_forest2 <- function(
    disproportionality_df,
    index = "IC",
    row = "substance",
    levs_row = NA,
    nested = FALSE,
    show_legend = TRUE,
    transformation = "identity",
    custom_threshold = NA,
    text_size_legend = 15,
    text_size_axis = 15,
    text_size_strip = 15,
    dodge = .3,
    nested_colors = NA,
    facet_v = NA,
    facet_h = NA,
    legend_position = "right",
    main = "Forest plot for DPA"
) {
  
  disproportionality_df <- data.table::copy(disproportionality_df)
  
  disproportionality_df[, median := .SD[[1]], .SDcols = paste0(index, "_median")]
  disproportionality_df[, lower  := .SD[[1]], .SDcols = paste0(index, "_lower")]
  disproportionality_df[, upper  := .SD[[1]], .SDcols = paste0(index, "_upper")]
  
  threshold <- if (index == "IC") 0 else 1
  
  if (!is.na(custom_threshold)) {
    threshold <- custom_threshold
  }
  
  if (!identical(nested, FALSE)) {
    
    vals <- disproportionality_df[[nested]]
    
    disproportionality_df[, nested_plot :=
                            factor(vals, levels = unique(vals))
    ]
    
    disproportionality_df[, shape_group :=
                            data.table::fifelse(
                              grepl("HF", as.character(nested_plot)),
                              "HF",
                              "Raw"
                            )
    ]
    
  } else {
    
    disproportionality_df[, nested_plot := factor("Overall")]
    disproportionality_df[, shape_group := "Raw"]
    
  }
  
  if (!"Bonferroni" %in% names(disproportionality_df)) {
    
    disproportionality_df[, Bonferroni := lower > threshold]
    
  }
  
  disproportionality_df[, Significance :=
                          factor(
                            data.table::fifelse(
                              Bonferroni,
                              "Significance",
                              "Not significant"
                            ),
                            levels = c("Significance", "Not significant")
                          )
  ]
  
  nested_levels <- levels(disproportionality_df$nested_plot)
  
  legend_shapes <- ifelse(
    grepl("HF", nested_levels),
    24,
    21
  )
  
  p <- ggplot2::ggplot(
    disproportionality_df,
    ggplot2::aes(
      x = median,
      xmin = lower,
      xmax = upper,
      y = .data[[row]]
    )
  ) +
    
    ggplot2::geom_linerange(
      ggplot2::aes(
        group = nested_plot,
        linetype = nested_plot,
        colour = nested_plot
      ),
      position = ggplot2::position_dodge(0.5),
      linewidth = 1,
      show.legend = FALSE
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        group = nested_plot,
        shape = shape_group,
        fill = nested_plot,
        size = D_E,
        stroke = Significance
      ),
      colour = "black",
      position = ggplot2::position_dodge(0.5)
    ) +
    
    ggplot2::geom_vline(
      xintercept = threshold,
      linetype = "dashed",
      linewidth = 1
    ) +
    
    ggplot2::scale_shape_manual(
      values = c(
        "Raw" = 21,
        "HF" = 24
      )
    ) +
    
    ggplot2::scale_discrete_manual(
      aesthetics = "stroke",
      values = c(
        "Significance" = 1.5,
        "Not significant" = 0
      ),
      breaks = "Significance",
      labels = "Significance"
    ) +
    
    ggplot2::scale_size_continuous(
      range = c(2.5, 7),
      guide = "none"
    ) +
    
    ggplot2::labs(
      x = index,
      y = "",
      colour = NULL,
      fill = NULL,
      shape = NULL,
      linetype = NULL,
      stroke = NULL,
      title = main
    ) +
    
    ggplot2::guides(
      colour = "none",
      shape = "none",
      fill = ggplot2::guide_legend(
        override.aes = list(
          shape = legend_shapes,
          size = 5,
          colour = "black",
          stroke = 0
        )
      ),
      stroke = "none"
    ) +
    
    ggplot2::theme_bw(base_size = 15) +
    
    ggplot2::theme(
      text = ggplot2::element_text(size = 17, face = "bold"),
      axis.text = ggplot2::element_text(size = 17, face = "bold"),
      axis.title = ggplot2::element_text(size = 17, face = "bold"),
      plot.title = ggplot2::element_text(size = 17, face = "bold", hjust = 0.5),
      strip.text = ggplot2::element_text(size = 17, face = "bold"),
      legend.title = ggplot2::element_text(size = 17, face = "bold"),
      legend.text = ggplot2::element_text(size = 15, face = "bold"),
      legend.position = legend_position
    )
  
  if (!is.na(facet_v)) {
    p <- p +
      ggplot2::facet_wrap(
        stats::as.formula(paste("~", facet_v)),
        scales = "fixed"
      )
  }
  
  return(p)
}

# ACE
dpa_df_ace <- run_dpa_by_availability(
  drug_selected_ace, reac_selected, Drug_temp, Reac_temp, restrictions,
  availability_years = c("Diuretics" = 1966, "ACEi" = 1981), demo = Demo)

#path <- paste0(here::here(), "/projects/Diuretici paper/")
write_xlsx(dpa_df_ace, paste0("DPA_ace.xlsx"))
p_ace_mono <- render_forest2(dpa_df_ace, index = "IC", facet_v = "event", nested = "nested")

# ACE + Diur
dpa_df_ace_diur <- run_dpa_by_availability(
  drug_selected_ace_diur, reac_selected, Drug_temp, Reac_temp, restrictions_c,
  availability_years = c("Diuretics + ACEi" = 1981), demo = Demo)

write_xlsx(dpa_df_ace_diur, paste0("DPA_ace_diur.xlsx"))
p_ace_diur <- render_forest2(dpa_df_ace_diur, index = "IC", facet_v = "event", nested = "nested", main = "")
p_ace <- p_ace_mono / p_ace_diur

# ARNi
dpa_df_arni <- run_dpa_by_availability(
  drug_selected_arni, reac_selected, Drug_temp, Reac_temp, restrictions,
  availability_years = c("Diuretics" = 1966, "ARNI" = 2015), demo = Demo)

write_xlsx(dpa_df_arni, paste0("DPA_arni.xlsx"))
p_arni_mono <- render_forest2(dpa_df_arni, index = "IC", facet_v = "event", nested = "nested")

# ARNi + Diur
dpa_df_arni_diur <- run_dpa_by_availability(
  drug_selected_arni_diur, reac_selected, Drug_temp, Reac_temp, restrictions_c,
  availability_years = c("Diuretics + ARNI" = 2015), demo = Demo)

write_xlsx(dpa_df_arni_diur, paste0("PA_arni_diur.xlsx"))
p_arni_diur <- render_forest2(dpa_df_arni_diur, index = "IC", facet_v = "event", nested = "nested", main = "")
p_arni <- p_arni_mono / p_arni_diur

# ARB
dpa_df_arb <- run_dpa_by_availability(
  drug_selected_arb, reac_selected, Drug_temp, Reac_temp, restrictions,
  availability_years = c("Diuretics" = 1966, "ARB" = 1996), demo = Demo)

write_xlsx(dpa_df_arb, paste0("DPA_arb.xlsx"))
p_arb_mono <- render_forest2(dpa_df_arb, index = "IC", facet_v = "event", nested = "nested")

# ARB + Diur
dpa_df_arb_diur <- run_dpa_by_availability(
  drug_selected_arb_diur, reac_selected, Drug_temp, Reac_temp, restrictions_c,
  availability_years = c("Diuretics + ARB" = 1996), demo = Demo)

write_xlsx(dpa_df_arb_diur, paste0("DPA_arb_diur.xlsx"))
p_arb_diur <- render_forest2(dpa_df_arb_diur, index = "IC", facet_v = "event", nested = "nested", main = "")
p_arb <- p_arb_mono / p_arb_diur

# Beta-blockers
dpa_df_beta <- run_dpa_by_availability(
  drug_selected_beta, reac_selected, Drug_temp, Reac_temp, restrictions,
  availability_years = c("Diuretics" = 1966, "Beta-blockers" = 1978), demo = Demo)

write_xlsx(dpa_df_beta, paste0("DPA_beta.xlsx"))
p_beta_mono <- render_forest2(dpa_df_beta, index = "IC", facet_v = "event", nested = "nested")

# Beta-blockers + Diur
dpa_df_beta_diur <- run_dpa_by_availability(
  drug_selected_beta_diur, reac_selected, Drug_temp, Reac_temp, restrictions_c,
  availability_years = c("Diuretics + Beta-blockers" = 1978), demo = Demo)

write_xlsx(dpa_df_beta_diur, paste0("DPA_beta_diur.xlsx"))
p_beta_diur <- render_forest2(dpa_df_beta_diur, index = "IC", facet_v = "event", nested = "nested", main = "")
p_beta <- p_beta_mono / p_beta_diur

# MRA
dpa_df_mra <- run_dpa_by_availability(
  drug_selected_mra, reac_selected, Drug_temp, Reac_temp, restrictions,
  availability_years = c("Diuretics" = 1966, "MRA" = 1960), demo = Demo)

write_xlsx(dpa_df_mra, paste0("DPA_mra.xlsx"))
p_mra_mono <- render_forest2(dpa_df_mra, index = "IC", facet_v = "event", nested = "nested")

# MRA + Diur
dpa_df_mra_diur <- run_dpa_by_availability(
  drug_selected_mra_diur, reac_selected, Drug_temp, Reac_temp, restrictions_c,
  availability_years = c("Diuretics + MRA" = 1966), demo = Demo)

write_xlsx(dpa_df_mra_diur, paste0("DPA_mra_diur.xlsx"))
p_mra_diur <- render_forest2(dpa_df_mra_diur, index = "IC", facet_v = "event", nested = "nested", main = "")
p_mra <- p_mra_mono / p_mra_diur

# SGLT-2
dpa_df_sglt <- run_dpa_by_availability(
  drug_selected_sglt, reac_selected, Drug_temp, Reac_temp, restrictions,
  availability_years = c("Diuretics" = 1966, "SGLT-2" = 2013), demo = Demo)

write_xlsx(dpa_df_sglt, paste0("DPA_sglt.xlsx"))
p_sglt_mono <- render_forest2(dpa_df_sglt, index = "IC", facet_v = "event", nested = "nested")

# SGLT-2 + Diur
dpa_df_sglt_diur <- run_dpa_by_availability(
  drug_selected_sglt_diur, reac_selected, Drug_temp, Reac_temp, restrictions_c,
  availability_years = c("Diuretics + SGLT-2" = 2013), demo = Demo)

write_xlsx(dpa_df_sglt_diur, paste0("DPA_sglt_diur.xlsx"))
p_sglt_diur <- render_forest2(dpa_df_sglt_diur, index = "IC", facet_v = "event", nested = "nested", main = "")
p_sglt <- p_sglt_mono / p_sglt_diur

# All df in one .xlsx file, one per sheet
all_dpa_dfs <- list("ACEi mono" = dpa_df_ace, "ACEi + Diuretics" = dpa_df_ace_diur,
                    "ARB mono" = dpa_df_arb, "ARB + Diuretics" = dpa_df_arb_diur,
                    "ARNI mono" = dpa_df_arni, "ARNI + Diuretics" = dpa_df_arni_diur,
                    "Beta-blockers mono" = dpa_df_beta, "Beta-blockers + Diuretics" = dpa_df_beta_diur,
                    "MRA mono" = dpa_df_mra, "MRA + Diuretics" = dpa_df_mra_diur,
                    "SGLT-2 mono" = dpa_df_sglt, "SGLT-2 + Diuretics" = dpa_df_sglt_diur)

write_xlsx(all_dpa_dfs, paste0("DPA_all.xlsx"))

## Save plots -----
dir.create("DPAs", recursive = TRUE, showWarnings = FALSE)

plots <- list(ace = p_ace, arb = p_arb, arni = p_arni, beta = p_beta, mra = p_mra, sglt = p_sglt)
invisible(mapply(function(p, n) ggsave(paste0("DPA_", n, ".png"), 
                                       p, width = 12, height = 10, dpi = 300), plots, names(plots)))
