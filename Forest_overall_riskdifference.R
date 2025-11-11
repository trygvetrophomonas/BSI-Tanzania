library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)

df <- read_excel("~/Documents/R/Project Tanzania/master.xlsx")

# Recode species
df$org <- NA
df$org[df$orgsimp == "01.klebs"] <- "Klebsiella pneumoniae"
df$org[df$orgsimp == "02.salm"] <- "Salmonella enterica"
df$org[df$orgsimp == "03.ecoli"] <- "Escherichia coli"
df$org[df$orgsimp %in% c("04.entales", "08.non-ent")] <- "Other Gram-negatives"
df$org[df$orgsimp == "07.pseudo"] <- "Pseudomonas spp."
df$org[df$orgsimp == "06.acinet"] <- "Acinetobacter spp."
df$org[df$orgsimp == "10.saureus"] <- "Staphylococcus aureus"
df$org[df$orgsimp == "11.entcoc"] <- "Enterococcus spp."
df$org[df$orgsimp == "12.strep"] <- "Streptococcus spp."

# per patient ----
df_patient <- df %>%
  group_by(id) %>%
  summarise(
    org    = ifelse(n_distinct(as.character(org)) > 1, "Polymicrobial", as.character(first(org))),
    AMP_CN  = any(AMP == 0 & CN == 0),
    AMP_AK  = any(AMP == 0 & AK == 0),
    CR      = any(CR == 0)
  ) %>%
  mutate(
    reg_ampigenta = as.integer(AMP_CN),
    reg_ampiamik  = as.integer(AMP_AK),
    reg_cro       = as.integer(CR)
  )

# Function forestplot-data ----
make_results <- function(data, var1, var2) {
  results <- data %>%
    group_by(org) %>%
    summarise(
      n = n(),
      n1 = sum(.data[[var1]] == 1, na.rm = TRUE),
      n2 = sum(.data[[var2]] == 1, na.rm = TRUE),
      p1 = mean(.data[[var1]], na.rm = TRUE),
      p2 = mean(.data[[var2]], na.rm = TRUE),
      mean = p1 - p2,
      SE = sqrt((p1 * (1 - p1) / n) + (p2 * (1 - p2) / n)),
      lower = mean - 1.96 * SE,
      upper = mean + 1.96 * SE,
      .groups = "drop"
    ) %>%
    mutate(
      pct1 = paste0(round(p1 * 100), "% (n = ", n1, ")"),
      pct2 = paste0(round(p2 * 100), "% (n = ", n2, ")"),
      diff95ci = paste0(round(mean * 100), "% (", round(lower * 100), "% to ", round(upper * 100), "%)")
    )
  
  summary <- data %>%
    summarise(
      org = "Overall",
      n = n(),
      n1 = sum(.data[[var1]] == 1, na.rm = TRUE),
      n2 = sum(.data[[var2]] == 1, na.rm = TRUE),
      p1 = mean(.data[[var1]], na.rm = TRUE),
      p2 = mean(.data[[var2]], na.rm = TRUE),
      mean = p1 - p2,
      SE = sqrt((p1 * (1 - p1) / n) + (p2 * (1 - p2) / n)),
      lower = mean - 1.96 * SE,
      upper = mean + 1.96 * SE
    ) %>%
    mutate(
      pct1 = paste0(round(p1 * 100), "% (n = ", n1, ")"),
      pct2 = paste0(round(p2 * 100), "% (n = ", n2, ")"),
      diff95ci = paste0(round(mean * 100), "% (", round(lower * 100), "% to ", round(upper * 100), "%)")
    )
  
  bind_rows(results, summary) %>%
    mutate(org = ifelse(org == "Overall", paste0("Overall (n = ", n, ")"), paste0(org, " (n = ", n, ")")))
}

# Create dataset both regimens
results_amik_genta <- make_results(df_patient, "reg_ampiamik", "reg_ampigenta") %>%
  mutate(regime = "Amikacin vs Gentamicin")
results_amik_cro   <- make_results(df_patient, "reg_ampiamik", "reg_cro") %>%
  mutate(regime = "Amikacin vs Ceftriaxone")

combined <- bind_rows(results_amik_genta, results_amik_cro)

# polygon-romb
make_diamond <- function(df) {
  df <- df %>% filter(grepl("^Overall", org))
  if (nrow(df) == 0) return(NULL)
  
  y <- as.numeric(factor(df$org, levels = levels(df$org)))
  data.frame(
    x = c(df$lower, df$mean, df$upper, df$mean) * 100,
    y = c(y, y + 0.3, y, y - 0.3),
    group = df$org
  )
}

# Function for panel per regimen
make_panel <- function(data, regime_title) {
  data <- data %>% filter(regime == regime_title)
  
  desired_order <- c(
    "Klebsiella pneumoniae", "Escherichia coli", "Salmonella enterica", "Other Gram-negatives", 
    "Pseudomonas spp.", "Acinetobacter spp.", "Staphylococcus aureus", 
    "Enterococcus spp.", "Streptococcus spp.", "Polymicrobial", "Overall"
  )
  
  actual_levels <- map_chr(desired_order, function(sp) {
    match <- grep(paste0("^", sp, " \\(n = "), data$org, value = TRUE)
    if (length(match) == 0) NA_character_ else match
  })
  
  full_levels <- rev(na.omit(actual_levels))
  data$org <- factor(data$org, levels = full_levels)
  
  diamond_df <- make_diamond(data)
  
  table_plot <- ggplot(data, aes(y = org)) +
    geom_text(aes(x = 0.4, label = org), hjust = 0, size = 4, fontface = "bold") +
    geom_text(aes(x = 2.1, label = pct2), hjust = 0, size = 4) +  
    geom_text(aes(x = 3.2, label = pct1), hjust = 0, size = 4) +
    geom_text(aes(x = 4.2, label = diff95ci), hjust = 0, size = 4) +
    scale_x_continuous(limits = c(0.3, 5.2), breaks = NULL) +
    theme_void() +
    theme(plot.margin = margin(5, 10, 5, 5))
  
  graph_plot <- ggplot(data, aes(x = mean * 100, y = org)) +
    geom_point(aes(size = n, alpha = ifelse(grepl("^Overall", org), 0, 1)),
               shape = 15, color = "royalblue", fill = "royalblue") +
    geom_errorbarh(aes(xmin = lower * 100, xmax = upper * 100,
                       alpha = ifelse(grepl("^Overall", org), 0, 1)),
                   height = 0.2, color = "royalblue") +
    geom_polygon(data = diamond_df, aes(x = x, y = y, group = group),
                 fill = "royalblue", color = "royalblue", inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    scale_shape_identity() +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    scale_alpha_identity() +
    scale_x_continuous(breaks = c(-100, -75, -50, -25, 0, 25), limits = c(-100, 40)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(5, 5, 5, 10),
      legend.position = "none"
    )
  
  table_plot + graph_plot + plot_layout(widths = c(3, 2))
}

# Panel both regimens
panel1 <- make_panel(combined, "Amikacin vs Gentamicin")
panel2 <- make_panel(combined, "Amikacin vs Ceftriaxone")

# headers
make_headers <- function(label1, label2) {
  format_label <- function(label) {
    if (label == "Ceftriaxone") label else paste0(label, "-Ampicillin")
  }
  
  ggplot() +
    annotate("text", x = 0.3, y = 1.0, label = "Species (n)", hjust = 0, size = 4, fontface = "bold") +
    annotate("text", x = 1.25, y = 1.0, label = paste0(format_label(label2), "\nResistance (%)"), hjust = 0, size = 4, fontface = "bold") +
    annotate("text", x = 1.9, y = 1.0, label = paste0(format_label(label1), "\nResistance (%)"), hjust = 0, size = 4, fontface = "bold") +
    annotate("text", x = 2.5, y = 1.0, label = "Resistance difference\n(95% CI)", hjust = 0, size = 4, fontface = "bold") +
    xlim(0.3, 5.2) + ylim(0.5, 1.8) +
    theme_void()
  
}

# separate headers
headers1 <- make_headers("Amikacin", "Gentamicin")
headers2 <- make_headers("Amikacin", "Ceftriaxone")

# Footer
make_footer <- function() {
  ggplot() +
    annotate("text", x = 0.8, y = 0.5, label = "Favours Amikacin–Ampicillin", hjust = 0.5, size = 4, fontface = "bold") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

footer1 <- make_footer()
footer2 <- make_footer()

# Combine everything
final_plot <- (
  headers1 / panel1 / footer1 /
    headers2 / panel2 / footer2
) + plot_layout(heights = c(0.12, 1, 0.05, 0.12, 1, 0.05))

print(final_plot)

