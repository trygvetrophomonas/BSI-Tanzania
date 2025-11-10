library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

klebsiella <- read_excel("~/Documents/R/Project Tanzania/klebsiella.xlsx")
ecoli <- read_excel("~/Documents/R/Project Tanzania/ecoli.xlsx")
salmonella <- read_excel("~/Documents/R/Project Tanzania/salmonella.xlsx")

ecoli$organism <- "E.coli"
klebsiella$organism <- "Klebsiella"
salmonella$organism <- "Salmonella"

all <- bind_rows(ecoli, salmonella, klebsiella) %>%
  select(ID, organism, everything()) %>%
  mutate(study = ifelse(str_detect(ID, "^BSI"), "Study 1", "Study 2"))

has_gene <- function(x, pattern) {
  if (all(is.na(x))) return(FALSE)
  genes <- unlist(str_split(paste(na.omit(x), collapse = " "), "\\s|,"))
  genes <- genes[genes != ""]
  any(str_detect(genes, pattern))
}

gen_cols <- all %>%
  select(-ID, -organism, -study, -Penicillinase, -ESBL) %>%
  select(where(is.character)) %>%
  names()

# Counts per cathegory
all <- all %>%
  rowwise() %>%
  mutate(
    Penicillinase = !is.na(Penicillinase),
    ESBL = !is.na(ESBL),
    `aac(6')-Ib-cr` = has_gene(c_across(all_of(gen_cols)), "aac\\(6'\\)-Ib-cr"),
    `aac(3)-II` = has_gene(c_across(all_of(gen_cols)), "aac\\(3\\)-II"),
    QRDR = has_gene(c_across(all_of(gen_cols)), "^(gyrA|gyrB|parC|parE|marR)")
  ) %>%
  ungroup()

# === Function for dumbbell-plot ===
make_dumbbell_plot <- function(df, organism_name) {
  summary <- df %>%
    group_by(study) %>%
    summarise(across(c(Penicillinase, ESBL, `aac(6')-Ib-cr`, `aac(3)-II`, QRDR),
                     ~ mean(.x) * 100), .groups = "drop")
  
  long <- summary %>%
    pivot_longer(cols = -study, names_to = "category", values_to = "percent")
  
  # Unicode for β
  long <- long %>%
    mutate(category = case_when(
      category == "Penicillinase" ~ "Narrow spectrum\n\u03B2-lactamase",
      TRUE ~ category
    ))
  
  
  category_order <- c("Narrow spectrum\n\u03B2-lactamase", "ESBL", "aac(3)-II", "aac(6')-Ib-cr", "QRDR")
  
  wide <- long %>%
    pivot_wider(names_from = study, values_from = percent) %>%
    mutate(category = factor(category, levels = rev(category_order))) %>%
    arrange(desc(`Study 2`))
  
  plotdata <- wide %>%
    pivot_longer(cols = c(`Study 1`, `Study 2`), names_to = "study", values_to = "percent")
  
  ggplot(plotdata, aes(x = percent, y = category, color = study)) +
    geom_segment(data = wide,
                 aes(x = `Study 1`, xend = `Study 2`, y = category, yend = category),
                 color = "gray70", size = 1.5, inherit.aes = FALSE) +
    geom_point(size = 4) +
    scale_color_manual(values = c("Study 1" = "#0e668b", "Study 2" = "#c64737")) +
    labs(
      title = organism_name,
      x = "Percentage of Isolates (%)",
      y = "Gene Category",
      color = "Study"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

p_ecoli <- make_dumbbell_plot(filter(all, organism == "E.coli"), "Escherichia coli")
p_kleb <- make_dumbbell_plot(filter(all, organism == "Klebsiella"), "Klebsiella pneumoniae")
p_salm <- make_dumbbell_plot(filter(all, organism == "Salmonella"), "Salmonella enterica")
p_all  <- make_dumbbell_plot(all, "All organisms")

combined_plot <- (p_ecoli | p_kleb) / 
  plot_spacer() / 
  (p_salm | p_all) +
  plot_layout(heights = c(1, 0.1, 1)) +
  plot_annotation(
    title = "Resistance Genes",
    theme = theme(plot.title = element_text(hjust = 0.5))
  )

print(combined_plot)

