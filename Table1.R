library(tidyverse)
library(readxl)
library(forestplot)
library(grid)
library(dplyr)
library(tidyr)
library(gtsummary)
library(gt)
library(stringr)

df <- read_excel("~/Documents/R/Project Tanzania/master.xlsx")

# Recode Species
df$org <- NA
df$org[df$orgsimp == "01.klebs"]       <- "Klebsiella"
df$org[df$orgsimp == "02.salm"]        <- "Salmonella"
df$org[df$orgsimp == "03.ecoli"]       <- "E. coli"
df$org[df$orgsimp == "04.entales"]     <- "Other Gram negatives*"
df$org[df$orgsimp == "08.non-ent"]     <- "Other Gram negatives*"
df$org[df$orgsimp == "07.pseudo"]      <- "Pseudomonas"
df$org[df$orgsimp == "06.acinet"]      <- "Acinetobacter"
df$org[df$orgsimp == "10.saureus"]     <- "S. aureus"
df$org[df$orgsimp == "11.entcoc"]      <- "Enterococcus"
df$org[df$orgsimp == "12.strep"]       <- "Streptococcus"
df$org <- factor(df$org, levels = c("Klebsiella", "Salmonella", "E. coli", "Other Gram negatives*", "Pseudomonas", "Acinetobacter", "S. aureus", "Enterococcus", "Streptococcus"))

# Patient level
df_patient <- df %>%
  group_by(id) %>%
  summarise(
    ageadm = first(ageadm),
    neo    = first(neo),
    died   = first(died),       
    outcome = first(outcome),  
    
    org    = ifelse(n_distinct(as.character(org)) > 1, "Polymicrobial", as.character(first(org))),
    
    # Combination regimens: resistant if one bacteria is resistant to all drugs in regimen
    reg_genta_amp         = any(CN == 0 & AMP == 0),
    reg_genta_amp_cloxa   = any(CN == 0 & AMP == 0 & Cloxa == 0),
    reg_amik_amp          = any(AK == 0 & AMP == 0),
    reg_plazo_amp         = any(PLZ == 0 & AMP == 0),
    reg_amik_cro          = any(AK == 0 & CR == 0),
    reg_cro              = any(CR == 0),
    reg_cip              = any(CIP == 0)
  ) %>%
  mutate(across(starts_with("reg_"), ~if_else(.x, 1, 0)))  # 1 = resistent, 0 = følsom

df_patient$org <- factor(df_patient$org, levels = c("Klebsiella", "Salmonella", "E. coli", "Other Gram negatives*", "Pseudomonas", "Acinetobacter", "S. aureus", "Enterococcus", "Streptococcus", "Polymicrobial"))

# Add N counts and dynamic labels
org_counts <- df_patient %>% filter(!is.na(org)) %>% count(org, name = "N")

data_labeled <- df_patient %>%
  left_join(org_counts, by = "org") %>%
  mutate(
    org_labeled = paste0(org, " (N=", N, ")"),
    org_labeled = factor(org_labeled, levels = c(
      paste0("Klebsiella (N=",     org_counts$N[org_counts$org == "Klebsiella"], ")"),
      paste0("E. coli (N=",        org_counts$N[org_counts$org == "E. coli"], ")"),
      paste0("Salmonella (N=",     org_counts$N[org_counts$org == "Salmonella"], ")"),
      paste0("Other Gram negatives* (N=", org_counts$N[org_counts$org == "Other Gram negatives*"], ")"),
      paste0("Pseudomonas (N=",    org_counts$N[org_counts$org == "Pseudomonas"], ")"),
      paste0("Acinetobacter (N=",  org_counts$N[org_counts$org == "Acinetobacter"], ")"),
      paste0("S. aureus (N=",      org_counts$N[org_counts$org == "S. aureus"], ")"),
      paste0("Enterococcus (N=",   org_counts$N[org_counts$org == "Enterococcus"], ")"),
      paste0("Streptococcus (N=",  org_counts$N[org_counts$org == "Streptococcus"], ")"),
      paste0("Polymicrobial (N=",  org_counts$N[org_counts$org == "Polymicrobial"], ")")
    ))
  )

# Summarize resistance
res_summary <- data_labeled %>%
  select(
    org = org_labeled,
    reg_genta_amp, reg_genta_amp_cloxa,
    reg_cro, reg_cip,
    reg_amik_amp, reg_plazo_amp, reg_amik_cro
  ) %>%
  pivot_longer(
    cols = starts_with("reg_"),
    names_to = "regimen",
    values_to = "resistant"
  ) %>%
  mutate(
    regimen = recode(regimen,
                     reg_genta_amp        = "Gentamicin + Ampicillin",
                     reg_genta_amp_cloxa  = "Gentamicin + Ampicillin + Cloxacillin",
                     reg_cro              = "Ceftriaxone",
                     reg_cip              = "Ciprofloxacin",
                     reg_amik_amp         = "Amikacin + Ampicillin",
                     reg_plazo_amp        = "Plazomicin + Ampicillin",
                     reg_amik_cro         = "Amikacin + Ceftriaxone"
    )
  ) %>%
  filter(!is.na(resistant)) %>%
  group_by(org, regimen) %>%
  summarise(
    resistance_pct = mean(resistant == 1) * 100,
    resistance_n = sum(resistant == 1),
    .groups = "drop"
  ) %>%
  mutate(
    resistance_label = sprintf("%d%% (n=%d)", round(resistance_pct), resistance_n)
  ) %>%
  select(-resistance_pct, -resistance_n) %>%
  pivot_wider(names_from = regimen, values_from = resistance_label) %>%
  mutate(group = case_when(
    str_detect(as.character(org), "Klebsiella|E. coli|Salmonella|Other Gram negatives*|Pseudomonas|Acinetobacter") ~ "Gram-negative bacteria",
    str_detect(as.character(org), "S. aureus|Enterococcus|Streptococcus") ~ "Gram-positive bacteria",
    TRUE ~ "Polymicrobial"
  ))

# Create total row across all organisms
total_N <- nrow(data_labeled)

res_totals <- data_labeled %>%
  pivot_longer(cols = starts_with("reg_"),
               names_to = "regimen", values_to = "resistant",
               values_drop_na = TRUE) %>%
  mutate(
    regimen = recode(regimen,
                     reg_genta_amp        = "Gentamicin + Ampicillin",
                     reg_genta_amp_cloxa  = "Gentamicin + Ampicillin + Cloxacillin",
                     reg_cro              = "Ceftriaxone",
                     reg_cip              = "Ciprofloxacin",
                     reg_amik_amp         = "Amikacin + Ampicillin",
                     reg_plazo_amp        = "Plazomicin + Ampicillin",
                     reg_amik_cro         = "Amikacin + Ceftriaxone"
    )
  ) %>%
  group_by(regimen) %>%
  summarise(
    resistance_pct = mean(resistant == 1) * 100,
    resistance_n = sum(resistant == 1),
    .groups = "drop"
  ) %>%
  mutate(
    resistance_label = sprintf("%d%% (n=%d)", round(resistance_pct), resistance_n)
  ) %>%
  select(-resistance_pct, -resistance_n) %>%
  pivot_wider(names_from = regimen, values_from = resistance_label) %>%
  mutate(org = paste0("Overall (N=", total_N, ")"), group = " ")

# Combine and display table with regimen order
regimen_order <- c(
  "Gentamicin + Ampicillin",
  "Gentamicin + Ampicillin + Cloxacillin",
  "Ceftriaxone",
  "Ciprofloxacin",
  "Amikacin + Ampicillin",
  "Plazomicin + Ampicillin",
  "Amikacin + Ceftriaxone"
)

res_summary_combined <- bind_rows(res_summary, res_totals) %>%
  mutate(across(all_of(regimen_order), ~replace_na(.x, "-"))) %>%
  select(group, org, all_of(regimen_order)) %>%
  rename(" " = org)

final_table <- res_summary_combined %>%
  gt(groupname_col = "group") %>%

  opt_table_lines(extent = "default")

final_table

gtsave(final_table, "Tabel_all.docx")
