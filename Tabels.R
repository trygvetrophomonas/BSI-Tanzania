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
    
    # Combination regimens: resistant if one bacteria is resistant to both antibiotics
    reg_ampigenta = any(AMP == 0 & CN == 0),
    reg_ampiamik  = any(AMP == 0 & AK == 0),
    reg_ampiplazo = any(AMP == 0 & PLZ == 0),
    reg_croamik   = any(CR == 0 & AK == 0),
    
    reg_cro       = any(CR == 0),
    reg_cip       = any(CIP == 0)
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
    reg_ampigenta, reg_ampiamik, reg_ampiplazo,
    reg_cro, reg_croamik, reg_cip
  ) %>%
  pivot_longer(
    cols = starts_with("reg_"),
    names_to = "regimen",
    values_to = "resistant"
  ) %>%
  mutate(
    regimen = recode(regimen,
                     reg_ampigenta  = "Ampicillin + Gentamicin",
                     reg_ampiamik   = "Ampicillin + Amikacin",
                     reg_ampiplazo  = "Ampicillin + Plazomicin",
                     reg_cro        = "Ceftriaxone",
                     reg_croamik    = "Ceftriaxone + Amikacin",
                     reg_cip        = "Ciprofloxacin"
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
    str_detect(as.character(org), "Klebsiella|E. coli|Salmonella|Other Gram negatives*|Pseudomonas|Acinetobacter") ~ "Gram negative",
    str_detect(as.character(org), "S. aureus|Enterococcus|Streptococcus") ~ "Gram positive",
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
                     reg_ampigenta  = "Ampicillin + Gentamicin",
                     reg_ampiamik   = "Ampicillin + Amikacin",
                     reg_ampiplazo  = "Ampicillin + Plazomicin",
                     reg_cro        = "Ceftriaxone",
                     reg_croamik    = "Ceftriaxone + Amikacin",
                     reg_cip        = "Ciprofloxacin"
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
  mutate(org = paste0("Total (all organisms) (N=", total_N, ")"), group = " ")

# Combine and display table
regimen_order <- c(
  "Ampicillin + Gentamicin",
  "Ceftriaxone",
  "Ampicillin + Amikacin",
  "Ciprofloxacin",
  "Ampicillin + Plazomicin",
  "Ceftriaxone + Amikacin"
)

res_summary_combined <- bind_rows(res_summary, res_totals) %>%
  mutate(across(all_of(regimen_order), ~replace_na(.x, "-"))) %>%
  select(group, org, all_of(regimen_order)) %>%
  rename(" " = org)

final_table <- res_summary_combined %>%
  gt(groupname_col = "group") %>%
  tab_header(title = "Percentage predicted failure because of resistance") %>%
  opt_table_lines(extent = "default")

final_table

gtsave(final_table, "Tabel_all.docx")

##################################################################

# Split into Study 1 and Study 2
df_patient_study1 <- df_patient %>% filter(str_detect(id, "^BSI"))
df_patient_study2 <- df_patient %>% filter(str_detect(id, "^[0-9]+$"))

create_res_table <- function(df_subset) {
  total_N <- nrow(df_subset)
  
  org_counts <- df_subset %>% filter(!is.na(org)) %>% count(org, name = "N")
  
  data_labeled <- df_subset %>%
    left_join(org_counts, by = "org") %>%
    mutate(
      org_labeled = paste0(org, " (N=", N, ")"),
      org_labeled = factor(org_labeled, levels = paste0(levels(df_patient$org), " (N=", org_counts$N[match(levels(df_patient$org), org_counts$org)], ")"))
    )
  
  res_summary <- data_labeled %>%
    select(
      org = org_labeled,
      reg_ampigenta, reg_ampiamik, reg_ampiplazo,
      reg_cro, reg_croamik, reg_cip
    ) %>%
    pivot_longer(
      cols = starts_with("reg_"),
      names_to = "regimen",
      values_to = "resistant"
    ) %>%
    mutate(
      regimen = recode(regimen,
                       reg_ampigenta  = "Ampicillin + Gentamicin",
                       reg_ampiamik   = "Ampicillin + Amikacin",
                       reg_ampiplazo  = "Ampicillin + Plazomicin",
                       reg_cro        = "Ceftriaxone",
                       reg_croamik    = "Ceftriaxone + Amikacin",
                       reg_cip        = "Ciprofloxacin"
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
      str_detect(as.character(org), "Klebsiella|E. coli|Salmonella|Other Gram negatives*|Pseudomonas|Acinetobacter") ~ "Gram negative",
      str_detect(as.character(org), "S. aureus|Enterococcus|Streptococcus") ~ "Gram positive",
      TRUE ~ "Polymicrobial"
    ))
  
  res_totals <- data_labeled %>%
    pivot_longer(cols = starts_with("reg_"),
                 names_to = "regimen", values_to = "resistant",
                 values_drop_na = TRUE) %>%
    mutate(
      regimen = recode(regimen,
                       reg_ampigenta  = "Ampicillin + Gentamicin",
                       reg_ampiamik   = "Ampicillin + Amikacin",
                       reg_ampiplazo  = "Ampicillin + Plazomicin",
                       reg_cro        = "Ceftriaxone",
                       reg_croamik    = "Ceftriaxone + Amikacin",
                       reg_cip        = "Ciprofloxacin"
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
    mutate(org = paste0("Total (all organisms) (N=", total_N, ")"), group = " ")
  
  regimen_order <- c(
    "Ampicillin + Gentamicin",
    "Ceftriaxone",
    "Ampicillin + Amikacin",
    "Ciprofloxacin",
    "Ampicillin + Plazomicin",
    "Ceftriaxone + Amikacin"
  )
  
  # Combine and make table
  res_summary_combined <- bind_rows(res_summary, res_totals) %>%
    mutate(across(all_of(regimen_order), ~replace_na(.x, "-"))) %>%
    select(group, org, all_of(regimen_order)) %>%
    rename(" " = org)
  
  gt(res_summary_combined, groupname_col = "group") %>%
    tab_header(title = "Percentage predicted failure because of resistance") %>%
    opt_table_lines(extent = "default")
}


table_study1 <- create_res_table(df_patient_study1)
table_study2 <- create_res_table(df_patient_study2)

table_study1
table_study2

gtsave(table_study1, "Tabel_study1.docx")
gtsave(table_study2, "Tabel_study2.docx")

#####################################################################

df_patient <- df_patient %>%
  mutate(
    ageadm = as.numeric(ageadm),
    neo = as.character(neo),
    age_group = case_when(
      (!is.na(neo) & neo == "1") | (!is.na(ageadm) & ageadm <= 28) ~ "neonate",
      (!is.na(neo) & neo == "0") | (!is.na(ageadm) & ageadm > 28) ~ "older_child",
      TRUE ~ NA_character_
    )
  )

table(df_patient$age_group, useNA = "ifany")

df_patient_neonates <- df_patient %>% filter(age_group == "neonate")
df_patient_older_children <- df_patient %>% filter(age_group == "older_child")

table_neonates <- create_res_table(df_patient_neonates)
table_older_children <- create_res_table(df_patient_older_children)

table_neonates
table_older_children

gtsave(table_neonates, "Table_neonates.docx")
gtsave(table_older_children, "Table_older_children.docx")

#######################################
#Died/survived

df_patient <- df_patient %>%
  mutate(
    mortality = case_when(
      died == 1 ~ "Died",
      died == 0 ~ "Alive",
      is.na(died) & outcome == "Died" ~ "Died",
      is.na(died) & outcome == "Improved and discharged" ~ "Alive",
      TRUE ~ NA_character_
    )
  )

df_patient_alive <- df_patient %>% filter(mortality == "Alive")
df_patient_died  <- df_patient %>% filter(mortality == "Died")

table_alive <- create_res_table(df_patient_alive)
table_died <- create_res_table(df_patient_died)

table_alive
table_died

gtsave(table_alive, "Table_alive.docx")
gtsave(table_died, "Table_died.docx")

