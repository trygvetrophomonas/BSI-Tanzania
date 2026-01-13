library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

df <- read_excel("~/Documents/pHd/Dataset_AST.xlsx")

df_patient <- df %>%
  group_by(id) %>%
  summarise(
    ageadm = first(ageadm),
    neo    = first(neo),
    
    # Combination regimens: resistent if 1 bacteria resistent to both antibiotics
    AMP_CN  = any(AMP == 0 & CN == 0),
    AMP_AK  = any(AMP == 0 & AK == 0),
    AMP_PLZ = any(AMP == 0 & PLZ == 0),
    CR_AK   = any(CR == 0 & AK == 0),
    
    CR      = any(CR == 0),
    CIP     = any(CIP == 0)
  ) %>%
  mutate(study = ifelse(grepl("^BSI", id), "Study 1", "Study 2")) %>%
  mutate(
    reg_ampigenta  = ifelse(AMP_CN, 1, 0),
    reg_ampiamik   = ifelse(AMP_AK, 1, 0),
    reg_ampiplazo  = ifelse(AMP_PLZ, 1, 0),
    reg_croamik    = ifelse(CR_AK, 1, 0),
    reg_cro        = ifelse(CR, 1, 0),
    reg_cip        = ifelse(CIP, 1, 0)
  )

# Summarize resistance
res_summary_patient <- df_patient %>%
  pivot_longer(cols = starts_with("reg_"),
               names_to = "regimen", values_to = "resistant") %>%
  filter(!is.na(resistant)) %>%
  group_by(study, regimen) %>%
  summarise(resistance_rate = mean(resistant == 1) * 100, .groups = "drop") %>%
  mutate(regimen = recode(regimen,
                          reg_ampigenta  = "Gentamicin-Ampicillin",
                          reg_cro        = "Ceftriaxone",
                          reg_cip        = "Ciprofloxacin",
                          reg_ampiamik   = "Amikacin-Ampicillin",
                          reg_ampiplazo  = "Plazomicin-Ampicillin",
                          reg_croamik    = "Amikacin-Ceftriaxone"))

# Format for dumbbell-plot
res_wide <- res_summary_patient %>%
  pivot_wider(names_from = study, values_from = resistance_rate) %>%
  arrange(desc(`Study 2`)) %>%
  mutate(regimen = factor(regimen, levels = rev(regimen)))

res_wide_long <- res_wide %>%
  pivot_longer(cols = c(`Study 1`, `Study 2`), names_to = "study", values_to = "rate")

# === Dumbbell-plot ===
p <- ggplot(res_wide_long, aes(x = rate, y = regimen, color = study)) +
  geom_segment(data = res_wide, aes(x = `Study 1`, xend = `Study 2`, y = regimen, yend = regimen),
               color = "gray70", size = 1.5, inherit.aes = FALSE) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Study 1" = "#0e668b", "Study 2" = "#c64737")) +
  labs(
    title = "",
    x = "Resistance Rate (%)",
    y = "Antibiotic Regimen",
    color = "Study"
  ) +
  theme_minimal()

print(p)

ggsave("dumbbell_patient_resistance.pdf", p, width = 8, height = 6)
