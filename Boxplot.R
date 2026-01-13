library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)

klebsiella <- read_excel("~/Documents/pHd/Klebsiella_ARG.xlsx")
ecoli <- read_excel("~/Documents/pHd/Ecoli_ARG.xlsx")

ecoli <- ecoli %>% mutate(Species = "Escherichia coli")
klebsiella <- klebsiella %>% mutate(Species = "Klebsiella pneumoniae")

combined <- bind_rows(ecoli, klebsiella)

combined <- combined %>%
  mutate(`aac(6')-Ib-cr gene` = case_when(
    Amikacin == "aac(6')-Ib-cr" ~ "Positive",
    is.na(Amikacin) ~ "Negative"
  ))

# Create dataset for Amikacin
amikacin_data <- combined %>%
  select(Species, `aac(6')-Ib-cr gene`, Amikacin_MIC) %>%
  rename(MIC = Amikacin_MIC) %>%
  mutate(
    Antibiotic = "Amikacin",
    log2_MIC = log2(MIC)
  )

# Create dataset for Plazomicin
plazomicin_data <- combined %>%
  select(Species, `aac(6')-Ib-cr gene`, Plazomicin_MIC) %>%
  rename(MIC = Plazomicin_MIC) %>%
  mutate(
    Antibiotic = "Plazomicin",
    log2_MIC = log2(MIC)
  )

# Combine
plot_data <- bind_rows(amikacin_data, plazomicin_data)

# Wilcoxon p-verdier
wilcoxon_p <- plot_data %>%
  group_by(Species, Antibiotic) %>%
  summarise(
    w_p = wilcox.test(log2_MIC ~ `aac(6')-Ib-cr gene`)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    w_label = ifelse(w_p < 0.001, "<0.001", sprintf("%.3f", w_p)),
    label = paste0("p = ", w_label)
  )

# Count isolates
counts <- plot_data %>%
  group_by(Species, Antibiotic, `aac(6')-Ib-cr gene`) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(`aac(6')-Ib-cr gene`, " (n = ", n, ")"))

# Breakpoint-labels
breakpoints <- data.frame(
  Antibiotic = c("Amikacin", "Plazomicin"),
  Species = c("Escherichia coli", "Escherichia coli"),
  x = c(1.0, 1.0),
  y = c(3.4, 2.4),
  label = c("Breakpoint (MIC = 8)", "Breakpoint (MIC = 4)")
)

# Make Amikacin-plot
plot_ami <- ggplot(filter(plot_data, Antibiotic == "Amikacin"),
                   aes(x = `aac(6')-Ib-cr gene`, y = log2_MIC, color = `aac(6')-Ib-cr gene`)) +
  geom_jitter(width = 0.05, height = 0.05, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  facet_wrap(~Species) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
  geom_text(data = filter(breakpoints, Antibiotic == "Amikacin"), aes(x = x, y = y, label = label),
            inherit.aes = FALSE, color = "red", size = 4, hjust = 0.5) +
  geom_text(data = filter(wilcoxon_p, Antibiotic == "Amikacin"), aes(x = 1.5, y = 4.8, label = label),
            inherit.aes = FALSE, size = 4, hjust = 0.5) +
  geom_text(data = filter(counts, Antibiotic == "Amikacin"),
            aes(x = `aac(6')-Ib-cr gene`, y = -0.5, label = label),
            inherit.aes = FALSE, size = 3.5, color = "black") +
  labs(
    title = expression(""),
    x = "Presence of aac(6')-Ib-cr gene",
    y = expression("Log"[2]*" MIC value Amikacin")
  ) +
  scale_color_manual(values = c("Positive" = "black", "Negative" = "darkgrey")) +
  scale_y_continuous(breaks = 1:5, labels = 1:5, limits = c(-1, 5.5)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

# Make Plazomicin-plot
plot_plaz <- ggplot(filter(plot_data, Antibiotic == "Plazomicin"),
                    aes(x = `aac(6')-Ib-cr gene`, y = log2_MIC, color = `aac(6')-Ib-cr gene`)) +
  geom_jitter(width = 0.05, height = 0.05, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  facet_wrap(~Species) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_text(data = filter(breakpoints, Antibiotic == "Plazomicin"), aes(x = x, y = y, label = label),
            inherit.aes = FALSE, color = "red", size = 4, hjust = 0.5) +
  geom_text(data = filter(wilcoxon_p, Antibiotic == "Plazomicin"), aes(x = 1.5, y = 4.8, label = label),
            inherit.aes = FALSE, size = 4, hjust = 0.5) +
  geom_text(data = filter(counts, Antibiotic == "Plazomicin"),
            aes(x = `aac(6')-Ib-cr gene`, y = -0.5, label = label),
            inherit.aes = FALSE, size = 3.5, color = "black") +
  labs(
    title = expression(""),
    x = "Presence of aac(6')-Ib-cr gene",
    y = expression("Log"[2]*" MIC value Plazomicin")
  ) +
  scale_color_manual(values = c("Positive" = "black", "Negative" = "darkgrey")) +
  scale_y_continuous(breaks = 1:5, labels = 1:5, limits = c(-1, 5.5)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

# Combine with patchwork
plot_ami / plot_plaz + plot_layout(guides = "collect")
