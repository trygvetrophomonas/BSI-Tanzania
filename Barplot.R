library(readxl)
library(tidyverse)
library(tidytext)

klebsiella <- read_excel("~/Documents/pHd/Klebsiella_ARG.xlsx")
ecoli <- read_excel("~/Documents/pHd/Ecoli_ARG.xlsx")
salmonella <- read_excel("~/Documents/pHd/salmonella_ARG.xlsx")

ecoli$organism <- "Escherichia coli"
klebsiella$organism <- "Klebsiella pneumoniae"
salmonella$organism <- "Salmonella enterica"

all <- bind_rows(ecoli, salmonella, klebsiella) %>%
  select(ID, organism, everything())

# Long format
data_long <- all %>%
  pivot_longer(
    cols = c(ESBL, Amikacin, Gentamicin, Ciprofloxacin, Penicillinase),
    names_to = "Type",
    values_to = "Gen"
  ) %>%
  filter(!is.na(Gen)) %>%
  mutate(
    Gen = str_split(Gen, ", |/"),
    Type = case_when(
      Type == "Penicillinase" ~ "Ampicillin",
      Type == "ESBL" ~ "Ceftriaxone",
      TRUE ~ Type
    )
  ) %>%
  unnest(Gen) %>%
  mutate(
    Gen = str_trim(Gen),
    Gen = case_when(
      str_starts(Gen, "parE") ~ "parE",
      str_starts(Gen, "parC") ~ "parC",
      str_starts(Gen, "gyrA") ~ "gyrA",
      str_starts(Gen, "gyrB") ~ "gyrB",
      str_starts(Gen, "marR") ~ "marR",
      TRUE ~ Gen
    ),
    ceftriaxone_flag = Type == "Ceftriaxone"
  ) %>%
  filter(!(organism == "Klebsiella pneumoniae" & Type == "Ampicillin"))

# Duplicate Ceftriaxone-genes to Ampicillin, exept Klebsiella
ceftriaxone_genes <- data_long %>%
  filter(Type == "Ceftriaxone" & organism != "Klebsiella pneumoniae") %>%
  mutate(Type = "Ampicillin")

data_long <- bind_rows(data_long, ceftriaxone_genes)

gene_counts_label <- data_long %>%
  group_by(organism, Type, Gen, ceftriaxone_flag) %>%
  summarise(Antall = n(), .groups = "drop") %>%
  group_by(organism, Type) %>%
  arrange(desc(ceftriaxone_flag), Antall) %>%
  mutate(
    Gen_clean = Gen,
    Gen = reorder_within(Gen, Antall, interaction(organism, Type)),
    pos = cumsum(Antall) - Antall / 2
  )

gene_counts_label <- gene_counts_label %>%
  mutate(Type = factor(Type, levels = c(
    "Ampicillin", "Ceftriaxone", "Gentamicin", "Amikacin", "Ciprofloxacin"
  )))

label_positions <- gene_counts_label %>%
  group_by(organism, Type) %>%
  summarise(Antall = sum(Antall), .groups = "drop")

all_types <- expand.grid(
  organism = unique(gene_counts_label$organism),
  Type = levels(gene_counts_label$Type)
)

label_positions <- full_join(label_positions, all_types, by = c("organism", "Type")) %>%
  mutate(Antall = replace_na(Antall, 0)) %>%
  filter(!(organism == "Klebsiella pneumoniae" & Type == "Ampicillin"))

multi_type_genes <- data_long %>%
  distinct(Gen, Type) %>%
  count(Gen) %>%
  filter(n > 1) %>%
  pull(Gen)

gene_colors <- setNames(
  ifelse(
    str_remove(levels(gene_counts_label$Gen), "__.+$") == "aac(6')-Ib-cr", "skyblue",
    ifelse(str_remove(levels(gene_counts_label$Gen), "__.+$") %in% multi_type_genes, "red", "grey70")
  ),
  levels(gene_counts_label$Gen)
)

# Plot
plot <- ggplot(gene_counts_label, aes(
  x = Type,
  y = Antall,
  fill = Gen
)) +
  geom_text(data = label_positions,
            aes(x = Type, y = 0, label = Type),
            inherit.aes = FALSE,
            vjust = 1.3,
            fontface = "bold",
            size = 3) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Gen_clean,
                vjust = ifelse(Antall == 1, -0.5, 0.5)),
            size = 2, color = "black",
            position = position_stack(vjust = 0.5)) +
  facet_wrap(~ organism, ncol = 2) +
  labs(
    x = "",
    y = "Frequencies"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  scale_fill_manual(values = gene_colors) +
  guides(fill = "none")

print(plot)

