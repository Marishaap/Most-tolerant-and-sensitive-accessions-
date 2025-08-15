setwd("C:/Users/Gebruiker/Downloads/Thesis data")

Size_data <- read.table(
  file="size_GWAS.csv",
  header = TRUE,
  sep = ","
)

library(ggplot2)
library(dplyr)
library(ggpmisc)

Size_data <- Size_data %>%
  filter(Round.Order == 95)

groen <- c("304", "216", "286", "347", "346")
rood <- c("268", "355", "322", "277", "287")

# Voeg een nieuwe kolom toe met labels
Size_data <- Size_data %>%
  mutate(kleurcode = case_when(
    LK_ID %in% groen ~ "High performance",
    LK_ID %in% rood ~ "Low performance",
    TRUE ~ "Anders"
  ))

plot1 <- ggplot(Size_data, aes(x = Size_Control, y = Size_LowP)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text(
    data = subset(Size_data, kleurcode != "Anders"),
    aes(label = LK_ID),
    hjust = -0.1, vjust = 0.5, size = 3, check_overlap = FALSE
  ) +
  
  # Eerst de grijze punten onderop
  geom_point(data = subset(Size_data, kleurcode == "Anders"), aes(color = kleurcode), size = 3, alpha = 0.8) +
  
  # Daarna de groene en rode punten erbovenop
  geom_point(data = subset(Size_data, kleurcode != "Anders"), aes(color = kleurcode), size = 3, alpha = 0.8) +
  
  scale_color_manual(values = c("High performance" = "darkgreen", "Low performance" = "red")) +
  labs(x = "PLA Control", y = "PLA Low P", color = "Accessiongroup") +
  theme_minimal(base_size = 14)

size_test <- cor.test(Size_data$Size_Control, Size_data$Size_LowP, method = "spearman")

PSII_df <- read.table(
  file="PSII_allrounds.csv",
  header = TRUE,
  sep = ","
)

PSII_df <- PSII_df %>%
  mutate(P_status = case_when(
    TABLE %in% c(1, 2) ~ "LowP",
    TABLE == 3 ~ "C"
  ))

PSII_df <- PSII_df %>%
  filter(Round.Order == 95)

groen <- c("304", "216", "286", "347", "346")
rood <- c("268", "355", "322", "277", "287")

# Voeg een nieuwe kolom toe met labels
PSII_df <- PSII_df %>%
  mutate(kleurcode = case_when(
    LK_ID %in% groen ~ "High performance",
    LK_ID %in% rood ~ "Low performance",
    TRUE ~ "Anders"
  ))

PSII_df <- PSII_df %>%
  group_by(LK_ID, P_status) %>%
  summarise(PSII = mean(PSII, na.rm = TRUE), .groups = "drop")

PSII_df <- PSII_df %>%
  select(LK_ID, PSII, P_status) %>%
  pivot_wider(names_from = P_status, values_from = PSII, names_prefix = "PSII_")



plot2 <-  ggplot(PSII_df, aes(x = PSII_C, y = PSII_LowP)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text(
    data = subset(PSII_df, kleurcode != "Anders"),
    aes(label = LK_ID),
    hjust = -0.1, vjust = 0.5, size = 3, check_overlap = FALSE
  ) +
  
  # Eerst de grijze punten onderop
  geom_point(data = subset(PSII_df, kleurcode == "Anders"), aes(color = kleurcode), size = 3, alpha = 0.8) +
  
  # Daarna de groene en rode punten erbovenop
  geom_point(data = subset(PSII_df, kleurcode != "Anders"), aes(color = kleurcode), size = 3, alpha = 0.8) +
  
  scale_color_manual(values = c("High performance" = "darkgreen", "Low performance" = "red")) +
  labs(x = "ΦPSII Control", y = "ΦPSII Low P", color = "Accessiongroup") +
  theme_minimal(base_size = 14)

psii_test <- cor.test(PSII_df$PSII_C, PSII_df$PSII_LowP, method = "spearman")

FvFm_data <- read.table(
  file="FvFm_GWAS.csv",
  header = TRUE,
  sep = ","
)

FvFm_data <- FvFm_data %>%
  filter(Round.Order == 97)

groen <- c("304", "216", "286", "347", "346")
rood <- c("268", "355", "322", "277", "287")

# Voeg een nieuwe kolom toe met labels
FvFm_data <- FvFm_data %>%
  mutate(kleurcode = case_when(
    LK_ID %in% groen ~ "High performance",
    LK_ID %in% rood ~ "Low performance",
    TRUE ~ "Anders"
  ))

FvFm_data <- FvFm_data %>%
  mutate(LK_ID = as.character(LK_ID))

plot3 <- ggplot(FvFm_data, aes(x = FvFm_Control, y = FvFm_LowP)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text(
    data = subset(FvFm_data, kleurcode != "Anders"),
    aes(label = LK_ID),
    hjust = -0.1, vjust = 0.5, size = 3, check_overlap = FALSE
  ) +
  
  # Eerst de grijze punten onderop
  geom_point(data = subset(FvFm_data, kleurcode == "Anders"), aes(color = kleurcode), size = 3, alpha = 0.8) +
  
  # Daarna de groene en rode punten erbovenop
  geom_point(data = subset(FvFm_data, kleurcode != "Anders"), aes(color = kleurcode), size = 3, alpha = 0.8) +
  
  scale_color_manual(values = c("High performance" = "darkgreen", "Low performance" = "red")) +
  labs(x = "Fv/Fm Control", y = "Fv/Fm Low P", color = "Accessiongroup") +
  theme_minimal(base_size = 14)

fvfm_test <- cor.test(FvFm_data$FvFm_Control, FvFm_data$FvFm_LowP, method = "spearman")

NPQ_data <- read.table(
  file="NPQ_GWAS.csv",
  header = TRUE,
  sep = ","
)

NPQ_data <- NPQ_data %>%
  filter(Round.Order == 95)

groen <- c("304", "216", "286", "347", "346")
rood <- c("268", "355", "322", "277", "287")

# Voeg een nieuwe kolom toe met labels
NPQ_data <- NPQ_data %>%
  mutate(kleurcode = case_when(
    LK_ID %in% groen ~ "High performance",
    LK_ID %in% rood ~ "Low performance",
    TRUE ~ "Anders"
  ))


NPQ_data$Ratio_Control_LowP <- NPQ_data$NPQ_Control / NPQ_data$NPQ_LowP

npq_test <- cor.test(NPQ_data$NPQ_Control, NPQ_data$NPQ_LowP, method = "spearman")

plot4 <- ggplot(NPQ_data, aes(x = NPQ_Control, y = NPQ_LowP)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text(
    data = subset(NPQ_data, kleurcode != "Anders"),
    aes(label = LK_ID),
    hjust = -0.1, vjust = 0.5, size = 3, check_overlap = FALSE
  ) +
  
  # Eerst de grijze punten onderop
  geom_point(data = subset(NPQ_data, kleurcode == "Anders"), aes(color = kleurcode), size = 3, alpha = 0.8) +
  
  # Daarna de groene en rode punten erbovenop
  geom_point(data = subset(NPQ_data, kleurcode != "Anders"), aes(color = kleurcode), size = 3, alpha = 0.8) +
  
  scale_color_manual(values = c("High performance" = "darkgreen", "Low performance" = "red")) +
  labs(x = "NPQ Control", y = "NPQ Low P", color = "Accessiongroup") +
  theme_minimal(base_size = 14)

library(patchwork)


(
  (plot1 + plot2 + plot3 + plot4) + 
    plot_layout(guides = "collect") +
    plot_annotation(
      tag_levels = 'A',
      theme = theme(
        plot.tag = element_text(size = 10)  # Pas dit getal aan voor kleinere/grotere letters
      )
    )
) &
  theme(
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 5)
  )
    





    
