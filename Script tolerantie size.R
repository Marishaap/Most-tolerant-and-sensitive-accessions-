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




    