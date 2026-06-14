library(ggplot2)
library(tidyr)

df <- read.csv("Hall_GPP.csv")

# Sort by mean_GPP ascending
df$stream <- factor(df$stream, levels = df$stream[order(df$mean_GPP)])

# Pivot to long format for grouped bars
df_long <- pivot_longer(df, cols = c(mean_GPP, mean_ER), names_to = "metric", values_to = "value")

ggplot(df_long, aes(x = stream, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(mean_GPP = "#4CAF50", mean_ER = "#E07B54"),
                    labels = c(mean_GPP = "GPP", mean_ER = "ER")) +
  labs(x = "Stream", y = expression(paste("g O"[2], " m"^{-2}, " d"^{-1})),
       fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Hall_GPP_plot.png", width = 7, height = 4, dpi = 300)
