library(arrow)
library(auk)
library(tidyverse)

ebd <- read_parquet("data/erd_checklists_1km_tw.paquet") |>
  filter(year > 2022)

# daily checklist counts
ebd |>
  filter(year == 2024) |>
  count(day_of_year) |>
  write_csv("~/Desktop/tw-daily-lists.csv")

ggplot(ebd) +
  aes(x = effort_hours) +
  geom_density(linewidth = 1, color = "#457999") +
  labs(x = "Checklist duration [hours]",
       y = "Proportion of checklists") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  coord_cartesian(xlim = c(0, 8),
                  ylim = c(0, NA),
                  expand = FALSE) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(10, 10, 5, 5))
ggsave("~/Desktop/tw-duration.png", width = 5, height = 5)

ggplot(ebd) +
  aes(x = effort_distance_km) +
  geom_density(linewidth = 1, color = "#457999") +
  labs(x = "Distance traveled [km]",
       y = "Proportion of checklists") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(0, 10),
                  ylim = c(0, NA),
                  expand = FALSE) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(10, 10, 5, 5))
ggsave("~/Desktop/tw-distance.png", width = 5, height = 5)
