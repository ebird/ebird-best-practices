library(auk)
library(tidyverse)

ebd <- read_sampling("data-raw/ebd_CO_bahtan1_smp_relAug-2025_sampling.txt") |>
  filter(year(observation_date) > 2022)

# daily checklist counts
ebd |>
  filter(year(observation_date) == 2024) |>
  count(month = month(observation_date),
        day = yday(observation_date)) |>
  write_csv("~/Desktop/co-daily-lists.csv")

#
ggplot(ebd) +
  aes(x = duration_minutes / 60) +
  geom_density(linewidth = 1, color = "#457999") +
  labs(x = "Checklist duration [hours]",
       y = "Proportion of checklists") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(0, 10),
                  ylim = c(0, NA),
                  expand = FALSE) +
  theme_classic() +
  theme(#axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(10, 10, 5, 5))
ggsave("~/Desktop/co-duration.png", width = 5, height = 5)

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
  theme(#axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 10, 5, 5))
ggsave("~/Desktop/co-distance.png", width = 5, height = 5)
