# /usr/bin/env R

library(ggplot2)
library(ggpubr)


text_sizes <- theme(
  axis.text.x = element_text(size = 5, colour = "black"),
  axis.text.y = element_text(size = 5, colour = "black"),
  axis.title.y = element_text(size = 5, colour = "black", margin = margin(t = 10, l = 5, r = 5, b = 10, unit = "pt")),
  axis.title.x = element_text(size = 5, colour = "black", margin = margin(t = 5, l = 10, r = 10, b = 5, unit = "pt")),
  legend.text = element_text(size = 5, colour = "black"),
  legend.title = element_text(size = 5, colour = "black", margin = margin(t = 5, l = 10, r = 10, b = 5, unit = "pt")),
  strip.text.x = element_text(size = 5, color = "black", face = "bold", angle = 0),
  strip.text.y = element_text(size = 5, color = "black", face = "bold", angle = 0, vjust = 0.5, hjust = 0),
  axis.ticks = element_line(color = "black", size = 0.2),
  axis.line = element_line(colour = "black", size = 0.2, linetype = "solid"),
  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
  plot.title = element_text(size = 7, face = "bold", colour = "black", margin = margin(t = 5, l = 10, r = 10, b = 5, unit = "pt"))
)


common_minimal <- text_sizes + theme(
  plot.background = element_rect(fill = NA, colour = NA),
  strip.background = element_rect(fill = NA, colour = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_rect(fill = NA, colour = "black")
)



# commonly used, x axis text 45 degree
common_45x <- common_minimal + theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1))

# commonly used,  x axis text 90 degree
common_90x <- common_minimal + theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1))

# for better display in jupyter
common_jupyter <- theme(
  axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = 1, hjust = 1),
  axis.text.y = element_text(size = 12, colour = "black"),
  axis.title.y = element_text(size = 15, colour = "black", margin = margin(t = 10, l = 5, r = 5, b = 10, unit = "pt")),
  axis.title.x = element_text(size = 15, colour = "black", margin = margin(t = 5, l = 10, r = 10, b = 5, unit = "pt")),
  legend.text = element_text(size = 12, colour = "black"),
  legend.title = element_text(size = 15, colour = "black", margin = margin(t = 5, l = 10, r = 10, b = 5, unit = "pt")),
  strip.text.x = element_text(size = 12, color = "black", face = "bold", angle = 0),
  strip.text.y = element_text(size = 12, color = "black", face = "bold", angle = 0, vjust = 0.5, hjust = 0),
  axis.ticks = element_line(color = "black", size = 0.2),
  axis.line = element_line(colour = "black", size = 0.2, linetype = "solid"),
  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
  plot.title = element_text(size = 12, face = "bold", colour = "black", margin = margin(t = 5, l = 10, r = 10, b = 5, unit = "pt"))
) +
  theme(
    plot.background = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(fill = NA, colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = NA, colour = "black")
  )
