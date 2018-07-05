library(tidyverse)

# Load data
prop <- read_csv("data/rosehip_ish_count.csv")
prop$`PDGFRA+` <- 100 * prop$`GAD/PDGFRA`/prop$GAD
prop$`PDGFRA+/TRPC3+` <- 100 * prop$`GAD/PDGFRA/TRPC3`/prop$GAD
prop.l <- prop %>% 
  gather(probe_label, gad_prop_pct, c("PDGFRA+", "PDGFRA+/TRPC3+"))
prop.stat <- data.frame(probe_label = c("PDGFRA+", "PDGFRA+/TRPC3+"),
                           freq_mean = apply(prop[, c("PDGFRA+", "PDGFRA+/TRPC3+")], 2, mean),
                           freq_sd = apply(prop[, c("PDGFRA+", "PDGFRA+/TRPC3+")], 2, sd))


g.prop <- ggplot(prop.l, aes(x = probe_label, y = gad_prop_pct)) +
  geom_bar(data = prop.stat, aes(x = probe_label, y = freq_mean), 
           stat = "identity", fill = "light blue", width = 0.5) +
  geom_errorbar(data = prop.stat, 
                aes(x = probe_label, y = freq_mean, 
                    ymin = freq_mean - freq_sd,
                    ymax = freq_mean + freq_sd), width = 0.05) +
  geom_jitter(width = 0.15, height = 0, size = 4, shape = 1) +
  ylim(c(0, 20)) +
  xlab("") +
  ylab("Percent Layer 1 GAD1+ Cells") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="italic"))
plot(g.prop)

ggsave(g.prop, filename = "output/rosehip_prop.pdf", width = 3.5, height = 3)
