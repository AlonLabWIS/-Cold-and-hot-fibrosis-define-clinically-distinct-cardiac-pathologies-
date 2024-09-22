# Generate proportion descriptions of Myeloid and myofibroblasts

library(tidyverse)
library(ggpubr)

liana_scores <- read_csv("./results/liana_scores.csv")[,-1]


liana_scores <- liana_scores %>%
  dplyr::filter(!major_labl %in% c("BZ")) %>%
  dplyr::mutate(major_labl = ifelse(major_labl %in% c("RZ", "CTRL"),
                                    "RZ_CTRL", major_labl)) %>%
  dplyr::mutate(major_labl = factor(major_labl,
                                    levels = c("RZ_CTRL", "IZ", "FZ")))

liana_scores %>%
  dplyr::select(hca_sample_id, major_labl) %>%
  unique() %>%
  group_by(major_labl) %>%
  summarise(n())


# Make plots:
my_comparisons <- list( c("RZ_CTRL", "IZ"),
                        c("FZ", "IZ"))

# Myeloid data:
pw_data <- liana_scores %>%
  group_by(interaction) %>%
  nest() %>%
  dplyr::mutate(comp_res = map(data, function(dat) {

    ggpubr::compare_means(mean ~ major_labl,
                          comparisons = my_comparisons,
                          data = dat) %>%
      dplyr::select(group1, group2, p , p.adj)

  })) %>%
  dplyr::select(-data)


TIMP1_LRP1_plot <- liana_scores %>%
  dplyr::filter(interaction == "TIMP1^LRP1") %>%
  ggplot(aes(x = major_labl, y = mean)) +
  geom_boxplot() +
  geom_point() +
  ggpubr::stat_pvalue_manual((pw_data %>%
                               dplyr::filter(interaction == "TIMP1^LRP1") %>%
                               pull("comp_res"))[[1]], label = "p.adj",
                             y.position = 1,
                             step.increase = 0.1,
                             tip.length = 0.01,size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =11)) +
  xlab("") +
  ylab("TIMP1^LRP1 avg. cosine")

TIMP1_POSTN_plot <- liana_scores %>%
  dplyr::filter(interaction == "TIMP1^POSTN") %>%
  ggplot(aes(x = major_labl, y = mean)) +
  geom_boxplot() +
  geom_point() +
  ggpubr::stat_pvalue_manual((pw_data %>%
                                dplyr::filter(interaction == "TIMP1^POSTN") %>%
                                pull("comp_res"))[[1]], label = "p.adj",
                             y.position = 1,
                             step.increase = 0.1,
                             tip.length = 0.01,size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =11)) +
  xlab("") +
  ylab("TIMP1^POSTN avg. cosine")

# Plot collective plot

spatial_plt <- cowplot::plot_grid(TIMP1_LRP1_plot,
                               TIMP1_POSTN_plot,
                               nrow = 1,
                               align = "hv")


pdf("./results/mean_cos_spat.pdf", height = 3.5, width = 3.5)

plot(spatial_plt)

dev.off()
