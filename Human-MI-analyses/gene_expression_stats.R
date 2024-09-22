library(tidyverse)

gene_expression <- read_csv("results/gene_expression.csv")[,-1] %>%
  dplyr::filter(! major_labl %in% c("BZ")) %>%
  dplyr::mutate(major_labl = ifelse(major_labl %in% c("RZ", "CTRL"),
                                    "RZ_CTRL", major_labl)) %>%
  group_by(patient_region_id, major_labl) %>%
  summarise(mean_LRP1 = mean(LRP1),
            mean_TIMP1 = mean(TIMP1),
            mean_POSTN = mean(POSTN))


my_comparisons <- list( c("RZ_CTRL", "IZ"),
                        c("FZ", "IZ"))

#LRP1

pw_data_LRP1 <- ggpubr::compare_means(mean_LRP1 ~ major_labl,
                                      comparisons = my_comparisons,
                                      data = gene_expression) %>%
  dplyr::select(group1, group2, p , p.adj)

LRP1_plt <- gene_expression %>%
  dplyr::mutate(major_labl = factor(major_labl,
                                    levels = c("RZ_CTRL", "IZ", "FZ"))) %>%
  ggplot(aes(x = major_labl, y = mean_LRP1)) +
  geom_boxplot() +
  geom_point(alpha =0.5) +
  ggpubr::stat_pvalue_manual(pw_data_LRP1, label = "p.adj",
                             y.position = max(gene_expression$mean_LRP1),
                             step.increase = 0.1,
                             tip.length = 0.01,size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =11)) +
  ylab("LRP1 mean expression")

#TIMP1

pw_data_TIMP1 <- ggpubr::compare_means(mean_TIMP1 ~ major_labl,
                                      comparisons = my_comparisons,
                                      data = gene_expression) %>%
  dplyr::select(group1, group2, p , p.adj)

TIMP1_plt <- gene_expression %>%
  dplyr::mutate(major_labl = factor(major_labl,
                                    levels = c("RZ_CTRL", "IZ", "FZ"))) %>%
  ggplot(aes(x = major_labl, y = mean_TIMP1)) +
  geom_boxplot() +
  geom_point(alpha =0.5) +
  ggpubr::stat_pvalue_manual(pw_data_TIMP1, label = "p.adj",
                             y.position = max(gene_expression$mean_TIMP1),
                             step.increase = 0.1,
                             tip.length = 0.01,size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =11)) +
  ylab("TIMP1 mean expression")

#POSTN

pw_data_POSTN <- ggpubr::compare_means(mean_POSTN ~ major_labl,
                                       comparisons = my_comparisons,
                                       data = gene_expression) %>%
  dplyr::select(group1, group2, p , p.adj)

POSTN_plt <- gene_expression %>%
  dplyr::mutate(major_labl = factor(major_labl,
                                    levels = c("RZ_CTRL", "IZ", "FZ"))) %>%
  ggplot(aes(x = major_labl, y = mean_POSTN)) +
  geom_boxplot() +
  geom_point(alpha =0.5) +
  ggpubr::stat_pvalue_manual(pw_data_POSTN, label = "p.adj",
                             y.position = max(gene_expression$mean_POSTN),
                             step.increase = 0.1,
                             tip.length = 0.01,size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =11)) +
  ylab("POSTN mean expression")


gene_plts <- cowplot::plot_grid(LRP1_plt, TIMP1_plt, POSTN_plt, ncol = 3, align = "hv")

pdf("./results/expr_stats.pdf", height = 3., width = 4.5)

plot(gene_plts)

dev.off()

write_csv(gene_expression,
          "results/mean_gene_expression.csv")

