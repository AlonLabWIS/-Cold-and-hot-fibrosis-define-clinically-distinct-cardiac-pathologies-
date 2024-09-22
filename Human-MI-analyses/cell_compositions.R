# Generate proportion descriptions of Myeloid and myofibroblasts

library(HDF5Array)
library(tidyverse)
library(compositions)
library(SummarizedExperiment)

mi_pb <- loadHDF5SummarizedExperiment("data/scell/hca_rnasamples_sce/")

obs_dat <- SummarizedExperiment::colData(mi_pb) %>%
  as.data.frame() %>%
  dplyr::select(patient_region_id, cell_type, cell_state, major_labl)

rm(mi_pb)

pat_meta <- obs_dat  %>%
  dplyr::select(patient_region_id, major_labl) %>%
  unique()

# Spatial data only

momics_compositions <- read_tsv("./data/scell/spatial_compositions.txt")

momics_props <- momics_compositions %>%
  dplyr::select(-sp_n_cells) %>%
  dplyr::rename("patient_region_id" = patient_id) %>%
  pivot_wider(names_from = cell_type, values_from = sp_prop_cells,values_fill = 0) %>%
  dplyr::select(patient_region_id, Fib, Myeloid) %>%
  left_join(pat_meta, by = "patient_region_id") %>%
  dplyr::filter(! major_labl %in% c("BZ")) %>%
  dplyr::mutate(major_labl = ifelse(major_labl %in% c("RZ", "CTRL"), "RZ_CTRL", major_labl)) %>%
  dplyr::mutate(major_labl = factor(major_labl,
                                    levels = c("RZ_CTRL", "IZ", "FZ")))

my_comparisons <- list( c("RZ_CTRL", "IZ"),
                        c("FZ", "IZ"))


pw_data_myel <- ggpubr::compare_means(Myeloid ~ major_labl,
                                      comparisons = my_comparisons,
                                      data = momics_props) %>%
  dplyr::select(group1, group2, p , p.adj)

myeloid_plot <- momics_props %>%
  ggplot(aes(x = major_labl, y = Myeloid)) +
  geom_boxplot() +
  geom_point(alpha =0.5) +
  ggpubr::stat_pvalue_manual(pw_data_myel, label = "p.adj",
                             y.position = 0.49,
                             step.increase = 0.1,
                             tip.length = 0.01,size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =11)) +
  xlab("")

pw_data_fib <- ggpubr::compare_means(Fib ~ major_labl,
                                     comparisons = my_comparisons,
                                     data = momics_props) %>%
  dplyr::select(group1, group2, p , p.adj)

fib_plot <- momics_props %>%
  ggplot(aes(x = major_labl, y = Fib)) +
  geom_boxplot() +
  geom_point(alpha =0.5) +
  ggpubr::stat_pvalue_manual(pw_data_fib, label = "p.adj",
                             y.position = 0.55,
                             step.increase = 0.1,
                             tip.length = 0.01,size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =11)) +
  xlab("")

spat_plt <- cowplot::plot_grid(myeloid_plot,
                               fib_plot,
                               nrow = 1,
                               align = "hv")


pdf("./results/props_spatial.pdf", height = 3.5, width = 3.25)

plot(spat_plt)

dev.off()




