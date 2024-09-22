library(Seurat)
library(tidyverse)

#This generates a mask for non-informative spots based on compositions
# by_prop: TRUE
filter_states <- function(slide,
                          by_prop = F,
                          prop_thrsh = 0.1,
                          assay = "cell_states") {

  # Get the cell-types associated with the states
  # As a name they have ct-state_id
  state_cms <- rownames(GetAssayData(slide, assay = assay)) %>%
    strsplit(., "-") %>%
    map_chr(., ~.x[[1]]) %>%
    unique() %>%
    set_names()

  state_mats <- map(state_cms, function(ct) {
    print(ct)

    prop_vector <- GetAssayData(slide, assay = "c2l_props")[ct,]

    #If you don't want to consider pure proportions use:
    if(!by_prop){
      prop_vector <- ifelse(prop_vector >= prop_thrsh, 1, 0)
    }

    state_mat <- GetAssayData(slide, assay = assay)

    ix <-  rownames(state_mat)[grepl(ct, rownames(state_mat))]
    state_mat <- state_mat[ix, ]

    apply(state_mat, 1, function(x) {
      x * prop_vector
    }) %>%
      t()

  })

  state_mats <- purrr::reduce(state_mats, rbind)

  state_mats <- state_mats %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("barcode") %>%
    pivot_longer(-barcode,names_to = "state", values_to = "module_score") %>%
    dplyr::filter(state %in% c("Fib-Fib-SCARA5", "Fib-Myofib"))

  return(state_mats)
}

# MAIN ##########################################################################

visium_meta <- read_csv("./data/visium_meta.csv") %>%
  dplyr::rename("slide_id" = "slide_name")

all_slides <- tibble("slide_id" = list.files("/Volumes/Rico_2024/MI_slides/objects/")) %>%
  dplyr::mutate(slide_file = paste0("/Volumes/Rico_2024/MI_slides/objects/", slide_id)) %>%
  dplyr::mutate(slide_id = gsub("[.]rds", "", slide_id)) %>%
  dplyr::mutate(state_info = map(slide_file, function(sf) {
    print("processing")
    print(sf)
    filter_states(readRDS(sf))
  })) %>%
  dplyr::select(-slide_file)

all_slides <- all_slides %>%
  unnest() %>%
  left_join(visium_meta[, c("slide_id", "major_labl")], by = "slide_id") %>%
  dplyr::filter(! major_labl %in% c("BZ")) %>%
  dplyr::mutate(major_labl = ifelse(major_labl %in% c("RZ", "CTRL"),
                                    "RZ_CTRL", major_labl))

mean_info <- all_slides %>%
  dplyr::filter(state == "Fib-Myofib") %>%
  group_by(slide_id, state, major_labl) %>%
  summarise(mean_score = mean(module_score))

mean_roi <- all_slides %>%
  dplyr::filter(module_score != 0) %>%
  dplyr::filter(state == "Fib-Myofib") %>%
  group_by(slide_id, state, major_labl) %>%
  summarise(mean_score = mean(module_score))

my_comparisons <- list( c("RZ_CTRL", "IZ"),
                        c("FZ", "IZ"))

pw_data_myof <- ggpubr::compare_means(mean_score ~ major_labl,
                                      comparisons = my_comparisons,
                                      data = mean_roi) %>%
  dplyr::select(group1, group2, p , p.adj)

sptl_plt <- mean_roi %>%
  dplyr::mutate(major_labl = factor(major_labl,
                                    levels = c("RZ_CTRL", "IZ", "FZ"))) %>%
  ggplot(aes(x = major_labl, y = mean_score)) +
  geom_boxplot() +
  geom_point(alpha =0.5) +
  ggpubr::stat_pvalue_manual(pw_data_myof, label = "p.adj",
                             y.position = 25,
                             step.increase = 0.1,
                             tip.length = 0.01,size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
        axis.text = element_text(size =11)) +
  ylab("Myofibroblast score")

write_csv(mean_roi, "./results/myofib_sptl_meanscores.csv")

pdf("./results/myofib_sptl_meanscores.pdf", height = 3.5, width = 2)

plot(sptl_plt)

dev.off()

all_slides %>%
  dplyr::filter(state == "Fib-Myofib") %>%
  write_csv("./results/myofib_sptl_spotscores.csv")







