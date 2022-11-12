p_tsne_temp=readRDS("data-raw/p_tsne.rds")

p_umap_temp=readRDS("data-raw/p_umap.rds")

obj_234_HTSeq=readRDS("data-raw/obj_234_HTSeqCountMini.rds")
usethis::use_data(obj_234_HTSeq)

obj_2042_HTSeq=readRDS("data-raw/obj_2042_HTSeqCountMini.rds")
usethis::use_data(obj_2042_HTSeq)
