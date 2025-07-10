example_ids <- prot.cor %>% select(Gene, Loadings.PC1) %>% arrange(desc(abs(Loadings.PC1))) %>% slice(1:100) %>% pull(Gene) %>% str_c(collapse=", ")

df.ann <- metadata %>%
  select(Model.detailed,Cell.line) %>%
  distinct() %>%
  mutate(Group = colnames(prot.means))

temp <- trimws(unlist(strsplit(example_ids, ","))) %>%
  c("KRT18")

mat.prot <- prot.means %>% 
  filter(rownames(prot.means) %in% temp)

mat.rna <- rna.means %>%
  filter(rownames(rna.means) %in% temp)

mat = mat.rna

mat.L <- mat %>% 
  scale() %>%
  as.data.frame() %>%
  rownames_to_column("Genes") %>%
  pivot_longer(2:16,
               names_to = "Group",
               values_to = "Mean.Intensity")

p_hm_prot <- (ggplot(mat.L, aes(x=Group, y= Genes, fill=Mean.Intensity)) + 
                geom_tile()+
                scale_fill_distiller(palette = "RdBu")) %>%
  ggplotly()

p_hm_rna <- (ggplot(mat.L, aes(x=Group, y= Genes, fill=Mean.Intensity)) + 
               geom_tile()+
               scale_fill_distiller(palette = "RdBu")) %>%
  ggplotly()

prot.subplot <- subplot(p_ann2, p_ann1, p_hm_prot, heights = c(ann_prop, ann_prop, hm_prop), nrows=3, shareX = TRUE, margin=0)
rna.subplot <- subplot(p_ann2, p_ann1, p_hm_rna, heights = c(ann_prop, ann_prop, hm_prop), nrows=3, shareX = TRUE, margin=0)

subplot(prot.subplot, rna.subplot)
