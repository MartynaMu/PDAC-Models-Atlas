# GLOBAL==========================================================================
ds = as.data.frame(qnorm)
## RNA
metadata <- data.frame("CDKN2A.state" = c(rep("Deletion", 30), rep("WT",15)),
                       "SMAD4.state" = c(rep("WT",30), rep("Deletion",15)),
                       "CFTR.state" = c(rep("WT",30), rep("Mutation",15)),
                       "Living model" = c(rep(c(rep("NO",9), rep("YES",6)),3)),
                       "Spatial model" = rep(c(rep("No", 3), rep("Yes", 12)),3),
                       "Model detailed" = as.factor(
                         c(panc1.l1, miapaca.l1, cfpac.l1)),
                       "Model general" = as.factor(
                         c(panc1.l2, miapaca.l2, cfpac.l2)),
                       "Cell line" =  as.factor(celllines))
## PROT
metadata <- data.frame("CDKN2A.state" = c(rep("Deletion", 33), rep("WT",15)),
                       "SMAD4.state" = c(rep("WT",33), rep("Deletion",15)),
                       "CFTR.state" = c(rep("WT",33), rep("Mutation",15)),
                       "Living model" = c(rep("No",12), rep("Yes",6), rep(c(rep("No",9),rep("Yes",6)),2)),
                       "Spatial model" = c(rep("No",4), rep("Yes",14), rep(c(rep("No",3), rep("Yes",12)),2)),
                       "Model detailed" = as.factor(
                         c(panc1.l1, miapaca.l1, cfpac.l1)),
                       "Model general" = as.factor(
                         c(panc1.l2, miapaca.l2, cfpac.l2)),
                       "Cell line" =  as.factor(celllines))

rownames(metadata) <- colnames(ds)
p <- ds %>% drop_na() %>% pca(metadata=metadata)
metadata_colors <- list(Model.detailed = rev(RColorBrewer::brewer.pal(5,"Paired")))

biplot(p,
       lab = NULL,
       colkey = metadata_colors$Model.detailed,
       colby = "Model.detailed",
       shape = "Cell.line",
       title = "Global proteome PCA, no filter",
       hline = 0,
       vline = 0,
       x = "PC1",
       y = "PC2")
ggsave("figures/cross-analysis/PCA/pca-global-biplot.svg", device = "svg", scale = .8, dpi = 100, width = 5.5, height = 6)

screeplot(p)
ggsave("figures/cross-analysis/PCA/pca-global-screeplot.svg", device = "svg", scale = .8, dpi = 100, width = 12, height = 6)

pairsplot(p,
          colkey = metadata_colors$Model.detailed,
          colby = "Model.detailed",
          shape="Cell.line",
          pointSize=3,
          hline = 0,
          vline = 0,
          plotaxes = FALSE)
ggsave("figures/cross-analysis/PCA/pca-global-pairsplot.svg", device = "svg", scale = 1.2, dpi = 100, width = 6, height = 6)

loadings <- plotloadings(p,
                         labSize = 3, 
                         title = "Global PCA loadings, no filter",
                         legendPosition = "none",
                         rangeRetain = 1, #retain all 100% of variables
                         components = paste0("PC", c(1:45))) 
loadings <- loadings$data
plotloadings(p,
             labSize = 3, 
             title = "Global proteome PCA loadings",
             subtitle = "top 10%, PC 1-5",
             legendPosition = "none",
             rangeRetain = .1, #retain TOP 10% variables
             components = paste0("PC", c(1:5)),
             shapeSizeRange = 5,
             col = c("blue", "white", "red")) 


ggsave("figures/cross-analysis/PCA/pca-global-loadings-PC1-5-TOP10.svg", device = "svg", scale = .8, dpi = 100, width = 15, height = 6)

eigen <- eigencorplot(p, 
                      metavars = colnames(metadata),
                      rotLabX = 45,
                      col = c('blue2', 'black', 'red2'),
                      colCorval = 'white',
                      main = "Global proteome PC correlation to sample classification",
                      posColKey = "top",
                      colFrame = "white",
                      cexMain = 1.3,
                      cexCorval = .8,
                      components = paste0("PC", c(1:10)))

svg("figures/cross-analysis/PCA/pca-global-eigencor.svg", width = 12, height = 4, pointsize = 10)
eigen
dev.off()

#
gghistogram(loadings %>% filter(PC == "PC1") %>% pull(Loading), bins = length(loadings %>% filter(PC == "PC1") %>% pull(Loading)))+
  scale_x_continuous(breaks = seq.int(-0.25,0.25,0.01))
pc1 <- loadings %>% filter(PC == "PC1", (between(abs(Loading),0.02,0.25))) %>% arrange(desc(Loading)) %>% pull(var) 
pc2 <- loadings %>% filter(PC == "PC2", (between(abs(Loading),0.01,0.25))) %>% arrange(desc(Loading)) %>% pull(var) 
pc3 <- loadings %>% filter(PC == "PC3", (between(abs(Loading),0.01,0.25))) %>% arrange(desc(Loading)) %>% pull(var) 
pc4 <- loadings %>% filter(PC == "PC4", (between(abs(Loading),0.01,0.25))) %>% arrange(desc(Loading)) %>% pull(var) 
pc5 <- loadings %>% filter(PC == "PC5", (between(abs(Loading),0.01,0.25))) %>% arrange(desc(Loading)) %>% pull(var) 
pc6 <- loadings %>% filter(PC == "PC6", (between(abs(Loading),0.01,0.25))) %>% arrange(desc(Loading)) %>% pull(var) 

# new.order <- colnames(zscored) %>% str_extract(pattern = "(?<=_)[:digit:]D_.+_[:digit:]|(?<=_)[:digit:]D_XENO_[:digit:]|(?<=_)2D_[:digit:]") %>% str_order()
# zscored.reloc <- relocate(zscored, new.order)
# zscored.reloc <- relocate(zscored.reloc, c(colnames(zscored.reloc) %>% str_which(pattern = "XENO",negate = TRUE),colnames(zscored.reloc) %>% str_which(pattern = "XENO")))
# colnames(zscored)
# colnames(zscored.reloc)

#or different order
# patterns <- c("2D_(?=[:digit:])", "3D_(?=YOUNG_[:digit:])", "3D_(?=OLD_[:digit:])", "2D_(?=XENO)", "3D_(?=XENO)")
# new.order <- c()
# for (i in patterns) {
#   temp <- colnames(zscored) %>% str_extract(pattern = i) %>% str_order(na_last = NA)
#   new.order <- c(new.order, temp)
# }
# zscored.reloc <- relocate(zscored, new.order)
# colnames(zscored.reloc)
# # anovad.reloc <- relocate(anovad, new.order)

# SOTA clustering of relevant PC---------------------------------------------------
pc = "pc1"
temp.sota <- means %>% filter(rownames(means) %in% get(pc)) %>% as.matrix() %>% sota(3)
plot(temp.sota)
assign(paste0("sota.", pc), temp.sota)

svg(paste0("figures/cross-analysis/PCA/sota/sota-", pc,".svg"), width = 10)
plot(temp.sota)
dev.off()

# zscored %>% 
#   filter(rownames(zscored) %in% pc4) %>% 
#   slice(which(sota.global.pc4$clust == 1)) %>% 
#   rownames() %>% 
#   clipr::write_clip()

# geneList <- zscored %>% 
#   filter(rownames(zscored) %in% pc4) %>% 
#   slice(which(sota.global.pc4$clust == 1)) %>% 
#   rownames()

prot.ref <- zscored %>% drop_na() %>% rownames()
prot.ref %>% writeLines("data/output/PROT/prot-ref.txt")
