

sample_data(dataSilva)<-sample_data(dataSilva) %>% 
  data.frame() %>% 
  mutate(Groups =str_c(sample_type, life, sep = "/"))

###########################################################################################
#                                                                                         #
#               Spearman correlation with BH Down syndrome - Children                     #
#                                                                                         #
###########################################################################################
DownChildren <- subset_samples(dataSilva, sample_data(dataSilva)$Groups =="Down syndrome/Children")
data_children_filter <- phyloseq_filter_prevalence(DownChildren,prev.trh = 0.25)
childrenSD<- tax_glom(data_children_filter, taxrank = "Genus")

otu_child_sd <- otu_table(childrenSD) %>% 
  data.frame

otu_childSD <- otu_table(childrenSD) %>% 
  data.frame() %>% 
  log1p()
tax_childSD <- tax_table(childrenSD) %>% 
  data.frame()
rownames(otu_childSD) <- tax_childSD$Genus
rownames(otu_child_sd) <- tax_childSD$Genus

#Normalization Z score
scaledownchildren <- sample_data(childrenSD)%>% 
  data.frame()

scaledownchildren <- scaledownchildren %>% 
  select(BMI,Total_Cholesterol, Triglicerydes,HDL.c,
         VLDL.c, LDL.c, Glucose) %>%
  dplyr::rename(
    `Total_Cholesterol` = Total_Cholesterol,
    `HDL-c` =  HDL.c,
    `VLDL-c`= VLDL.c,
    `LDL-c` = LDL.c) %>%
  scale()

#Correlación spearman
correlationstableDownchildren <- associate(t(otu_childSD), scaledownchildren, method = "spearman", mode = "table", 
                               p.adj.threshold = 0.05, p.adj.method = "BH", n.signif = 0)

write.xlsx(correlationstableDownchildren,"Spearman_correlation_DownChildren.xlsx")


grupos_interes <- c("Control/Children", "Control/Adolescents-Adults",
                    "Down syndrome/Children", "Down syndrome/Adolescents-Adults")

rel.abu_children <- transform_sample_counts(data_children_filter, function(x) x/sum(x) * 100)

Abundancia_Genus_children <- classic_taxa(rel.abu_children, Genus, Groups)

AbuGenus_children <- Abundancia_Genus_children %>%
  filter(Groups %in% grupos_interes) %>%
  dplyr::select(Sample, Genus, mean, Abundance, Groups,Glucose)

AbuGenus_all <- AbuGenus_children %>%
  dplyr::filter(Groups == "Down syndrome/Children", Genus == "Fusicatenibacter") %>%
  dplyr::select(Sample,Genus, mean, Abundance, Groups,Glucose)

AbuGenus_filtrado <- AbuGenus_all %>%
  dplyr::filter(Abundance > 0)

spearman_test <- cor.test(AbuGenus_all$Abundance, AbuGenus_all$Glucose, method = "spearman")
print(spearman_test)

(Fusi_nofilter<-ggplot(AbuGenus_all, aes(x = Abundance, y = Glucose)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "spearman", label.x = 4, label.y = max(AbuGenus_filtrado$Glucose)*1) +
  labs(title = "Spearman correlation: Fusicatenibacter/Glucose",
       x = "Fusicatenibacter abundance",
       y = "Glucose") +
  theme_bw())

ggsave("nueva data/SpearmanGlucose_fusicatenibacter.png", Fusi_nofilter,height = 7,width = 7,dpi = 300)

(Fusi_filter<-ggplot(AbuGenus_filtrado, aes(x = Abundance, y = Glucose)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    stat_cor(method = "spearman", label.x = 4, label.y = max(AbuGenus_filtrado$Glucose)*1) +
    labs(title = "Spearman correlation: Fusicatenibacter/Glucose",
         x = "Fusicatenibacter abundance",
         y = "Glucose") +
    theme_bw())

ggsave("nueva data/SpearmanGlucose_fusicatenibacter_filter.png", Fusi_filter,height = 7,width = 7,dpi = 300)




library(coin)
zeros_idx <- factor(AbuGenus_all$Abundance == 0, labels = c("non_zero", "zero"))
table(zeros_idx)

perm_p <- independence_test(Abundance ~ Glucose | zeros_idx, 
                            data = AbuGenus_all, 
                            distribution = approximate(nresample = 1000),
                            teststat = "quad")
print(perm_p)


###########################################################################################
#                                                                                         #
#               Spearman correlation with BH Down syndrome - Adolescents-Adults           #
#                                                                                         #
###########################################################################################
DownAdults <- subset_samples(dataSilva, sample_data(dataSilva)$Groups =="Down syndrome/Adolescents-Adults")
data_Adults_filter <- phyloseq_filter_prevalence(DownAdults,prev.trh = 0.25)
AdultsSD<- tax_glom(data_Adults_filter, taxrank = "Genus")

otu_adultSD <- otu_table(AdultsSD) %>% 
  data.frame() %>% 
  log1p()
tax_adultSD <- tax_table(AdultsSD) %>% 
  data.frame()
rownames(otu_adultSD) <- tax_adultSD$Genus

#Normalization Z score
scaledownadult <- sample_data(AdultsSD)%>% 
  data.frame()

scaledownadult <- scaledownadult %>% 
  select(IMC,Colesterol_total, Triglicéridos, C_HDL,
         C_VLDL, C_LDL, Glucosa) %>%
  dplyr::rename(
    `BMI` = IMC,
    `Total Cholesterol` = Colesterol_total,
    Triglycerides = Triglicéridos,
    `HDL-c` =  C_HDL,
    `VLDL-c`=C_VLDL,
    `LDL-c` = C_LDL,
    Glucose = Glucosa
  ) %>%
  scale()

#Correlación spearman
correlationstableDownadult <- associate(t(otu_adultSD), scaledownadult, method = "spearman", mode = "table", 
                               p.adj.threshold = 0.05, p.adj.method = "BH", n.signif = 0)

write.xlsx(correlationstableDownadult,"Spearman_correlation_DownAdults.xlsx")


###########################################################################################
#                                                                                         #
#               Spearman correlation with BH Control - Children                           #
#                                                                                         #
###########################################################################################
ControlChildren <- subset_samples(dataSilva, sample_data(dataSilva)$Groups =="Control/Children")
data_Control_filter <- phyloseq_filter_prevalence(ControlChildren,prev.trh = 0.25)
Controlchild<- tax_glom(data_Control_filter, taxrank = "Genus")

otu_controlchild <- otu_table(Controlchild) %>% 
  data.frame() %>% 
  log1p()
tax_controlchild <- tax_table(Controlchild) %>% 
  data.frame()
rownames(otu_controlchild) <- tax_controlchild$Genus

#Normalization Z score
scalecontrolchildren <- sample_data(Controlchild)%>% 
  data.frame()

scalecontrolchildren <- scalecontrolchildren %>% 
  select(IMC,Colesterol_total, Triglicéridos, C_HDL,
         C_VLDL, C_LDL, Glucosa) %>%
  dplyr::rename(
    `BMI` = IMC,
    `Total Cholesterol` = Colesterol_total,
    Triglycerides = Triglicéridos,
    `HDL-c` =  C_HDL,
    `VLDL-c`=C_VLDL,
    `LDL-c` = C_LDL,
    Glucose = Glucosa
  ) %>%
  scale()

#Correlación spearman
correlationstablecontrolchildren <- associate(t(otu_controlchild), scalecontrolchildren, method = "spearman", mode = "table", 
                                           p.adj.threshold = 0.05, p.adj.method = "BH", n.signif = 0)

write.xlsx(correlationstablecontrolchildren,"Spearman_correlation_controlChildren.xlsx")

###########################################################################################
#                                                                                         #
#               Spearman correlation with BH Control - Adolescents-Adults                 #
#                                                                                         #
###########################################################################################
Controladults <- subset_samples(dataSilva, sample_data(dataSilva)$Groups =="Control/Adolescents-Adults")
data_Control_filter <- phyloseq_filter_prevalence(Controladults,prev.trh = 0.30)
controladults<- tax_glom(data_Control_filter, taxrank = "Genus")

otu_controladults <- otu_table(controladults) %>% 
  data.frame() %>% 
  log1p()
tax_controladults <- tax_table(controladults) %>% 
  data.frame()
rownames(otu_controladults) <- tax_controladults$Genus

#Normalization Z score
scalecontroladults <- sample_data(controladults)%>% 
  data.frame()

scalecontroladults <- scalecontroladults %>% 
  select(IMC,Colesterol_total, Triglicéridos, C_HDL,
         C_VLDL, C_LDL, Glucosa) %>%
  dplyr::rename(
    `BMI` = IMC,
    `Total Cholesterol` = Colesterol_total,
    Triglycerides = Triglicéridos,
    `HDL-c` =  C_HDL,
    `VLDL-c`=C_VLDL,
    `LDL-c` = C_LDL,
    Glucose = Glucosa
  ) %>%
  scale()

#Correlación spearman
correlationstablecontroladults <- associate(t(otu_controladults), scalecontroladults, method = "spearman", mode = "table", 
                                              p.adj.threshold = 0.05, p.adj.method = "BH", n.signif = 0)

write.xlsx(correlationstablecontroladults,"Spearman_correlation_controladults.xlsx")


#########################################################################################
#
#                               HEATMAP CON MERGE SAMPLES 
#
##########################################################################################

#Transformación a una prevalencia del 20%
sample_data(dataSilva)<-sample_data(dataSilva) %>% 
  data.frame() %>% 
  mutate(Groups =str_c(sample_type, life, sep = "/"))

core.Genus <- transform_sample_counts(dataSilva, function(x) x/sum(x) * 100) %>%
  core(prevalence = 0.20, detection = 1/100) %>% 
  tax_glom(taxrank = 'Genus')

# merge samples
core.Genus <- merge_samples(core.Genus, sample_data(core.Genus)$Groups, fun = mean)
taxa.matrix <- otu_table(core.Genus) %>% 
  t()

taxa.labels<- tax_table(core.Genus)[,"Genus"] %>% 
  data.frame() %>% 
  mutate(Genus = gsub("(.*)(_.*)","\\1",Genus))

#Anotacion de los filos a la columna de row para identificar los generos
ann_row <- tax_table(core.Genus)[,"Phylum"] %>% 
  data.frame()

#Anotación de columnas correspondientes a observar como grupo,IMC y Edad
ann_col <- sample_data(core.Genus) %>% 
  data.frame() %>% 
  select(Groups)
ann_col$Groups <- rownames(ann_col)

orden_grupos <- c(
  "Control/Children", 
  "Down syndrome/Children", 
  "Control/Adolescents-Adults", 
  "Down syndrome/Adolescents-Adults"
)

# Convierte la columna 'Grupos' en un factor ordenado
ann_col$Groups <- factor(ann_col$Groups, levels = orden_grupos)

#Asignación de colores especificos para cada anotacion 

colorphylum <- viridis::plasma(length(levels(factor(ann_row$Phylum))))
names(colorphylum) <- levels(factor(ann_row$Phylum))

Colorsample <- c("#c18b41","orange","deepskyblue","aquamarine")
names(Colorsample) <- levels(factor(ann_col$Groups))

#Variable con la lista de variables con los colores asignados para las anotaciones
ann_colors <- list(Phylum = colorphylum,Grupos=Colorsample)

# Plot core microbiota heatmap

(core.heatmap <- ComplexHeatmap::pheatmap(taxa.matrix, scale = 'row',
                                          annotation_row = ann_row,
                                          annotation_colors = ann_colors,
                                          annotation_col = ann_col,
                                          cluster_cols = T,
                                          labels_row = taxa.labels$Genus,
                                          border_color = "black",
                                          labels_col = NULL,
                                          show_colnames = F,
                                          heatmap_legend_param = list(title = 'Scale'),
                                          col = circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red")),
                                          treeheight_row = 100))

(coreabundance<-pheatmap(taxa.matrix, scale = "row",
                         cluster_rows = T, 
                         cluster_cols = T,
                         annotation_col = ann_col,
                         annotation_row = ann_row,
                         border_color = "black",
                         #annotation_colors = ann_colors,
                         #clustering_distance_rows = "manhattan",
                         #clustering_distance_cols = "manhattan", 
                         #clustering_method = "average",
                         #cutree_rows = 7,
                         #cutree_cols = 4,
                         labels_row = taxa.labels$Genus,
                         show_colnames = F))

coreheatmapaper <- as.ggplot(coreabundance)

ggsave("nueva data/coreheatmap.png", coreheatmapaper, width = 10, height = 8, dpi = 300)

