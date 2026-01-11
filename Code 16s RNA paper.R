if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
install.packages("dunn.test")
install.packages("pwr")
install.packages("effectsize")
install.packages("effsize")

library(pairwiseAdonis)
library(dunn.test)
library(rstatix)
library(pwr)
library(effectsize)
library(effsize)
library(car)


#Importar a un archivo phyloseq directo de archivos qza
dataSilva <- qza_to_phyloseq(features = "filtered_table_138.2.qza",
                        taxonomy = "taxonomy_silva_138.2.qza",
                        metadata = "Metadata16s_paper.tsv",
                        tree = "rooted-dada2.qza")

#Filtrado de asignaciones

tax_table(dataSilva) <- tax_table(dataSilva) %>%
  data.frame() %>% 
  filter(!grepl("Unassigned", Kingdom) & 
           !grepl("Unassigned", Phylum) & 
           Kingdom != 'd__Eukaryota' & 
           Kingdom != 'd__Archaea' &
           Genus != 'Chloroplast' &
           Genus != 'uncultured') %>%
  mutate(Genus = ifelse(Genus == 'uncultured', paste('u', Family, sep = '_'), Genus)) %>%
  as.matrix()


sample_data(dataSilva)<-sample_data(dataSilva) %>% 
  data.frame() %>% 
  mutate(Groups =str_c(sample_type, life, sep = "<br>"))

data_silva <- sample_data(dataSilva) %>% 
  data.frame()

#Transformar counts en abundancia relativa
rel.abu <- transform_sample_counts(dataSilva, function(x) x/sum(x) * 100)
#head(otu_table(rel.abu), 10)


#función para clasificar como "Others" a los taxa 
#menores al promedio de abundancia
classic_taxa <- function(physeq, taxrank, groupvar, abundance_col = "Abundance", perc_threshold = 1) {
  df <- tax_glom(physeq = physeq, taxrank = deparse(substitute(taxrank)), NArm = T)  %>%
    psmelt() %>%
    group_by({{groupvar}} ,
             {{taxrank}} ) %>%
    mutate(mean = mean(Abundance),
           Classification = {{taxrank}},
           {{taxrank}} := ifelse(mean < perc_threshold, "Others", {{taxrank}}))
  return(df)
}


#Aplicar classic_taxa a niveles de phylu, class, genus

Abundancia_phylum <- classic_taxa(rel.abu, Phylum, Groups)
Abundancia_class <- classic_taxa(rel.abu, Class, Groups)
Abundancia_Genus <- classic_taxa(rel.abu, Genus, Groups)



#################################################################################
#
#                     ABUNDANCIA RELATIVA A NIVEL DE FILO
#
##################################################################################


coloresPhylum <-palette(c("slateblue","deepskyblue", "blue","deeppink2",
                          "green2", "black"))



Abundancia_phylum$Phylum <- factor(Abundancia_phylum$Phylum,
                                   levels = c("Bacillota",
                                              "Bacteroidota",
                                              "Actinomycetota",
                                              "Pseudomonadota",
                                              "Verrucomicrobiota",
                                              "Others"))

(G.Abu.phylum <- ggplot(Abundancia_phylum, aes(x =Groups, y = Abundance, fill =Phylum, color=Phylum))+
    geom_bar(position = "fill", stat = "identity", width= 0.6) +
    theme_classic() +
    scale_fill_manual(values = coloresPhylum) +
    scale_color_manual(values = coloresPhylum) +
    scale_y_continuous(label = percent) +
    labs(x=NULL, y="Relative abundance (%)") +
    theme(legend.key.size = unit(10, "pt"),
          legend.text = element_text(size = 12, color = "black")) +
    theme(axis.text.y = element_text(color = "black"),
          strip.text = element_text(face = "bold"),
          panel.background = element_blank(),
          axis.text.x = element_markdown(color = "black", size = rel(1.5))) +
    scale_x_discrete(limits = c("Control<br>Children",
                                "Down syndrome<br>Children",
                                "Control<br>Adolescents-Adults",
                                "Down syndrome<br>Adolescents-Adults")))


ggsave("nueva data/AbundanciaphylumSILVA.png", G.Abu.phylum, width = 10, height = 8, dpi = 300)


coloresgenus <-palette(c("slateblue","deepskyblue", "blue","deeppink2",
                         "green2","darkorchid2", "goldenrod2","#3c92a8",
                         "darkblue","#c18b41","#7ea342","#c9586b","red","skyblue","orange",
                         "#e13966","#158a2c","#d1b524","#7841de","#416dab","#8b98e1","#d7882b",
                         "#b5b736", "#866429","#d545c1","#9b3fdc","#4e8ec7","#6ac042","#2ddae4",
                         "#045c31","#b0c089","pink","grey","darkblue","darkgreen","white","yellow",
                         "#7ced31","#fe5900","#b8f39a","#ffb452", "#cc7744","#e70049","#01f174","#b31a60"))


(G.Abu.Genus <- ggplot(Abundancia_Genus, aes(x = Groups, y = Abundance, fill = Genus, color = Genus)) +
    geom_bar(position = "fill", stat = "identity", width= 0.6) +
    theme_classic() +
    scale_fill_manual(values = coloresgenus) +
    scale_color_manual(values = coloresgenus) +
    scale_y_continuous(label = percent) +
    labs(x=NULL, y="Abundancia relativa (%)") +
    theme(legend.key.size = unit(15, "pt"),
          legend.text = element_text(face = "italic", size = 9),
          legend.title = element_text(size = 18),
          axis.text.y = element_text(color = "black",size = 14),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text = element_text(face = "bold"),
          panel.background = element_blank(),
          axis.text.x = element_markdown(color = "black",size = rel(1.5))) +
    guides(fill = guide_legend(ncol = 1), 
           color = guide_legend(ncol = 1)) +
    scale_x_discrete(limits = c("Control<br>Children",
                                "Down syndrome<br>Children",
                                "Control<br>Adolescents-Adults",
                                "Down syndrome<br>Adolescents-Adults")))



ggsave("nueva data/AbundanciagenusSILVA.png", G.Abu.Genus, width = 13, height = 9, dpi = 300)

#################################################################################
#
#                                  Rarefaction
#
##################################################################################
sample_data(dataSilva)<-sample_data(dataSilva) %>% 
  data.frame() %>% 
  mutate(Groups =str_c(sample_type, life, sep = "/"))


set.seed(88)
raredata <- rarefy_even_depth(dataSilva, rngseed = 88, replace = FALSE)


mat <- t(otu_table(dataSilva))
class(mat) <- "matrix"

raremax <- min(rowSums(mat))
raremax

rare2 <- rarecurve(mat, step = 20, sample = raremax, col = "blue", cex = 0.6, label = FALSE, tidy = TRUE)

rarelabels <- rare2 %>%
  merge(sample_data(raredata), by.x = "Site", by.y = 0) %>%
  group_by(Site, Groups) %>%
  summarise(Y = max(Species), X = max(Sample))

# Graficar curvas con etiquetas y línea vertical indicando rarefacción aplicada
(Rarefaction <- rare2 %>%
  merge(sample_data(raredata), by.x = "Site", by.y = 0) %>%
  ggplot(aes(x = Sample, y = Species, color = Groups, group = Site)) +
  geom_line() +
  geom_label(data = rarelabels, aes(x = X, y = Y, label = Site), size = 1.5, nudge_x = 500) +
  geom_vline(xintercept = raremax) +
  theme_bw() +
  labs(color = "Groups", 
       x = "Sequences", 
       y = "Species") +
  theme(legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(color = "black"),
        strip.text = element_text(face = "bold"),
        panel.background = element_blank(),
        legend.position = "bottom"))


ggsave("nueva data/Rarefactioncurve.png",Rarefaction, width = 12, height = 6, dpi = 300)


datararefac <- dataSilva

otu_table(datararefac)<-raredata
###############################################################################
#
#                           Alpha diversity metrics
#
##############################################################################
(richness_plot<- plot_richness(datararefac, x = "Groups", color = "Groups",
                               measures = c ("Shannon", "Simpson","Chao1"))+
   geom_boxplot(width=0.6) +
   theme_bw() +
   labs(x=NULL, color="")+
   theme(axis.title.x = element_text(size = 12, face = "bold"),
         axis.text.x = element_blank(),
         legend.position = "bottom",
         axis.title.y = element_text(size = 12, face = "bold"),
         axis.text.y = element_text(size = 12)))

richness_plot$data$variable <- factor(richness_plot$data$variable, levels = c("Shannon", "Simpson", "Chao1"))
richness_plot$data$Groups <- factor(richness_plot$data$Groups,
                                    levels = c("Control/Children", 
                                               "Down syndrome/Children",
                                               "Control/Adolescents-Adults", 
                                               "Down syndrome/Adolescents-Adults"))
richness_plot

ggsave("nueva data/AlphaDiversity.png", width = 9, height = 6, dpi = 300)


richness_df <- estimate_richness(datararefac, measures = c("Shannon", "Simpson","Chao1"))
richness_df <- richness_df %>% 
  cbind(Groups = sample_data(dataSilva)$Groups)

richness_summary <- richness_df %>%
  group_by(Groups) %>%
  summarise(across(c(Shannon, Simpson, Chao1), 
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE),
                        median = ~median(.x, na.rm = TRUE),
                        min = ~min(.x, na.rm = TRUE),
                        max = ~max(.x, na.rm = TRUE)),
                   .names = "{col}_{fn}"))


################################################################################
#                   
#                  Pairwise wilcoxon test / Benjamini-Hochberg 
#                  
################################################################################

metrics <- c("Shannon", "Simpson", "Chao1")
wilcox_results <- list()

for (metric in metrics) {
  kw <- kruskal.test(as.formula(paste(metric, "~ Groups")), data = richness_df)
  
  if (kw$p.value < 0.05) {
    wilcox_res <- richness_df %>%
      rstatix::pairwise_wilcox_test(as.formula(paste(metric, "~ Groups")), p.adjust.method = "BH", exact = FALSE) %>%
      mutate(Metric = metric)
    wilcox_results[[metric]] <- wilcox_res
  } else {
    # Si no hay diferencia global, correr pairwise para Observed igual para obtener p
    if (metric == "Chao1") {
      wilcox_res <- richness_df %>%
        rstatix::pairwise_wilcox_test(as.formula(paste(metric, "~ Groups")), p.adjust.method = "BH", exact = FALSE) %>%
        mutate(Metric = metric)
      wilcox_results[[metric]] <- wilcox_res
    } else {
      wilcox_results[[metric]] <- tibble(
        .y. = metric,
        group1 = NA_character_,
        group2 = NA_character_,
        p = NA_real_,
        p.adj = NA_real_,
        p.adj.signif = NA_character_,
        Metric = metric
      )
    }
  }
}

final_wilcox_df <- bind_rows(wilcox_results)
final_wilcox_df

#######################################################################
#                    ADD statistics values
#                  
######################################################################
richness_summary_long <- richness_summary %>%
  pivot_longer(
    cols = -Groups,
    names_to = c("Metric", "Statistic"),
    names_sep = "_"
  ) %>%
  pivot_wider(names_from = Statistic, values_from = value)


# Unión para agregar resumen a cada comparación para grupo1
final_wilcox_with_summary <- final_wilcox_df %>%
  left_join(richness_summary_long, by = c("group1" = "Groups", "Metric" = "Metric")) %>%
  rename_at(vars(mean: max), ~ paste0("group1_", .)) %>%
  # Unión para grupo2
  left_join(richness_summary_long, by = c("group2" = "Groups", "Metric" = "Metric")) %>%
  rename_at(vars(mean:max), ~ paste0("group2_", .))

final_wilcox_with_summary


write.xlsx(final_wilcox_with_summary, "wilcox_BH_summary.xlsx")

#######################################################################
#                     Kruskal.Walis test / Dunn post-hoc
#                  
######################################################################
metrics <- c("Shannon", "Simpson", "Observed")
dunn_results <- list()

for (metric in metrics) {
  kw <- kruskal.test(as.formula(paste(metric, "~ Groups")), data = richness_df)
  
  if (kw$p.value < 0.05) {
    dunn_res <- richness_df %>%
      dunn_test(as.formula(paste(metric, "~ Groups")), p.adjust.method = "bonferroni") %>%
      mutate(Metric = metric)
    dunn_results[[metric]] <- dunn_res
  } else {
    dunn_results[[metric]] <- tibble(
      .y. = metric,
      group1 = NA_character_,
      group2 = NA_character_,
      p = NA_real_,
      p.adj = NA_real_,
      p.adj.signif = NA_character_,
      Metric = metric
    )
  }
}

# Combinar todos los resultados en un solo data frame
final_dunn_df <- bind_rows(dunn_results)

#######################################################################
#                    ADD statistics values
#                  
######################################################################
richness_long <- richness_summary %>%
  pivot_longer(-Groups, names_to = c("Metric", "Stat"), names_sep = "_") %>%
  pivot_wider(names_from = Stat, values_from = value)


final_dunn_with_summary <- final_dunn_df %>%
  left_join(richness_long, by = c("group1" = "Groups", "Metric" = "Metric")) %>%
  rename_with(~ paste0("group1_", .), c(mean, sd, median, min, max)) %>%
  left_join(richness_long, by = c("group2" = "Groups", "Metric" = "Metric")) %>%
  rename_with(~ paste0("group2_", .), c(mean, sd, median, min, max))

final_dunn_with_summary

write.xlsx(final_dunn_with_summary, "final_dunn_bonferroni_summary.xlsx")


###############################################################################
#
#                              SIZE EFFECT
#
##############################################################################
richalfa_df <- richness_df %>% 
  filter(Groups %in% c("Control/Children", "Down syndrome/Children"))


efectoshannon <- richalfa_df %>%
  wilcox_effsize(Shannon ~ Groups)
print(efectoshannon)

modeloshannon <- aov(Shannon ~ Groups, data = richalfa_df)
eta_squared(modeloshannon, partial = TRUE)

shannon_df <- data.frame(
  Index      = "Shannon",
  Wilcoxon_r   = efectoshannon$effsize, 
  Eta2_partial = eta_squared(modeloshannon, partial = TRUE)$Eta2,
  Grupo1       = efectoshannon$group1,
  Grupo2       = efectoshannon$group2,
  Magnitude    = efectoshannon$magnitude,
  n1           = efectoshannon$n1,
  n2           = efectoshannon$n2
)

efectosimpson <- richalfa_df %>%
  wilcox_effsize(Simpson ~ Groups)
print(efectosimpson)

modelosimpson <- aov(Simpson ~ Groups, data = richalfa_df)
eta_squared(modelosimpson, partial = TRUE)

simpson_df <- data.frame(
  Index      = "Simpson",
  Wilcoxon_r   = efectosimpson$effsize, 
  Eta2_partial = eta_squared(modelosimpson, partial = TRUE)$Eta2,
  Grupo1       = efectoshannon$group1,
  Grupo2       = efectoshannon$group2,
  Magnitude    = efectoshannon$magnitude,
  n1           = efectoshannon$n1,
  n2           = efectoshannon$n2)


efectobserved <- richalfa_df %>%
  wilcox_effsize(Observed ~ Groups)
print(efectobserved)

modelobserved <- aov(Observed ~ Groups, data = richalfa_df)
eta_squared(modelobserved, partial = TRUE)

observed_df <- data.frame(
  Index      = "Observed",
  Wilcoxon_r   = efectobserved$effsize, 
  Eta2_partial = eta_squared(modelobserved, partial = TRUE)$Eta2,
  Grupo1       = efectoshannon$group1,
  Grupo2       = efectoshannon$group2,
  Magnitude    = efectoshannon$magnitude,
  n1           = efectoshannon$n1,
  n2           = efectoshannon$n2
)

final_data <- rbind(shannon_df,simpson_df,observed_df)
write_xlsx(final_data, "final_data_effect_sizes.xlsx")

################################################################################
#
#                          COHEN´S D TEST
#
################################################################################
shapiro.test(richness_df$Shannon[richness_df$Groups == "Control/Children"])
shapiro.test(richness_df$Shannon[richness_df$Groups == "Down syndrome/Children"])

leveneTest(Shannon ~ Groups, data = richalfa_df)

cohen.d(richness_df$Shannon[richness_df$Groups == "Control/Children"], 
        richness_df$Shannon[richness_df$Groups == "Down syndrome/Children"],
        ci = T)

cohen_d_result_shannon <- cohen.d(richness_df$Shannon[richness_df$Groups == "Control/Children"], 
                          richness_df$Shannon[richness_df$Groups == "Down syndrome/Children"],
                          hedges.correction = TRUE,
                          ci = TRUE)
print(cohen_d_result_shannon)

cohen_d_result_simpson <- cohen.d(richness_df$Simpson[richness_df$Groups == "Control/Children"], 
                                  richness_df$Simpson[richness_df$Groups == "Down syndrome/Children"],
                                  hedges.correction = TRUE,
                                  ci = TRUE)
print(cohen_d_result_simpson)



################################################################################
#                 Beta diversity metrics NMDS Unifrac Distances
#                               TODOS LOS GRUPOS
################################################################################
dist = phyloseq::distance(datararefac, method="unifrac")
set.seed(999)
unifracor = ordinate(datararefac, method="NMDS", distance=dist)
stressplot(unifracor)
(NMDS_unifrac<- plot_ordination(datararefac, unifracor, color="Groups") + 
    geom_point() +
    stat_ellipse(aes(fill = Groups),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Unweighted Unifrac Distance") + # Añadir un título al gráfico
    labs(color = "Groups", fill = "Groups"))

ggsave("nueva data/BetaDiversidad_unifrac.png",NMDS_unifrac, height = 8,width = 10,dpi = 300)

dist_all_unifrac <- phyloseq::distance(datararefac, method = "unifrac")
metadata_unifrac <- as(sample_data(datararefac), "data.frame")

# Ejecutar PERMANOVA con adonis2
set.seed(999)
permanova_res <- adonis2(dist_all_unifrac ~ Groups, data = metadata_unifrac, 
                       permutations = 999)
permanova_res
beta_results <-as.data.frame(permanova_res)
write.xlsx(beta_results,"BetaResults.xlsx")


#FUNCTION PAIRWISE PERMANOVA

pairwise_permanova <- function(dist_mat, groups, perm = 999) {
  groups <- factor(groups)
  lvls <- levels(groups)
  res_df <- data.frame(Comparison = character(), 
                       F = numeric(), 
                       R2 = numeric(), 
                       p_value = numeric(), 
                       stringsAsFactors = FALSE)
  
  for (i in seq_along(lvls)[-length(lvls)]) {
    for (j in (i + 1):length(lvls)) {
      sel <- groups %in% c(lvls[i], lvls[j])
      dist_sub <- as.dist(as.matrix(dist_mat)[sel, sel])
      groups_sub <- droplevels(groups[sel])
      set.seed(999)
      ad <- adonis2(dist_sub ~ groups_sub, permutations = perm)
      res_df <- rbind(res_df, data.frame(
        Comparison = paste0(lvls[i], "_vs_", lvls[j]),
        F = ad$F[1],
        R2 = ad$R2[1],
        p_value = ad$`Pr(>F)`[1]
      ))
    }
  }
  
  # Corrección BH (FDR)
  res_df$p_adjusted <- p.adjust(res_df$p_value, method = "BH")
  res_df <- res_df[order(res_df$p_adjusted), ]
  return(res_df)
}

set.seed(999)
results_adonis2 <- pairwise_permanova(dist_all_unifrac, metadata_unifrac$Groups, 999)
results_adonis2

write.xlsx(results_adonis2,"Betadiversidad_padjust.xlsx")

################################################################################
#               Beta diversity metrics NMDS Weighted Unifrac Distances
#                                 TODOS LOS GRUPOS
################################################################################
dist_Wunifrac = phyloseq::distance(datararefac, method="Wunifrac")
set.seed(999)
unifracor = ordinate(datararefac, method="NMDS", distance=dist_Wunifrac)
stressplot(unifracor)
(NMDS_Wunifrac<- plot_ordination(datararefac, unifracor, color="Groups") + 
    geom_point() +
    stat_ellipse(aes(fill = Groups),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Weighted Unifrac Distance") + # Añadir un título al gráfico
    labs(color = "Groups", fill = "Groups"))

ggsave("nueva data/BetaDiversidad_Wunifrac.png",NMDS_Wunifrac, height = 8,width = 10,dpi = 300)

metadata_all <- as(sample_data(datararefac), "data.frame")

# Ejecutar PERMANOVA con adonis2
set.seed(999)
permanova_res <- adonis2(dist_Wunifrac ~ Groups, data = metadata_all, 
                         permutations = 999)
permanova_res
beta_results <-as.data.frame(permanova_res)
write.xlsx(beta_results,"BetaResults_Wunifrac.xlsx")


# Ejecutar pairwise PERMANOVA
set.seed(999)
resultado_pairwise <- pairwise_permanova(dist_Wunifrac, metadata_all$Groups, 999)
print(resultado_pairwise)
write.xlsx(resultado_pairwise,"Beta_Wunifrac.xlsx")

################################################################################
#               Beta diversity metrics NMDS Unifrac Distances
#                                      CHILDREN
################################################################################
data_children <- subset_samples(datararefac, 
                             !life %in% c("Adolescents-Adults"))

sample_data(data_children)<-sample_data(data_children) %>% 
  data.frame() %>% 
  mutate(Group =str_c(Groups, BMI2, sep = "/"))


dist_unifrac_children = phyloseq::distance(data_children, method="unifrac")
set.seed(999)
unifracor = ordinate(data_children, method="NMDS", distance=dist_unifrac_children)
stressplot(unifracor)
(NMDS_unifrac_children<- plot_ordination(data_children, unifracor, color="Groups") + 
    geom_point() +
    stat_ellipse(aes(fill = Groups),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Unweighted Unifrac Distance") + # Añadir un título al gráfico
    labs(color = "Groups", fill = "Groups"))

ggsave("nueva data/BetaDiversidad_unifrac_children.png",NMDS_unifrac_children, height = 8,width = 10,dpi = 300)

metadata_children_unifrac <- as(sample_data(data_children), "data.frame")

set.seed(999)
permanova_res <- adonis2(dist_unifrac_children ~ Groups, 
                         data = metadata_children_unifrac,
                         permutations = 999)
permanova_res
beta_results_children <-as.data.frame(permanova_res)
write.xlsx(beta_results_children,"Permanova_results_children.xlsx")


################################################################################
#          Beta diversity NMDS Unifrac CHILDREN (Groups + BMI2)
#                                      
################################################################################

dist = phyloseq::distance(data_children, method="unifrac")
set.seed(999)
unifracor = ordinate(data_children, method="NMDS", distance=dist)
stressplot(unifracor)
(NMDS_unifrac_children_IMC<- plot_ordination(data_children, unifracor, color="Group") + 
    geom_point() +
    stat_ellipse(aes(fill = Group),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Unweighted Unifrac Distance")) # Añadir un título al gráfico


ggsave("nueva data/BetaDiversidad_noponderado_children_BMI.png",NMDS_unifrac_children_IMC, height = 8,width = 10,dpi = 300)

set.seed(999)
pairwise_results_children <- adonis2(dist_unifrac_children ~ Groups + BMI2, 
                                     data = metadata_children_unifrac, 
                                     perm = 999)
print(pairwise_results_children)
results_df <- as.data.frame(pairwise_results_children)
write.xlsx(results_df,"Permanova_pairwise_children_BMI.xlsx")

set.seed(999)
resultado_pairwise_children <- pairwise_permanova(dist_unifrac_children, 
                                                  metadata_children_unifrac$Group, 999)
print(resultado_pairwise_children)
write.xlsx(resultado_pairwise_children,"Permanova_pairwise_Children_multiple.xlsx")

################################################################################
#               Beta diversity metrics NMDS Weigthed Unifrac Distances
#                                      CHILDREN
################################################################################
dist_Wunifrac_children = phyloseq::distance(data_children, method="Wunifrac")
set.seed(999)
unifracor = ordinate(data_children, method="NMDS", distance=dist_Wunifrac_children)
stressplot(unifracor)
(NMDS_Wunifrac_children<- plot_ordination(data_children, unifracor, color="Groups") + 
    geom_point() +
    stat_ellipse(aes(fill = Groups),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Weighted Unifrac Distance") + # Añadir un título al gráfico
    labs(color = "Groups", fill = "Groups"))

ggsave("nueva data/BetaDiversidad_Wunifrac_children.png",NMDS_Wunifrac_children, height = 8,width = 10,dpi = 300)

metadata_children_Wunifrac <- as(sample_data(data_children), "data.frame")


# Ejecutar PERMANOVA con adonis2
set.seed(999)
permanova_res <- adonis2(dist_Wunifrac_children ~ Groups, 
                         data = metadata_children_Wunifrac,
                         permutations = 999)
permanova_res
beta_results_children <-as.data.frame(permanova_res)
write.xlsx(beta_results_children,"Adonis2_children_Wunifrac.xlsx")

################################################################################
#          Beta diversity NMDS Weighted Unifrac CHILDREN (Groups + BMI2)
#                                      
################################################################################
set.seed(999)
unifracor = ordinate(data_children, method="NMDS", distance=dist_Wunifrac_children)
stressplot(unifracor)
(NMDS_Wunifrac_children<- plot_ordination(data_children, unifracor, color="top_var") + 
    geom_point(size = 2) +
    stat_ellipse(aes(fill = top_var),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Weighted Unifrac Distance")) # Añadir un título al gráfico

ggsave("nueva data/BetaDiversidad_Wunifrac_children_BMI.png",NMDS_Wunifrac_children, height = 8,width = 10,dpi = 300)

set.seed(999)
ponderado_children <- adonis2(dist_Wunifrac_children ~ Groups + BMI2,
                                       data = metadata_children_Wunifrac,
                                       permutations = 999)
ponderado_children
results_wuni_children <- as.data.frame(ponderado_children)
write.xlsx(results_wuni_children,"Adonis_children_BMI.xlsx")

set.seed(999)
wuni_children <- pairwise_permanova(dist_Wunifrac_children, metadata_children_Wunifrac$top_var, 999)
wuni_children
write.xlsx(wuni_children,"Adonis_children_BMI_multiple.xlsx")


################################################################################
#               Beta diversity metrics NMDS Unifrac Distances
#                               Adolescents-Adults
################################################################################

data_adolescents_adults <- subset_samples(datararefac, 
                                          !life %in% c("Children"))

sample_data(data_adolescents_adults)<-sample_data(data_adolescents_adults) %>% 
  data.frame() %>% 
  mutate(Group =str_c(Groups, BMI2, sep = "/"))


dist_uni_adults = phyloseq::distance(data_adolescents_adults, method="unifrac")
set.seed(999)
unifracor = ordinate(data_adolescents_adults, method="NMDS", distance=dist_uni_adults)
stressplot(unifracor)
(NMDS_unifrac_Adoles_adults<- plot_ordination(data_adolescents_adults, unifracor, color="Groups") + 
    geom_point() +
    stat_ellipse(aes(fill = Groups),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Unweighted Unifrac Distance")) # Añadir un título al gráfico

ggsave("nueva data/BetaDiversidad_unifrac_Adults.png",NMDS_unifrac_Adoles_adults, height = 8,width = 10,dpi = 300)

metadata_adults_unifrac <- as(sample_data(data_adolescents_adults), "data.frame")

# Ejecutar PERMANOVA con adonis2
set.seed(999)
permanova_res <- adonis2(dist_uni_adults ~ Groups, 
                         data = metadata_adults_unifrac,
                         permutations = 999)
permanova_res
beta_results_adults <-as.data.frame(permanova_res)
write.xlsx(beta_results_adults,"Adults_Adonis2_beta.xlsx")


beta_unifrac_all <- ggarrange(NMDS_unifrac_children,NMDS_unifrac_Adoles_adults,
                              ncol = 2, nrow = 1, labels = c("B","C"))
beta_unifrac_all

ggsave("nueva data/BetaDiversidad_unifrac.png",beta_unifrac_all, height = 6,width = 12,dpi = 300)

merge_alfa_beta <- ggarrange(richness_plot, beta_unifrac_all,
                             ncol = 1,nrow = 2, labels = "A")
merge_alfa_beta

ggsave("nueva data/alfa_beta_diversidad.png",merge_alfa_beta, height = 11,width = 14,dpi = 300)
################################################################################
#          Beta diversity NMDS Unifrac Adults (Groups + BMI2)
#                                      
################################################################################
dist = phyloseq::distance(data_adolescents_adults, method="unifrac")
set.seed(999)
unifracor = ordinate(data_adolescents_adults, method="NMDS", distance=dist)
stressplot(unifracor)
(NMDS_unifrac_adults_IMC<- plot_ordination(data_adolescents_adults, unifracor, color="Group") + 
    geom_point() +
    stat_ellipse(aes(fill = Group),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Unweighted Unifrac Distance")) # Añadir un título al gráfico

ggsave("nueva data/BetaDiversidad_unifrac_Adults_BMI.png",NMDS_unifrac_adults_IMC, height = 8,width = 10,dpi = 300)


set.seed(999)
no_ponderado_adults <- adonis2(dist ~ Groups + BMI2,
                              data = metadata_adults_unifrac,
                              permutations = 999)
no_ponderado_adults
results_uni_adults <- as.data.frame(no_ponderado_adults)
write.xlsx(results_uni_adults,"Adonis_adults_BMI.xlsx")

set.seed(999)
uni_adults <- pairwise_permanova(dist, metadata_adults_unifrac$Group, 999)
uni_adults
write.xlsx(uni_adults,"Adonis_adults_BMI_multiple.xlsx")



beta_IMC <- ggarrange(NMDS_unifrac_children_IMC,NMDS_unifrac_adults_IMC,
                      ncol = 1,nrow = 2,labels = c("A","B"))
beta_IMC



ggsave("nueva data/Beta_IMC_noponderado.png",beta_IMC, height = 11,width = 15,dpi = 300)

################################################################################
#               Beta diversity metrics NMDS Weighted Unifrac Distances
#                               Adolescents-Adults
################################################################################
dist_wunifrac_adults = phyloseq::distance(data_adolescents_adults, method="Wunifrac")
set.seed(999)
unifracor = ordinate(data_adolescents_adults, method="NMDS", distance=dist_wunifrac_adults)
stressplot(unifracor)
(NMDS_Wunifrac_adole_adults<- plot_ordination(data_adolescents_adults, unifracor, color="Groups") + 
    geom_point() +
    stat_ellipse(aes(fill = Groups),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Weighted Unifrac Distance") + # Añadir un título al gráfico
    labs(color = "Groups", fill = "Groups"))

ggsave("nueva data/BetaDiversidad_Wunifrac_Adults.png",NMDS_Wunifrac_adole_adults, height = 8,width = 10,dpi = 300)

metadata_wuni_adults <- as(sample_data(data_adolescents_adults), "data.frame")

# Ejecutar PERMANOVA con adonis2
set.seed(999)  # para reproducibilidad
permanova_res <- adonis2(dist_wunifrac_adults ~ Groups, 
                         data = metadata_wuni_adults, 
                         permutations = 999)
permanova_res
per_res_Wuni_adults <- as.data.frame(permanova_res)
write.xlsx(per_res_Wuni_adults,"Adonis_Wuni_adults.xlsx")



################################################################################
#          Beta diversity NMDS Weighted Unifrac Adults (Groups + BMI2)
#                                      
################################################################################
dist = phyloseq::distance(data_adolescents_adults, method="Wunifrac")
set.seed(999)
unifracor = ordinate(data_adolescents_adults, method="NMDS", distance=dist)
stressplot(unifracor)
(NMDS_Wunifrac_adults_IMC<- plot_ordination(data_adolescents_adults, unifracor, color="top_var") + 
    geom_point() +
    stat_ellipse(aes(fill = top_var),geom = "polygon", alpha = 0.02) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "bottom", # Posición de la leyenda
          legend.title = element_text(face = "bold"), # Hacer el título de la leyenda en negrita
          axis.text = element_text(size = 12), # Ajustar el tamaño del texto de los ejes
          plot.title = element_text(hjust = 0.5)) + # Centrar el título del gráfico
    ggtitle("Weighted Unifrac Distance")) # Añadir un título al gráfico

ggsave("nueva data/BetaDiversidad_Wunifrac_Adults_BMI.png",NMDS_Wunifrac_adults_IMC, height = 8,width = 10,dpi = 300)


set.seed(999)
ponderado_adults <- adonis2(dist ~ Groups + BMI2,
                               data = metadata_wuni_adults,
                               permutations = 999)
ponderado_adults
results_wuni_adults <- as.data.frame(ponderado_adults)
write.xlsx(results_wuni_adults,"Adonis_Wunifrac_adults_BMI.xlsx")

set.seed(999)
wuni_adults <- pairwise_permanova(dist, metadata_wuni_adults$top_var, 999)
wuni_adults
write.xlsx(wuni_adults,"Adonis_Wunifrac_adults_BMI_multiple.xlsx")



################################################################################
#
#                                DESeq2 CHILDREN
#
################################################################################
DESq2_children <- subset_samples(datararefac, 
                                !life %in% c("Adolescents-Adults"))

deseq<- phyloseq_to_deseq2(DESq2_children, ~ Groups)
deseq$Groups <- relevel(deseq$Groups, ref = "Control/Children")

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geomean<-apply(counts(deseq), 1, gm_mean)

des<- estimateSizeFactors(deseq, geoMeans = geomean)

des<-DESeq(des, test= "Wald",fitType="local")

res_complete<- results(des) %>% 
  data.frame() %>% 
  merge(data.frame(tax_table(DESq2_children)[,"Genus"]), by = 0)

res<-results(des) %>% 
  data.frame() %>% 
  filter(padj <0.05) %>% 
  merge(
    data.frame(tax_table(DESq2_children)[,"Genus"]), by = 0
  ) %>% 
  mutate(Groups=ifelse(log2FoldChange <0, #condicion
                       "Control/Children", # si se cumple
                       "Down syndrome/Children")) # si no se cumple

write.xlsx(res,"DESq2_children.xlsx")

(Deseq_foldchangeEarly <- ggplot(res, aes(x = log2FoldChange, y = Genus, fill = Groups)) +
    geom_col(width = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = NULL, fill = "Groups",
         title = "Differential abundance") + 
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(face = "bold", size = 12),   
      axis.title.y = element_text(face = "bold", size = 12),   
      axis.text.x = element_text(size = 10),                    
      axis.text.y = element_text(size = 10, face = "italic"),
      legend.text = element_text(size = 12),
      #panel.grid.major = element_blank(),                        
      #panel.grid.minor = element_blank(),                        
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    ))

ggsave("nueva data/DESeq2_children.png", Deseq_foldchangeEarly, height = 6, width = 8,dpi = 300)

################################################################################
#
#                                DESeq2 ADULTS
#
################################################################################
DESq2_adults <- subset_samples(datararefac, 
                                          !life %in% c("Children"))
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

deseq_adult<- phyloseq_to_deseq2(DESq2_adults, ~ Groups)
deseq_adult$Groups <- relevel(deseq_adult$Groups, ref = "Control/Adolescents-Adults")

geomean<-apply(counts(deseq_adult), 1, gm_mean)

des_adult<- estimateSizeFactors(deseq_adult, geoMeans = geomean)
des_adult<-DESeq(des_adult, test= "Wald",fitType="local")

res_adult<-results(des_adult) %>% 
  data.frame() %>% 
  filter(padj <0.05) %>% 
  merge(
    data.frame(tax_table(DESq2_adults)[,"Genus"]), by = 0
  ) %>% 
  mutate(Groups=ifelse(log2FoldChange <0,
                       "Control/Adolescents-Adults",
                       "Down syndrome/Adolescents-Adults"))

write.xlsx(res_adult,"DESq2_adults.xlsx")

(Deseq_foldchangeAdult <- ggplot(res_adult, aes(x = log2FoldChange, y = Genus, fill = Groups)) +
    geom_col(width = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = NULL, fill = "Groups",
         title = "Differential abundance") + 
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(face = "bold", size = 12),   
      axis.title.y = element_text(face = "bold", size = 12),   
      axis.text.x = element_text(size = 10),                    
      axis.text.y = element_text(size = 10, face = "italic"),
      legend.text = element_text(size = 12),
      #panel.grid.major = element_blank(),                        
      #panel.grid.minor = element_blank(),                        
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    ))

ggsave("nueva data/DESeq2_adults.png", Deseq_foldchangeAdult, height = 6, width = 8,dpi = 300)



DESeq2_all <- ggarrange(Deseq_foldchangeEarly,Deseq_foldchangeAdult,
                        ncol = 1,nrow = 2,labels = c("A","B"))
ggsave("nueva data/DESeq2_All.png", DESeq2_all, height = 7, width = 9,dpi = 300)
################################################################################
#
#                       ANCOM-BC2 DIFERENTIAL ABUNDANCE CHILDREN
#
################################################################################
library(ANCOMBC)

childrenSilva <- subset_samples(dataSilva, 
                                !life %in% c("Adolescents-Adults"))

Ancom_children <- phyloseq(otu_table(childrenSilva,taxa_are_rows = T), 
                           tax_table(childrenSilva),
                          sample_data(childrenSilva))


counts <- as.data.frame(otu_table(Ancom_children))
var_taxa <- apply(counts, 1, var)
keep_taxa <- names(var_taxa[var_taxa > 0])

Ancom_children_filt <- prune_taxa(keep_taxa, Ancom_children)

Ancom_res<- ancombc2(Ancom_children_filt, fix_formula = "Groups",
                       struc_zero = TRUE,
                       group = "Groups",
                       pseudo_sens = TRUE,
                       neg_lb = TRUE,
                       p_adj_method = "BH",
                       prv_cut = 0.21,
                       tax_level = "Genus")

Ancomres <- Ancom_res$res


write.xlsx(Ancomres, "ANCOMBC_resultsChildren.xlsx")


################################################################################
#
#                       ANCOM-BC2 DIFERENTIAL ABUNDANCE ADULTS
#
################################################################################
AdultsSilva <- subset_samples(dataSilva, 
                                !life %in% c("Children"))

Ancom_adults <- phyloseq(otu_table(AdultsSilva,taxa_are_rows = T), 
                           tax_table(AdultsSilva),
                           sample_data(AdultsSilva))


counts <- as.data.frame(otu_table(Ancom_adults))
var_taxa <- apply(counts, 1, var)
keep_taxa <- names(var_taxa[var_taxa > 0])

Ancom_adults_filt <- prune_taxa(keep_taxa, Ancom_adults)

Ancom_res_adults<- ancombc2(Ancom_adults_filt, fix_formula = "Groups",
                     struc_zero = TRUE,
                     group = "Groups",
                     pseudo_sens = TRUE,
                     neg_lb = TRUE,
                     p_adj_method = "BH",
                     prv_cut = 0.25,
                     tax_level = "Genus")

Ancomres_adults <- Ancom_res_adults$res

write.xlsx(Ancomres_adults, "ANCOMBC_resultsAdults.xlsx")

