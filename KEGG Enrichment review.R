library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(biomformat)
library(KEGGREST)
library(readr)
library(writexl)
library(openxlsx)
library(tidyverse)
library(jsonlite)
###############################################################################
#
#                     FUNTIONAL ANALYSIS (PICRUST2)
#
##############################################################################

sample_data(dataSilva)<-sample_data(dataSilva) %>% 
  data.frame() %>% 
  mutate(Groups =str_c(sample_type, life, sep = "/"))



data_children <- subset_samples(dataSilva, 
                                !life %in% c("Adolescents-Adults"))

df_children <- sample_data(data_children) %>% 
  data.frame()

data_adults <- subset_samples(dataSilva, 
                                !life %in% c("Children"))

df_adults <- sample_data(data_adults) %>% 
  data.frame()


ko_abundance <- read_tsv("pred_metagenome_unstrat.tsv.gz")
names(ko_abundance)[names(ko_abundance) == "function"] <- "KO_ID"
ko_abundance$KO_ID <- sub("^ko:", "", ko_abundance$KO_ID)

write.xlsx(ko_abundance,"KO_ABUNDANCE.xlsx")

KO_ALL <-read.xlsx("KO_ABUNDANCE.xlsx",rowNames = T)

nsti_weighted <- read_tsv("weighted_nsti.tsv.gz")

data_silva <- data_silva %>% 
  rownames_to_column(var = "sample")

nsti_weighted <- nsti_weighted %>%
  left_join(data_silva %>% select(sample, Groups), by = "sample")

write.xlsx(nsti_weighted, "NSTI_weighted.xlsx")

nsti_values <-nsti_weighted %>%
  group_by(Groups) %>%
  summarise(
    median_NSTI = median(weighted_NSTI, na.rm = TRUE),
    mean_NSTI = mean(weighted_NSTI, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop")

write.xlsx(nsti_values, "NSTI_values.xlsx")

(nsti_plot<-nsti_weighted %>%
  mutate(Groups = fct_relevel(Groups, 
                              "Control/Children", 
                              "Down syndrome/Children", 
                              "Control/Adolescents-Adults", 
                              "Down syndrome/Adolescents-Adults")) %>%
  ggplot(aes(x = Groups, y = weighted_NSTI, fill = Groups)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 2) +
  labs(x = "", y = "Weighted NSTI", 
       title = "NSTI distribution") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2"))


ggsave("comments papper/nsti_plot.png",nsti_plot, width = 10, height = 6, dpi = 300)

################################################################################################
#
#                              DATA PREPARATION ENRICHMENT CHILDREN
#
#################################################################################################
children <- rownames(df_children)

KO_niños <- KO_ALL[, colnames(KO_ALL) %in% children]

KO_childrenSD <- KO_niños[, grep("^SD", colnames(KO_niños))]
KO_childrencontrol <- KO_niños[, grep("^CT", colnames(KO_niños))]
################################################################################################
#
#                       ENRICHMENT DOWN SYNDROME-CHILDREN 40% ABUNDANCE OF KOs
#
#################################################################################################
umbralSDchildren<- ceiling(0.4 * ncol(KO_childrenSD)) 

umbral_koSDchildren <- rowSums(KO_childrenSD > 0) 


# Identifica los KO con baja abundancia (presentes en menos del umbral)
ko_lowSDchildren <- names(umbral_koSDchildren[umbral_koSDchildren < umbralSDchildren])

koSDchildren <- KO_childrenSD[!(rownames(KO_childrenSD) %in% ko_lowSDchildren), ]



enrichChildrenSD <- enrichKEGG(gene = rownames(koSDchildren), 
                               organism = "ko", 
                               keyType = "kegg",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               use_internal_data = FALSE)

Enrichmentsdchildren <- as.data.frame(enrichChildrenSD)

Enrichmentsdchildren$GeneFactor <- sapply(strsplit(Enrichmentsdchildren$GeneRatio, "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})

write.xlsx(Enrichmentsdchildren,"KEGG_enrichment_SD_children.xlsx",rownames = T)


(EnrichmentDownchildren<-ggplot(Enrichmentsdchildren, aes(x = reorder(Description,GeneFactor), 
                                                             y = GeneFactor, fill = -log10(p.adjust))) +
    geom_col(alpha = 0.8) +
    theme_bw() +
    labs(title = "KEGG Enrichment ",
         subtitle = "Down syndrome/Children",
         y = "Genes Ratio",
         x = NULL) +
    theme(axis.text.x = element_text()) +
    scale_fill_gradient(
      low = "blue", 
      high = "red", 
      breaks = c(10,20,30,40,50,60)) + # Personaliza aquí los valores de la barra
    scale_size(range = c(10, 20)) +
    guides(fill = guide_colorbar(title.position = "top")) +
    coord_flip())

ggsave("comments papper/EnrichmentDownchildren.png", EnrichmentDownchildren, width = 11, height = 12, dpi = 300)


################################################################################################
#
#                       DATA PREPARATION ENRICHMENT CONTROL-CHILDREN
#
#################################################################################################
umbralCTRLchildren<- ceiling(0.4 * ncol(KO_childrencontrol)) 

umbral_koCTRLchildren <- rowSums(KO_childrencontrol > 0) 

# Identifica los KO con baja abundancia (presentes en menos del umbral)
ko_lowCTRLchildren <- names(umbral_koCTRLchildren[umbral_koCTRLchildren < umbralCTRLchildren])

koCTRLchildren <- KO_childrencontrol[!(rownames(KO_childrencontrol) %in% ko_lowCTRLchildren), ]


enrichChildrenCTRL <- enrichKEGG(gene = rownames(koCTRLchildren), 
                               organism = "ko", 
                               keyType = "kegg",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               use_internal_data = FALSE)


EnrichmentchildrenCTRL <- as.data.frame(enrichChildrenCTRL)

EnrichmentchildrenCTRL$GeneFactor <- sapply(strsplit(EnrichmentchildrenCTRL$GeneRatio, "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})

write.xlsx(EnrichmentchildrenCTRL,"KEGG_enrichment_ChildCTRL.xlsx")


(EnrichmentchildCTRL<-ggplot(EnrichmentchildrenCTRL, aes(x = reorder(Description,GeneFactor), 
                                                          y = GeneFactor, fill = -log10(p.adjust))) +
    geom_col(alpha = 0.8) +
    theme_bw() +
    labs(title = "KEGG Enrichment ",
         subtitle = "Control/Children",
         y = "Genes ratio",
         x = NULL) +
    theme(axis.text.x = element_text()) +
    scale_fill_gradient(
      low = "blue", 
      high = "red", 
      breaks = c(10,20,30,40,50,60)) + # Personaliza aquí los valores de la barra
    scale_size(range = c(10, 20)) +
    guides(fill = guide_colorbar(title.position = "top")) +
    coord_flip())

ggsave("comments papper/EnrichmentchildrenControl.png", EnrichmentchildCTRL, width = 11, height = 12, dpi = 300)


################################################################################################
#
#                              DATA PREPARATION ENRICHMENT ADULTS
#
#################################################################################################
adults <- rownames(df_adults)

KO_adults <- KO_ALL[, colnames(KO_ALL) %in% adults]

KO_adultsSD <- KO_adults[, grep("^SD", colnames(KO_adults))]
KO_adultscontrol <- KO_adults[, grep("^CT", colnames(KO_adults))]

################################################################################################
#
#                   DATA PREPARATION ENRICHMENT Down syndrome-Adolescents Adults
#
#################################################################################################
umbralSDadults<- ceiling(0.4 * ncol(KO_adultsSD)) 

umbral_koSDadults <- rowSums(KO_adultsSD > 0) 

# Identifica los KO con baja abundancia (presentes en menos del umbral)
ko_lowSDadults <- names(umbral_koSDadults[umbral_koSDadults < umbralSDadults])

koSDadults <- KO_adultsSD[!(rownames(KO_adultsSD) %in% ko_lowSDadults), ]

enrichAdults_SD <- enrichKEGG(gene = rownames(koSDadults), 
                                 organism = "ko", 
                                 keyType = "kegg",
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 use_internal_data = FALSE)


EnrichmentAdults_SD <- as.data.frame(enrichAdults_SD)

EnrichmentAdults_SD$GeneFactor <- sapply(strsplit(EnrichmentAdults_SD$GeneRatio, "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})

write.xlsx(EnrichmentAdults_SD,"KEGG_enrichment_AdultSD.xlsx")


(EnrichmentAdultsSD<-ggplot(EnrichmentAdults_SD, aes(x = reorder(Description,GeneFactor), 
                                                         y = GeneFactor, fill = -log10(p.adjust))) +
    geom_col(alpha = 0.8) +
    theme_bw() +
    labs(title = "KEGG Enrichment ",
         subtitle = "Down syndrome/Adolescents-Adults",
         y = "Genes ratio",
         x = NULL) +
    theme(axis.text.x = element_text()) +
    scale_fill_gradient(
      low = "blue", 
      high = "red", 
      breaks = c(10,20,30,40,50,60)) + 
    scale_size(range = c(10, 20)) +
    guides(fill = guide_colorbar(title.position = "top")) +
    coord_flip())

ggsave("comments papper/EnrichmentAdultsDown.png",EnrichmentAdultsSD, width = 11, height = 12, dpi = 300)

################################################################################################
#
#                     DATA PREPARATION ENRICHMENT Control-Adolescents Adults
#
#################################################################################################
umbralCTRLadults<- ceiling(0.4 * ncol(KO_adultscontrol)) 

umbral_koCTRLadults <- rowSums(KO_adultscontrol > 0) 

# Identifica los KO con baja abundancia (presentes en menos del umbral)
ko_lowCTRLadults <- names(umbral_koCTRLadults[umbral_koSDadults < umbralCTRLadults])

koCTRLadults <- KO_adultscontrol[!(rownames(KO_adultscontrol) %in% ko_lowCTRLadults), ]

enrichAdults_CTRL <- enrichKEGG(gene = rownames(koCTRLadults), 
                              organism = "ko", 
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              use_internal_data = FALSE)


EnrichmentAdults_CTRL <- as.data.frame(enrichAdults_CTRL)

EnrichmentAdults_CTRL$GeneFactor <- sapply(strsplit(EnrichmentAdults_CTRL$GeneRatio, "/"), function(x) {
  as.numeric(x[1]) / as.numeric(x[2])
})

write.xlsx(EnrichmentAdults_CTRL,"KEGG_enrichment_AdultCTRL.xlsx")


(EnrichmentAdultsCTRL<-ggplot(EnrichmentAdults_CTRL, aes(x = reorder(Description,GeneFactor), 
                                                     y = GeneFactor, fill = -log10(p.adjust))) +
    geom_col(alpha = 0.8) +
    theme_bw() +
    labs(title = "KEGG Enrichment ",
         subtitle = "Control/Adolescents-Adults",
         y = "Genes ratio",
         x = NULL) +
    theme(axis.text.x = element_text()) +
    scale_fill_gradient(
      low = "blue", 
      high = "red", 
      breaks = c(10,20,30,40,50,60)) + 
    scale_size(range = c(10, 20)) +
    guides(fill = guide_colorbar(title.position = "top")) +
    coord_flip())

ggsave("comments papper/EnrichmentAdultsControl.png",EnrichmentAdultsCTRL, width = 11, height = 12, dpi = 300)


all_picrust <- ggarrange(EnrichmentDownchildren,EnrichmentchildCTRL,
                         EnrichmentAdultsSD,EnrichmentAdultsCTRL,
                         ncol = 2,nrow = 2,labels = c("A","B","C","D"))
all_picrust

ggsave("comments papper/Enrichment_Allgroups.png",all_picrust, width = 16, height = 14, dpi = 300)

################################################################################################
#
#                    IDENTIFICATION OF UNIQUE KEGG PATHWAYS CHILDREN
#
###############################################################################################
Downchildren <- unique(Enrichmentsdchildren$Description)
ControlChildren <- unique(EnrichmentchildrenCTRL$Description)

# Vías compartidas
compartidas <- intersect(ControlChildren, Downchildren)

children_shared_pathway <- Enrichmentsdchildren %>%
  filter(Description %in% compartidas) %>%
  select(Description, p.adjust_Down = p.adjust, Count_Down = Count, GeneRatio_Down = GeneRatio) %>%
  left_join(
    EnrichmentchildrenCTRL %>%
      filter(Description %in% compartidas) %>%
      select(Description, p.adjust_CTRL = p.adjust, Count_CTRL = Count, GeneRatio_CTRL = GeneRatio),
    by = "Description"
  ) %>%
  arrange(p.adjust_Down)

write.xlsx(children_shared_pathway,"KEGG_shared_pathways_children.xlsx")


solo_ControlChildren <- setdiff(ControlChildren, Downchildren)
solo_DownChildren <- setdiff(Downchildren, ControlChildren)
solo_DownChildren

KEGG_Children <- data.frame(
  Grupo = c("Compartidas", "Solo ControlChildren", "Solo DownChildren"),
  N_vias = c(length(compartidas), length(solo_ControlChildren), length(solo_DownChildren)),
  Porcentaje = c(
    round(length(compartidas)/length(union(ControlChildren, Downchildren))*100, 1),
    round(length(solo_ControlChildren)/length(union(ControlChildren, Downchildren))*100, 1),
    round(length(solo_DownChildren)/length(union(ControlChildren, Downchildren))*100, 1)))



################################################################################################
#
#                 IDENTIFICATION OF UNIQUE KEGG PATHWAYS CHILDREN
#
###############################################################################################

DownAdults <- unique(EnrichmentAdults_SD$Description)
ControlAdults <- unique(EnrichmentAdults_CTRL$Description)

# Vías compartidas
compartidas <- intersect(ControlAdults, DownAdults)

Adults_shared_pathways <- EnrichmentAdults_SD %>%
  filter(Description %in% compartidas) %>%
  select(Description, p.adjust_Down = p.adjust, Count_Down = Count, GeneRatio_Down = GeneRatio) %>%
  left_join(
    EnrichmentAdults_CTRL %>%
      filter(Description %in% compartidas) %>%
      select(Description, p.adjust_CTRL = p.adjust, Count_CTRL = Count, GeneRatio_CTRL = GeneRatio),
    by = "Description"
  ) %>%
  arrange(p.adjust_Down)

write.xlsx(Adults_shared_pathways,"KEGG_shared_pathways_Adults.xlsx")


solo_ControlAdults <- setdiff(ControlAdults, DownAdults)
solo_ControlAdults
solo_DownAdults <- setdiff(DownAdults, ControlAdults)
solo_DownAdults

KEGG_Children <- data.frame(
  Grupo = c("Compartidas", "Solo ControlChildren", "Solo DownChildren"),
  N_vias = c(length(compartidas), length(solo_ControlChildren), length(solo_DownChildren)),
  Porcentaje = c(
    round(length(compartidas)/length(union(ControlChildren, Downchildren))*100, 1),
    round(length(solo_ControlChildren)/length(union(ControlChildren, Downchildren))*100, 1),
    round(length(solo_DownChildren)/length(union(ControlChildren, Downchildren))*100, 1)))







