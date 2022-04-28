#### packages ####
library(ComplexHeatmap)
library(RColorBrewer)
library(ggpubr)
library(ggdendro)
library(dendextend)
library(tools)
library(tidyverse)

#### hclust ####
allFiles <- list.files("PPCG_vcfs_for_Sidi/", recursive=F, full.names=T)
allData <- tibble(file = allFiles) %>% 
  mutate(sampleName = str_extract(file, "vs_+[A-Za-z0-9]+_")  %>% 
           str_sub(4,-3)) %>% 
  mutate(extensionName = file_ext(file)) %>% 
  filter(extensionName == "gz") %>% 
  mutate(fileContent = map(file, ~read.table(gzfile(.x)))) %>% 
  unnest(fileContent) %>%  
  select(-file, -extensionName, -V3, -V6, -V7, -V8, -V9, -V10, -V11) %>% 
  rename(chr = V1, pos = V2, ref = V4, alt = V5) 

mutationType <- read_csv("mutationType2.csv", show_col_types = FALSE)

allData <- allData %>% 
  left_join(mutationType) %>% 
  group_by(sampleName) %>% 
  mutate(count1 = n()) %>% 
  group_by(sampleName, type) %>% 
  mutate(count2 = n()) %>% 
  mutate(propotion = count2 / count1) %>% 
  ungroup()

#### make sure there are over than 100 mutations in one chromosome
allData %>% group_by(sampleName) %>% 
  summarise(count = n()) %>% arrange(desc(count)) %>% print(n=100)

propData <- allData %>% 
  ungroup() %>%
  distinct(type, sampleName, .keep_all = T) %>% 
  select(-pos, -count1, -count2, -ref, -alt, -chr) %>% 
  pivot_wider(names_from = c(type), values_from = propotion, values_fill = 0) %>% 
  filter(sampleName %in% final$ppcg_donor_id) %>% 
  inner_join(final %>% 
               filter(type == "Clonal mutation") %>%
               select(ppcg_donor_id, new_metastatic_biology_indicator, relapse_ind, path_t_stage),
             by = c("sampleName" = "ppcg_donor_id")) %>% 
  filter(new_metastatic_biology_indicator != "not applicable") %>% 
  filter(relapse_ind != "not applicable") %>% 
  filter(relapse_ind != "missing") %>% 
  filter(!is.na(path_t_stage))


tibble(a = 1:3, b = 1:3) %>% left_join(tibble(a = 1:3, b = 1:3), by = "a")

hCluct <- hclust(dist(propData %>% column_to_rownames("sampleName")), "ave") %>%
  as.dendrogram %>%
  set("branches_k_color", k=26)

plot(hCluct)

propData2 <- propData
propData2 <-
  cutree(hCluct, k = 26) %>%
  tibble(sampleName = names(.), group = .) %>%
  inner_join(propData2, by = "sampleName") %>%
  inner_join(
    final %>%
      filter(type == "Clonal mutation") %>%
      select(ppcg_donor_id, relapse_ind, new_metastatic_biology_indicator, type, prop),
    by = c("sampleName" = "ppcg_donor_id")
  ) %>% 
  filter(group %in% c("1", "2")) %>% 
  filter(type == "Clonal mutation")

# group the patients into different clusters
cutree(hCluct, k = 26) %>%
  table()
#Statistics chi square
chisq_test(propData2$group, propData2$relapse_ind.x)
table(propData2$group, propData2$relapse_ind.x)



col = list(status = c("no_mets_biol" = "#65B33A", "mets_biol" = "#E74219"),
           relapse = c("relapsed" = "#E74219", "no relapse" = "#65B33A")
)

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  status = propData$new_metastatic_biology_indicator, 
  relapse = propData$relapse_ind, 
  col = col
)

# Combine the heatmap and the annotation

heatmapData <- as.data.frame(t(propData %>% 
                                 select(-new_metastatic_biology_indicator, -relapse_ind, -path_t_stage) %>% 
                                 column_to_rownames(var = ("sampleName"))))



library(InteractiveComplexHeatmap)
heat <- Heatmap(heatmapData %>% apply(1, scale) %>% t, name = "mutation",
                top_annotation = ha)  # %>% draw %>% htShiny
ggexport(heat, filename = "export/heatMap_of_mutation_type.pdf", 
         width = 8, height = 6)









