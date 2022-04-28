#### packages ####
library(survival)
library(survminer)
library(ranger)
library(ggfortify)
library(rstatix)
library(ggpubr)
library(tidyverse)

#### read table and calculate the clonality of mutations ####
clinicalData <- read_csv("PPCG_donors_clinical.20210831.shared.combimets.david.csv") #the clinical data is not public 
files <- list.files("PPCG_data/Cluster_CCFs", recursive=F, full.names=T)
geneDataRaw <- tibble(file = files) %>% 
  mutate(sampleName = str_extract(files, "s/[A-Za-z0-9]+_") %>% 
           str_sub(3,-3)) %>% 
  distinct(sampleName, .keep_all = TRUE) %>% 
  mutate(fileContent = map(file, ~read_table(.x, col_types = cols()))) %>% 
  unnest(fileContent) %>% 
  select(-file, -cluster.no) %>% 
  mutate(type = ifelse(location > 0.95 & location < 1.05, "Clonal mutation", "Subclonal mutation"))

geneFake <- expand.grid(sampleName= unique(geneDataRaw$sampleName), 
                        type = c("Subclonal mutation", "Clonal mutation"), 
                        prop1 = 0) %>% arrange(sampleName) %>% as_tibble

geneData1 <- geneDataRaw %>% 
  group_by(sampleName, type) %>% 
  summarise(total.mutation = sum(no.of.mutations)) %>% 
  group_by(sampleName) %>% 
  mutate(rate = total.mutation / sum(total.mutation))

geneData <- 
  geneFake %>% 
  left_join(geneData1, by = c("sampleName", "type")) %>% 
  rowwise() %>% 
  mutate(prop = sum(prop1, rate, na.rm = T)) %>% 
  select(-rate, -prop1)

#### join the table (SNV) #### 
stageData <- read_csv("Sample_Donor_Tissue_Origin.csv", show_col_types = FALSE) %>% 
  select(-PPCG_Sample_ID) %>% 
  distinct(PPCG_Donor_ID, .keep_all = TRUE) %>% 
  mutate(Tissue_Origin = Tissue_Origin %>% 
           recode("Normal" = "missing"))

final <- 
  clinicalData %>% inner_join(geneData, by = c("ppcg_donor_id" = "sampleName")) %>% 
  left_join(stageData, by = c("ppcg_donor_id" = "PPCG_Donor_ID")) %>% 
  mutate(path_t_stage = path_t_stage %>% 
           na_if("not applicable") %>% 
           na_if("missing")) %>% 
  mutate(path_t_stage = path_t_stage %>% 
           str_sub(2,3)) %>% 
  mutate(path_t_stage = toupper(path_t_stage)) %>% 
  mutate(gleason_predominantly = gleason_predominantly %>% 
           na_if("not applicable") %>% 
           na_if("missing")) %>%
  mutate(gleason_less_predominantly = gleason_less_predominantly %>% 
           na_if("not applicable") %>% 
           na_if("missing")) %>% 
  mutate(gleason = as.numeric(gleason_predominantly) + as.numeric(gleason_less_predominantly)) %>% 
  mutate(gleason = gleason %>% 
           na_if(5) %>% 
           na_if(10)) %>% 
  mutate(gleason = as.factor(gleason))

#### (SNV) draw boxplots add p value ####
#1 alive & dead
a1 <- ggboxplot(final %>%
                  filter(type == "Clonal mutation") %>% 
                  na_if("not applicable") %>% 
                  na_if("missing") %>% 
                  drop_na(death_ind), 
                x = "death_ind", y = "prop", facet.by = "type",
                color = "death_ind") + stat_compare_means(label = "p.format", label.y = 1.1) + theme_gray() +
  geom_jitter(aes(color = death_ind), width = 0.2, size = 0.7)
ggexport(a1, filename = "export/death_ind_bp.pdf", 
         width = 8, height = 6)
# we did not use this figure in the report



#2 path_t_stage, ns
my_comparisons1 <- list( c("T4", "T3"), c("T3", "T2"), c("T4", "T2") )
a2 <- final %>%
  filter(type == "Clonal mutation") %>% 
  na_if("not applicable") %>% 
  na_if("missing") %>% 
  drop_na(path_t_stage) %>% 
  ggboxplot(x = "path_t_stage", y = "prop", facet.by = "type",
            color = "path_t_stage") + stat_compare_means(comparisons = my_comparisons1) + theme_gray() +
  stat_compare_means(label.y = 1.5) + geom_jitter(aes(color = path_t_stage), width = 0.2, size = 0.7)
ggexport(a2, filename = "export/path_t_stage.pdf", 
         width = 8, height = 6)




#3 gleason
my_comparisons2 <- list( c("6", "7"), c("7", "8"), c("8", "9"), c("6", "8"), c("7", "9"), c("6", "9"))
a3 <- final %>%
  filter(type == "Clonal mutation") %>% 
  na_if("not applicable") %>% 
  na_if("missing") %>% 
  drop_na(gleason) %>% 
  ggboxplot(x = "gleason", y = "prop", facet.by = "type",
            color = "gleason") + stat_compare_means(comparisons = my_comparisons2) + theme_gray() +
  stat_compare_means(label.y = 1.8) + geom_jitter(aes(color = gleason), width = 0.2, size = 0.7)

ggexport(a3, filename = "export/gleason.pdf", 
         width = 8, height = 6)


#4 Risk & primary
my_comparisons3 <- list( c("high", "intermediate"), c("intermediate", "low"), c("high", "low"))

a4 <- final %>%
  # group_by(ppcg_donor_id) %>%
  # filter(sum(total.mutation, na.rm = T) > 500) %>%
  # ungroup %>%
  # filter(!is.na(total.mutation)) %>% 
  filter(type == "Clonal mutation") %>%
  # filter(total.mutation > ) %>% 
  filter(new_metastatic_biology_indicator == "no_mets_biol") %>% 
  na_if("not applicable") %>% 
  na_if("missing") %>% 
  drop_na(new_ppcg_risk_category) %>% 
  ggboxplot(x = "new_ppcg_risk_category", y = "prop", facet.by = c("type"),
            color = "new_ppcg_risk_category") + stat_compare_means(comparisons = my_comparisons3) + theme_gray() +
  stat_compare_means(label.y = 1.4) + geom_jitter(aes(color = new_ppcg_risk_category), width = 0.2, size = 0.7)
ggexport(a4, filename = "export/risk1.pdf", 
         width = 8, height = 6)

# we did not use this figure in the report

#5 add low + intermediate
a5 <-
  final %>%
  filter(type == "Clonal mutation") %>%
  mutate(new_ppcg_risk_category = new_ppcg_risk_category %>% 
           recode("low" = "low-intermediate", "intermediate" = "low-intermediate")) %>% 
  filter(new_metastatic_biology_indicator == "no_mets_biol") %>% 
  na_if("not applicable") %>% 
  na_if("missing") %>% 
  drop_na(new_ppcg_risk_category) %>% 
  ggboxplot(x = "new_ppcg_risk_category", y = "prop", facet.by = c("type"),
            color = "new_ppcg_risk_category") + stat_compare_means(label.y = 1.1) + theme_gray() +
  geom_jitter(aes(color = new_ppcg_risk_category), width = 0.2, size = 0.7)
ggexport(a5, filename = "export/risk2.pdf", 
         width = 8, height = 6)


#6 primary & meta in different tissue
change <- 
  final %>% 
  select(new_metastatic_biology_indicator, Tissue_Origin, prop, type) %>% 
  mutate(new_metastatic_biology_indicator  = new_metastatic_biology_indicator %>%
           recode("mets_biol" = "mets", 	
                  "no_mets_biol" = "no_mets")) %>% 
  mutate(state_tissue_origin = paste0(new_metastatic_biology_indicator, "_", Tissue_Origin)) %>% 
  mutate(state_tissue_origin = state_tissue_origin %>% 
           recode("mets_Metastasis" = "Metastasis", "mets_Primary" = "Primary (met)", "no_mets_Primary" = "Primary (no met)"))

my_comparisons4 <- list( c("Metastasis", "Primary (met)"), c("Primary (met)", "Primary (no met)"), c("Metastasis", "Primary (no met)"))


a6 <- change %>%
  filter(type == "Clonal mutation") %>% 
  na_if("not applicable_missing") %>% 
  na_if("missing") %>% 
  mutate(state_tissue_origin = state_tissue_origin %>% factor(levels = c("Primary (no met)", "Primary (met)", "Metastasis"))) %>% 
  drop_na(state_tissue_origin) %>%
  ggboxplot(x = "state_tissue_origin", y = "prop", facet.by = "type",
            color = "state_tissue_origin") + stat_compare_means(comparisons = my_comparisons4) + theme_gray() +
  stat_compare_means(label.y = 1.5) + geom_jitter(aes(color = state_tissue_origin), width = 0.2, size = 0.7)


ggexport(a6, filename = "export/primary & meta.pdf", 
         width = 8, height = 6)










#### (SNV) relapse & time survival analysis ####
primaryFinal <- 
  final %>% 
  filter(new_metastatic_biology_indicator =="no_mets_biol") %>% 
  filter(relapse_ind != "missing") %>% 
  filter(relapse_ind != "not applicable") %>% 
  filter(donor_relapse_interval != "missing") %>% 
  filter(donor_relapse_interval != "not applicable") %>% 
  mutate(time = as.numeric(donor_relapse_interval)) %>% 
  mutate(status = relapse_ind %>% 
           recode("no relapse" = 0, "relapsed" = 1)) %>%
  filter(type == "Clonal mutation")

primaryFinal <- primaryFinal %>% 
  mutate(med_of_cm_prop = ifelse(prop < median(prop), 1, 2))
km_med_fit <- survfit(Surv(time, status) ~ med_of_cm_prop, data = primaryFinal)
ggp <- ggsurvplot(
  km_med_fit,
  conf.int = FALSE,
  surv.median.line = c('hv'), 
  data = primaryFinal, 
  pval = TRUE,
  pval.method = TRUE,
  risk.table = FALSE) +
  ylab("Proportion of patients without relapse")
ggexport(ggp$plot, filename = "export/surv_on_clonal_mutation.pdf", 
         width = 8, height = 6)

#### data processing (CNV) ####
copyNumberData <- read_csv("Clinical_Genomic_data.csv")

addData <- 
  copyNumberData %>% 
  select(Patient, PGA, Clonal_PGA, Subclonal_PGA, Is_CombiMets, Tissue_origin) %>% 
  inner_join(final %>% 
               filter(type == "Clonal mutation") %>% 
               select(ppcg_donor_id, relapse_ind, donor_relapse_interval, death_ind, gleason, new_ppcg_risk_category, path_t_stage), by = c("Patient" = "ppcg_donor_id")) %>% 
  filter(Is_CombiMets == "no") %>% 
  filter(relapse_ind != "missing") %>% 
  filter(relapse_ind != "not applicable") %>% 
  filter(donor_relapse_interval != "missing") %>% 
  filter(donor_relapse_interval != "not applicable") %>% 
  mutate(time = as.numeric(donor_relapse_interval)) %>% 
  mutate(status = relapse_ind %>% 
           recode("no relapse" = 0, "relapsed" = 1)) %>% 
  filter(!is.na(PGA)) %>% 
  mutate(med_of_PGA = ifelse(PGA < median(PGA), 1, 2)) %>% 
  mutate(med_of_subclonal_PGA = ifelse(Subclonal_PGA < median(Subclonal_PGA), 1, 2)) %>% 
  mutate(med_of_clonal_PGA = ifelse(Clonal_PGA < median(Clonal_PGA), 1, 2)) 

#### (CNA) draw boxplots add p value ####
# wide to long

data_long <- pivot_longer(data = addData,     
                          cols = c("PGA", "Clonal_PGA", "Subclonal_PGA")) %>% 
  mutate(
    path_t_stage = path_t_stage %>% factor(levels = c("T2", "T3", "T4"))
  )

#1 alive & dead
b1 <- ggboxplot(data_long %>%
                  na_if("not applicable") %>% 
                  na_if("missing") %>% 
                  drop_na(death_ind), 
                x = "death_ind", y = "value", facet.by = "name",
                color = "death_ind") + stat_compare_means(label = "p.format", label.y = 1) + theme_gray() +
  geom_jitter(aes(color = death_ind), width = 0.2, size = 0.7)

ggexport(b1, filename = "export/death_ind_bp2.pdf", 
         width = 8, height = 6)
# we did not use this graph in the report

#2 path_t_stage
my_comparisons1 <- list( c("T4", "T3"), c("T3", "T2"), c("T4", "T2") )
b2 <- data_long %>%
  na_if("not applicable") %>% 
  na_if("missing") %>% 
  drop_na(path_t_stage) %>% 
  ggboxplot(x = "path_t_stage", y = "value", facet.by = "name",
            color = "path_t_stage") + stat_compare_means(comparisons = my_comparisons1) + theme_gray() +
  stat_compare_means(label.y = 1.5) + geom_jitter(aes(color = path_t_stage), width = 0.2, size = 0.7)
ggexport(b2, filename = "export/path_t_stage2.pdf", 
         width = 8, height = 6)


#3 gleason
my_comparisons2 <- list( c("6", "7"), c("7", "8"), c("8", "9"), c("6", "8"), c("7", "9"), c("6", "9"))
b3 <- data_long %>%
  na_if("not applicable") %>% 
  na_if("missing") %>% 
  drop_na(gleason) %>% 
  ggboxplot(x = "gleason", y = "value", facet.by = "name",
            color = "gleason") + stat_compare_means(comparisons = my_comparisons2) + theme_gray() +
  stat_compare_means(label.y = 1.8) + geom_jitter(aes(color = gleason), width = 0.2, size = 0.7)

ggexport(b3, filename = "export/gleason2.pdf", 
         width = 8, height = 6)


#5 add low + intermediate
b5 <-
  data_long %>%
  mutate(new_ppcg_risk_category = new_ppcg_risk_category %>% 
           recode("low" = "low-intermediate", "intermediate" = "low-intermediate")) %>% 
  mutate(
    new_ppcg_risk_category = new_ppcg_risk_category %>% factor(levels = c("low-intermediate", "high"))
  ) %>% 
  # filter(Is_CombiMets == "no") %>% 
  na_if("not applicable") %>% 
  na_if("missing") %>% 
  drop_na(new_ppcg_risk_category) %>% 
  ggboxplot(x = "new_ppcg_risk_category", y = "value", facet.by = c("name"),
            color = "new_ppcg_risk_category") + stat_compare_means(label.y = 1.1) + theme_gray() +
  geom_jitter(aes(color = new_ppcg_risk_category), width = 0.2, size = 0.7)
ggexport(b5, filename = "export/risk22.pdf", 
         width = 8, height = 6)



#6 primary & meta in different tissue
data_long2 <- pivot_longer(data = copyNumberData,     
                           cols = c("PGA", "Clonal_PGA", "Subclonal_PGA"))

my_comparisons4 <- list( c("Metastasis", "Primary (met)"), c("Primary (met)", "Primary (no met)"), c("Metastasis", "Primary (no met)"))

b6 <- data_long2 %>%
  na_if("BPH") %>% 
  na_if("Normal") %>% 
  na_if("Recurrence") %>% 
  drop_na(Tissue_origin) %>% 
  drop_na(value) %>% 
  mutate(Tissue_origin = Tissue_origin %>% factor(levels = c("Primary (no met)", "Primary (met)", "Metastasis"))) %>% 
  ggboxplot(x = "Tissue_origin", y = "value", facet.by = "name",
            color = "Tissue_origin") + 
  stat_compare_means(comparisons = my_comparisons4) +
  theme_gray() +
  stat_compare_means(label.y = 1.5) + geom_jitter(aes(color = Tissue_origin), width = 0.2, size = 0.7)


ggexport(b6, filename = "export/primary & meta2.pdf", 
         width = 11, height = 6)
#### (CNA) relapse & time survival analysis ####
km_med_fit2 <- survfit(Surv(time, status) ~ med_of_PGA, data = addData)
km_med_fit3 <- survfit(Surv(time, status) ~ med_of_subclonal_PGA, data = addData)
km_med_fit4 <- survfit(Surv(time, status) ~ med_of_clonal_PGA, data = addData)
survPGA1 <- ggsurvplot(
  km_med_fit2,
  conf.int = FALSE,
  surv.median.line = c('hv'), 
  data = addData, 
  pval = TRUE,
  pval.method = TRUE,
  risk.table = FALSE) +
  ylab("Proportion of patients without relapse")
ggexport(survPGA1$plot, filename = "export/surv_on_PGA.pdf", 
         width = 8, height = 6)

survPGA2 <- ggsurvplot(
  km_med_fit3,
  conf.int = FALSE,
  surv.median.line = c('hv'), 
  data = addData, 
  pval = TRUE,
  pval.method = TRUE,
  risk.table = FALSE) +
  ylab("Proportion of patients without relapse")
ggexport(survPGA2$plot, filename = "export/surv_on_sub_PGA.pdf", 
         width = 8, height = 6)

survPGA3 <- ggsurvplot(
  km_med_fit4,
  conf.int = FALSE,
  surv.median.line = c('hv'), 
  data = addData, 
  pval = TRUE,
  pval.method = TRUE,
  risk.table = FALSE) +
  ylab("Proportion of patients without relapse")
ggexport(survPGA3$plot, filename = "export/surv_on_clonal_PGA.pdf", 
         width = 8, height = 6)