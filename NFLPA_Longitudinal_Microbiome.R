###LIBRA Mndsource Microbiome 
library(vegan)
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ANCOMBC)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(vegan)
library(knitr)
library(forcats)
library(tidyverse)
library(DT)
library(outliers)
library(beeswarm)
library(WGCNA)
library(flextable)
library(magrittr)
library(dplyr)
library(tidyheatmaps)
library(devtools)
library(ggvegan)
library(janitor)
library(stringr)
library(gtsummary)
#microbiomestat object
library(microbiome)
library(MicrobiomeStat)
library(GUniFrac)

#save and load data
load("/Path/To/NFLPA_Longidutinal_Microbiome.RData")
my_colors = c("#A77903", 
              "#1B3558",
              "#3DC8E0",
              "#5A595A")
options(
  ggplot2.discrete.fill = my_colors,
  ggplot2.discrete.colour = my_colors
)
theme_set(theme_bw())

#Set theme bw
theme_set(theme_bw())
#ste up path for later use
path = '/Out/Path/'

#Loads in meta data for all timepoints data
meta_mindsource

#Adding LIBRA scores to meta data
#LIBRA Data
libra_meta 
libra_meta_bl
#merge Libra information together
libra_meta_bl2 = libra_meta_bl %>%
  mutate(Ever_Homeless = case_when(homelessn >= 1 ~ "Yes",
                                   homelessn == 0 ~ "No",
                                   T~"Unknown") ) %>%
  mutate(LIBRA_Half = case_when(Total_LIBRA_Score >= median(Total_LIBRA_Score) ~"Upper",
                                Total_LIBRA_Score < median(Total_LIBRA_Score) ~"Lower")) %>%
  dplyr::rename(Participant = participant) 

libra_meta_bl2$Education = factor(libra_meta_bl2$Education, levels=c("High school diploma or equivalent", "Some college, no degree", "Associate's degree", 
                                                                     "Bachelor's degree", "Master's degree"))

libra_meta_bl2$Ever_Homeless = factor(libra_meta_bl2$Ever_Homeless, levels=c("Yes", "No", "Unknown"))
libra_meta_bl2$sex = factor(libra_meta_bl2$sex)

#Filters mindsource meta data to only those in the study
participant_list_mindsource

meta_mindsource_filter = meta_mindsource %>% filter(Participant %in% participant_list_mindsource$Participants)
meta_mindsource_filter_update_ID = meta_mindsource_filter %>% 
  mutate(Sample_ID = str_replace_all(SampleID, fixed("."), "_")) %>%
  mutate(Sample_ID = str_replace_all(Sample_ID, fixed("-"), "_"))


meta_and_libra = meta_mindsource_filter_update_ID %>% mutate(Time_num = case_when(Time_Point == "BL" ~ 0,
                                                                                  Time_Point == "48h" ~ 0.07,
                                                                                  Time_Point == "2m" ~ 2,
                                                                                  Time_Point == "3m" ~ 3,
                                                                                  Time_Point == "4m" ~ 4,
                                                                                  Time_Point == "5m" ~ 5,
                                                                                  Time_Point == "6m" ~ 6,
                                                                                  Time_Point == "7m" ~ 7,
                                                                                  Time_Point == "8m" ~ 8,
                                                                                  Time_Point == "9m" ~ 9,
                                                                                  Time_Point == "10m" ~ 10,
                                                                                  Time_Point == "11m" ~ 11,
                                                                                  Time_Point == "12m" ~ 12))


meta_and_libra$Time_Point = factor(meta_and_libra$Time_Point, levels = c("BL", "48h","2m", "3m", "4m", "5m", "6m", "7m","8m","9m","10m","11m","12m" ))
meta_and_libra_final = full_join(meta_and_libra, libra_meta_bl2, by="Participant") %>% filter(!Sample_ID %in% c('19_1307_P24_8m_g', '19_1307_P34_2m_g'))

meta_and_libra_final %>% group_by(Participant) %>% summarize(count=n())
######Table 1 Demographics#####
#Creates table1 and saves to the working directory
table1 <-
  tbl_summary(
    libra_meta_bl2,
    include = c(Total_LIBRA_Score, age, sex, race, Ever_Homeless, BMI_Current, Education),
    by = LIBRA_Half, # split table by group
    missing = "ifany", # don't list missing data separately
    statistic = all_continuous() ~ "{mean} (± {sd})" # Custom format for mean and SD
  ) |> 
  add_n() |> # add column with total number of non-missing observations
  add_p() |> # test for a difference between groups
  add_overall() %>%
  modify_header(label = "**Variable**", stat_1 = "**Lower Brain Health Risk Group**  \nN = {n}", stat_2 = "**Higher Brain Health Risk Group**  \nN = {n}") |> # update the column header
  bold_labels() |>
  modify_spanning_header( ~ "**Table 1 Summary of Demographics for Mindsource Microbiome Study**") |>
  as_flex_table() |> 
  flextable::save_as_docx(path = paste0(path, "Table_1_Summary_Mindsource.docx"))

libra_meta_bl2$Healthy_Diet =factor(libra_meta_bl2$Healthy_Diet)
libra_meta_bl2$Low_Alcohol =factor(libra_meta_bl2$Low_Alcohol)
libra_meta_bl2$Physical_Inactivity =factor(libra_meta_bl2$Physical_Inactivity)
libra_meta_bl2$Obesity =factor(libra_meta_bl2$Obesity)
libra_meta_bl2$Depression =factor(libra_meta_bl2$Depression)
libra_meta_bl2$Heart_Disease =factor(libra_meta_bl2$Heart_Disease)
libra_meta_bl2$CKD =factor(libra_meta_bl2$CKD)
libra_meta_bl2$Diabetes =factor(libra_meta_bl2$Diabetes)
libra_meta_bl2$Hypercholesterolemia =factor(libra_meta_bl2$Hypercholesterolemia)
libra_meta_bl2$Smoking =factor(libra_meta_bl2$Smoking)
libra_meta_bl2$Hypertension =factor(libra_meta_bl2$Hypertension)

#Creates table 2 for LIBRA Score characteristics
table2 <-
  tbl_summary(
    libra_meta_bl2,
    include = c(Total_LIBRA_Score, Healthy_Diet, Low_Alcohol, Physical_Inactivity, Obesity, Depression,Heart_Disease,CKD,Diabetes,Hypercholesterolemia,Smoking,Hypertension),
    by = LIBRA_Half, # split table by group
    missing = "ifany", # don't list missing data separately
    statistic = all_continuous() ~ "{mean} (± {sd})" # Custom format for mean and SD
  ) |> 
  add_n() |> # add column with total number of non-missing observations
  add_p() |> # test for a difference between groups
  add_overall() %>%
  modify_header(label = "**Variable**", stat_1 = "**Lower Brain Health Risk Group**  \nN = {n}", stat_2 = "**Higher Brain Health Risk Group**  \nN = {n}") |> # update the column header
  bold_labels() |>
  modify_spanning_header( ~ "**Table 2 Summary of Characteristics for Mindsource Microbiome Study**") |>
  as_flex_table() |> 
  flextable::save_as_docx(path = paste0(path, "Table_2_Summary_Mindsource.docx"))
######


####build out the phyloseq object
taxa_table_1 = read.table('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0001_DeepSeq_Jul_2023_Diversigen/filtered taxonomy tables/taxatables-aggregated-absolute/taxatable-species-absolute.tsv', header =T)
trans_taxa_table_1 = t(data.frame(taxa_table_1))
trans_taxa_table_1 = trans_taxa_table_1 %>% row_to_names(row_number = 1)
trans_taxa_table_1 = data.frame(trans_taxa_table_1)
trans_taxa_table_1$Sample_ID = rownames(trans_taxa_table_1)
rownames(trans_taxa_table_1) = NULL
trans_taxa_table_1 = trans_taxa_table_1 %>% dplyr::select(Sample_ID, everything())


taxa_table_2 = read.table('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0003_DeepSeq_May_2024_Diversigen/filtered taxonomy tables/taxatables-by-level/taxatable-species-absolute.tsv', header=T)
trans_taxa_table_2 = t(data.frame(taxa_table_2))
trans_taxa_table_2 = trans_taxa_table_2 %>% row_to_names(row_number = 1)
trans_taxa_table_2 = data.frame(trans_taxa_table_2)
trans_taxa_table_2$Sample_ID = rownames(trans_taxa_table_2)
rownames(trans_taxa_table_2) = NULL
trans_taxa_table_2 = trans_taxa_table_2 %>% dplyr::select(Sample_ID, everything())


taxa_table_all = bind_rows(trans_taxa_table_1,trans_taxa_table_2)
taxa_table_all[is.na(taxa_table_all)] <- 0
taxa_table_all_ID_Update = taxa_table_all %>% mutate(Sample_ID = str_replace(Sample_ID, "X", ""))
taxa_table_all_ID_Update$Sample_ID
taxa_table_all_filtered = taxa_table_all_ID_Update %>% dplyr::filter(Sample_ID %in% meta_mindsource_filter_update_ID$Sample_ID & (!Sample_ID %in% c('19_1307_P24_8m_g', '19_1307_P34_2m_g')))
rownames(taxa_table_all_filtered) = taxa_table_all_filtered$Sample_ID
taxa_table_all_filtered_mat = taxa_table_all_filtered %>% select(-Sample_ID) %>% select(where(~ any(. != 0)))
taxa_table_all_filtered_mat[] = lapply(taxa_table_all_filtered_mat, as.numeric)
taxa_table_all_filtered_mat = as.matrix(taxa_table_all_filtered_mat)

#otu table
otu_table = t(taxa_table_all_filtered_mat)
#taxa names
tax_names = rownames(otu_table)
#Update
tax_names_2 = data.frame(tax_names) %>%
  separate(tax_names, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\.")

rownames(tax_names_2) = tax_names

tax_names_2_mat = as.matrix(tax_names_2)

OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(tax_names_2_mat)
sample = data.frame(meta_and_libra_final)
rownames(sample) = sample$Sample_ID
sampledata = sample_data(sample)

physeq = phyloseq(OTU, TAX, sampledata)
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq_object_mindsource = phyloseq(OTU, TAX, sampledata, random_tree)
physeq_object_mindsource_bac <- subset_taxa(physeq_object_mindsource, Kingdom == "k__Bacteria")

random_tree2 = rtree(ntaxa(physeq_object_mindsource_bac), rooted=TRUE, tip.label=taxa_names(physeq_object_mindsource_bac))
physeq_object_mindsource_final = phyloseq(OTU, TAX, sampledata, random_tree2)


####alpha diversity
alpha_meta = estimate_richness(physeq_object_mindsource_final,measures=c("Shannon","Observed")) %>%
  rownames_to_column(var = "Sample_ID") 
alpha_meta2 = alpha_meta %>% mutate(Sample_ID = str_replace(Sample_ID, "X", ""))
alpha_meta2$Sample_ID

alpha_meta_final = left_join(alpha_meta2, as.data.frame(sample_data(physeq_object_mindsource_final)), by = "Sample_ID") 



#colvec <- c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A")

#Spaghetti plots
shannon_line = ggplot(alpha_meta_final, aes(x=Time_num, y=Shannon,group=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = alpha_meta_final)+
  labs(tag = "A", x="Time (Months)", y="Shannon Index", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #scale_y_continuous(limits = c(0,max(alpha_meta_final$Observed)), breaks =c(0,1, 2, 3, max(alpha_meta_final$Observed)))+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  ) 


shannon_linetest = ggplot(alpha_meta_final, aes(x=Time_num, y=Shannon,group=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = alpha_meta_final)+
  labs(tag = "A", x="Time (Months)", y="Shannon Index", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 14, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  facet_wrap(~Participant)+
  #scale_y_continuous(limits = c(0,max(alpha_meta_final$Observed)), breaks =c(0,1, 2, 3, max(alpha_meta_final$Observed)))+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  ) 


alpha_meta_final %>% 
  filter((Time_Point == "BL" | Time_Point == "48H")  & LIBRA_Half == "Lower") %>% 
  summarise(value = mean(Shannon))
alpha_meta_final %>% 
  filter((Time_Point == "BL" | Time_Point == "48H")  & LIBRA_Half == "Lower") %>% 
  summarise(value = sd(Shannon))

t.test(Shannon ~ LIBRA_Half, data = alpha_meta_final %>%  filter((Time_Point == "BL" | Time_Point == "48H") ))

t.test(Shannon ~ LIBRA_Half, data = alpha_meta_final %>%  filter((Time_num == 12 |Time_num == 11)  ))


alpha_meta_final %>% 
  filter((Time_Point == "BL" | Time_Point == "48H")  & LIBRA_Half == "Upper") %>% 
  summarise(value = mean(Shannon))
alpha_meta_final %>% 
  filter((Time_Point == "BL" | Time_Point == "48H")  & LIBRA_Half == "Upper") %>% 
  summarise(value = sd(Shannon))


alpha_meta_final %>% 
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Lower") %>% 
  summarise(value = mean(Shannon))
alpha_meta_final %>% 
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Lower") %>% 
  summarise(value = sd(Shannon))


alpha_meta_final %>% 
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Upper") %>% 
  summarise(value = mean(Shannon))
alpha_meta_final %>% 
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Upper") %>% 
  summarise(value = sd(Shannon))


observed_line = ggplot(alpha_meta_final, aes(x=Time_num, y=Observed,group=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = alpha_meta_final)+
  labs(tag = "A", x="Time (Months)", y="Observed OTUs", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 12, face="bold"),
        plot.tag = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #scale_y_continuous(limits = c(0,max(alpha_meta_final$Observed)), breaks =c(0, 50, 100, max(alpha_meta_final$Observed)))+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Better Score", "Worse Score")
  )



spaghetti_plot = ggarrange(observed_line, shannon_line)

ggsave(paste0(path, "alpha_spaghetti_plot_", Sys.Date(), ".jpeg"), spaghetti_plot, height = 8, width = 15, units="in",device="jpeg")


###LIBRA Score
#meta_and_alpha$Time_bin = factor(meta_and_alpha$Time_bin, levels = c("Base", "2-4m", "5-8m","9-12m"))
shannon_libra_bin = ggplot(alpha_meta_final, aes(x=reorder(Participant, Total_LIBRA_Score), y=Shannon, color=LIBRA_Half))+
  geom_boxplot()+
  geom_point()+
  labs(tag = "B", x="Participants \n(Ordered by LIBRA Score)", y="Shannon Index", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_x_reverse()+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )


observed_libra_bin = ggplot(alpha_meta_final, aes(x=reorder(Participant, Total_LIBRA_Score), y=Observed, color=LIBRA_Half))+
  geom_boxplot()+
  geom_point()+
  labs(tag = "B", x="Participants \n(Ordered by LIBRA Score)", y="Observed OTUs", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 12, face="bold"),
        plot.tag = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  #scale_x_reverse()+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Better Score", "Worse Score")
  )

alpha_libra_bin_plot = ggarrange(shannon_libra_bin, observed_libra_bin)

ggsave(paste0(path, "alpha_libra_bin_plot_", Sys.Date(), ".jpeg"), alpha_libra_bin_plot, height = 6, width = 10, units="in",device="jpeg")

shannon_merged_spaghetti_box = ggarrange(shannon_line, shannon_libra_bin)
otu_merged_spaghetti_box = ggarrange(observed_line, observed_libra_bin)

ggsave(paste0(path, "shannon_merged_spaghetti_box_", Sys.Date(), ".jpeg"), shannon_merged_spaghetti_box, height = 6, width = 14, units="in",device="jpeg")
ggsave(paste0(path, "otu_merged_spaghetti_box_", Sys.Date(), ".jpeg"), otu_merged_spaghetti_box, height = 8, width = 14, units="in",device="jpeg")


#alpha lmer
library(broom.mixed)
library(modelsummary)
library(lme4)
library(lmerTest)
model_score_shan <- lmer(Shannon ~ Time_num*Total_LIBRA_Score + age+ as.factor(sex) + BMI_Current +as.factor(Ever_Homeless)+ (1 | Participant), data = alpha_meta_final)
#model_half_shan <- lmer(Shannon ~ Time_num*LIBRA_Half + age+sex + Ever_Homeless+ (1 | Participant), data = alpha_meta_final)
summary(model_score_shan)
anova(model_score_shan)
#summary(model_half_shan)
modelsummary(list(model_score_shan, model_half_shan), output = paste0(path, "shannon_model_summary.docx"))

sup_table_shan_model = flextable(tidy(model_score_shan))
sup_table_shan_model = set_caption(sup_table_shan_model, caption = "Supplementary Table x. Mixed Effects Model Shannon Index") %>% autofit()
save_as_docx(sup_table_shan_model, path = paste0(path, "Supplementary_Table_shan_mod._", Sys.Date(), ".docx"))



library(lme4)
library(lmerTest)
model_score_ob <- lmer(Observed ~ Time_num*Total_LIBRA_Score + age+sex +BMI_Current+ Ever_Homeless+ (1 | Participant), data = alpha_meta_final)
model_half_ob <- lmer(Observed ~ Time_num*LIBRA_Half + age+sex + Ever_Homeless+ (1 | Participant), data = alpha_meta_final)
summary(model_score_ob)
summary(model_half_ob)
modelsummary(list(model_score_shan, model_half_shan), output = paste0(path, "obs_model_summary.docx"))

sup_table_obs_model = flextable(tidy(model_score_ob))
sup_table_obs_model = set_caption(sup_table_obs_model, caption = "Supplementary Table x. Mixed Effects Model Observed OTUs") %>% autofit()
save_as_docx(sup_table_obs_model, path = paste0(path, "Supplementary_Table_obs_mod._", Sys.Date(), ".docx"))

calculate_cv <- function(x, na.rm = TRUE) {
  (sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)) * 100
}


cv_by_libra_time <- alpha_meta_final %>%
  group_by(Time_num, LIBRA_Half) %>%
  summarize(
    Mean_Value = mean(Shannon, na.rm = TRUE),
    SD_Value = sd(Shannon, na.rm = TRUE),
    CV_Percentage = calculate_cv(Shannon)
  )
cv_by_libra_time  %>%group_by(LIBRA_Half) %>% summarize(mean(CV_Percentage))

cv_by_libra_half <- alpha_meta_final %>%
  group_by(LIBRA_Half) %>%
  summarize(
    Mean_Value = mean(Shannon, na.rm = TRUE),
    SD_Value = sd(Shannon, na.rm = TRUE),
    CV_Percentage = calculate_cv(Shannon)
  )

# View the results
print(cv_by_group)


#####Form microbiomestat object for volatility tests
data.obj <- mStat_convert_phyloseq_to_data_obj(physeq_object_mindsource_final)

alpha_volatility_test_results <- generate_alpha_volatility_test_long(
  data.obj = data.obj,
  alpha.obj = NULL,
  alpha.name = c("shannon","observed_species"),
  time.var = "Time_num",
  subject.var = "Participant",
  group.var = "LIBRA_Half",
  adj.vars ="Total_LIBRA_Score"
)


alpha_shannon_vol = flextable(alpha_volatility_test_results$shannon)
save_as_docx(alpha_shannon_vol, path = paste0(path,"Supplementary_Table_alpha_shannon_vol_results.docx"))
alpha_obs_vol = flextable(alpha_volatility_test_results$observed_species)
save_as_docx(alpha_obs_vol, path = paste0(path,"Supplementary_Table_alpha_obs_vol_results.docx"))


#beta
#####Beta Diversity from Phyloseq
#Function to create beta diversity measure and figure from phyloseq object
#inputs include phyloseq object, distance measure, weighted true or false, the coloring group, title for plot, outpath, and the adonis p-value
weighted_fun = function (data_set,da_method, weight_T_F, color_group, title_name,path_out){
  weight_ordination = ordinate(data_set, da_method, "unifrac", weighted=weight_T_F)
  p_weight <- plot_ordination(data_set,weight_ordination, type="samples", color=color_group, title=title_name)+ theme_bw()
  p_weight2 = p_weight + 
    stat_ellipse(type = "norm", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 16, face="bold"))+
    labs(color="Brain Health Group \n(LIBRA Score)") +
    #theme(plot.tag.position = c(.885, .98), plot.tag = element_text(size=7.25))+
    scale_colour_manual(
      values = c( "#A77903", 
                  "#1B3558"),
      labels= c("Lower Risk Group", "Higher Risk Group")
    )
  ggsave(path_out, p_weight2,device="jpeg", height = 6, width = 8, units = "in")
  return(p_weight2)
}

?ordinate
weighted_fun_part = function (data_set,da_method, weight_T_F, color_group, title_name,path_out){
  weight_ordination = ordinate(data_set, da_method, "unifrac", weighted=weight_T_F)
  p_weight <- plot_ordination(data_set,weight_ordination, type="samples", color=color_group, title=title_name)+ theme_bw()
  p_weight2 = p_weight + 
    stat_ellipse(type = "norm", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 12, face="bold"))+
    labs(color="Participant")
  #theme(plot.tag.position = c(.885, .98), plot.tag = element_text(size=7.25))+
  ggsave(path_out, p_weight2,device="jpeg", height = 6, width = 8, units = "in")
  return(p_weight2)
}


ps.meta_mindsource = meta(physeq_object_mindsource_final)
#changes data to abundance counts and moves 0 to a low value
phylo_mindsource_beta_filter_tranform = transform_sample_counts(physeq_object_mindsource_final, function(x) x/sum(x))

ps_transform = transform_sample_counts(physeq_object_mindsource_final, function(x) 1E6 * x/sum(x))
ps.meta_mindsource = meta(physeq_object_mindsource_final)
weighted.unifrac.dist = UniFrac(ps_transform ,  weighted = T, normalized = TRUE, parallel = FALSE,fast = TRUE)
Adonis_group_control_one_category <- adonis2(weighted.unifrac.dist~Time_num*Total_LIBRA_Score*BMI_Current*age*as.factor(sex)*as.factor(Ever_Homeless)*Participant, data = ps.meta_mindsource, permutations = 9999)
adonis = Adonis_group_control_one_category$`Pr(>F)`[1]
adonis

unweighted.unifrac.dist = UniFrac(ps_transform ,  weighted = F, normalized = TRUE, parallel = FALSE,fast = TRUE)
unAdonis_group_control_one_category <- adonis2(unweighted.unifrac.dist~Time_num*Total_LIBRA_Score*BMI_Current*age*as.factor(sex)*as.factor(Ever_Homeless)*Participant, data = ps.meta_mindsource, permutations = 9999)
unadonis = unAdonis_group_control_one_category$`Pr(>F)`[1]
unadonis



#Functions for beta, heatmap, and stack bar plot
p_weight_group = weighted_fun(phylo_mindsource_beta_filter_tranform,"PCoA", T, "LIBRA_Half",paste0("Weighted UniFrac LIBRA Groups") ,paste0(path, "/weighted_mindsource_libra_half_",Sys.Date(),".jpeg"))
p_unweight_group = weighted_fun(phylo_mindsource_beta_filter_tranform,"PCoA", F, "LIBRA_Half",paste0("Unweighted UniFrac LIBRA Groups") ,paste0(path, "/unweighted_mindsource_libra_half_",Sys.Date(),".jpeg"))



p_weight_part = weighted_fun_part(phylo_mindsource_beta_filter_tranform,"PCoA", T, "Participant",paste0("Weighted UniFrac Participant") ,paste0(path, "/weighted_mindsource_Participant_",Sys.Date(),".jpeg"))
p_unweight_part = weighted_fun_part(phylo_mindsource_beta_filter_tranform,"PCoA", F, "Participant",paste0("Unweighted UniFrac Participant") ,paste0(path, "/unweighted_mindsource_Participant",Sys.Date(),".jpeg"))

coord_beta_w_part = p_weight_part + labs(tag="B") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
coord_beta_uw_part = p_unweight_part + labs(tag="A") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())

ggsave(paste0(path, "supp1B_beta_w_part_figure_", Sys.Date(), ".jpeg"), coord_beta_w_part, height = 6, width = 8, units="in",device="jpeg")
ggsave(paste0(path, "supp1A_beta_uw_part_figure_", Sys.Date(), ".jpeg"), coord_beta_uw_part, height = 6, width = 8, units="in",device="jpeg")


#will need to edit the time point data to have 48 go to zero
meta_and_alpha_libra_beta_long = meta_and_libra_final %>% mutate(Time_num = case_when(Time_Point == "BL" ~ 0,
                                                                                      Time_Point == "48h" ~ 0,
                                                                                      Time_Point == "2m" ~ 2,
                                                                                      Time_Point == "3m" ~ 3,
                                                                                      Time_Point == "4m" ~ 4,
                                                                                      Time_Point == "5m" ~ 5,
                                                                                      Time_Point == "6m" ~ 6,
                                                                                      Time_Point == "7m" ~ 7,
                                                                                      Time_Point == "8m" ~ 8,
                                                                                      Time_Point == "9m" ~ 9,
                                                                                      Time_Point == "10m" ~ 10,
                                                                                      Time_Point == "11m" ~ 11,
                                                                                      Time_Point == "12m" ~ 12))



sample = data.frame(meta_and_alpha_libra_beta_long)
rownames(sample) = sample$Sample_ID
sampledata = sample_data(sample)
physeq3 = phyloseq(OTU, TAX, sampledata, random_tree2)


data.obj <- mStat_convert_phyloseq_to_data_obj(physeq3)
dist.obj_use = mStat_calculate_beta_diversity(data.obj,
                                              dist.name=c("Unifrac","WUnifrac"))



#beta long
beta_long = generate_beta_change_spaghettiplot_long(
  data.obj = data.obj,
  dist.obj = NULL,
  subject.var = "Participant",
  time.var = "Time_num", 
  t0.level = sort(unique(data.obj$meta.dat$Time_num))[1],
  ts.levels = sort(unique(data.obj$meta.dat$Time_num))[-1],
  group.var = "LIBRA_Half",
  strata.var = NULL,
  dist.name = c("UniFrac","WUniFrac"),
  base.size = 20,
  theme.choice = "bw",
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11, 
  pdf.hei = 8.5
)


ggsave(paste0(path, "beta_long_plot_", Sys.Date(), ".jpeg"), beta_long, height = 10, width = 12, units="in",device="jpeg")

uni_data_long = data.frame(layer_data(beta_long$UniFrac, 1)) %>% select(group, colour, x,y) %>% mutate(LIBRA_Half = case_when(colour=="#E31A1C" ~ "Lower",
                                                                                                                              T~"Upper"))

uni_distance = ggplot(uni_data_long, aes(x=x, y=y,group=LIBRA_Half))+
  geom_line(aes(group=group),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = uni_data_long)+
  labs(tag = "A", x="Time (Months)", y="Unweighted Unifrac\n Distance from Baseline", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c(  "#A77903", 
                 "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )


#'#0047AB' '#006666' '#002D72' '#FFB81C' "#1F3A5F"


#beta lmer
library(lme4)
library(lmerTest)
uni_model <- lmer(y ~ LIBRA_Half*x+ (1 | group), data = uni_data_long)
summary(uni_model)
(anova(uni_model))

sup_table_uni_dist_model = flextable(tidy(uni_model))
sup_table_uni_dist_model = set_caption(sup_table_uni_dist_model, caption = "Supplementary Table x. Mixed Effects Model Unweighted Unifrac Distance") %>% autofit()
save_as_docx(sup_table_uni_dist_model, path = paste0(path, "Supplementary_Table_uni_mod._", Sys.Date(), ".docx"))





wuni_data_long = data.frame(layer_data(beta_long$WUniFrac, 1)) %>% select(group, colour, x,y) %>% mutate(LIBRA_Half = case_when(colour=="#E31A1C" ~ "Lower",
                                                                                                                                T~"Upper"))

wuni_distance = ggplot(wuni_data_long, aes(x=x, y=y,group=LIBRA_Half))+
  geom_line(aes(group=group),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = wuni_data_long)+
  labs(tag = "A", x="Time (Months)", y="Weighted Unifrac\n Distance from Baseline", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c(  "#A77903", 
                 "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )



#beta lmer
library(lme4)
library(lmerTest)
wuni_model <- lmer(y ~ LIBRA_Half*x+ (1 | group), data = wuni_data_long)
summary(wuni_model)
(anova(wuni_model))

sup_table_wuni_dist_model = flextable(tidy(wuni_model))
sup_table_wuni_dist_model = set_caption(sup_table_wuni_dist_model, caption = "Supplementary Table x. Mixed Effects Model Weighted Unifrac Distance") %>% autofit()
save_as_docx(sup_table_wuni_dist_model, path = paste0(path, "Supplementary_Table_wuni_mod._", Sys.Date(), ".docx"))

generate_beta_volatility_test_long(
  data.obj = data.obj,
  dist.obj = NULL,
  subject.var = "Participant",
  time.var = "Time_num", 
  group.var = "LIBRA_Half",
  dist.name = c("UniFrac","WUniFrac"),
)


######Create figure 1 for all diversity 

#alpha spaghetti, alpha box, w coor, w distance, uw coord, uw distane
coord_beta_w_fig1 = p_weight_group + labs(tag="C") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())
coord_beta_uw_fig1 = p_unweight_group + labs(tag="E") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())
dist_beta_w_fig1 = wuni_distance + labs(tag="D") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())
dist_beta_uw_fig1 = uni_distance+ labs(tag="F") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())

beta_dista = ggarrange(wuni_distance,uni_distance, ncol=2,nrow=1)

figure_1 = ggarrange(shannon_line, coord_beta_w_fig1,coord_beta_uw_fig1, shannon_libra_bin,dist_beta_w_fig1,dist_beta_uw_fig1, ncol=3,nrow=2)
ggsave(paste0(path, "figure_1_merged_", Sys.Date(), ".jpeg"), figure_1, height = 10, width = 22, units="in",device="jpeg")
ggsave(paste0(path, "beta_dist_figure_pres_merged_", Sys.Date(), ".jpeg"), beta_dista, height = 6, width = 12, units="in",device="jpeg")





########Differential abundance work
#maaslin
meta_and_libra_final = as.data.frame(meta_and_libra_final)
meta_and_libra_final$Ever_Homeless = factor(meta_and_libra_final$Ever_Homeless, levels=c("No", "Unknown", "Yes"))
rownames(meta_and_libra_final) = meta_and_libra_final$Sample_ID

ps_genus <- tax_glom(physeq_object_mindsource_final, "Genus")
#pull out taxa table into separate data frame
otu_ps <- otu_table(physeq_object_mindsource_final)
otu_ps <- t(otu_ps)
df_maslin <- as.data.frame(otu_ps)
colnames(df_maslin) <- as.data.frame(tax_table(physeq_object_mindsource_final))$Species



de_graph_fun = function(data_set, title, meta){
  de_data = data_set %>% 
    filter(qval_individual<0.3 & name == meta)%>%
    mutate(Up_Down = case_when(coef > 0 ~ "Up",
                               coef<  0  ~ "Down")) %>%
    arrange((qval_individual))%>% 
    mutate(feature = str_remove(feature, ".*?(?=s__)")) %>% 
    mutate(feature = str_replace(feature, "s__", "")) %>% 
    mutate(feature = str_replace(feature, "_", " "))
  
  de_plot = ggplot(de_data, aes(x=reorder(feature,coef), y=coef, fill=Up_Down)) +
    geom_bar(stat='identity') +
    theme_bw()+
    theme(legend.position = "none",
          text = element_text(size = 16, face="bold"))+
    #scale_x_discrete(position = "top")+
    theme(axis.text.y = element_text(face = "bold.italic")) + # Assign specific colors to groups
    coord_flip()+
    labs(#title="Species LFC by LIBRA Score", 
      x="Species", 
      y="Log Fold Change",
      title= title) +
    scale_fill_manual(values = c("Up" = "#1B3558", "Down" = "#A77903"))
  
  #if(length(de_data$name > 1)){
  # de_plot = de_plot + facet_wrap(~name, scales="free", nrow=3)
  #}else
  #{
  #  de_plot = de_plot
  #}
  
  return(de_plot)
}

library(maaslin3)
fit_data = maaslin3(
  input_data = df_maslin,
  input_metadata = meta_and_libra_final,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("/Path/to/Output",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~Time_num*Total_LIBRA_Score + Age+Sex + Ever_Homeless + BMI_Current+(1|Participant)",
  min_prevalence = 0.25,
  min_abundance = 10E-6,
  plot_associations = F,
  small_random_effects=F)


fit_data_abund = maaslin3(
  input_data = taxa_table_all_filtered_mat,
  input_metadata = meta_and_libra_final,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("/path/to/output", Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~Time_num*LIBRA_Score_Total + Age+Sex+ any_homeless +(1|Participant)",
  min_prevalence = 0.25,
  min_abundance = 0,
  plot_associations = F,
  small_random_effects=F,                      
  plot_summary_plot = F,
  correction="BH",
  transform = "PLOG",
  evaluate_only = "abundance",
  warn_prevalence = F,
  zero_threshold = -1,
  median_comparison_abundance = F)


de_abund_full_time = de_graph_fun(fit_data$fit_data_abundance$results, "Differential Expression by Long Time for Mindsource, Abundance", "Time_num")
de_prevfull_time = de_graph_fun(fit_data$fit_data_prevalence$results, "Differential Expression by Long Time for Mindsource, Prevelence", "Time_num")
de_abund_full_libra = de_graph_fun(fit_data$fit_data_abundance$results, "Differential Expression by Long Time for Mindsource, Abundance", "Total_LIBRA_Score")
de_prevfull_libra = de_graph_fun(fit_data$fit_data_prevalence$results, "Differential Expression by Long Time for Mindsource, Prevelence", "Total_LIBRA_Score")
de_abund_full_time_libra = de_graph_fun(fit_data$fit_data_abundance$results, "Differential Expression of Time x LIBRA Score, Abundance", "Time_num:Total_LIBRA_Score")
de_prevfull_time_libra = de_graph_fun(fit_data$fit_data_prevalence$results, "Differential Expression of Time x LIBRA Score, Prevelence", "Time_num:Total_LIBRA_Score")


ggsave(paste0(path, "DE_long_mindsource_time_abundance_full_plot_", Sys.Date(), ".jpeg"), de_abund_full_time, height = 10, width = 10, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_time_prev_full_plot_", Sys.Date(), ".jpeg"), de_prevfull_time, height = 10, width = 10, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_libra_abundance_full_plot_", Sys.Date(), ".jpeg"), de_abund_full_libra, height = 10, width = 10, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_libra_prev_full_plot_", Sys.Date(), ".jpeg"), de_prevfull_libra, height = 10, width = 10, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_time_libra_abundance_full_plot_", Sys.Date(), ".jpeg"), de_abund_full_time_libra, height = 6, width = 6, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_time_libra_prev_full_plot_", Sys.Date(), ".jpeg"), de_prevfull_time_libra, height = 6, width = 6, units="in",device="jpeg")


abund_list = c("Clostridium clostridioforme", 'Akkermansia muciniphila', "Blautia hydrogenotrophica", "Sutterella wadsworthensis")

prev_list = c('Bilophila wadsworthia', 'Roseburia intestinalis', 'Collinsella aerofaciens', 'Eubacterium eligens')

scfa_list = c("Faecalibacterium prausnitzii", "Roseburia intestinalis", "Eubacterium rectale", "Eubacterium hallii",
              "Subdoligranulum variabile", "Roseburia hominis", "Anaerobutyricum soehngenii", "Coprococcus comes",
              "Coprococcus catus", "Butyricicoccus pullicaecorum", "Agathobacter rectalis", "Bacteroides thetaiotaomicron",
              "Bacteroides vulgatus", "Bacteroides uniformis", "Bacteroides fragilis", "Parabacteroides distasonis",
              "Megasphaera elsdenii", "Veillonella parvula", "Veillonella dispar", "Ruminococcus obeum", "Blautia obeum",
              "Bifidobacterium longum", "Bifidobacterium adolescentis", "Bifidobacterium breve", "Lactobacillus rhamnosus",
              "Lactobacillus plantarum", "Akkermansia muciniphila", "Ruminococcus gnavus")

# Create a regex pattern using the pipe '|' (which means "OR" in regex)
scfa_pattern <- paste(scfa_list, collapse = "|")

#Abundance plots
libra_filter_abd = fit_data$fit_data_abundance$results %>%
  filter(qval_individual<0.3 & name == "Time_num:Total_LIBRA_Score")%>%
  mutate(Up_Down = case_when(coef > 0 ~ "Up",
                             coef<  0  ~ "Down")) %>%
  arrange((qval_individual))%>%
  mutate(feature = str_remove(feature, ".*?(?=s__)")) %>%
  mutate(feature = str_replace(feature, "s__", "")) %>%
  mutate(feature = str_replace(feature, "_", " "))

libra_filter_prev = fit_data$fit_data_prevalence$results %>%
  filter(qval_individual<0.3 & name == "Time_num:Total_LIBRA_Score")%>%
  mutate(Up_Down = case_when(coef > 0 ~ "Up",
                             coef<  0  ~ "Down")) %>%
  arrange((qval_individual))%>%
  mutate(feature = str_remove(feature, ".*?(?=s__)")) %>%
  mutate(feature = str_replace(feature, "s__", "")) %>%
  mutate(feature = str_replace(feature, "_", " "))


sup_table_mal_abd = flextable(libra_filter_abd %>% dplyr::select(feature, name, coef, stderr,pval_individual, qval_individual))
sup_table_mal_abd = set_caption(sup_table_mal_abd, caption = "Supplementary Table x. Significant Species Abundance Associated with LIBRA x Time") %>% autofit()
save_as_docx(sup_table_mal_abd, path = paste0(path, "Supplementary_Table_libra_time_abnd._", Sys.Date(), ".docx"))


sup_table_mal_prev = flextable(libra_filter_prev %>% dplyr::select(feature, name, coef, stderr,pval_individual, qval_individual))
sup_table_mal_prev = set_caption(sup_table_mal_prev, caption = "Supplementary Table x. Significant Species Prevalence Associated with LIBRA x Time") %>% autofit()
save_as_docx(sup_table_mal_prev, path = paste0(path, "Supplementary_Table_libra_time_prev._", Sys.Date(), ".docx"))



libra_filter_prev$feature
#####creaate relative abundance 
physeq_object_mindsource_final_relative <- transform_sample_counts(physeq_object_mindsource_final, function(x) x / sum(x))
OTU_df <- as.data.frame(as(otu_table(physeq_object_mindsource_final_relative), "matrix"))


OTU_df$feature = rownames(OTU_df)
taxa_table_all_df2 = OTU_df %>%
  mutate(feature = str_remove(feature, ".*?(?=s__)")) %>%
  mutate(feature = str_replace(feature, "s__", "")) %>%
  mutate(feature = str_replace(feature, "_", " ")) %>%
  mutate(feature = str_remove(feature, "\\..*"))
#taxa_table_all_df_filter = taxa_table_all_df2 %>% filter(feature %in% libra_filter_abd$feature)
#taxa_table_all_df_filter = taxa_table_all_df2 %>% filter(str_detect(feature, scfa_pattern))


orgo_list = data.frame(taxa_table_all_df2$feature)
#rownames(taxa_table_all_df_filter) = taxa_table_all_df_filter$feature
relative_taxa_adbund_filtered_update2 <- data.frame(t(taxa_table_all_df2[-140]))
#colnames(relative_taxa_adbund_filtered_update2) <- taxa_table_all_df_filter[, 1]
relative_taxa_adbund_filtered_update3 <- relative_taxa_adbund_filtered_update2 %>%
  rownames_to_column(var = "Sample_ID")
relative_taxa_adbund_filtered_update3_melt = reshape2::melt(relative_taxa_adbund_filtered_update3, id.var="Sample_ID", variable.name="Species", value.name = "Abundance")
relative_taxa_adbund_filtered_update3_melt2 = relative_taxa_adbund_filtered_update3_melt %>%
  mutate(Species = str_remove(Species, ".*?(?=s__)")) %>%
  mutate(Species = str_replace(Species, "s__", "")) %>%
  mutate(Species = str_replace(Species, "_", " ")) %>%
  mutate(Species = str_remove(Species, "\\..*")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "X", "")) %>%
  dplyr::filter(Sample_ID %in% meta_and_libra_final$Sample_ID)

relative_taxa_adbund_filtered_update3_melt3 = relative_taxa_adbund_filtered_update3_melt2 %>%
  group_by(Sample_ID, Species) %>%
  mutate(total_abd = sum(as.numeric(Abundance))) %>%
  ungroup() %>%
  select(-Abundance) %>%
  distinct()
all_taxa_meta_data = merge(relative_taxa_adbund_filtered_update3_melt3, meta_and_libra_final, by="Sample_ID")
all_taxa_meta_data_update = all_taxa_meta_data %>% mutate(Prev = case_when(total_abd >0 ~ 1,
                                                                           T~0))

length(unique(meta_and_libra_final$Sample_ID))


sig_abnd_df = all_taxa_meta_data_update %>% filter(Species %in% abund_list)

abund_fig = ggplot(sig_abnd_df, aes(y=total_abd,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, method="lm",data = sig_abnd_df)+
  labs(tag = "A", x="Time (Months)", y="Relative Abundance", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  facet_wrap(~Species, scales="free") +
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)

sig_abnd_df_supp = all_taxa_meta_data_update %>% filter(Species %in% libra_filter_abd$feature)

abund_fig_supp = ggplot(sig_abnd_df_supp, aes(y=total_abd,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, method="lm",data = sig_abnd_df_supp)+
  labs( x="Time (Months)", y="Relative Abundance", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  facet_wrap(~Species, scales="free") +
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)


sig_prev_df = all_taxa_meta_data_update %>% filter(Species %in% prev_list)

prev_plot = ggplot(sig_prev_df, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  #scale_size_continuous(limits=c(0, 1), breaks=seq(0, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  labs(tag = "B", x="Time (Months)", y="Participant", color="Brain Health Group \n(LIBRA Score)", size="Prevalence") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  scale_size(range = c(1, 3))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(group=LIBRA_Half)) +
  theme(strip.text.x = element_text(face = "italic"))+ # For x-axis facet titles)
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )


sig_prev_df_supp = all_taxa_meta_data_update %>% filter(Species %in% libra_filter_prev$feature)

prev_plot_supp = ggplot(sig_prev_df_supp, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  #scale_size_continuous(limits=c(0, 1), breaks=seq(0, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  labs( x="Time (Months)", y="Participant", color="Brain Health Group \n(LIBRA Score)", size="Prevalence") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  scale_size(range = c(1, 3))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(group=LIBRA_Half)) +
  theme(strip.text.x = element_text(face = "italic"))+ # For x-axis facet titles)
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )



####Figure 3
de_abund_fig3 = de_abund_full_time_libra + labs(tag="A") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank(),    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
de_prev_fig3 = de_prevfull_time_libra + labs(tag="C") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank(),    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
species_abund_fig3 = abund_fig + labs(tag="B") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
prev_abund_fig3 = prev_plot+ labs(tag="D") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
ggsave(paste0(path, "figure_3b_abund_", Sys.Date(), ".jpeg"), abund_fig, height = 4, width = 6, units="in",device="jpeg")
ggsave(paste0(path, "figure_3d_prev_", Sys.Date(), ".jpeg"), prev_plot, height = 6, width = 10, units="in",device="jpeg")


species_abund_supp = abund_fig_supp +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
prev_abund_supp = prev_plot_supp +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
ggsave(paste0(path, "figure_supp_abund_", Sys.Date(), ".jpeg"), species_abund_supp, height = 10, width = 18, units="in",device="jpeg")
ggsave(paste0(path, "figure_supp_prev_", Sys.Date(), ".jpeg"), prev_abund_supp, height = 12, width = 14, units="in",device="jpeg")


figure_3 = ggarrange(de_abund_fig3,species_abund_fig3, de_prev_fig3, prev_abund_fig3, ncol=2,nrow=2)
ggsave(paste0(path, "figure_3_merged_", Sys.Date(), ".jpeg"), figure_3, height = 14, width = 18, units="in",device="jpeg")


####supplemental figure for time and libra seperate
####Figure 3
de_abund_time_supp = de_abund_full_time + labs(tag="A") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
de_prev_time_supp = de_prevfull_time + labs(tag="B") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
de_abund_libra_supp = de_abund_full_libra + labs(tag="C") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
de_prev_libra_supp = de_prevfull_libra+ labs(tag="D") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())



figure_supp_time_and_libra = ggarrange(de_abund_time_supp, de_prev_time_supp,de_abund_libra_supp, de_prev_libra_supp, ncol=2,nrow=2)
ggsave(paste0(path, "figure_supplemental_time_and_libra_merged_", Sys.Date(), ".jpeg"), figure_supp_time_and_libra, height = 10, width = 10, units="in",device="jpeg")




#####Start of SCFA Bacteria abundance figures
all_taxa_meta_data_update_scfa = all_taxa_meta_data_update %>% filter(Species %in% scfa_list) %>% group_by(Participant, Time_num) %>% mutate(scfa_abund = sum(total_abd))

all_taxa_meta_data_update_scfa_2 = all_taxa_meta_data_update_scfa %>% select(Participant, age, sex, Ever_Homeless, Time_num, scfa_abund, LIBRA_Half,Total_LIBRA_Score, BMI_Current) %>% distinct()
all_taxa_meta_data_update_scfa_2 %>%
  filter(Time_num == 12 | Time_num == 11) %>%
  group_by(LIBRA_Half) %>%
  summarise(sd_scfa_abund = sd(scfa_abund, na.rm = TRUE))

scfa_abund_fig = ggplot(all_taxa_meta_data_update_scfa_2, aes(y=scfa_abund,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  #geom_smooth(aes(color=LIBRA_Half),se=F, method="lm",data = all_taxa_meta_data_update)+
  geom_smooth(aes(color=LIBRA_Half),se=F,data = all_taxa_meta_data_update_scfa_2)+
  labs(tag = "A", x="Time (Months)", y="Summed Relative Abundance of \nSCFA Associated Bacteria", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #facet_wrap(~Species, scales="free") +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
ggsave(paste0(path, "scfa_spaghetti_plot_", Sys.Date(), ".jpeg"), scfa_abund_fig, height = 8, width = 15, units="in",device="jpeg")

scfa_abund_fig2 = ggplot(all_taxa_meta_data_update_scfa_2, aes(y=scfa_abund,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  #geom_point()+
  facet_wrap(~Participant)+
  #geom_smooth(aes(color=LIBRA_Half),se=F, method="lm",data = all_taxa_meta_data_update)+
  geom_smooth(aes(color=LIBRA_Half),se=F,data = all_taxa_meta_data_update_scfa_2)+
  labs(tag = "A", x="Time (Months)", y="Summed Relative Abundance of \nSCFA Associated Bacteria", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 14, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #facet_wrap(~Species, scales="free") +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)





scfa_abund_bin = ggplot(all_taxa_meta_data_update_scfa_2, aes(x=reorder(Participant, Total_LIBRA_Score), y=scfa_abund, color=LIBRA_Half))+
  geom_boxplot()+
  geom_point()+
  labs(tag = "B", x="Participants \n(Ordered by LIBRA Score)", y="Summed Relative Abundance of \nSCFA Associated Bacteria", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  #scale_x_reverse()+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )

ggsave(paste0(path, "scfa_box_plot_", Sys.Date(), ".jpeg"), scfa_abund_bin, height = 8, width = 15, units="in",device="jpeg")

figure_2= ggarrange(scfa_abund_fig, scfa_abund_bin)
ggsave(paste0(path, "figure_2_merged_", Sys.Date(), ".jpeg"), figure_2, height = 6, width = 14, units="in",device="jpeg")

library(lme4)
library(lmerTest)

model_score_scfa <- lmer(scfa_abund ~ Time_num*Total_LIBRA_Score + age+sex +BMI_Current+ Ever_Homeless+ (1 | Participant), data = all_taxa_meta_data_update_scfa_2)
summary(model_score_scfa)
anova(model_score_scfa)


all_taxa_meta_data_update_scfa_2_Test = all_taxa_meta_data_update_scfa_2 %>% dplyr::filter(!(Participant=="P24" & Time_num == 12))
model_score_scfa <- lmer(scfa_abund ~ Time_num*Total_LIBRA_Score + age+sex +BMI_Current+ Ever_Homeless+ (1 | Participant), data = all_taxa_meta_data_update_scfa_2)
summary(model_score_scfa)
anova(model_score_scfa)


sup_table_scfab_model = flextable(tidy(model_score_scfa))
sup_table_scfab_model = set_caption(sup_table_scfab_model, caption = "Supplementary Table x. Mixed Effects Model SCFA-Bacteria Abundance Association LIBRA x Time") %>% autofit()
save_as_docx(sup_table_scfab_model, path = paste0(path, "Supplementary_Table_wuni_mod._", Sys.Date(), ".docx"))




model_half_scfa <- lmer(scfa_abund ~ Time_num*LIBRA_Half + age+sex + Ever_Homeless+ (1 | Participant), data = all_taxa_meta_data_update_scfa_2)
summary(model_half_scfa)
modelsummary(model_score_scfa, output = paste0(path, "scfa_model_summary.docx"))


cv_by_libra_time <- all_taxa_meta_data_update_scfa_2 %>%
  ungroup() %>%
  group_by(Time_num, LIBRA_Half) %>%
  summarize(
    Mean_Value = mean(scfa_abund, na.rm = TRUE),
    SD_Value = sd(scfa_abund, na.rm = TRUE),
    CV_Percentage = calculate_cv(scfa_abund)
  )

cv_by_libra_time  %>%group_by(LIBRA_Half) %>% summarize(mean(CV_Percentage))

cv_by_libra_half <- all_taxa_meta_data_update_scfa_2 %>%
  group_by(LIBRA_Half) %>%
  summarize(
    Mean_Value = mean(scfa_abund, na.rm = TRUE),
    SD_Value = sd(scfa_abund, na.rm = TRUE),
    CV_Percentage = calculate_cv(scfa_abund)
  )


all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 0| Time_num == 0.07 )  & LIBRA_Half == "Lower") %>% 
  summarise(value = mean(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 0| Time_num == 0.07 )  & LIBRA_Half == "Lower") %>% 
  summarise(value = sd(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 0 | Time_num == 0.07)  & LIBRA_Half == "Upper") %>% 
  summarise(value = mean(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 0 | Time_num == 0.07)  & LIBRA_Half == "Upper") %>% 
  summarise(value = sd(scfa_abund))


all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Lower") %>% 
  summarise(value = mean(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Lower") %>% 
  summarise(value = sd(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Upper") %>% 
  summarise(value = mean(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Upper") %>% 
  summarise(value = sd(scfa_abund))




abund_fig_test = ggplot(all_taxa_meta_data_update_scfa_2 , aes(y=scfa_abund,x=BMI_Current))+
  #geom_line(aes(group=Participant),color="grey")+
  geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(se=F, method="lm",data = all_taxa_meta_data_update_scfa_2) +#%>% filter(Time_num < 1  ))+
  labs(tag = "A", x="LIBRA Score", y="Relative Abundance")+#, color="Brain Health Group \n(LIBRA Score)") +
  theme(plot.tag = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  stat_cor(label.y.npc="top", label.x = 30, size = 3)

#facet_wrap(~Species, scales="free") +
#scale_colour_manual(
# values = c( "#4F94CD",  "#CD3333"),
#labels= c("Better Score", "Worse Score")
#)+
#theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)


all_taxa_meta_data_update_scfa_2_test = all_taxa_meta_data_update_scfa_2 %>% filter(Time_num > 6 & Total_LIBRA_Score < 10)
model_score_scfa_test <- lm(scfa_abund ~ Total_LIBRA_Score + age+sex +BMI_Current+ Ever_Homeless, data = all_taxa_meta_data_update_scfa_2)
summary(model_score_scfa_test)
(anova(model_score_scfa_test))



#####Pathways
#taxa_table_1 = read.table('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0001_DeepSeq_Jul_2023_Diversigen/filtered functional tables/uniref-pathway-filtered-absolute.tsv', header =T, sep="\t", check.names = F)
path_table_1
path_table_1 <- path_table_1[-c(1, 2), ]
path_table_2
path_table_2 <- path_table_2[-c(1, 2), ]
path_table_all = merge(path_table_1,path_table_2, by=0, all=T)
path_table_all[is.na(path_table_all)] <- 0
rownames(path_table_all) = path_table_all$Row.names
path_table_all = path_table_all %>% dplyr::select(-Row.names)
# Remove the first two rows
#taxa_table_all <- taxa_table_all[-c(1, 2), ]
path_table_all_df = data.frame(path_table_all)
names(path_table_all_df) <- gsub("X", "", names(path_table_all_df))
path_table_all_df_filter <- path_table_all_df[, names(path_table_all_df) %in% meta_and_libra_final$Sample_ID]
#meta_and_alpha %>% filter(SampleID == "19-1307.P24.8m-g")
#dat[-1] <-sweep(dat[-1], 2, colSums(dat[,-1]), `/`) * 10000
#Round data so that all the data are integers


BiocManager::install("edgeR")
library(edgeR)

?pData
library(DESeq2)
library(splines)




pathway_df = data.frame(path_table_all_df_filter)
dat_pathway <- pathway_df %>%
  dplyr::mutate_if(is.numeric, round)
names(dat_pathway) <- gsub("X", "", names(dat_pathway))
keep_genes <- colSums( dat_pathway > 100 ) >= 70
dat_pathway <- dat_pathway[keep_genes,]
dat_pathway_2 <- as.matrix(dat_pathway)
d0 <- DGEList(dat_pathway)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
meta_and_libra_final

X <- ns(meta_and_libra_final$Time_num, df=3)
X2 = ns(meta_and_libra_final$Total_LIBRA_Score, df=3)
design <- model.matrix(~Time_num*Total_LIBRA_Score+age+as.factor(sex)+BMI_Current + as.factor(Ever_Homeless), (meta_and_libra_final))

y <- voom(d, design, plot = T)

y <- voom(d, design)
corfit <- duplicateCorrelation(y, design, block = meta_and_libra_final$Participant)
y <- voom(d, design, block = meta_and_libra_final$Participant, correlation =
            corfit$consensus, plot=T)

fit <- lmFit(y, design,block = meta_and_libra_final$Participant, correlation =
               corfit$consensus )
fit <- eBayes(fit)
top.table = topTable(fit, coef= "Time_num:Total_LIBRA_Score", n=Inf)
top.table2 = top.table %>% filter(adj.P.Val<0.3)
which(top.table$adj.P.Val < 0.3)


BiocManager::install("variancePartition")

library("variancePartition")
detach("package:lmer", unload=T,character.only = TRUE )
library("BiocParallel")
library(lmerTest)
param <- SnowParam(4, "SOCK", progressbar = TRUE)
# The variable to be tested must be a fixed effect
form <- (~ age + (1 | Participant))

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(d, form, meta_and_libra_final)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation
fitmm <- dream(vobjDream, form, meta_and_libra_final)
fitmm <- eBayes(fitmm)

# Examine design matrix
head(fitmm$design, 3)



#update_dat_pathway = dat_pathway_2 + 1
#align data
ind <- match(x = colnames(dat_pathway_2), table = rownames(meta_and_libra_final))
meta1 <- meta_and_libra_final[ind,]
meta1 = data.frame(meta1)
rownames(meta1) = meta1$Sample_ID
meta1[1:10, 1:10]
#Check that data are aligned
all(colnames(dat_pathway_2) == rownames(meta1))
#make deseq object
meta1
library(caret)
#scale and center data
meta1_scale_factors = subset(meta1, select = c( Time_num, Total_LIBRA_Score, age, sex,BMI_Current))
meta1_scaled  = scale(meta1_scale_factors, scale = T, center = T)
meta1_scaled = data.frame(meta1_scaled)
meta1_scaled$Ever_Homeless = as.factor(meta1$Ever_Homeless)
meta1_scaled$Participant = (meta1$Participant)

#Create deseq object
dds <- DESeqDataSetFromMatrix(countData = dat_pathway_2, colData = meta1_scaled, design = ~ Participant + age+ sex+ BMI_Current + Ever_Homeless + Time_num*Total_LIBRA_Score)


#generate normalized counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
#QC for DE analysis using DESeq2
#If input file is very large > 800 MB then, use the below command
rld <- varianceStabilizingTransformation(dds, blind=T)
#differential expression analysis
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds)
#shrinkage factor to normalize dataa
dds_shrink = lfcShrink(dds = dds, res = res, coef="Time_num.Total_LIBRA_Score")
summary(res)
summary(dds_shrink)
summary(dds)
plotMA(res, ylim=c(-2,2))
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(top.table))
# Extract significant results
signif_res <- dds_shrink[dds_shrink$padj < 0.3 & dds_shrink$pvalue < 0.05 & !is.na(dds_shrink$padj), ]
plotMA(signif_res, ylim=c(-4,4))
#This will be empty if no genes were significant
signif_genes <- as.character(rownames(top.table2))
#Below is old code to export sig-taxa
sigtab = cbind(as(top.table, "data.frame"), all_genes)


##Function to create volcano plot
volcplot <- function(data, padj_threshold = 0.3, fc = 1, plot_title = 'Volcano Plot', plot_subtitle = NULL, genelist_vector = genes$all_genes, genelist_filter = FALSE) {
  # Set the fold-change thresholds
  neg_log2fc <- -log2(fc)
  pos_log2fc <- log2(fc)
  # Make a dataset for plotting, add the status as a new column
  plot_ready_data <- data %>%
    mutate_at('adj.P.Val', ~replace(.x, is.na(.x), 1)) %>%
    mutate_at('logFC', ~replace(.x, is.na(.x), 0)) %>%
    mutate(
      log2fc_threshold = ifelse(logFC >= pos_log2fc & adj.P.Val <= padj_threshold, 'up',
                                ifelse(logFC <= neg_log2fc & adj.P.Val <= padj_threshold, 'down', 'ns')
      )
    ) %>%
    mutate(all_genes = replace_na(all_genes, 'none'))
  if (genelist_filter) {
    plot_ready_data <- plot_ready_data %>% dplyr::filter(all_genes %in% genelist_vector)
  }
  if(!is.null(genelist_vector)) {
    plot_ready_data <- plot_ready_data %>% mutate(all_genes = ifelse(all_genes %in% genelist_vector & adj.P.Val <= padj_threshold & log2fc_threshold != 'ns', all_genes, ''))
  }
  # Get the number of up, down, and unchanged genes
  up_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'up') %>% nrow()
  down_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'down') %>% nrow()
  unchanged_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'ns') %>% nrow()
  # Make the labels for the legend
  legend_labels <- c(
    str_c('Up: ', up_genes),
    str_c('NS: ', unchanged_genes),
    str_c('Down: ', down_genes)
  )
  # Set the x axis limits, rounded to the next even number
  x_axis_limits <- DescTools::RoundTo(
    max(abs(plot_ready_data$logFC)),
    0.1,
    ceiling
  )
  
  # Set the plot colors
  plot_colors <- c(
    'up' = "#1B3558",
    'ns' = 'gray',
    'down' = "#A77903"
  )
  # Make the plot, these options are a reasonable strting point
  plot <- ggplot(plot_ready_data) +
    geom_point(
      alpha = 0.25,
      size = 1.5
    ) +
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      color = log2fc_threshold,
      label = all_genes
    ) +
    geom_vline(
      xintercept = c(neg_log2fc, pos_log2fc),
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = -log10(padj_threshold),
      linetype = 'dashed'
    ) +
    scale_x_continuous(
      'log2(FC)',
      limits = c(-x_axis_limits, x_axis_limits)
    ) +
    scale_y_continuous('-log10(padj)',
                       limits = c(-0.25,1.5)
    )+
    scale_color_manual(
      values = plot_colors#,
      #labels = legend_labels
    ) +
    #labs(
    #  color = str_c(fc, '-fold, padj ≤', padj_threshold),
    #  title = plot_title,
    #  subtitle = plot_subtitle
    #) +
    theme_bw(base_size = 24) +
    theme(text = element_text(size = 14, face="bold"),
          legend.position = "none",
          aspect.ratio = 1,
          axis.text = element_text(color = 'black'),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),  # Reduces dead area around legend
          legend.spacing.x = unit(0.2, 'cm'),
          legend.title = element_blank()
    )
  # Add gene labels if needed
  if (!is.null(genelist_vector)) {
    plot <- plot +
      geom_label_repel(
        size = 2.5,
        force = 1,
        max.overlaps = 100000,
        nudge_x = 0,
        nudge_y = 0,
        segment.color = 'black',
        min.segment.length = 0,
        show.legend = FALSE
      )
  }
  plot
}
?geom_label_repel
top.table

sigtab2 = na.omit(sigtab)
rownames(sigtab2) <- NULL
#pulls out genes for labeling
genes = sigtab2 %>% dplyr::filter(logFC > log2(1) | logFC< -log2(1)) %>%
  dplyr::filter(adj.P.Val < 0.3 & P.Value < 0.05) %>%
  arrange(desc(abs(logFC))) %>% slice_head(n=6)
out_table = sigtab2 %>% dplyr::filter(adj.P.Val < 0.3 & P.Value < 0.05) %>%
  arrange(desc(logFC))
#creates volcano plot and prints
library(ggrepel)
sigtab2_udpate = sigtab2 %>% rename(padj = adj.P.Val, log2FoldChange = logFC)
volc = volcplot(sigtab2)
volc

ggsave(paste0(path, "figure_4_merged_", Sys.Date(), ".jpeg"), volc, height = 7.5, width = 7.5, units="in",device="jpeg")

#volc = volc + scale_y_continuous(limits=c(0))
out_table_pathways = genes %>% dplyr::select(all_genes, logFC, P.Value, adj.P.Val) %>% arrange(desc(logFC))
out_table_pathways

out_table_pathways %>% filter(log2FoldChange < 0)
genes
sup_table_3 = flextable(out_table_pathways)
sup_table_3 = set_caption(sup_table_3, caption = "Supplementary Table 3. Metabolic Pathways associated with LIBRA Scores") %>% autofit()
save_as_docx(sup_table_3, path = paste0(path, "Supplementary_Table_3._", Sys.Date(), ".docx"))




#MaAslin Extra 
## Supplemental
######Extra for MaAsLin
#MaAsLin analysis for each libra factor
mal_fun= function(effect,ref){
  fit = Maaslin2(input_data = input_data, input_metadata = meta, min_prevalence = .5, 
                 min_abundance = .01, normalization = "TSS", 
                 output = paste0("MaAsLin_Out_species-level_",effect,"_",Sys.Date()), 
                 fixed_effects = effect,
                 reference = ref,
                 plot_scatter = F,
                 plot_heatmap = F,                          
                 correction="BH",
                 analysis_method = "LM",
                 transform = "LOG",
                 max_significance = 0.3)#,
  #reference = "Current_alcohol,No")}
  return(fit)
}

meta$Physical.inactivity_lib.weight = as.factor(meta$Physical.inactivity_lib.weight)
#head(meta$Healthy_diet)
depressed_mal = mal_fun("Depressed", "Depressed,No")

obese_mal = mal_fun("Obese", "Obese,No")

alc_mal = mal_fun("Current_alcohol","Current_alcohol,No")

chd_mal = mal_fun("Coronary.heart.disease", "Coronary.heart.disease,No")

cdd_mal = mal_fun("Chronic_kidney_disease", "Chronic_kidney_disease,No")

diabetes_mal = mal_fun("Diabetes","Diabetes,No")

Hypercholesterolemia_mal = mal_fun("Hypercholesterolemia", "Hypercholesterolemia,No")

Lifetime_smoking_mal = mal_fun("Lifetime_smoking", "Lifetime_smoking,No")

Hypertension_mal = mal_fun("Hypertension","Hypertension,No")

Physical.inactivity_lib.weight_mal = mal_fun("Physical.inactivity_lib.weight", "Physical.inactivity_lib.weight,0")

Healthy_diet_mal = mal_fun("Healthy_diet", "Healthy_diet,Yes")




abundance_plot = ggplot(all_taxa_meta_data_update, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  scale_size_continuous(limits=c(1, 1), breaks=seq(1, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(color=Total_LIBRA_Score)) +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
abundance_plot
all_taxa_meta_data_update

abundance_plot = ggplot(all_taxa_meta_data_update, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  scale_size_continuous(limits=c(1, 1), breaks=seq(1, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(color=Total_LIBRA_Score)) +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
abundance_plot

abundance_plot = ggplot(all_taxa_meta_data_update, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  scale_size_continuous(limits=c(0, 1), breaks=seq(0, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(color=Total_LIBRA_Score)) +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
abundance_plot

prev_plot = ggplot(all_taxa_meta_data_update, aes(y=total_abd, x=Time_num,group=LIBRA_Half, color=LIBRA_Half)) +
  #geom_point(aes(size=Prev)) +
  stat_summary(fun = "mean", geom = "point")+
  #scale_size_continuous(limits=c(10E-10, 0.2), breaks=seq(10E-10, 0.2, by=0.01))+
  facet_wrap(~Species, scales="free") +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(group=LIBRA_Half)) +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
prev_plot

