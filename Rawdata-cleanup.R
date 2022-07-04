## PLACENTA PROJECT - ANALYSIS OF ANNOTATION DATA FROM IHC 

###    1. GENERAL SETUP ----
#### 1.1. Library ----
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(ggplot2)
library(VennDiagram)
library(nVennR)
library(grImport2)
library(rsvg)
library(UpSetR)
library(tidyverse)
library(ComplexUpset)

#### 1.2 Work directory ----
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

#### 1.3 Input Data ----
# Annotation Raw Data
Annotation_Rawdata <- read_excel("Placenta_Rawdata/Annotation_Rawdata.xlsx")%>%
  rename(identifiers = ...1)  #rename column 1
# Meta Data with Patient Information
meta_data <- read_excel("Placenta_Rawdata/patient_info.xlsx")
# List of Trophoblast-specific Proteins (Manually Annotated due to lack of Annotation Data)
tropho_man <- read_excel("Placenta_Rawdata/troph_spec_manual.xlsx") 
# List of Gene Names (corresponding to Antibodies)
gene_names <- read_excel("Placenta_Rawdata/Genes.xlsx")

###    2. DATA PREPERATION ----

#### 2.1 Formatting ----
#change character "NA" to NA
Annotation_Rawdata[ Annotation_Rawdata == "NA" ] <- NA                
cleaned_data <- Annotation_Rawdata  %>%
  #separate first column into four columns
  separate(identifiers,
           c("Observation", "Celltype", "ID", "Time"),
           sep = "_",
           remove = TRUE) %>% 
  #reshape "AB" columns to rows
  gather(Antibody, value, 5:ncol(.)) %>%    
  #spread "Observation" into 3 separate columns
  spread(Observation, value) %>%  
  #change chr to num in columns + change format of Intensity
  mutate(Qt = as.numeric(Qt), 
         Int = ifelse(as.numeric(Int)<=1, "Undistinct", "Distinct"))

#### 2.2 Figure TMA ----
meta_data %>%
  ggplot(aes(x = xcor, y = ycor)) +
  scale_x_continuous(expand = expansion(mult = 0.0005), limits = c(-1, 4)) +
  scale_y_continuous(expand = expansion(mult = 0.1),limits = c(-1, 8)) +
  scale_color_manual(breaks =c("Decidua", "Decidua+EVT", "Villi"),
                     values=c("#CADBBB", "#FFCC99","#9AABEE")) +
  geom_point(size = 13, aes(col = tissue)) + 
  geom_text(aes(label = ID), size =2) +
  labs(x = "", y = "") + theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) 

###    3. QC TISSUE ----
#### 3.1 Find Outliers ----
#count different Celltypes
total_cell <- cleaned_data %>%
  count(Celltype) %>% nrow()
#count number of Antobodies
total_ab <- cleaned_data %>%
  count(Antibody) %>% nrow()

tissue_QC <- cleaned_data %>%
  #summarize patients in a group
  group_by(ID, Antibody) %>%  
  #count cores, where no cells (NA=4) could be found
  summarise(NA_count = sum(is.na(Int))) %>% ungroup() %>%   
  filter(NA_count == total_cell) %>%
  count(ID) %>%  
  #create new column with percentage of empty cores
  mutate(Per_Missing = n/total_ab*100) %>%   
  #filter for those with more than 25% empty (treshold)
  filter(Per_Missing > 25)
#export table
write.table(tissue_QC, file = "Output_tables/Outlier_tissue.txt", sep = "\t", row.names = FALSE)

#### 3.2 Figure QC Tissue ----
TMA_Missing <- cleaned_data %>%
  #summarize patients in a group
  group_by(ID, Antibody) %>%    
  #count cores, where no cells (NA=4) could be found
  summarise(NA_count = sum(is.na(Int))) %>% ungroup() %>%   
  filter(NA_count == total_cell) %>%
  count(ID) %>%  
  #create new column with percentage of empty cores
  full_join(meta_data, by = "ID") %>%
  mutate(Per_Missing = n/total_ab*100) %>%                  
  replace_na(list(n = 0, Per_Missing = 0)) 

TMA_Missing %>%
  ggplot(aes(x = xcor, y = ycor)) +
  scale_x_continuous(expand = expansion(mult = 0.05), limits = c(-1, 4)) +
  scale_y_continuous(expand = expansion(mult = 0.05),limits = c(-1, 8)) +
  geom_point(size = 13, aes(col = Per_Missing)) + 
  geom_text(aes(label = ID), size =2) +
  labs(x = "", y = "") + theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_colour_gradient(low = "#90CE90", high = "#FF99CC",limits = c(0, 100) ) 


#### 3.3 Remove outlier(s) from table ----
filtered_tissue <- cleaned_data %>%                                 
  filter(!ID %in% tissue_QC$ID) 

###    4. QC AB ----
#### 4.1 Find Outliers ----
ab_QC_outliers <- filtered_tissue %>%
  #summarize antibodies into a group
  group_by(Antibody,Int) %>%  
  #count intensity values (NA, undist., distinct) for each AB
  summarise(count = n()) %>% 
  #reshape table, 1 line per antibody
  spread(Int, count) %>% 
  rename(NAInt = "<NA>") %>%
  # replace to numerical values for calculation
  replace_na(list(Distinct = 0, Undistinct = 0, NAInt = 0)) %>%
  #treshold for AB to pass QC: at least 2 distinct stainings, define outliers
  filter(Distinct <2 )

#### 4.2 Figure QC AB ----
ab_QC <- filtered_tissue %>% 
  group_by(Antibody, Int) %>%  
  summarise(count = n()) %>% 
  spread(Int, count) %>% 
  rename(NAInt = "<NA>") %>%
  replace_na(list(Distinct = 0, Undistinct = 0, NAInt = 0)) %>%
  gather(AB_qual, "count", 2:4)
#export table
write.table(ab_QC, file = "Output_tables/Outlier_ab.txt", sep = "\t", row.names = FALSE)
  
#Order Antibodys based on number of distinct staining
AB_order <- ab_QC %>%
  filter(AB_qual == "Distinct") %>% 
  arrange(count) %>%
  pull(Antibody) 

#figure
ab_QC %>%
  mutate(Antibody = factor(Antibody, levels = AB_order)) %>%
  ggplot(aes(x="", y=count, fill=AB_qual)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = c("#90CE90", "#C0C0C0","#FFCCE5")) +
  theme(axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid  = element_blank(),
      text = element_text(size = 3.5),
      ) +
  facet_wrap(~ Antibody) 
 
#### 4.3 Remove Outliers ----
filter_tiss_ab <- filtered_tissue %>% 
  filter(!Antibody %in% ab_QC_outliers$Antibody)
#export final table with data cleaned for tissue & AB
write.table(filter_tiss_ab, file = "Output_tables/superclean.txt", sep = "\t", row.names = FALSE)

### 5. TROPHO SPEC ----
#### 5.1 Find Outliers ----
# remove antibodies, that are not spec. for trophoblasts (= stain in DC)
tropho_unspec <- filter_tiss_ab %>%
  filter(Celltype == "DC") %>%
  count(Antibody, Int) %>% 
  filter(Int == "Distinct") 
#export table
write.table(tropho_unspec, file = "tropho_unspec.txt", sep = "\t")

#### 5.2 Remove Outliers ----
tropho_spec <- filter_tiss_ab %>% 
  filter(!Antibody %in% tropho_unspec$Antibody) %>%
  distinct(Antibody) 
#export final table with data cleaned for tissue & AB
write.table(tropho_spec, file = "Output_tables/tropho_spec.txt", sep = "\t", row.names = FALSE)

#### 5.4 Incorporate Manual Annotation Data ----
tspec_raw <- filter_tiss_ab %>%
# keep only antibodies, that are tropho-specific in both lists
inner_join(tropho_man, by = "Antibody")

### 6. FIGURES: TROPHOBLAST MARKER  ----
#### 6.1 Venn Diagram - Trophoblast specific marker ----
#Filtering for trophoblast-spec antibodies 
AB_Venn <- tspec_raw %>% 
  filter(Int == "Distinct") %>%
  group_by(Antibody, Celltype) %>%
  count(Int)
# create lists Venn
set_SCT <- AB_Venn %>% filter(Celltype == "SCT") %>% pull(Antibody)
set_CT <- AB_Venn %>% filter(Celltype == "CT") %>% pull(Antibody)
set_EVT <- AB_Venn %>% filter(Celltype == "EVT") %>% pull(Antibody)
# Venn plot
grid.newpage()   
myV <- plotVenn(list(SCT = set_SCT, CT = set_CT, EVT= set_EVT), nCycles = 2000 )
showSVG(myV, opacity = 0.3, borderWidth = 2 , outFile = "Output_figures/Venn.tiff", labelRegions = FALSE,
        setColor = c("#99CCFF", "#CC99FF","#FFCC99"))

#### 6.2 Upset Plot - Trophoblast specific marker ----
# create lists 
ElevatedGenes <- list("SCT" = set_SCT, "CT" = set_CT, "EVT"= set_EVT)
# Upset Plot
upset(
  fromList(ElevatedGenes),
  c("SCT","CT", "EVT"),
  queries=list(
    upset_query(set="EVT", fill='#FFCC99'),
    upset_query(set="CT", fill='#CC99FF'),
    upset_query(set="SCT", fill='#99CCFF')),
  base_annotations=list(
    'trophoblast genes'=(
      intersection_size(
        # show all numbers on top of bars
        bar_number_threshold=1, 
        # reduce width of the bars
        width=0.5,)
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')))),
  stripes=upset_stripes(
    geom=geom_segment(size=12),  
    colors=c('grey95', 'white')),
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.3)),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.6))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line())),
  sort_sets='FALSE')

#### 6.3 Lineplot - Timepoint specific marker ----
tspec_raw <- cleaned_data %>%
  inner_join(tropho_man, by = "Antibody") %>%
  filter(ID != 420.1, 
         Celltype%in% c("SCT", "CT"), is.na(Int)==FALSE) %>%
  left_join(gene_names %>%
              select(Antibody, Gene), by = "Antibody")

tspec_raw %>%
  ggplot(aes(x=Time, y=Int, color=Celltype, group=Celltype)) +
  scale_color_manual(breaks =c("SCT", "CT"),
                     values=c("#99CCFF", "#CC99FF")) +
  geom_point() +
  geom_smooth(se=FALSE) +
  facet_wrap(~Gene.x)
