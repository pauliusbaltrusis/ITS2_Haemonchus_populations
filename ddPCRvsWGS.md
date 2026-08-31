# Importing data
```r
library(tidyverse)
library(readr)
library(patchwork)

setwd("C:/Users/pauli/Desktop/SLU/data/ddPCR data pops 1-17")

# Generate a list metafile for all ddPCR .csvs
files_ddPCR<-list.files(pattern = "_Amplitude\\.csv")

# Iterate through the .csvs and save each ggplot under a variable name
## Even better - map_dfr!

combined_ddPCR <- files_ddPCR %>%
  set_names(~ str_remove(.x, "_Amplitude\\.csv$")) %>% ## every vector is a vector with it's own name
  map_dfr(
    ~ readr::read_csv(.x,skip=3, show_col_types=FALSE) %>%
            select(-starts_with("...")),
    .id = "Sample" ) %>% ## looping through every file applying a function to bind rows of each into a big table (_dfr data frame) ## adds a column .id = "Sample" to id where that file came from
  mutate(Sample = gsub("Amlan_Nancy_", "",Sample)) # edit the column source file

combined_ddPCR<- combined_ddPCR %>% group_by(Sample) %>%
  slice_sample(prop = 1) %>%  # Randomize the order of droplets
  mutate(Number=row_number()) %>%
  ungroup()
```
# Generating 1D ddPCR plots
```r
# Need a longer table format with amplitudes and Ch1 and Ch2 column in order to plot FAM and HEX together as subplots
library(ggh4x)

combined_long <- combined_ddPCR %>% 
pivot_longer(cols= c(Ch1Amplitude, Ch2Amplitude),
             names_to="Channel",values_to="Amplitude") %>%
  mutate(Status=case_when(
    Channel=="Ch1Amplitude" & Hc == 1 ~ "Ch1_pos", # When the case is that Channel==Ch1Amplitude and Hc == 1, that row is then named CH1_pos under column Status
    Channel=="Ch1Amplitude" & Hc == 0 ~ "Ch1_neg", 
    Channel=="Ch2Amplitude" & Univ == 1 ~ "Ch2_pos",
    Channel=="Ch2Amplitude" & Univ == 0 ~ "Ch2_neg"
  ))

ddPCR_plots_17<-ggplot(combined_long, aes(x=Number, y=Amplitude, colour=factor(Status))) +
  geom_point(show.legend = FALSE,size = 1, alpha = 0.7, shape = 16)+
  scale_colour_manual(values = c(
   "Ch1_pos" = "#000666", "Ch1_neg" = "#333333",
    "Ch2_pos" = "#006600", "Ch2_neg" = "#333333"
  ), name = "") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Ch1 Amplitude")+
 theme_minimal(base_size = 11)+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        strip.text = element_text(face="bold"))+
  facet_nested_wrap(vars(Sample, Channel), ncol=8, nest_line=element_line(),strip.position = "top")
ddPCR_plots_17

```
# Calculating k (Separability score) for each population
```r
#Calculating k separation score Mean pos amplitude / Mean neg amplitude

sep_scores_fam<- combined_ddPCR %>% group_by(Sample) %>% summarize(mean_pos= mean(Ch1Amplitude[Hc==1], na.rm=TRUE),
                                                            mean_neg= mean(Ch1Amplitude[Hc==0], na.rm=TRUE),
                                                            sd_pos=sd(Ch1Amplitude[Hc==1], na.rm=TRUE),
                                                            sd_neg=sd(Ch1Amplitude[Hc==0], na.rm=TRUE),
                                                            n_pos = sum(Hc==1),
                                                            n_neg = sum(Hc==0),
                                                            k= (mean_pos-mean_neg)/(sd_pos+sd_neg),  # Gerdes et al 2016
                                                            .groups = "drop")

sep_scores_hex<- combined_ddPCR %>% group_by(Sample) %>% summarize(mean_pos= mean(Ch2Amplitude[Univ==1], na.rm=TRUE),
                                                            mean_neg= mean(Ch2Amplitude[Univ==0], na.rm=TRUE),
                                                            sd_pos=sd(Ch2Amplitude[Univ==1], na.rm=TRUE),
                                                            sd_neg=sd(Ch2Amplitude[Univ==0], na.rm=TRUE),
                                                            n_pos = sum(Univ==1),
                                                            n_neg = sum(Univ==0),
                                                            k= (mean_pos-mean_neg)/(sd_pos+sd_neg), # Gerdes et al 2016
                                                            .groups = "drop")



# FAM #

# empirical percentile test leave-one-out
rank_test_fam<- sep_scores_fam$k

target_valJ<-sep_scores_fam$k[sep_scores_fam$Sample=="J"]
target_valM<-sep_scores_fam$k[sep_scores_fam$Sample=="M"]

withoutJ_rank_test_fam<- rank_test_fam[rank_test_fam!=target_valJ]
withoutM_rank_test_fam<- rank_test_fam[rank_test_fam!=target_valM]

# empirical percentile of target relative to others

## sample J
rank_below_J<- sum(withoutJ_rank_test_fam<target_valJ) # how many samples are smaller than the 16th 

p_percentile_J <- (rank_below_J +0.5)/(length(withoutJ_rank_test_fam) + 1) # sum of all values below/less than sample J value + 0.5 / all values without J value (n) + 1.
  # What fraction of the ref distribution (without J sample) does the sample J sit at? At which percentile! In this case its: J sample is at ~15st percentile of the distribution without it. Not extreme
p_percentile_J # 15th percentile
# Two tailed p-value

p_value_J <- 2* min(p_percentile_J, 1 - p_percentile_J) # looking at both tails (whichever is smaller p or 1-p) and doubles it 

p_value_J # p=0.29

## sample M
rank_below_M<- sum(withoutM_rank_test_fam<target_valM)

p_percentile_M <- (rank_below_M +0.5)/(length(withoutM_rank_test_fam) + 1) # Sample M value is at 3rd percentile of the distribution without it

p_percentile_M # 3rd percentile
# Two tailed p-value

p_value_M <- 2* min(p_percentile_M, 1 - p_percentile_M)

p_value_M # p=0.058

# HEX #

rank_test_hex<- sep_scores_hex$k

target_hex_valJ<-sep_scores_hex$k[sep_scores_hex$Sample=="J"]
target_hex_valM<-sep_scores_hex$k[sep_scores_hex$Sample=="M"]

withoutJ_rank_test_hex<- rank_test_hex[rank_test_hex!=target_hex_valJ]
withoutM_rank_test_hex<- rank_test_hex[rank_test_hex!=target_hex_valM]

# J sample
rank_below_J_hex<- sum(withoutJ_rank_test_hex<target_hex_valJ)

p_percentile_J_hex <- (rank_below_J_hex +0.5)/(length(withoutJ_rank_test_hex) + 1) # Sample J value is at 38th percentile

p_percentile_J_hex # 38th percentile

p_value_J_hex <- 2* min(p_percentile_J_hex, 1 - p_percentile_J_hex) 

p_value_J_hex # p=0.76

# M sample
rank_below_M_hex<- sum(withoutM_rank_test_hex<target_hex_valM)

p_percentile_M_hex <- (rank_below_M_hex +0.5)/(length(withoutM_rank_test_hex) + 1) # Sample M value is at 3rd percentile

p_percentile_M_hex

p_value_M_hex <- 2* min(p_percentile_M_hex, 1 - p_percentile_M_hex) 

p_value_M_hex # p=0.06

```
# Plotting WGS HC % composition in samples against ddPCR results 
```r
# Bar plot % HC
WGS_Hc<- bracken_prct_with_unclassified %>% filter(species=="Haemonchus contortus") %>% select(sample,prct) %>% rename(., WGS=prct)

ddPCR_Hc<- read.csv("C:/Users/pauli/Desktop/SLU/data/ddPCR data pops 1-17/Amlan_Nancy_Hc_XJ4126_1_1000_20260326_122548_985.csv") %>% 
  select(Sample.description.1, Target, Copies.20µLWell) %>% 
  rename(., sample=Sample.description.1) %>%
  mutate(sample=gsub(" ", "-", sample)) %>%
  mutate(sample=gsub("XJ*", "Sample_XJ-", sample)) %>%
  pivot_wider(names_from = Target, values_from = Copies.20µLWell ) %>%
  mutate(ddPCR=(Hc/Univ)*100)


wgs_ddpcr_combined_short<-inner_join(by="sample", ddPCR_Hc, WGS_Hc) %>% select(sample, ddPCR, WGS)

wgs_ddpcr_combined<-inner_join(by="sample", ddPCR_Hc, WGS_Hc) %>% select(sample, ddPCR, WGS) %>% 
  pivot_longer(cols=c(ddPCR, WGS ), names_to = "method", values_to = "fraction")

wgs_ddpcr_combined$sample<- factor(wgs_ddpcr_combined$sample, levels = true_order)


bar_plot_ddPCR_WGS<-ggplot(wgs_ddpcr_combined, aes(x=sample, y= fraction, fill = method)) +
  geom_col(position=position_dodge(width=0.8), width = 0.7) +
  labs(x= "sample", y="Haemonchus, %", fill="Method")+
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(1,5)])+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=45, hjust=1))

bar_plot_ddPCR_WGS
```
# Bland-Altmann plot
```r
# Bland-Altmann instead

baltmann<- wgs_ddpcr_combined_short %>% mutate(mean_val= (ddPCR+WGS)/2,
                                               diff_val= ddPCR-WGS)


baltmann$sample<-real_names[match(baltmann$sample, temp_names)]


mean_diff<- mean(baltmann$diff_val)
sd_diff<- sd(baltmann$diff_val)

upper_limit<-mean_diff+ 1.96*sd_diff
lower_limit<-mean_diff-1.96*sd_diff


baltmann_plot<-ggplot(baltmann, aes(x=mean_val, y=diff_val))+
  geom_point(size=3, color="steelblue")+
  geom_hline(yintercept = mean_diff, colour = "black", linetype = "solid")+
  geom_hline(yintercept = upper_limit, colour = "red", linetype = "dashed")+
  geom_hline(yintercept = lower_limit, colour = "red", linetype = "dashed")+
  annotate("text", x=max(baltmann$mean_val), y=mean_diff,
           label= paste0("Mean diff = ", round(mean_diff, 1)),
           vjust=-0.5, hjust=1, size=3.5)+
   annotate("text", x=max(baltmann$mean_val), y=upper_limit,
           label= paste0("+1.96 SD = ", round(upper_limit, 1)),
           vjust=-0.5, hjust=1, size=3.5, color="red")+
   annotate("text", x=max(baltmann$mean_val), y=lower_limit,
           label= paste0("-1.96 SD = ", round(lower_limit, 1)),
           vjust=-0.5, hjust=1, size=3.5, color="red")+
  ggrepel::geom_text_repel(aes(label=sample), size=3)+
  labs(x = "(ddPCR + WGS) / 2",
       y = "ddPCR - WGS") + theme_minimal()

baltmann_plot
```
