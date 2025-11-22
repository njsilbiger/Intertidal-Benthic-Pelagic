### Code for Analysis on Benthic Pelagic Coupling in the Rocky Intertidal###
## Data Collected in Coastal Oregon ##
### Code developed by Nyssa Silbiger ###

### change umol to mmol and mg to gram and fix the fmol C 

### Load libraries   

library(here)
library(tidyverse)
library(janitor)
library(ggtext)
library(patchwork)
library(calecopal)
library(lubridate)
library(lme4)
library(lmerTest)
library(vegan)
library(broom)
library(broom.mixed)
library(ggh4x)
library(brms)
library(scales)
library(emmeans)


### read in data #########
# Benthic Data
BenthicData<-read_csv(here("Data","CommunityData.csv"))

# Tide pool meta data
MetaData<-read_csv(here("Data","TidePoolDescriptions.csv"))

## All Chemistry data
data_all<-read_csv(here("Data","JoinedChemData.csv"))

#Read In CUTI
CT<-read_csv(here("Data","CUTI.csv"))
#m2 s-1
# https://oceanview.pfeg.noaa.gov/products/upwelling/cutibeuti

#CUTI from
#Jacox, M. G., C. A. Edwards, E. L. Hazen, and S. J. Bograd (2018) Coastal upwelling revisited: Ekman, Bakun, and improved upwelling indices for the U.S. west coast. Journal of Geophysical Research, doi:10.1029/2018JC014187. 
#Coastal Upwelling Transport Index


# Plot the CUTI data for our sampling dates

# function to reverse the axes
reverse2_trans <- function() {
  trans_new(
    "reverse2",
    function(x) -1 * as.numeric(x), # Force values to be numeric for Date objects
    function(x) -1 * as.numeric(x)
  )
}

CT %>%
  filter(latitude == 45) %>%
  mutate(date = as_date(date)) %>%
  ggplot(aes(x = date, y = CUTI))+
  geom_hline(aes(yintercept = 0), linetype = 2)+
  geom_ribbon(aes(xmin = ymd("2019-07-08"), 
                  xmax =ymd("2019-07-09") ), alpha = 0.2, fill = "red")+
  geom_ribbon(aes(xmin = ymd("2019-08-05"), 
                  xmax =ymd("2019-08-06") ), alpha = 0.2, fill = "red")+
  geom_point()+
  geom_point(data = CT %>%
               mutate(date = as_date(date)) %>%
               filter(latitude == 45,
                      date >= ymd("2019-07-08")&
                        date <= ymd("2019-07-09") |
                        date >= ymd("2019-08-05")&
                        date <= ymd("2019-08-06")
                        ),
             color = "red")+
  geom_path()+
  labs(x = "",
       y = "Coastal Upwelling Transport Index <br> (m<sup>2</sup> s<sup>-1</sup>)")+
  coord_flip()+
  scale_x_continuous(trans = c("date", "reverse2"))+
  theme_bw()+
  theme(axis.title.x = element_markdown(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank())
 
ggsave(here("Output","Figure_1b.pdf"), height = 8, width = 6)



### show difference from the ocean
data_all<-data_all %>%
  mutate(nh4_umol_l = ifelse(nh4_umol_l>60, NA, nh4_umol_l)) %>% # drop the 2 crazy outliers
  mutate(time_point = factor(time_point, levels = c("start","end")),
         do_umol_l = do_mg_l/0.032) #%>% # covert to micromol using 32 g/mol conversion
  #  filter(sampling_day !=ymd("2019-07-11"), time_point != 5)%>%
  #group_by(before_after, time_point, day_night, sampling_group)%>%
  #mutate(do_change = do_mg_l - do_mg_l[foundation_spp == "Ocean"],
  #       temp_change = temp_pool - temp_pool[foundation_spp == "Ocean"],
  #       nn_change = nn_umol_l - nn_umol_l[foundation_spp == "Ocean"],
  #       nh4_change = nh4_umol_l - nh4_umol_l[foundation_spp == "Ocean"],
  #       hetero_change = heterotrophic_bacterioplankton_m_l - heterotrophic_bacterioplankton_m_l[foundation_spp == "Ocean"],
  #       mc_change = m_c - m_c[foundation_spp == "Ocean"],
  #       bix_change = bix - bix[foundation_spp == "Ocean"],
  #       hix_change = hix - hix[foundation_spp == "Ocean"],
  #       fi_change = fi - fi[foundation_spp == "Ocean"],
  #       uv_humic_change = ultra_violet_humic_like - ultra_violet_humic_like[foundation_spp == "Ocean"],
  #       marine_change = marine_humic_like - marine_humic_like[foundation_spp == "Ocean"],
  #       visible_change = visible_humic_like - visible_humic_like[foundation_spp == "Ocean"],
  #       trypto_change = tryptophan_like - tryptophan_like[foundation_spp == "Ocean"],
  #       tyrosine_change = tyrosine_like - tyrosine_like[foundation_spp == "Ocean"],
  #       phenyl_change = phenylalanine_like - phenylalanine_like[foundation_spp == "Ocean"]
  #) 


## add in the metadata
data_all<-data_all %>%
  left_join(MetaData %>%
              clean_names() %>%
              mutate(pool_id = as.character(pool_id)))%>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)"))) 

# make it long  and select just the parameters I am working with
data_long<-data_all %>%
  ungroup()%>%
  mutate(together = paste(before_after, removal_control),
         manipulated = ifelse(together == "After Removal","Manipulated", "Not Manipulated"))%>%
  
  #filter(removal_control != "Removal") %>%
  select(month, pool_id,day_night, time_point,removal_control, manipulated,foundation_spp,do_umol_l,heterotrophic_bacterioplankton_m_l,tyrosine_like,tryptophan_like, phenylalanine_like, ultra_violet_humic_like, visible_humic_like, marine_humic_like, bix, hix, m_c, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = c(do_umol_l:nh4_umol_l)) %>%
  group_by(month, name, foundation_spp,time_point, day_night, manipulated)%>%
  summarise(mean_val = mean(value, na.rm = TRUE), # calculate means and SE
            se_val= sd(value, na.rm = TRUE)/sqrt(n())) %>%
  mutate(time_point_clean = ifelse(time_point == "start","Early", "Late"))%>%
  mutate(nicenames = case_when(  ## make pretty names for plotting
  #  name == "do_mg_l" ~ "DO <br> (mg L<sup>-1</sup>)",
    name == "do_umol_l" ~ "DO <br> (&mu;mol L<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic <br> (cells &mu;L<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
    name == "nn_umol_l" ~ "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
    name == "hix"~"HIX",
    name == "bix"~"BIX",
    name == "m_c"~"M:C",
    name =="tyrosine_like" ~"Tyrosine <br> (Raman units)",
    name == "tryptophan_like" ~ "Tryptophan <br> (Raman units)",
    name == "phenylalanine_like" ~"Phenylalanine <br> (Raman units)",
    name == "ultra_violet_humic_like" ~"UV Humic <br> (Raman units)",
    name == "visible_humic_like" ~"Visible Humic <br> (Raman units)",
    name == "marine_humic_like" ~ "Marine Humic <br> (Raman units)"
    )
  ) %>%
  mutate(nicenames = factor(nicenames, levels = c("DO <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Tyrosine <br> (Raman units)",
                                                  "Tryptophan <br> (Raman units)",
                                                  "Phenylalanine <br> (Raman units)",
                                                  "UV Humic <br> (Raman units)",
                                                  "Visible Humic <br> (Raman units)",
                                                  "Marine Humic <br> (Raman units)",
                                                  "HIX","BIX","M:C",
                                                  "Heterotrophic <br> (cells &mu;L<sup>-1</sup>)"
                                                  
                                                  
  )))




### create a function that calculated difference between start and end values for each tide pool

change_val<-function(x, time){
  val<-x[time == "end"]-x[time == "start"]
  val<-as.numeric(val)
  return(val)
}

# first get the times 
## calculate the change in time (hours) for each tidepool to join
times<-data_all%>%
  ungroup()%>%
  mutate(sampling_datetime = mdy_hms(paste(sampling_day, sampling_time))) %>% 
  select(pool_id,before_after,day_night, time_point, sampling_datetime,sampling_group) %>%
  group_by(pool_id, before_after, day_night,sampling_group) %>%
  reframe(diff_time = change_val(sampling_datetime, time_point)) 

# now calculate all the rates per hour
Rates<-data_all%>%
  ungroup()%>% 
  select(pool_id:removal_control, time_point, do_umol_l, nn_umol_l:nh4_umol_l, heterotrophic_bacterioplankton_m_l:fi, tyrosine_like,tryptophan_like, phenylalanine_like,
         ultra_violet_humic_like,visible_humic_like,marine_humic_like) %>%
  mutate(heterotrophic_bacterioplankton_m_l = heterotrophic_bacterioplankton_m_l *1e6 # convert to per L        
  )%>%
  mutate_at(vars(do_umol_l:nh4_umol_l), .funs = function(x){x/1000}) %>% #convert to mmol and g for chem
  pivot_longer(cols = do_umol_l:marine_humic_like)%>%
  group_by(pool_id, before_after, removal_control,foundation_spp, sampling_group, name) %>%
  reframe(change = change_val(value, time_point)) %>% # calculate the difference between start and end
  left_join(times)%>%# join with the times
  left_join(MetaData %>%
              clean_names() %>%  # clean the names
              mutate(pool_id = as.character(pool_id)) %>%
              select(pool_id, before_after, surface_area,vol))%>% # add in the tide pool info
  mutate(rate_hr = change/diff_time,# difference in value per hour
         rate_m2_hr = rate_hr*vol/surface_area  # normalize it to surface area for a flux (mmol m-2 hr-1)
  ) 

### run models of rates with a BACI design
mods_BACI<-Rates %>%
  filter(foundation_spp != "Ocean",
         removal_control != "Ocean")%>%
  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "m_c","bix","do_umol_l") )%>%
  mutate(rate_m2_hr_sqrt = sign(rate_m2_hr)*sqrt(abs(rate_m2_hr)), # sqrt transform the data
         #rate_m2_hr_log = sign(rate_m2_hr)*log(abs(rate_m2_hr))
         )%>%
  group_by(foundation_spp, name)%>%
  mutate(rate_diff_scale = as.numeric(scale(rate_m2_hr_sqrt, scale = TRUE,center = TRUE)) # z-score the data
        ) %>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) { # run an lmer for every model
                       lmer(rate_diff_scale ~ removal_control*before_after +(1|pool_id),
                          data = df)
                     })) %>%
  mutate(
    tidy = map(model, function(x)tidy(x,effects="fixed")), # extract the fixed effects
    glance = map(model, glance)
  ) %>%
  mutate(nicenames = case_when(
    name == "do_umol_l" ~ "&Delta;Dissolved Oxygen",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria",
    name == "nh4_umol_l" ~ "&Delta;Ammonium",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite",
    name == "bix"~"&Delta;BIX" ,
    name == "m_c"~"&Delta;M:C")
   )%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;Dissolved Oxygen",
                                                  "&Delta;Ammonium",
                                                  "&Delta;Nitrate+Nitrite",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic Bacteria",
                                                  "&Delta;fDOM"
                                                  
  )))

# extract the coefs for the interaction term
r2<-mods_BACI%>%
  unnest(tidy)%>%
  filter(
      term =="removal_controlRemoval:before_afterBefore",
    ) %>%
  mutate(alpha = ifelse(p.value<= 0.055,1, 0.6)) # make the color grey or dark for significance

write_csv(mods_BACI%>%
            unnest(tidy), file = here("Output","Rates_models.csv"))

## Do the models with the average concentration 
values<-data_all %>%
  ungroup()%>%
  filter(
         foundation_spp != "Ocean"
  )  %>%
  select(month,pool_id, time_point,removal_control, foundation_spp,do_umol_l,heterotrophic_bacterioplankton_m_l,m_c,
         tyrosine_like, tryptophan_like, phenylalanine_like, 
         ultra_violet_humic_like, visible_humic_like, marine_humic_like,bix, hix,fi, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_umol_l:nh4_umol_l) %>%
  group_by(foundation_spp, pool_id,removal_control,name, month) %>%
  summarise(mean_val = mean(value, na.rm = TRUE))%>%
  mutate(mean_val_sqrt = sqrt(mean_val),
         mean_val_log = log(mean_val))%>%
  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "do_umol_l","m_c","bix","hix", "tyrosine_like", "tryptophan_like", "phenylalanine_like", 
                     "ultra_violet_humic_like", "visible_humic_like", "marine_humic_like") )%>%
  group_by(foundation_spp, name)%>%
  mutate(value_scale = as.numeric(scale(mean_val_sqrt, scale = TRUE,center = TRUE))) %>%
  ungroup()%>%
  mutate(nicenames = case_when(
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria",
    name == "nh4_umol_l" ~ "&Delta;Ammonium",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite",
    name == "bix"~"&Delta;BIX" ,
    name == "hix"~"&Delta;HIX",
    name == "m_c"~"&Delta;M:C",
     name =="tyrosine_like" ~"&Delta;Tyrosine",
    name == "tryptophan_like" ~ "&Delta;Tryptophan",
    name == "phenylalanine_like" ~"&Delta;Phenylalanine",
    name == "ultra_violet_humic_like" ~"&Delta;UV Humic",
    name == "visible_humic_like" ~"&Delta;Visible Humic",
    name == "marine_humic_like" ~ "&Delta;Marine Humic",
    name == "do_umol_l" ~ "&Delta;Dissolved Oxygen"
    
  )
  )%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;Dissolved Oxygen",
                                                  "&Delta;Ammonium",
                                                  "&Delta;Nitrate+Nitrite",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Tyrosine",
                                                  "&Delta;Tryptophan",
                                                  "&Delta;Phenylalanine",
                                                  "&Delta;UV Humic",
                                                  "&Delta;Visible Humic",
                                                  "&Delta;Marine Humic",
                                                  "&Delta;Heterotrophic Bacteria",
                                                  "&Delta;fDOM"
                                                  
  )))

# run the models
mods2<-values %>%
  mutate(before_after = ifelse(month == "July", "Before","After"))%>%
  group_by(foundation_spp,nicenames)%>%
  nest() %>%
  mutate(model = map(data,
                     function(df) {
                       lmer(value_scale~ removal_control*before_after+(1|pool_id), data = df) #sqrt transformed 
                     })) %>%
  mutate(
    map(model, function(x)tidy(x,effects="fixed")),
    glance = map(model, glance)
  ) %>%
  rename(tidy = `map(model, function(x) tidy(x, effects = "fixed"))`)

# tidy the concentrations
conc2<-mods2%>%
  unnest(tidy)%>%
  filter(
    term =="removal_controlRemoval:before_afterBefore",
     ) %>%
  mutate(alpha = ifelse(p.value<= 0.055,1, 0.6))

write_csv(mods2%>%
            unnest(tidy), file = here("Output","Concentration_models.csv"))

### make a plot of the interaction terms
con_int<-conc2 %>%
  filter(nicenames %in%c("&Delta;Dissolved Oxygen","&Delta;Ammonium", "&Delta;Nitrate+Nitrite",
  "&Delta;BIX","&Delta;M:C","&Delta;Heterotrophic Bacteria"))%>%
ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha, color = foundation_spp))+
  scale_color_manual(values = c("black","#34c230"))+
  geom_point(size = 3, position = position_dodge(width = 0.9))+
  geom_errorbar(aes(xmin = estimate-std.error, xmax = estimate+std.error), width = 0.1, position = position_dodge(width = 0.9))+
  geom_vline(xintercept = 0)+
  labs(x = "Standardized effect size",
       y = "")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_markdown(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_markdown(size = 12),
        legend.position = "none")

r_int<-r2 %>%
  filter(name %in%c("do_umol_l","nh4_umol_l", "nn_umol_l","bix","m_c","heterotrophic_bacterioplankton_m_l"))%>%
  ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha, color = foundation_spp))+
  scale_color_manual(values = c("black","#34c230"))+
  geom_point(size = 3, position = position_dodge(width = 0.9))+
  geom_errorbar(aes(xmin = estimate-std.error, xmax = estimate+std.error), width = 0.1, position = position_dodge(width = 0.9))+
  geom_vline(xintercept=0)+
  labs(x = "Standardized effect size <br> (rates m<sup>-2</sup> hr<sup>-1</sup>)",
       y = "")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_markdown(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_markdown(size = 12),
        legend.position = "none")

BC_int<-r_int|con_int

# To get the upwelling effect we cal only look at the unmanipulated pools 
## do a two-way ANOVA between month and foundation species to see effect of upwelling and species identity

Unmamipulated_mean<-data_all %>%
  ungroup()%>%
  filter(foundation_spp != "Ocean")  %>%
  select(month, pool_id,before_after, removal_control, foundation_spp,do_umol_l,heterotrophic_bacterioplankton_m_l,m_c, bix,hix, nn_umol_l, nh4_umol_l) %>%
  pivot_longer(cols = do_umol_l:nh4_umol_l) %>%
  filter(name %in% c("nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "m_c","bix","do_umol_l") )%>%
  group_by(foundation_spp, pool_id,removal_control, before_after,name) %>%
  summarise(mean_val = mean(value, na.rm = TRUE)) %>% # get mean for the low tide
  mutate(together = paste(removal_control, before_after),
         mean_val_sqrt = sqrt(mean_val))%>%
  filter(together != "Removal After"
         
  )%>%
  group_by(name)%>%
  mutate(value_scale = as.numeric(scale(mean_val_sqrt, scale = TRUE,center = TRUE))) %>%
  ungroup()%>%
  mutate(nicenames = case_when(
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic Bacteria",
    name == "nh4_umol_l" ~ "Ammonium",
    name == "nn_umol_l" ~ "Nitrate+Nitrite",
    name == "bix"~"BIX" ,
    name == "m_c"~"M:C",
    name == "do_umol_l" ~ "Dissolved Oxygen"
      )
  )%>%
  mutate(nicenames = factor(nicenames, levels = c("Dissolved Oxygen",
                                                  "Ammonium",
                                                  "Nitrate+Nitrite",
                                                  "BIX" ,
                                                  "M:C",
                                                  "Heterotrophic Bacteria"
                                                  
  )))

# run the models
mods3<-Unmamipulated_mean %>%
  group_by(nicenames)%>%
  nest() %>%
  mutate(model = map(data,
                     function(df) {
                       lm(value_scale~ foundation_spp*before_after, data = df) #sqrt transformed 
                     })) %>%
  mutate(
    tidy = map(model, function(x)tidy(x,effects="fixed")),
    glance = map(model, glance)
  )

conc1<-mods3%>%
  unnest(tidy)%>%
  #filter(!nicenames %in% c("HIX",NA))%>%
  filter(
    term != "(Intercept)",
  ) %>%
  mutate(term = case_when(term == "before_afterBefore"~"Upwelling Effect",
                          term == "foundation_sppPhyllospadix"~"Foundation Species Effect",
                          term =="foundation_sppPhyllospadix:before_afterBefore"~"Interaction"))%>%
  mutate(term = factor(term,levels = c("Foundation Species Effect",
                                       "Upwelling Effect",
                                       "Interaction")))%>%
  mutate(alpha = ifelse(p.value <= 0.055, 1, 0.5))%>%
  ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha))+
  geom_point(size = 3)+
  geom_errorbar(aes(xmin = estimate-std.error, xmax = estimate+std.error), width = 0.1)+
  geom_vline(xintercept=0)+
  labs(x = "Standardized effect size <br> (concentration or density)",
       y = "")+
  facet_wrap(~term, scale = "free_y", ncol = 1, strip.position = "right")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_markdown(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_markdown(size = 14),
        legend.position = "none")


write_csv(mods3%>%
            unnest(tidy), here("Output","Concentration_controlonlymods.csv"))

#### do it again for the rates
mods<-Rates %>%
  filter(foundation_spp != "Ocean",
         removal_control != "Ocean")%>%
  mutate(together = paste(removal_control, before_after),
         rate_m2_hr_sqrt = sign(rate_m2_hr)*sqrt(abs(rate_m2_hr)))%>%
  filter(together != "Removal After"
        
  )%>%
  filter(name %in% c("m_c","bix", "nn_umol_l","heterotrophic_bacterioplankton_m_l","nh4_umol_l", "do_umol_l") )%>%
   group_by(name)%>%
  mutate(
      rate_hr_scale = as.numeric(scale(rate_m2_hr_sqrt, scale = TRUE,center = TRUE))) %>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) {
                       lm(rate_hr_scale  ~ foundation_spp*before_after, #sqrt transformed
                          data = df)
                     })) %>%
  
   mutate(
    tidy = map(model, function(x)tidy(x,effects="fixed")),
    glance = map(model, glance)
  ) %>%
  mutate(nicenames = case_when(
    name == "do_umol_l" ~ "&Delta;Dissolved Oxygen",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic bacteria",
     name == "nh4_umol_l" ~ "&Delta;Ammonium",
      name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite",
     name == "bix"~"&Delta;BIX" ,
     name == "m_c"~"&Delta;M:C"
  ))%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;Dissolved Oxygen",
                                                  "&Delta;Ammonium",
                                                  "&Delta;Nitrate+Nitrite",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic bacteria"
                                                  
  )))


r1<-mods%>%
  unnest(tidy)%>%
  filter(
    term != "(Intercept)",
   ) %>%
  mutate(term = case_when(term == "before_afterBefore"~"Upwelling Effect",
                          term == "foundation_sppPhyllospadix"~"Foundation Species Effect",
                          term =="foundation_sppPhyllospadix:before_afterBefore"~"Interaction"))%>%
  mutate(term = factor(term,levels = c("Foundation Species Effect",
                                       "Upwelling Effect",
                                       "Interaction")))%>%
  mutate(alpha = ifelse(p.value <= 0.055, 1, 0.5))%>%
  ggplot(aes(y = fct_rev(nicenames), x = estimate, alpha = alpha))+
  geom_point(size = 3)+
  geom_errorbar(aes(xmin = estimate-std.error, xmax = estimate+std.error), width = 0.1)+
  geom_vline(xintercept=0)+
  labs(x = "Standardized effect size <br> (rates m<sup>-2</sup> hr<sup>-1</sup>)",
       y = "")+
  facet_wrap(~term, scale = "free_y", ncol = 1, strip.position = "right")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_markdown(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_markdown(size = 14),
        legend.position = "none")

write_csv(mods%>%
            unnest(tidy), here("Output","Rates_controlonlymods.csv"))


### Make a long dataframe to get ready for models
Long_all<-Rates %>% 
  filter(
         foundation_spp != "Ocean",
         name %in% c("bix","do_umol_l","heterotrophic_bacterioplankton_m_l", "m_c", "nh4_umol_l", "nn_umol_l", "marine_humic_like", 
                     "phenylalanine_like", "tryptophan_like",
                     "tyrosine_like", "ultra_violet_humic_like", "visible_humic_like"))%>%
  mutate(together = paste(before_after, removal_control),
         manipulated = ifelse(together == "After Removal","Manipulated", "Not Manipulated"))%>%
  mutate(rate_sqrt = sign(rate_m2_hr)*sqrt(abs(rate_m2_hr)),
         rate_log = sign(rate_m2_hr)*log(abs(rate_m2_hr)))

# make a square root function with negatives
sqrt_fcn<-function(x){sign(x)*sqrt(abs(x))}
# transform the axes with this function
tn <- trans_new("sqrt_fcn",
                function(x){sign(x)*sqrt(abs(x))},
                function(y){sign(y)*abs(y)^2}
               )
# Plot showing relationship between DO and HBac

Rates_wide<- Rates %>%  
  filter(day_night == "Day") %>%
  filter(foundation_spp != "Ocean") %>%
  select(-c(change:vol, rate_hr))%>%
  pivot_wider(values_from = rate_m2_hr, names_from = name) %>%  
  mutate(month = factor(ifelse(before_after == "Before", "July", "August (Upwelling)"), levels = c("July", "August (Upwelling)")))

# Annotation data frame
annotations <- data.frame(
  removal = c("Foundation species removed","Unmanipulated"),
  label = c("b)", "a)"),
  x_pos = c(-1.875, -2.1875), # Adjust position as needed
  y_pos = c(.370, .370)  # Adjust position as needed
)%>%
  mutate( removal = factor(removal, levels = c("Foundation species removed","Unmanipulated")),)

Rates_wide %>%
  mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                             removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                             removal_control == "Removal"& month != "July" ~ "Foundation species removed"))%>%
  mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation species removed")))%>%
  mutate(foundation_spp = ifelse(foundation_spp == "Mytilus", "Mussel-Dominated","Surfgrass-Dominated"))%>%
  mutate(het_carbon = heterotrophic_bacterioplankton_m_l*20/1e12) %>% # 20 fmol C - millimol
  ggplot(aes(y = het_carbon, x = do_umol_l))+
  geom_hline(aes(yintercept  = 0), lty = 2, alpha = 0.5)+
  geom_vline(aes(xintercept  = 0), lty = 2, alpha = 0.5)+
  geom_point(aes(color = foundation_spp))+
  geom_smooth(method ="lm", color = "grey3", data = Rates_wide %>%  
                mutate(het_carbon = heterotrophic_bacterioplankton_m_l*20/1e12) %>%
                mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                                           removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                                           removal_control == "Removal"& month != "July" ~ "Foundation species removed"))%>%
                mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation species removed")))%>%
                filter(foundation_spp != "Ocean",
                       removal == "Unmanipulated"))+
  geom_text(data = annotations, aes(x = x_pos, y = y_pos, label = label), 
            hjust = 0, vjust = 1, size = 5)+
  scale_color_manual(values = c("black","#34c230"))+
  #scale_shape_manual(values = c(1,16))+
  labs(y = "Heterotrophic bacteria <br> (mmol  C  m<sup>-2</sup> hr<sup>-1</sup>)",
       x = "Dissolved oxygen flux <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
       color = " Foundation Species")+
  facet_wrap(~removal, scales = "free_x")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_markdown(size = 14),
        axis.title.y = element_markdown(size = 14),
        axis.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
  )
ggsave(here("Output","Figure_5.pdf"), width = 8, height = 4)


DO_het_mod<-lm(data = Rates_wide %>%  
                 mutate(het_carbon = heterotrophic_bacterioplankton_m_l*20/1e6) %>%
                   mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                                              removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                                              removal_control == "Removal"& month != "July" ~ "Foundation spp. removed"))%>%
                   mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation spp. removed")))%>%
                   filter(foundation_spp != "Ocean"), het_carbon~do_umol_l*removal)


anova(DO_het_mod)
summary(DO_het_mod)

# get the ocean temperatures
data_all %>%
  filter(removal_control == "Ocean") %>%
  select(before_after, temp_pool) %>%
  group_by(before_after)%>%
  summarise(mean_temp = mean(temp_pool, na.rm = TRUE),
            se_temp = sd(temp_pool)/sqrt(n()))

## run one-way t-tests
Rates_ttest <- Long_all %>%
  group_by(name,foundation_spp, manipulated, before_after)%>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) {
                       t.test(df$rate_sqrt)
                     })) %>%
  
  mutate(
    tidy = map(model, tidy),
    glance = map(model, glance)
  )  %>% 
  unnest(tidy) %>%
  select(-c(data, model, glance, method, alternative)) %>%
  ungroup()%>%
  mutate(before_after = factor(before_after, levels = c("Before","After"))) %>%
  mutate(significant = ifelse(p.value<0.055, 1, 0)) %>%
  arrange(foundation_spp, before_after) %>%
  select(before_after:manipulated, estimate, statistic, p.value, significant) 

# write the t-test output
  write_csv(Rates_ttest, here("Output","ttest_rates.csv"))

  # Look at all the individual fDOM values
Long_wfDOM<-Long_all %>% 
  mutate(man = ifelse(manipulated == "Manipulated", "Foundation spp removed", "Unmanipulated"),
         foundation_spp = ifelse(foundation_spp == "Mytilus", "Mussels", "Surfgrass")) %>%
  mutate(nicenames = case_when(
    name == "do_umol_l" ~ "&Delta;Dissolved Oxygen <br> (&mu;mol)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria <br> (counts)",
    name == "nh4_umol_l" ~ "&Delta; Ammonium <br> (&mu;mol)",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite <br> (&mu;mol)",
    name == "bix"~"&Delta;BIX" ,
    name == "hix"~"&Delta;HIX",
    name == "m_c"~"&Delta;M:C",
    name =="tyrosine_like" ~"&Delta;Tyrosine <br> (Raman units)",
    name == "tryptophan_like" ~ "&Delta;Tryptophan <br> (Raman units)",
    name == "phenylalanine_like" ~"&Delta;Phenylalanine <br> (Raman units)",
    name == "ultra_violet_humic_like" ~"&Delta;UV Humic <br> (Raman units)",
    name == "visible_humic_like" ~"&Delta;Visible Humic <br> (Raman units)",
    name == "marine_humic_like" ~ "&Delta;Marine Humic <br> (Raman units)"
  ))%>%
  mutate(nicenames = factor(nicenames, levels = c("&Delta;Dissolved Oxygen <br> (&mu;mol)",
                                                  "&Delta; Ammonium <br> (&mu;mol)",
                                                  "&Delta;Nitrate+Nitrite <br> (&mu;mol)",
                                                  "&Delta;BIX" ,
                                                  "&Delta;M:C",
                                                  "&Delta;HIX",
                                                  "&Delta;Heterotrophic Bacteria <br> (counts)",
                                                  "&Delta;Tyrosine <br> (Raman units)",
                                                  "&Delta;Tryptophan <br> (Raman units)",
                                                  "&Delta;Phenylalanine <br> (Raman units)",
                                                  "&Delta;UV Humic <br> (Raman units)",
                                                  "&Delta;Visible Humic <br> (Raman units)",
                                                  "&Delta;Marine Humic <br> (Raman units)"
  )))


##axes the same
scale_x <- Long_wfDOM %>%
  filter(name %in% c("do_umol_l","heterotrophic_bacterioplankton_m_l",
                     "nh4_umol_l","nn_umol_l","bix", "hix","m_c",
                     "tyrosine_like" , "tryptophan_like",
                     "phenylalanine_like","ultra_violet_humic_like",
                     "visible_humic_like", "marine_humic_like")) %>%
  ungroup()%>%
  group_by(foundation_spp, nicenames, before_after, manipulated)%>%
  summarise(mean_r = mean(rate_m2_hr, na.rm = TRUE),
            se_r = sd(rate_m2_hr, na.rm = TRUE)/sqrt(n())) %>%
  mutate(max = mean_r+se_r,
         min = mean_r - se_r) %>%
  select(foundation_spp, nicenames, max, min) %>%
  pivot_longer(cols = max:min) %>%
  rename(rate_m2_hr = value) %>%
  bind_rows(tibble(nicenames = "&Delta; Ammonium <br> (&mu;mol)", rate_m2_hr = 0)) %>%
  split(~nicenames) |>
  map(~range(.x$rate_m2_hr)) |> 
  imap(
    ~ scale_x_facet(
      nicenames == .y,
      limits = .x
    )
  )

# plot the rates
#Long_wfDOM %>%
#  filter(name %in% c("do_umol_l","heterotrophic_bacterioplankton_m_l",
#                     "nh4_umol_l","nn_umol_l","bix", "hix","m_c" )) %>%
#  left_join(Rates_ttest %>%
#              select(before_after, foundation_spp, name, manipulated, significant))%>%
#ggplot(aes(x= rate_m2_hr, y = before_after, shape = man, color = foundation_spp))+
#  geom_vline(xintercept = 0,  lty = 2)+
#  stat_summary(size = 0.8, fun.data = "mean_se")+
#  scale_shape_manual(values = c(21,19))+
#  scale_color_manual(values = c("black","#34c230"), guide = "none")+
#  labs(x = "Flux (m<sup>-2</sup> hr<sup>-1</sup>)", 
#       y = "",
#       shape = "")+
#  ggh4x::facet_grid2(nicenames~foundation_spp,scales = "free_x", independent = "x")+
#  theme_bw()+
#  theme(strip.background = element_blank(),
#        strip.placement = "outside", 
#        panel.grid.minor = element_blank(),
#        strip.text = element_markdown(size = 12),
#        axis.title.x = element_markdown(),
#        legend.position = "bottom", 
#        axis.text = element_text(size = 10),
#        axis.title = element_text(size=12))+
#  scale_x

## cobble params for the supplement

#Long_wfDOM %>%
#  filter(name %in% c("tyrosine_like" , "tryptophan_like",
#                     "phenylalanine_like","ultra_violet_humic_like",
#                     "visible_humic_like", "marine_humic_like" )) %>%
#  ggplot(aes(x= rate_m2_hr, y = before_after, shape = man, color = foundation_spp))+
#  geom_vline(xintercept = 0,  lty = 2)+
#  # geom_point(alpha = 0.1)+
#  stat_summary(size = 0.8, fun.data = "mean_se")+
#  scale_shape_manual(values = c(21,19))+
#  scale_color_manual(values = c("black","#34c230"), guide = "none")+
#  labs(x = "Flux (m<sup>-2</sup> hr<sup>-1</sup>)", 
#       y = "",
#       shape = "")+
#  ggh4x::facet_grid2(nicenames~foundation_spp,scales = "free_x", independent = "x")+
#  theme_bw()+
#  theme(strip.background = element_blank(),
#        strip.placement = "outside", 
#        panel.grid.minor = element_blank(),
#        strip.text.y.right = element_markdown(size = 12),
#        axis.title.x = element_markdown(),
#        legend.position = "bottom", 
#        axis.text = element_text(size = 10),
#        axis.title = element_text(size=12))+
#  scale_x

## same for the values

# get the ocean values
ocean <-data_all %>%
  ungroup()%>%
  filter(
         removal_control == "Ocean") %>%
  select(before_after, do_umol_l, heterotrophic_bacterioplankton_m_l, nn_umol_l, nh4_umol_l, bix, hix, m_c, tryptophan_like, tyrosine_like, phenylalanine_like, ultra_violet_humic_like, marine_humic_like, visible_humic_like) %>%
  group_by(before_after)%>%
  summarise_all(.funs = function(x){mean(x, na.rm = TRUE)}) %>%
  pivot_longer(cols = do_umol_l:visible_humic_like) %>%
  mutate(nicenames = case_when(
    name == "do_umol_l" ~ "Dissolved Oxygen <br> (&mu;mol L<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic Bacteria <br> (cells &mu;L<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
    name == "nn_umol_l" ~ "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
    name == "bix"~"BIX" ,
    name == "hix"~"HIX",
    name == "m_c"~"M:C",
    name =="tyrosine_like" ~"Tyrosine <br> (Raman units)",
    name == "tryptophan_like" ~ "Tryptophan <br> (Raman units)",
    name == "phenylalanine_like" ~"Phenylalanine <br> (Raman units)",
    name == "ultra_violet_humic_like" ~"UV Humic <br> (Raman units)",
    name == "visible_humic_like" ~"Visible Humic <br> (Raman units)",
    name == "marine_humic_like" ~ "Marine Humic <br> (Raman units)"))%>%
  mutate(nicenames = factor(nicenames, levels = c("Dissolved Oxygen <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
                                                  "BIX" ,
                                                  "M:C",
                                                  "HIX",
                                                  "Tyrosine <br> (Raman units)",
                                                  "Tryptophan <br> (Raman units)",
                                                  "Phenylalanine <br> (Raman units)",
                                                  "UV Humic <br> (Raman units)",
                                                  "Visible Humic <br> (Raman units)",
                                                  "Marine Humic <br> (Raman units)",
                                                  "Heterotrophic Bacteria <br> (cells &mu;L<sup>-1</sup>)"
  )))%>%
  mutate(removal = "Ocean")

## Plot all the values
value_plotdata<-values %>%
  mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                             removal_control == "Removal"& month == "July" ~ "Unmanipulated",
                             removal_control == "Removal"& month != "July" ~ "Foundation species removed"))%>%
  mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation species removed")))%>%
  mutate(foundation_spp = ifelse(foundation_spp == "Mytilus", "Mussels", "Surfgrass")) %>%
  mutate(before_after = ifelse(month == "July","Before","After"))%>%
  mutate(nicenames = case_when(
    name == "do_umol_l" ~ "Dissolved Oxygen <br> (&mu;mol L<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "Heterotrophic Bacteria <br> (cells &mu;L<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
    name == "nn_umol_l" ~ "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
    name == "bix"~"BIX" ,
    name == "m_c"~"M:C",
   ))%>%
  mutate(nicenames = factor(nicenames, levels = c("Dissolved Oxygen <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Ammonium <br> (&mu;mol L<sup>-1</sup>)",
                                                  "Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
                                                  "BIX" ,
                                                  "M:C",
                                                  "Heterotrophic Bacteria <br> (cells &mu;L<sup>-1</sup>)"
  )))


scale_x2 <- value_plotdata %>%
  ungroup()%>%
  group_by(foundation_spp, nicenames, before_after, removal)%>%
  summarise(mean_r = mean(mean_val , na.rm = TRUE),
            se_r = sd(mean_val , na.rm = TRUE)/sqrt(n())) %>%
  mutate(max = mean_r+se_r,
         min = mean_r - se_r) %>%
  select(foundation_spp, nicenames, max, min) %>%
  pivot_longer(cols = max:min) %>%
  rename(mean_val = value) %>%
  bind_rows(tibble(nicenames = "Ammonium <br> (&mu;mol L<sup>-1</sup>)", mean_val = 0)) %>%
  bind_rows(ocean %>% select(before_after, nicenames, mean_val = value)) %>%
  split(~nicenames) |>
  map(~range(.x$mean_val )) |> 
  imap(
    ~ scale_x_facet(
      nicenames == .y,
      limits = .x
    )
  )

#stocksplot<-value_plotdata %>%
#  filter(name %in% c("do_umol_l","heterotrophic_bacterioplankton_m_l",
#                     "nh4_umol_l","nn_umol_l","bix", "m_c" )) %>%
#  mutate(removal = factor(removal, levels = c("Foundation species removed", "Unmanipulated")),
#         foundation_spp = ifelse(foundation_spp == "Mussels", "Mussel-Dominated", "Surfgrass-Dominated"))%>%
#  ggplot(aes(x= mean_val, y = before_after, shape = removal))+
#  #geom_vline(xintercept = 0, color = "grey")+
#  geom_point(alpha = 0.1)+
#  stat_summary(size = 0.8, aes(color = foundation_spp))+
#  geom_point(data = ocean %>%
#               filter(name %in% c("do_umol_l","heterotrophic_bacterioplankton_m_l",
#                                  "nh4_umol_l","nn_umol_l","bix", "m_c" )), aes(x= value, y = before_after), 
#             color = "lightblue", size = 3)+
#  scale_shape_manual(values = c(21,19,18))+
#  scale_color_manual(values = c("black","#34c230"), guide = "none")+
#  labs(x = "Stock", 
#       y = "",
#       shape = "")+
#  ggh4x::facet_grid2(nicenames~foundation_spp, scales = "free_x", independent = "x")+
#  theme_bw()+
#  theme(strip.background = element_blank(),
#        strip.placement = "outside", 
#        panel.grid.minor = element_blank(),
#        strip.text = element_markdown(size = 12),
#        axis.title.x = element_markdown(),
#        legend.position = "bottom", 
#        axis.text = element_text(size = 10),
#        axis.title = element_text(size=12),
#        legend.text = element_text(size = 10))+
#  scale_x2


#######
## Make a plot of the benthic data

### order the pools to make it prettier
PoolOrder<-as_tibble(list(PoolOrder = c(1:32), 
                          PoolID = factor(c(4,1,20, 18, 5, 27,8, 29,
                                            21,26,6,3,2,28,22,7,13,31,14,19,30,9,25,
                                            12,32,15,11,10,17,16,23,24))))%>%
  mutate(PoolOrder2 = c(1:16,1:16))


BenthicData %>%
  # filter(Before_After == "After")%>%
  mutate(Before_After = factor(Before_After, levels = c("Before","After")))%>%
  select(PoolID, Before_After, Removal_Control, Foundation_spp, Macroalgae = macroalgae , Microphytobenthos = Diatoms, `Non-mussel Consumers` = consumers, `Crustose Coralline Algae` = allCCA, Mussels = AdjMusselCover, Surfgrass = AdjSurfgrassCover)%>%
  mutate(PoolID =  as.factor(PoolID)) %>%
  mutate(`Rock/Sand` = (100-(Macroalgae+`Non-mussel Consumers`+`Crustose Coralline Algae`+Mussels+Surfgrass)),
         Macroalgae = Macroalgae - Microphytobenthos,
         PoolID = fct_reorder2(PoolID,`Rock/Sand`, Mussels)) %>% # the macroalgae includes diatom cover and I want to see the difference
  left_join(PoolOrder)%>%
  pivot_longer(cols = Macroalgae:`Rock/Sand`) %>%
  mutate(name = factor(name, levels = c("Surfgrass","Macroalgae","Mussels","Non-mussel Consumers", "Microphytobenthos","Crustose Coralline Algae","Rock/Sand")),
         Foundation_spp = ifelse(Foundation_spp == "Mytilus","Mussel-dominated","Surfgrass-dominated"),
         Before_After = ifelse(Before_After == "Before","Before","After-Impact"),
         Before_After = factor(Before_After, levels = c("Before","After-Impact")))%>%
  mutate(Foundation_spp = factor(Foundation_spp, levels = c("Surfgrass-dominated","Mussel-dominated"))) %>%
  ggplot(aes(x = PoolOrder2, y = value, fill = name))+
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = cal_palette("tidepool", n = 7, type = "continuous"))+
  geom_vline(aes(xintercept = 8.5), linewidth = 1.25, color = "white")+
  labs(x = "",
       y = "% Cover",
       fill = "")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16)
  )+
  facet_grid(Before_After~Foundation_spp, scale = "free_x")

ggsave(here("Output","Figure_1c.pdf"), width = 10, height = 8, device = cairo_pdf)

BenthicData %>%
  # filter(Before_After == "After")%>%
  mutate(Before_After = factor(Before_After, levels = c("Before","After")))%>%
  select(PoolID, Before_After, Removal_Control, Foundation_spp, Macroalgae = macroalgae , Microphytobenthos = Diatoms, `Non-mussel Consumers` = consumers, `Crustose Coralline Algae` = allCCA, Mussels = AdjMusselCover, Surfgrass = AdjSurfgrassCover)%>%
  mutate(PoolID =  as.factor(PoolID)) %>%
  mutate(`Rock/Sand` = (100-(Macroalgae+`Non-mussel Consumers`+`Crustose Coralline Algae`+Mussels+Surfgrass)),
         Macroalgae = Macroalgae - Microphytobenthos,
         PoolID = fct_reorder2(PoolID,`Rock/Sand`, Mussels)) %>% # the macroalgae includes diatom cover and I want to see the difference
  pivot_longer(cols = Macroalgae:`Rock/Sand`) %>%
mutate(Removal_Control = ifelse(Before_After == "Before", "Control", Removal_Control)) %>%
 group_by(Foundation_spp, Removal_Control, Before_After, name) %>%
  summarise(mean_cover = mean(value, na.rm = TRUE),
            se_cover = sd(value, na.rm = TRUE)/sqrt(n()))


Value_rates<-values %>%
  select(foundation_spp:month, mean_val, name) %>%
  pivot_wider(names_from = name, values_from = mean_val, names_prefix = "value_")%>%
  mutate(value_total_fDOM = value_marine_humic_like+value_phenylalanine_like+value_tryptophan_like+value_tyrosine_like+value_ultra_violet_humic_like+value_visible_humic_like)%>%
  left_join(Rates_wide) %>%
  mutate(manipulated = ifelse(month == "July", "Control", removal_control))


#### Add in the paper figures for rates and stocks using the raw data 

scale_y2 <- Long_wfDOM %>%
  filter(name %in% c("heterotrophic_bacterioplankton_m_l",
                     "nh4_umol_l","nn_umol_l","bix", "m_c", "do_umol_l")) %>%
   mutate(nicenames2 = case_when(
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria <br> (cells m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "&Delta;Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "bix"~"&Delta;BIX <br> (m<sup>-2</sup> hr<sup>-1</sup>)" ,
    name == "m_c"~"&Delta;M:C <br> (m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "do_umol_l"~"&Delta;Dissolved Oxygen <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)"
  ))%>%
  mutate(nicenames2 = factor(nicenames2, levels = c(
    "&Delta;Dissolved Oxygen <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
    "&Delta;Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
    "&Delta;Nitrate+Nitrite <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
    "&Delta;BIX <br> (m<sup>-2</sup> hr<sup>-1</sup>)"  ,
    "&Delta;M:C <br> (m<sup>-2</sup> hr<sup>-1</sup>)" ,
    "&Delta;Heterotrophic Bacteria <br> (cells m<sup>-2</sup> hr<sup>-1</sup>)"
  ))) %>%
  ungroup()%>%
  group_by(foundation_spp, nicenames2, before_after, manipulated)%>%
  summarise(max = max(rate_m2_hr, na.rm = TRUE),
            min = min(rate_m2_hr, na.rm = TRUE)) %>%
  select(foundation_spp, nicenames2, max, min) %>%
  pivot_longer(cols = max:min) %>%
  rename(rate_m2_hr = value) %>%
  bind_rows(tibble(nicenames2 = "Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)", rate_m2_hr = 0)) %>%
  split(~nicenames2) |>
  map(~range(.x$rate_m2_hr)) |> 
  imap(
    ~ scale_y_facet(
      nicenames2 == .y,
      limits = .x
    )
  )

rate_mean_plot<-Long_wfDOM %>%
  mutate(nicenames2 = case_when(
    name == "do_umol_l" ~"&Delta;Dissolved Oxygen <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria <br> (cells m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nh4_umol_l" ~ "&Delta;Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
    name == "bix"~"&Delta;BIX <br> (m<sup>-2</sup> hr<sup>-1</sup>)" ,
    name == "m_c"~"&Delta;M:C <br> (m<sup>-2</sup> hr<sup>-1</sup>)",
    ))%>%
  mutate(nicenames2 = factor(nicenames2, levels = c("&Delta;Dissolved Oxygen <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
                                                  "&Delta;Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
                                                  "&Delta;Nitrate+Nitrite <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
                                                  "&Delta;BIX <br> (m<sup>-2</sup> hr<sup>-1</sup>)"  ,
                                                  "&Delta;M:C <br> (m<sup>-2</sup> hr<sup>-1</sup>)" ,
                                                 "&Delta;Heterotrophic Bacteria <br> (cells m<sup>-2</sup> hr<sup>-1</sup>)"
  ))) %>%
  
  filter(before_after == "After")%>%
  filter(name %in% c("do_umol_l","bix","heterotrophic_bacterioplankton_m_l","m_c","nh4_umol_l","nn_umol_l")) %>%
  ggplot(aes(shape = removal_control, y = rate_m2_hr, color = foundation_spp, x = foundation_spp))+
  geom_hline(yintercept = 0, lty = 2)+
 # geom_violin(aes(fill = foundation_spp), alpha = 0.3, color = NA)+
  stat_summary(size = 0.7)+
  scale_color_manual(values = c("black","#34c230"), guide = "none")+
  labs(x = "", 
      # y = expression("Rate (value m"^-2~"hr"^-1~")"),
       y = "",
       shape = "")+
  scale_shape_manual (values = c(16,1))+
  facet_wrap(~nicenames2, scales = "free_y",nrow = 1, strip.position = "left")+
  
  #ggh4x::facet_grid2(nicenames2~foundation_spp, scales = "free_y", independent = "y", switch = "y")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside", 
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_markdown(size = 12),
        axis.title.x = element_markdown(),
        legend.position = "bottom", 
        axis.text = element_text(size = 10),
        axis.title = element_text(size=12),
        legend.text = element_text(size = 10),
        aspect.ratio = 1)
 # scale_y2


### same with stocks for just the after period
scale_y4 <- value_plotdata %>%
  filter(before_after == "After")%>%
  filter(name %in% c("heterotrophic_bacterioplankton_m_l",
                     "nh4_umol_l","nn_umol_l","bix", "m_c", "do_umol_l")) %>%
   ungroup()%>%
  group_by(foundation_spp, nicenames, removal_control)%>%
  summarise(max = max(mean_val, na.rm = TRUE),
            min = min(mean_val, na.rm = TRUE)) %>%
  select(foundation_spp, nicenames, max, min) %>%
  pivot_longer(cols = max:min) %>%
  rename(mean_val = value) %>%
  split(~nicenames) |>
  map(~range(.x$mean_val)) |> 
  imap(
    ~ scale_y_facet(
      nicenames == .y,
      limits = .x
    )
  )

# get the ocean data
ocean_line <- ocean %>%
  filter(name %in% c("do_umol_l","bix","heterotrophic_bacterioplankton_m_l","m_c","nh4_umol_l","nn_umol_l")) %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August"), levels = c("July", "August"))) %>%
  mutate(mean_val = value)


value_mean_plot<-value_plotdata %>%
  filter(before_after == "After")%>%
  filter(name %in% c("do_umol_l","bix","heterotrophic_bacterioplankton_m_l","m_c","nh4_umol_l","nn_umol_l")) %>%
  ggplot(aes(x = foundation_spp, color = foundation_spp, y = mean_val, shape = removal_control))+
   geom_hline(data = ocean_line %>% filter(before_after == "After"), aes(yintercept = mean_val), color = "lightblue", linewidth = 1.5)+
  stat_summary(size = 0.7)+
  scale_color_manual(values = c("black","#34c230"), guide = "none")+
  scale_shape_manual (values = c(16,1))+
  labs(x = "", 
       # y = expression("Rate (value m"^-2~"hr"^-1~")"),
       y = "",
       shape = "")+
  facet_wrap(~nicenames, scales = "free_y",nrow = 1, strip.position = "left")+
   theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside", 
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_markdown(size = 12),
        axis.title.x = element_markdown(),
        legend.position = "bottom", 
        axis.text = element_text(size = 10),
        axis.title = element_text(size=12),
        legend.text = element_text(size = 10),
        aspect.ratio = 1)
  #scale_y4


Fig4<-((rate_mean_plot&theme(legend.position = "none", axis.text.x = element_blank()))+r_int+theme(aspect.ratio = 1, axis.title.x = element_blank()) +plot_layout(nrow = 1, widths = c(8,1)))/((value_mean_plot&theme(legend.position = "none"))+con_int+theme(aspect.ratio = 1, axis.title.x = element_blank()) +plot_layout(nrow = 1, widths = c(8,1)))

Fig4
ggsave(filename = here("Output","Figure_4.pdf"), height = 10, width = 20, device = cairo_pdf )


## same with the stocks but include both July and August control only

scale_y3 <- value_plotdata %>%
  mutate(month = factor(ifelse(before_after == "Before", "July", "August"), levels = c("July", "August")))%>%
  mutate(mean_val = ifelse(name == "m_c" & mean_val>2, NA, mean_val)) %>%
  mutate(mean_val = ifelse(name == "bix" & mean_val>2, NA, mean_val)) %>%
  filter(name %in% c("do_umol_l","heterotrophic_bacterioplankton_m_l",
                     "nh4_umol_l","nn_umol_l","bix", "m_c")) %>%
  ungroup()%>%
  group_by(foundation_spp, nicenames, before_after, month)%>%
  summarise(max = max(mean_val, na.rm = TRUE),
            min = min(mean_val, na.rm = TRUE)) %>%
  select(foundation_spp, nicenames, month, max, min) %>%
  pivot_longer(cols = max:min) %>%
  rename(mean_val = value) %>%
  split(~nicenames) |>
  map(~range(.x$mean_val)) |> 
  imap(
    ~ scale_y_facet(
      nicenames == .y,
      limits = .x
    )
  )

 mean_box<-value_plotdata %>%
   mutate(month = factor(ifelse(before_after == "Before", "July", "August"), levels = c("July", "August"))) %>%
   mutate(mean_val = ifelse(name == "m_c" & mean_val>2, NA, mean_val)) %>%
   mutate(mean_val = ifelse(name == "bix" & mean_val>2, NA, mean_val)) %>%
   filter(removal_control == "Control")%>%
   filter(name %in% c("do_umol_l","bix","heterotrophic_bacterioplankton_m_l","m_c","nh4_umol_l","nn_umol_l")) %>%
   ggplot(aes(x = foundation_spp, y = mean_val))+
   geom_hline(data = ocean_line, aes(yintercept = value), color = "lightblue", linewidth = 1.5)+
   geom_boxplot(aes(fill = foundation_spp), alpha = 0.3)+
    labs(x = "",
        y = " ")+
   scale_fill_manual(values = c("black","#34c230"), guide = "none")+
   
   ggh4x::facet_grid2(nicenames~month, scales = "free_y", independent = "y", switch = "y")+
   theme_bw()+
   theme(strip.background = element_blank(),
         strip.placement = "outside", 
         panel.grid.minor = element_blank(),
         strip.text.y.left = element_markdown(size = 12),
         axis.title.x = element_markdown(),
         legend.position = "bottom", 
         axis.text = element_text(size = 10),
         axis.title = element_text(size=12),
         legend.text = element_text(size = 10))+
   scale_y3
 
 # boxplots with the rates
 scale_y5 <- Long_wfDOM %>%
   mutate(month = factor(ifelse(before_after == "Before", "July", "August"), levels = c("July", "August"))) %>%
   filter(removal_control == "Control")%>%
   filter(name %in% c("do_umol_l","heterotrophic_bacterioplankton_m_l",
                      "nh4_umol_l","nn_umol_l","bix", "m_c")) %>%
   mutate(nicenames2 = case_when(
     name == "do_umol_l"~"&Delta;Dissolved Oxygen <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria <br> (cells m<sup>-2</sup> hr<sup>-1</sup>)",
     name == "nh4_umol_l" ~ "&Delta;Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     name == "bix"~"&Delta;BIX <br> (m<sup>-2</sup> hr<sup>-1</sup>)" ,
     name == "m_c"~"&Delta;M:C <br> (m<sup>-2</sup> hr<sup>-1</sup>)",
   ))%>%
   mutate(nicenames2 = factor(nicenames2, levels = c(
     "&Delta;Dissolved Oxygen <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     "&Delta;Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     "&Delta;Nitrate+Nitrite <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     "&Delta;BIX <br> (m<sup>-2</sup> hr<sup>-1</sup>)"  ,
     "&Delta;M:C <br> (m<sup>-2</sup> hr<sup>-1</sup>)" ,
     "&Delta;Heterotrophic Bacteria <br> (cells m<sup>-2</sup> hr<sup>-1</sup>)"
   ))) %>%
   ungroup()%>%
   group_by(foundation_spp, nicenames2, before_after, month)%>%
   summarise(max = max(rate_m2_hr, na.rm = TRUE),
             min = min(rate_m2_hr, na.rm = TRUE)) %>%
   select(foundation_spp, nicenames2, max, min) %>%
   pivot_longer(cols = max:min) %>%
   rename(rate_m2_hr = value) %>%
   bind_rows(tibble(nicenames2 = "Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)", rate_m2_hr = 0)) %>%
   split(~nicenames2) |>
   map(~range(.x$rate_m2_hr)) |> 
   imap(
     ~ scale_y_facet(
       nicenames2 == .y,
       limits = .x
     )
   )
 
  
 rate_box<-Long_wfDOM %>%
   mutate(nicenames2 = case_when( name == "do_umol_l"~"&Delta;Dissolved Oxygen <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     name == "heterotrophic_bacterioplankton_m_l" ~ "&Delta;Heterotrophic Bacteria <br> (cells m<sup>-2</sup> hr<sup>-1</sup>)",
     name == "nh4_umol_l" ~ "&Delta;Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     name == "nn_umol_l" ~ "&Delta;Nitrate+Nitrite <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     name == "bix"~"&Delta;BIX <br> (m<sup>-2</sup> hr<sup>-1</sup>)" ,
     name == "m_c"~"&Delta;M:C <br> (m<sup>-2</sup> hr<sup>-1</sup>)",
   ))%>%
   mutate(nicenames2 = factor(nicenames2, levels = c(
     "&Delta;Dissolved Oxygen <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     "&Delta;Ammonium <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     "&Delta;Nitrate+Nitrite <br> (mmol m<sup>-2</sup> hr<sup>-1</sup>)",
     "&Delta;BIX <br> (m<sup>-2</sup> hr<sup>-1</sup>)"  ,
     "&Delta;M:C <br> (m<sup>-2</sup> hr<sup>-1</sup>)" ,
     "&Delta;Heterotrophic Bacteria <br> (cells m<sup>-2</sup> hr<sup>-1</sup>)"
   ))) %>%
 mutate(month = factor(ifelse(before_after == "Before", "July", "August"), levels = c("July", "August"))) %>%
   filter(removal_control == "Control")%>%
   filter(name %in% c("do_umol_l","bix","heterotrophic_bacterioplankton_m_l","m_c","nh4_umol_l","nn_umol_l")) %>%
   ggplot(aes(x = foundation_spp, y = rate_m2_hr))+
   geom_hline(yintercept = 0, lty = 2)+
   geom_boxplot(aes(fill = foundation_spp), alpha = 0.3)+
   #  stat_summary(size = 0.7)+
   # geom_point(aes(fill = foundation_spp), shape = 23)+
   labs(x = "",
        y = " ")+
   scale_fill_manual(values = c("black","#34c230"), guide = "none")+
   
   ggh4x::facet_grid2(nicenames2~month, scales = "free_y", independent = "y", switch = "y")+
   theme_bw()+
   theme(strip.background = element_blank(),
         strip.placement = "outside", 
         panel.grid.minor = element_blank(),
         strip.text.y.left = element_markdown(size = 12),
         axis.title.x = element_markdown(),
         legend.position = "bottom", 
         axis.text = element_text(size = 10),
         axis.title = element_text(size=12),
         legend.text = element_text(size = 10))+
   scale_y5

 
rate_box|r1
ggsave(filename = here("Output","Figure_2.pdf"), height = 12, width = 12, device = cairo_pdf)

mean_box|conc1
ggsave(filename = here("Output","Figure_3.pdf"), height = 12, width = 12, device = cairo_pdf)


##### Make the supplemental Plot Reaction Norms
P_July<-data_long %>%
  filter(name %in% c("do_umol_l","nn_umol_l","nh4_umol_l","heterotrophic_bacterioplankton_m_l","bix","m_c"))%>%
  mutate(foundation_spp = case_when( foundation_spp == "Ocean"~ "Ocean",
                                     foundation_spp == "Mytilus" ~ "Mussels",
                                     foundation_spp == "Phyllospadix" ~"Surfgrass")) %>%
  mutate(linetype = ifelse(foundation_spp == "Ocean", "dashed","solid"))%>% # make different lines
    filter(month == "July") %>%
    ggplot(aes(x = time_point_clean, y = mean_val, color = foundation_spp, 
               group = foundation_spp, lty = linetype)
     )+
    geom_point(size = 3)+
    geom_line()+
    geom_errorbar(aes(x = time_point_clean, ymin = mean_val - se_val, ymax = mean_val+se_val), width = 0.1, linetype = "solid")+
  scale_linetype_identity(guide = NULL)+
    facet_wrap(~nicenames, scales = "free_y", strip.position = "left", nrow = 6)+
    facetted_pos_scales(
      y = rep(list(
        scale_y_continuous(limits=c(200, 800)),
        scale_y_continuous(limits=c(0, 30)),
        scale_y_continuous(limits=c(0, 16)),
        scale_y_continuous(limits = c(0.8, 1.2)),
        scale_y_continuous(limits=c(0.9, 1.4)),
      scale_y_continuous(limits = c(0, 1000))
        ), each = 1))+
    labs(x = "",
         y = "",
         color = "",
         title = "Before \n (July)"
         )+
    scale_color_manual(values = c("grey30","#79ACBD","#567d46"))+
    theme_bw()+
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_markdown(size = 14),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 14), 
          plot.title = element_text(hjust = 0.5, size = 14)
    )

P_Aug_control<-data_long %>%
  filter(name %in% c("do_umol_l","nn_umol_l","nh4_umol_l","heterotrophic_bacterioplankton_m_l","bix","m_c"))%>%
  mutate(foundation_spp = case_when( foundation_spp == "Ocean"~ "Ocean",
                                     foundation_spp == "Mytilus" ~ "Mussels",
                                     foundation_spp == "Phyllospadix" ~"Surfgrass")) %>%
  mutate(linetype = ifelse(foundation_spp == "Ocean", "dashed","solid"))%>% # make different lines
  filter(month != "July",
         manipulated == "Not Manipulated"
  ) %>%
  ggplot(aes(x = time_point_clean, y = mean_val, color = foundation_spp, 
             group = foundation_spp, lty = linetype)
  )+
  geom_point(size = 3)+
  geom_line()+
  geom_errorbar(aes(x = time_point_clean, ymin = mean_val - se_val, ymax = mean_val+se_val), width = 0.1, linetype = "solid")+
  scale_linetype_identity(guide = NULL)+
  facet_wrap(~nicenames, scales = "free_y", strip.position = "left", nrow = 6)+
  facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits=c(200, 800)),
      scale_y_continuous(limits=c(0, 30)),
      scale_y_continuous(limits=c(0, 16)),
      scale_y_continuous(limits = c(0.8, 1.2)),
      scale_y_continuous(limits=c(0.9, 1.4)),
      scale_y_continuous(limits = c(0, 1000))
    ), each = 1))+
  labs(x = "",
       y = "",
       color = "",
       title =  "After Control \n (August upwelling)"
  )+
  scale_color_manual(values = c("grey30","#79ACBD","#567d46"))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_markdown(size = 14),
        strip.text = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 14)
  )




P_Aug_impact<-data_long %>%
  filter(name %in% c("do_umol_l","nn_umol_l","nh4_umol_l","heterotrophic_bacterioplankton_m_l","bix","m_c"))%>%
  mutate(foundation_spp = case_when( foundation_spp == "Ocean"~ "Ocean",
                                     foundation_spp == "Mytilus" ~ "Mussels",
                                     foundation_spp == "Phyllospadix" ~"Surfgrass")) %>%
  filter(month != "July",
         manipulated != "Not Manipulated"
  )  %>% bind_rows( # add back the ocean
    data_long %>%
      filter(name %in% c("do_umol_l","nn_umol_l","nh4_umol_l","heterotrophic_bacterioplankton_m_l","bix","m_c")) %>%
      filter(foundation_spp == "Ocean",
             month != "July")
  ) %>%
  mutate(linetype = ifelse(foundation_spp == "Ocean", "dashed","solid"))%>% # make different lines%>%
  ggplot(aes(x = time_point_clean, y = mean_val, color = foundation_spp, 
             group = foundation_spp, lty = linetype)
  )+
  geom_point(size = 3)+
  geom_line()+
  geom_errorbar(aes(x = time_point_clean, ymin = mean_val - se_val, ymax = mean_val+se_val), width = 0.1, linetype = "solid")+
  scale_linetype_identity(guide = NULL)+
  facet_wrap(~nicenames, scales = "free_y", strip.position = "left", nrow = 6)+
  facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits=c(200, 800)),
      scale_y_continuous(limits=c(0, 30)),
      scale_y_continuous(limits=c(0, 16)),
      scale_y_continuous(limits = c(0.8, 1.2)),
      scale_y_continuous(limits=c(0.9, 1.4)),
      scale_y_continuous(limits = c(0, 1000))
    ), each = 1))+
  labs(x = "",
       y = "",
       color = "",
       title =  "After Impact \n (August upwelling)"
  )+
  scale_color_manual(values = c("grey30","#79ACBD","#567d46"))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_markdown(size = 14),
        strip.text = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 14)
  )

(P_July|P_Aug_control|P_Aug_impact)+plot_layout(guides = "collect")&theme(legend.position = "bottom")
ggsave(here("Output","Supp_Fig_1.pdf"), width = 10, height = 12)

### Does pool size affect any of the results
Long_all%>% 
  filter(foundation_spp != "Ocean")%>%
  mutate(removal = case_when(removal_control == "Control"~"Unmanipulated",
                             removal_control == "Removal"& before_after == "Before" ~ "Unmanipulated",
                             removal_control == "Removal"& before_after != "Before" ~ "Foundation spp. removed"))%>%
  mutate(removal = factor(removal, levels = c("Unmanipulated","Foundation spp. removed"))) %>%
  filter(name %in% c("bix","do_umol_l","heterotrophic_bacterioplankton_m_l","m_c","nh4_umol_l","nn_umol_l")) %>%
  group_by(name, foundation_spp, removal)%>%
  nest() %>%
  mutate(model = map(data, 
                     function(df) {
                       lm(rate_sqrt  ~ surface_area, #sqrt transformed
                          data = df)
                     })) %>%
  
  mutate(
    tidy = map(model, function(x)tidy(x,effects="fixed")),
    glance = map(model, glance)
  )  %>%
  unnest(tidy) %>%
  select(name, term, statistic,p.value) %>%
  ungroup() %>%
  filter(term == "surface_area") %>%
  filter(p.value < 0.05) # unmanipulated mussel M_C is the only significant


MetaData %>%
  filter(Before_After == "Before") %>%
  ggplot(aes(x = Foundation_spp, y = SurfaceArea))+
  geom_point(alpha = 0.25)+
  stat_summary(size =1)+
  labs(x = "",
       y = "Surface Area (m<sup>2)")+
  facet_wrap(~Removal_Control)+
  theme_minimal()+
  theme(axis.title.y = element_markdown(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14))
ggsave(here("Output","Supp_Fig2.pdf"), width = 5, height = 5)
