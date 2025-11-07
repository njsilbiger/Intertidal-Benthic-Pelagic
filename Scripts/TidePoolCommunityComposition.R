##Community Composition data for Sessile and Mobiles from Surfgrass and Mussel
##Oregon Tide Pools
##By: Jenn Fields
##Last updated: 4.26.2021


#load libraries
library(tidyverse)

#load data from community comp surveys
#source scripts
source("scripts/tidepoolphysicalparameters.R")
#data
Sessiles <- read_csv("Data/CommunityComposition/SessilesAll.csv")
Mobiles <- read_csv("Data/CommunityComposition/Mobiles.csv")
SessilesGroupings<-read_csv("Data/CommunityComposition/SessilesFun_withcitations.csv")
MobileGroupings<-read_csv("Data/CommunityComposition/MobilesFun_withcitations.csv")

#replace NA with 0
Mobiles[is.na(Mobiles)]<-0
#remove immediate time pt
Mobiles<-Mobiles %>%
  filter(Before_After !="Immediate")
#write.csv(Mobiles,file="Data/CommunityComposition/Mobilespplist.csv") #save spp list without mid point sampling

#######Sessiles#########
#convert characters to numeric in sessile sheet
Sessiles$Epiactis.prolifera<-as.numeric(Sessiles$Epiactis.prolifera)
Sessiles$Chaetomorpha.linum<-as.numeric(Sessiles$Chaetomorpha.linum)
Sessiles$Costaria.costata<-as.numeric(Sessiles$Costaria.costata)
Sessiles[is.na(Sessiles)]<-0 

# Make all the community data a relative percent
PercentSessile<-100*Sessiles[7:ncol(Sessiles)]/Sessiles$Squares #change to rock--end spp
#View(PercentSessile)
#normalize to the sum of the total cover (since it can be greater than 100%)
StandardizedSessile<- 100*PercentSessile/rowSums(PercentSessile)


Communitymetrics<-cbind(Sessiles$PoolID, Sessiles$Foundation_spp, Sessiles$Removal_Control, 
                        Sessiles$Before_After,StandardizedSessile,PercentSessile$Mytilus.californianus,PercentSessile$Phyllospadix.spp)

#adjusted mussel and surfgrass are made relative sum of the total cover of the pool whereas Surfgrass and Mussel cover are not by size
Communitymetrics <- Communitymetrics %>%
  rename(PoolID = "Sessiles$PoolID", Foundation_spp = "Sessiles$Foundation_spp",Removal_Control ="Sessiles$Removal_Control",
         Before_After ="Sessiles$Before_After",
         MusselCover = "PercentSessile$Mytilus.californianus",
         SurfgrassCover = "PercentSessile$Phyllospadix.spp",
         AdjSurfgrassCover = "Phyllospadix.spp", # these the standardized metrics,
         AdjMusselCover = "Mytilus.californianus" # standardized values
         ) %>% #mussel and surfgrass cover are not standardized to pool 
  dplyr::filter(Before_After != "Immediate")
#rename joined columns

#write.csv(Communitymetrics,file="Data/CommunityComposition/Sessilesspplist.csv") #save updated spp list without mid point sampling

#merge with functional groups
###phylum, class and functional group
SessilesFun<-Communitymetrics %>%
  dplyr::group_by(PoolID, Foundation_spp, Before_After, Removal_Control) %>%
  tidyr::pivot_longer(
    cols = AdjMusselCover:Stylantheca.spp, #leave out unadjusted surfgrass and mussel cover
    names_to = "Species", #creates column with species in longformat
    values_to = "Cover", #adds column for % cover
    values_drop_na = TRUE
  ) 

SessilesFun<-left_join(SessilesFun,SessilesGroupings) #combine with fun groups 
write.csv(SessilesFun,file="Data/CommunityComposition/Sessilefunctionalgroups.csv")

MobilesFun<-Mobiles %>%
  dplyr::filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID, Foundation_spp, Before_After, Removal_Control) %>%
  tidyr::pivot_longer(
    cols = Nuttalina.spp:Gunnel,
    names_to = "Species", #creates column with species in longformat
    values_to = "Count", #adds column for % cover
    values_drop_na = TRUE
  ) 

MobilesFun<-left_join(MobilesFun,MobileGroupings)

#save spp with functional group categories
#write.csv(Mobiles,file="Data/CommunityComposition/Mobilespp.csv")
#write.csv(MobilesFun,file="Data/CommunityComposition/Mobilefunctionalgroups.csv")


#summarized groups of sessile community shift between before and after
#This is TP Sessile Community Metrics Dataset in Data folder
Summarizedgroups<-Communitymetrics %>%
  dplyr::filter(Before_After != "Immediate") %>%
  dplyr:: mutate(allCCA = (Crustose.coralline + Bosiella.spp + Corallina.spp + Corallina.vancouveriensis +
                             Calliarthron.tuberculosum), #creates column for all CCA
                 Diatoms = Diatoms, # pull out the diatoms
                 macroalgae = (Diatoms + Algae.film + Turf.algae +Acrosiphonia.coalita + Analipus.japonicus + Chaetomorpha.linum +
                                 Centroceras.Ceramium + Cladophora.columbiana + Cryptosiphonia.wooddii + Costaria.costata +
                                 Cumagloia.andersonii + Erythrophyllum.delessrioides +	Endocladia.muricata	+
                                 Fucus.gardneri +	Halosaccion.glandiforme	+ Pyropia.spp	+ Polysiphonia.spp +
                                 Ptilota.spp +	Cryptopluera.spp + Odonthalia.floccosa +	Microcladia.borealis +	Mazzaella.splendens +
                                 Mazzaella.flaccida	+ Mazzaella.oregona +	Neorhodomela.larix + Odanthalia.washingtoniensis +Scytosiphon.lomentaria	
                               + Osmundea.spectabilis	+ Ulva.spp + Plocamium.pacificum + Smithura.naiadum + Mastocarpus	+ Farlowia.mollis	+
                                 Savoiea.robusta +	Palmaria.hecatensis	+ Melobesia.mediocris +	Leathesia.marina + Callithamnion.pikeanum	+
                                 Laminara.setchellii + Schizymenia.pacifica),
                 macrophytes = (AdjSurfgrassCover + macroalgae), #includes phyllospadix & macroalgae
                 macroCCA = (macroalgae + allCCA), #includes macroalgae and CCA
                 consumers = (Chthamalus +	Semibalanus.cariosus +	Balanus.nibulis	+ Balanus.glandula +
                                Pollicipes.polymerus +	tube.worm	+ Ophlitaspongia.pennata + Halichondria +	
                                Haliclona.permollis	+ Anthropluera.elegantissima	+ Anthropluera.xanthogrammica	+ 
                                Urticina.coriacea	+ Epiactis.prolifera +Anthopleura.artemisia	+ Stylantheca.spp),
                 #includes all consumers (no mytilus)
                 allconsumers = (consumers + AdjMusselCover), #consumers and mytilus
                 prodphyllodom = (macroalgae - (consumers + AdjMusselCover)), #producer dominance for phyllo model (so subtract mytilus too)
                 allproddom = (macrophytes - allconsumers)) %>%#prod dominance with foundation spp
  dplyr::select(PoolID,Foundation_spp,Removal_Control,Before_After,AdjMusselCover, AdjSurfgrassCover, MusselCover, SurfgrassCover,Diatoms, allCCA, macroalgae,macrophytes, macroCCA,
                consumers,allconsumers,  prodphyllodom, allproddom)

#write.csv(Summarizedgroups,file="Data/TPSessileCommunityMetrics.csv")
write_csv(Summarizedgroups,here("Data","Microbe_Clean","CommunityData.csv"))

#Change in foundation species cover
Funsppcover<- Communitymetrics%>%
  filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID,Foundation_spp, Removal_Control) %>%
  summarise(Mytilusdelta = -1*(MusselCover[Before_After == 'After'] - MusselCover[Before_After == 'Before']),
            Phyllodelta = -1*(SurfgrassCover[Before_After == 'After'] - SurfgrassCover[Before_After == 'Before']))

#average between time periods for physical parameters
PP<-TidePooldes %>%
  dplyr::group_by(PoolID,Removal_Control) %>%
  dplyr::summarise(SAVav = mean(SAtoV), #ave between before and after since SA/V changed with fspp removal
                   THav = mean(TideHeight),SAav=mean(SurfaceArea),Vav=mean(Vol),Depthav=mean(MaxDepth),
                   loggerdepth=mean(LoggerDepth)) #tide height didn't change 

PP$PoolID<-as.factor(PP$PoolID)
Funsppcover$PoolID<-as.factor(Funsppcover$PoolID)
Funsppandpp<-left_join(Funsppcover,PP)
#write.csv(Funsppandpp, file="Data/Fspplossandphysicalparameters.csv")
