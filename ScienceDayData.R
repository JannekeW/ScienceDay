## Make data frame Science Day


## We have 100 strains (1:100)
## There are 6 phenotypes
## For each phenotype there are 3 mutations of relevance = 18 in total
# ULTIMATELY 66 SHOULD BE THE BEST STRAIN

##PACKAGES
library('ggbeeswarm')
library(wakefield)

Sample <- 1:100
Draught <- c(round(runif(n=50, min=50, max=90), 0),round(runif(n=50, min=80, max=120), 0))
Locating <- round(runif(n=100, min=300, max=400), 0)
infection_even <- round(runif(n=50, min=2, max=20), 0)
infection_uneven <- round(runif(n=50, min=18, max=40), 0)
cold_low <- round(runif(n=50, min=-10, max=-2), 1)
cold_high <- round(runif(n=50, min=-3, max=5), 1)
mobility_low <- round(runif(n=91, min=50, max=100), 0)
mobility_high <- round(runif(n=9, min=120, max=140), 0)
heat_high <- round(runif(n=50, min=37, max=38.5),1)
heat_low <- round(runif(n=50, min=32, max=38),1)
pB1 = upper_factor(100, k = 2, x = c("T","C"), prob = NULL, name = "Upper")
pB2 = upper_factor(100, k = 2, x = c("T","G"), prob = NULL, name = "Upper")
pB3 = upper_factor(100, k = 2, x = c("A","C"), prob = NULL, name = "Upper")
pD2 = upper_factor(100, k = 2, x = c("T","G"), prob = NULL, name = "Upper")
pD3 = upper_factor(100, k = 2, x = c("G","C"), prob = NULL, name = "Upper")
pF2 = upper_factor(100, k = 2, x = c("A",  "C"), prob = NULL, name = "Upper")
pF3 = upper_factor(100, k = 2, x = c("C",  "T"), prob = NULL, name = "Upper")

Science.Day.Base <- data_frame(Sample, Draught, Locating)

Science.Day <- Science.Day.Base %>% 
    dplyr::mutate(Infection = ifelse(row_number()%%2, infection_even , infection_uneven) ,
                  Cold = ifelse(row_number() %in% c(1:25, 51:75) , cold_low,  cold_high),
                  Mobility = ifelse(row_number() %in% c(11,22,33,44,55,66,77,88,99) , mobility_high,  mobility_low),
                  Heat = ifelse(row_number() %in% c(1:10,21:30,41:50,61:70,81:90) , heat_high,  heat_low),
                  Position_1 =  ifelse(row_number() %in% c(1:50) , "A",  "G"), ## first number is if false, second if true
                  Position_2 =  ifelse(row_number() %in% c(1:40,51:60) , "C",  "T"),
                  Position_3 =  ifelse(row_number() %in% c(1:48,59:60) , "T",  "G"),
                  Position_7 =  pB1,
                  Position_8 =  pB2,
                  Position_9 =  pB3,
                  Position_16 =  ifelse(row_number()%%2, "G" , "T"),
                  Position_17 =  ifelse(row_number()%%2, "A" , "C"),
                  Position_18 =  ifelse(row_number()%%2, "C" , "G"),
                  Position_18 =  ifelse(row_number() %in% c(12,24,36,48,52,64,76,88) , "C",  Position_18),
                  Position_18 =  ifelse(row_number() %in% c(11,23,35,47,51,63,75,87) , "G",  Position_18),
                  Position_4 =  ifelse(row_number() %in% c(1:25,51:75), "A" , "C"),
                  Position_5 =  pD2,
                  Position_6 =  pD3,
                  Position_10 =  ifelse(row_number() %in% c(11,22,33,44,55,66,77,88,99) , "T", "A"),
                  Position_11 =  ifelse(row_number() %in% c(11,22,33,44,55,66,77,88,99) , "T", "C"),
                  Position_12 =  ifelse(row_number() %in% c(11,22,33,44,55,66,77,88,99) , "G", "T"),
                  Position_13 =  ifelse(row_number() %in% c(1:10,21:30,41:50,61:70,81:90) , "C",  "T"),
                  Position_14 =  pF2,
                  Position_15 =  pF3
  )
View(Science.Day[66,])
write.table(Science.Day, file='~/Desktop/ScienceDay/ScienceDayData.tsv', quote=FALSE, sep='\t', col.names = NA)

SAVE_PRES <- Science.Day %>% ggplot(aes(x=posF1, y=Heat, alpha = sample, fill = sample)) +
  #geom_beeswarm(priority='ascending', color = "blue") +
  geom_quasirandom(color = "orange", width = .4, size = 2.5, shape = 16)+
  scale_x_discrete(name = "Chromosome 7 - position 1") +
  theme(panel.background = element_rect(fill = NA, colour = "white", size = 0.5),
        panel.border = element_rect(color = "black",fill=NA),
        panel.grid.major.y = element_line(color = "grey", linetype= 3),
        legend.position="none", 
        axis.text = element_text(size = 10),
        axis.title.x = element_text( size=12, face = "bold"), 
        axis.title.y = element_text( size=12, face = "bold", vjust = 1.5), 
  )

#?geom_beeswarm()
SAVE_PRES
ggsave(paste0("~/Desktop/ScienceDay/Plots/PhenGen.png"), SAVE_PRES, height = 5, width = 5)

## Make students pick color and trait, as well as position (but link positions to trait)


############MANHATTAN


# Load the library
install.packages("qqman")
library(qqman)

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P" )
head(gwasResults)
snpsOfInterest


##MAKE 9 CHROMOSOME DATASET, replace the peak
Play <- gwasResults[gwasResults$CHR %in% c(1:9),]
NewVal <- runif(n=100, min=0, max=1.4)
Play <- Play %>%
  mutate(P = ifelse(SNP %in%  snpsOfInterest, NewVal, P))

## Original MAnhattan plot from package
manhattan(Play, chr="CHR", bp="BP", snp="SNP", p="P" )+
  theme(panel.background = element_rect(fill = NA, colour = "white", size = 0.5),
        panel.border = element_rect(color = "black",fill=NA),
        panel.grid.major.y = element_line(color = "grey", linetype= 3),
        legend.position="none", 
        axis.text = element_text(size = 10),
        axis.title.x = element_text( size=12, face = "bold"), 
        axis.title.y = element_text( size=12, face = "bold", vjust = 1.5), 
  )

## Change order of chromosomes ... not needed!
TypesMan <- Play %>%
  dplyr::mutate(CHR = case_when(CHR == 1 ~ 10,
                                CHR == 2 ~ 4,
                                CHR == 3 ~ 6,
                                CHR == 4 ~ 3,
                                CHR == 5 ~ 9,
                                CHR == 6 ~ 8,
                                CHR == 7 ~ 5,
                                CHR == 8 ~ 2,
                                CHR == 9 ~ 7)) %>%
  dplyr::mutate(CHR = case_when(CHR == 10 ~ 1,
                                CHR != 10 ~ CHR))
manhattan(TypesMan, chr="CHR", bp="BP", snp="SNP", p="P" )

max(Play$BP[Play$CHR == 1])
max(Play$BP[Play$CHR == 2])
max(Play$BP[Play$CHR == 3])
max(Play$BP[Play$CHR == 4])
max(Play$BP[Play$CHR == 5])
max(Play$BP[Play$CHR == 6])
max(Play$BP[Play$CHR == 7])
max(Play$BP[Play$CHR == 8])
max(Play$BP[Play$CHR == 9])


#Color scheme for chromosomes: black-grey repeat
ChromCol <- c("1" = "Black","2" = "Grey","3" = "Black", "4" = "Grey", '5' = "Black", "6" = "Grey", "7" = "Black", "8" = "Grey", "9" = "Black")


LocateMan <- Play

HeatMan <- Play %>%
  dplyr::mutate(CHR = case_when(CHR == 1 ~ 10,
                                CHR == 2 ~ 4,
                                CHR == 3 ~ 6,
                                CHR == 4 ~ 3,
                                CHR == 5 ~ 9,
                                CHR == 6 ~ 8,
                                CHR == 7 ~ 5,
                                CHR == 8 ~ 2,
                                CHR == 9 ~ 7)) %>%
  dplyr::mutate(CHR = case_when(CHR == 10 ~ 1,
                                CHR != 10 ~ CHR)) %>% 
  dplyr::mutate(CHR = as_factor(CHR)) %>%
  dplyr::arrange(CHR) %>%
  dplyr::mutate(ORDER = row_number()) %>%
  dplyr::filter(-log(P)>0,
                -log(P)<5) 

random_highs <- c(runif(n=50, min=0, max=0.02), runif(n=50, min=0, max=0.01),runif(n=100, min=0.02, max=0.3)  ) 
PEAK <- data_frame("P" = random_highs,CHR = as_factor(2), ORDER = min(HeatMan[HeatMan$CHR == 2,5]) + 100.5+ c(1:length(random_highs)))

HeatMan <- HeatMan %>%
  dplyr::bind_rows(PEAK)

HeatMan %>%
  ggplot(aes(x = ORDER, y = -log10(P), axes = FALSE, colour = CHR)) +
  geom_point(size = 0.5) +
  scale_x_continuous(name = "Chromosome", breaks=NULL, labels = NULL, expand =c(0.02,0.02)) +
  scale_colour_manual(values = ChromCol,  aesthetics = c("colour", "fill"))+
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.4), "cm")) +
  theme(panel.background = element_rect(fill = NA, colour = "white", size = 0.5),
        panel.border = element_rect(color = "black",fill=NA),
        panel.grid.major.y = element_line(color = "grey", linetype= 3),
        legend.position="none", 
        axis.text = element_text(size = 10),
        axis.title.x = element_text( size=12, face = "bold"), 
        axis.title.y = element_text( size=12, face = "bold", vjust = 1.5), 
  )
