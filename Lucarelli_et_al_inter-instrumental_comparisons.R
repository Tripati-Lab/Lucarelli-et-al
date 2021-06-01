# Inter-instrumental comparisons
library(tidyverse)
library(readxl)
library(nlme)
library(emmeans)
library(pwr)
library(rstatix)
library(pracma) # Note - the pracma implementation of standard error seems to be less finicky than the base function


# D48

sheets <- excel_sheets("Final D48 Data Config 2 All Standards.xlsx")

list_all <- lapply(sheets, function(x) read_excel("Final D48 Data Config 2 All Standards.xlsx", sheet = x))

str(list_all)

finalconfig2dat <- do.call("rbind", list_all)
str(finalconfig2dat)

sheets2 <- excel_sheets("Final D48 Data Config 3 All Standards.xlsx")

list_all2 <- lapply(sheets2, function(x) read_excel("Final D48 Data Config 3 All Standards.xlsx", sheet = x))

finalConfig3dat <- do.call("rbind", list_all2)


sheets3 <- excel_sheets("Final Data Config 1a D48.xlsx")

list_all3 <- lapply(sheets3, function(x) read_excel("Final Data Config 1a D48.xlsx", sheet = x))

finalConfig1adat <- do.call("rbind", list_all3)

allD48dat <- merge(finalConfig2dat, finalConfig3dat, all = TRUE)
allD48dat <- merge(allD48dat, finalConfig1adat, all = TRUE)

# D48 comparisons

allD48dat <- allD48dat %>%
  group_by(Mass_Spec, Standard) %>%
  mutate(D48CDESse = std.error(D48CDES_Final),
         D47CDESse = std.error(D47CDES_Final))

# Set up lme
lme_D48 <- lme(D48CDES_Final ~ Standard + Mass_Spec, 
               data=allD48dat, random = ~1|D48CDESse, method = "REML")
summary(lme_D48)

# Get pairwise contrasts
modpairwise <- emmeans(lme_D48, pairwise ~ Mass_Spec)



# All D47 - easier to pull from a truncated file of the supplementary tables

sheets <- excel_sheets("Supplementary Tables.xlsx")

list_all <- lapply(sheets, function(x) read_excel("Supplementary Tables.xlsx", sheet = x))

str(list_all)

finaldatall <- do.call(dplyr::bind_rows, list_all)
str(finaldatall)

names(finaldatall)

finaldatall <- finaldatall %>%
rename(
  D48.CDES.90  = `Δ48 CDES 90 (‰)`,
  D47.ICDES = `Δ47I-CDES (‰)`,
  D47.CDES.90 = `Δ47CDES 90 (‰)`,
  D47.CDES.70 = `Δ47CDES 70 (‰)`,
  MassSpec = `Mass Spectrometer`,
  SampleName = `Sample Name`
)

finaldatall <- finaldatall %>%
  mutate(Standard = case_when(
    !grepl("ETH", SampleName) ~ SampleName,
    grepl("ETH-1", SampleName) ~ "ETH-1",
    grepl("ETH-2", SampleName) ~ "ETH-2", 
    grepl("ETH-3", SampleName) ~ "ETH-3", 
    grepl("ETH-4", SampleName) ~ "ETH-4"
  )
  )

finaldatall <- finaldatall %>%
  group_by( MassSpec, Standard, Analysis) %>%
  mutate(D48CDES90se = pracma::std_err(D48.CDES.90),
         D47ICDESse = pracma::std_err(D47.ICDES),
         D47CDES90se = pracma::std_err(D47.CDES.90),
         D47CDES70se = pracma::std_err(D47.CDES.70)
  )

D47ICDES <- finaldatall[finaldatall$Analysis == "D47 I-CDES",]
D47CDES90 <- finaldatall[finaldatall$Analysis == "D47 CDES 90",]
D48CDES90 <- finaldatall[finaldatall$Analysis == "D48 CDES 90",]

unique(finaldatall$Analysis)

# Set up lme
lme_D47ICDES <- lme(D47.ICDES ~ MassSpec+Standard, 
               data=D47ICDES, random = ~1|D47ICDESse, method = "REML", na.action = na.omit)
summary(lme_D47ICDES)

modpairwise1 <- emmeans(lme_D47ICDES, pairwise ~ MassSpec)
modpairwise2 <- emmeans(lme_D47ICDES, pairwise ~ MassSpec|Standard)

plot(modpairwise1,comparisons = TRUE, adjust = "mvt", 
     horizontal = FALSE, colors = "darkgreen") + theme_bw() +
  xlab("Estimated marginal means") + ylab("")

ggsave("D47ICDES_emmeansplot.tiff", dpi = 800, compression = "lzw")

write.csv(modpairwise1, file = "D47ICDES_emmeans.csv", row.names = FALSE)

lme_D48CDES90 <- lme(D48.CDES.90 ~ MassSpec+Standard, 
                    data=D48CDES90, random = ~1|D48CDES90se, method = "REML", na.action = na.omit)
summary(lme_D48CDES90)

modpairwise3 <- emmeans(lme_D48CDES90, pairwise ~ MassSpec)
modpairwise4 <- emmeans(lme_D48CDES90, pairwise ~ Standard+MassSpec)

write.csv(modpairwise4, file = "D48CDES90_emmeans_bystandard.csv", row.names = FALSE)

plot(modpairwise3 ,comparisons = TRUE, adjust = "mvt", 
     horizontal = FALSE, colors = "darkgreen") + theme_bw()+
  xlab("Estimated marginal means") + ylab("")
ggsave("D48CDES90_emmeansplot.tiff", dpi = 800, compression = "lzw", scale = 1.25)

plot(modpairwise4 ,comparisons = TRUE, adjust = "mvt", 
     horizontal = FALSE, colors = "darkgreen") + theme_bw()+
  xlab("Estimated marginal means") + ylab("")
ggsave("D48CDES90_emmeansplot.tiff", dpi = 800, compression = "lzw", scale = 1.25)
# 
lme_D47CDES90 <- lme(D47.CDES.90 ~ Standard, 
                     data=D47CDES90, random = ~1|D47CDES90se, method = "REML", na.action = na.omit)
summary(lme_D48CDES90)

modpairwise5 <- emmeans(lme_D48CDES90, pairwise ~ Standard)


# Power analysis

# ETH-1
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-1"& D48CDES90$MassSpec == "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-1"& D48CDES90$MassSpec == "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-1"& D48CDES90$MassSpec == "Configuration 1a"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 40, length = 44

# ETH-2
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-2"& D48CDES90$MassSpec == "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-2"& D48CDES90$MassSpec == "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-2"& D48CDES90$MassSpec == "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 116, length = 38

# ETH-3
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-3"& D48CDES90$MassSpec == "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-3"& D48CDES90$MassSpec == "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-3"& D48CDES90$MassSpec == "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 4, length = 45

# ETH-4
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-4"& D48CDES90$MassSpec == "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-4"& D48CDES90$MassSpec == "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-4"& D48CDES90$MassSpec == "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 21, length = 45

# TV03
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "TV03"& D48CDES90$MassSpec == "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "TV03"& D48CDES90$MassSpec == "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "TV03"& D48CDES90$MassSpec == "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 4, length = 55

# Veinstrom
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Veinstrom"& D48CDES90$MassSpec == "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Veinstrom"& D48CDES90$MassSpec == "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Veinstrom"& D48CDES90$MassSpec == "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 7, length = 74

# Carmel Chalk
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carmel Chalk"& D48CDES90$MassSpec == "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carmel Chalk"& D48CDES90$MassSpec == "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carmel Chalk"& D48CDES90$MassSpec == "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 4, length = 71

# Carrara Marble
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carrara Marble"& D48CDES90$MassSpec == "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carrara Marble"& D48CDES90$MassSpec == "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carrara Marble"& D48CDES90$MassSpec == "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 29, length = 64




# Config 1b and 2


# ETH-1
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-1"& D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-1"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-1"& D48CDES90$MassSpec != "Configuration 1a"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 11, length = 652

# ETH-2
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-2"& D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-2"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-2"& D48CDES90$MassSpec != "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 17, length = 643

# ETH-3
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-3"& D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-3"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-3"& D48CDES90$MassSpec != "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 4, length = 381

# ETH-4
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-4"& D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-4"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-4"& D48CDES90$MassSpec != "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 6, length = 428

# TV03
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "TV03"& D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "TV03"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "TV03"& D48CDES90$MassSpec != "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 3, length = 90

# Veinstrom
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Veinstrom"& D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Veinstrom"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Veinstrom"& D48CDES90$MassSpec != "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 4, length = 629

# Carmel Chalk
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carmel Chalk"& D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carmel Chalk"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carmel Chalk"& D48CDES90$MassSpec != "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 4, length = 485

# Carrara Marble
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carrara Marble"& D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carrara Marble"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carrara Marble"& D48CDES90$MassSpec != "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 20, length 215

# CMTile
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "CMTile"& D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "CMTile"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "CMTile"& D48CDES90$MassSpec != "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 13, length = 453

# MERCK
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "MERCK" & D48CDES90$MassSpec != "Configuration 1a"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "MERCK"& D48CDES90$MassSpec != "Configuration 1a"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "MERCK"& D48CDES90$MassSpec != "Configuration 1a"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 4, length = 70


# Configuration 3


# ETH-1
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-1"& D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-1"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-1"& D48CDES90$MassSpec == "Configuration 3"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 250, length 188

# ETH-2
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-2"& D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-2"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-2"& D48CDES90$MassSpec == "Configuration 3"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 116, length = 204

# ETH-3
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-3"& D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-3"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-3"& D48CDES90$MassSpec == "Configuration 3"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 8, length = 145

# ETH-4
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-4"& D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-4"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "ETH-4"& D48CDES90$MassSpec == "Configuration 3"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 31, length = 171

# TV03
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "TV03"& D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "TV03"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "TV03"& D48CDES90$MassSpec == "Configuration 3"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 6, length = 32

# Veinstrom
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Veinstrom"& D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Veinstrom"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Veinstrom"& D48CDES90$MassSpec == "Configuration 3"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 7, length = 193

# Carmel Chalk
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carmel Chalk"& D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carmel Chalk"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carmel Chalk"& D48CDES90$MassSpec == "Configuration 3"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 9, length = 166

# Carrara Marble
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carrara Marble"& D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carrara Marble"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "Carrara Marble"& D48CDES90$MassSpec == "Configuration 3"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 3116, length = 80

# CMTile
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "CMTile"& D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "CMTile"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "CMTile"& D48CDES90$MassSpec == "Configuration 3"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 77, length = 144

# MERCK
pwr.t.test(d = (mean(D48CDES90$D48.CDES.90[D48CDES90$Standard == "MERCK" & D48CDES90$MassSpec == "Configuration 3"])-sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "MERCK"& D48CDES90$MassSpec == "Configuration 3"]))/sd(D48CDES90$D48.CDES.90[D48CDES90$Standard == "MERCK"& D48CDES90$MassSpec == "Configuration 3"]), 
           sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 31103.84, length = 11

# Normality of all final datasets

sheets <- excel_sheets("Supplementary Tables.xlsx")

list_all <- lapply(sheets, function(x) read_excel("Supplementary Tables.xlsx", sheet = x))

str(list_all)

finaldatall <- do.call(dplyr::bind_rows, list_all)
str(finaldatall)

names(finaldatall)

finaldatall <- finaldatall %>%
  mutate(Standard = case_when(
    !grepl("ETH", `Sample Name`) ~ `Sample Name`,
    grepl("ETH-1", `Sample Name`) ~ "ETH-1",
    grepl("ETH-2", `Sample Name`) ~ "ETH-2", 
    grepl("ETH-3", `Sample Name`) ~ "ETH-3", 
    grepl("ETH-4", `Sample Name`) ~ "ETH-4"
  )
  )


names(finaldatall)<-make.names(names(finaldatall),unique = TRUE) # Standardize names

# split by analysis

D48CDES90 <- finaldatall %>% filter(Analysis == "Δ48 CDES 90")
D47CDES90 <- finaldatall %>% filter(Analysis == "Δ47 CDES 90")
D47ICDES <- finaldatall %>% filter(Analysis == "Δ47 I-CDES")

unique(D48CDES90$Mass.Spectrometer)
unique(finaldatall$Standard)


# D48 CDES 90
D48CDES90_summary <- D48CDES90 %>%
  group_by(`Mass.Spectrometer`, Standard) %>%
  summarize(`Mean Δ48 CDES 90 (‰)` = round(mean(`Δ48.CDES.90....`), 3),
            `Min Δ48 CDES 90 (‰)` = min(`Δ48.CDES.90....`),
            `Max Δ48 CDES 90 (‰)` = max(`Δ48.CDES.90....`),
            `SD Δ48 CDES 90 (‰)` = round(sd(`Δ48.CDES.90....`), 3),
            `N Δ48 CDES 90 (‰)` = length(`Δ48.CDES.90....`),
                        na.action(na.omit)
            
            )

D48CDES90_shap <- D48CDES90 %>%
  group_by(`Mass.Spectrometer`, Standard) %>%
  filter(n() >=3) %>%
  shapiro_test(`Δ48.CDES.90....`)
D48CDES90_summary_final <- merge(D48CDES90_summary, D48CDES90_shap, all = TRUE)

write.csv(D48CDES90_summary_final, "D48CDES90_summary_final2.csv", row.names = FALSE)

D48CDES90 <- D48CDES90 %>%
  group_by(`Mass.Spectrometer`, Standard) %>%
  mutate(D48CDES90se = pracma::std_err(`Δ48.CDES.90....`), na.action(na.omit))

lme_D48CDES90 <- lme(`Δ48.CDES.90....` ~ Mass.Spectrometer + Standard, 
                    data=D48CDES90[D48CDES90$Mass.Spectrometer != "Multiple" & !is.na(D48CDES90$D48CDES90se),], random = ~1|D48CDES90se , method = "REML", na.action = na.omit)
summary(lme_D48CDES90)

modpairwise3 <- emmeans(lme_D48CDES90, pairwise ~ `Mass.Spectrometer`)


plot(modpairwise3,comparisons = TRUE, adjust = "mvt", 
     horizontal = FALSE, colors = "darkgreen") + theme_bw() +
  xlab("Estimated marginal means") + ylab("")

ggsave("D48CDES90_emmeansplot.tiff", dpi = 800, compression = "lzw", scale = 1.25)

write.csv(modpairwise3, file = "D48CDES90_emmeans.csv", row.names = FALSE)

# D47 CDES 90
D47CDES90_summary <- D47CDES90 %>%
  group_by(`Mass.Spectrometer`, Standard) %>%
  summarize(`Mean Δ47 CDES 90 (‰)` = round(mean(`Δ47CDES.90....`), 3),
            `Min Δ47 CDES 90 (‰)` = min(`Δ47CDES.90....`),
            `Max Δ47 CDES 90 (‰)` = max(`Δ47CDES.90....`),
            `SD Δ47 CDES 90 (‰)` = round(sd(`Δ47CDES.90....`), 3),
            `N Δ47 CDES 90 (‰)` = length(`Δ47CDES.90....`),
            na.action(na.omit)
            
  )

D47CDES90_shap <- D47CDES90 %>%
  group_by(`Mass.Spectrometer`, Standard) %>%
  filter(n() >=3) %>%
  shapiro_test(`Δ47CDES.90....`)

# Make a nice summary table
D47CDES90_summary_final <- merge(D47CDES90_summary, D47CDES90_shap, all = TRUE)

write.csv(D47CDES90_summary_final, "D47CDES90_summary_final.csv", row.names = FALSE)



# D47 ICDES
D47ICDES_summary <- D47ICDES %>%
  group_by(`Mass.Spectrometer`, Standard) %>%
  summarize(`Mean Δ47 ICDES (‰)` = round(mean(`Δ47I.CDES....`), 3),
            `Min Δ47 ICDES (‰)` = min(`Δ47I.CDES....`),
            `Max Δ47 ICDES (‰)` = max(`Δ47I.CDES....`),
            `SD Δ47 ICDES (‰)` = round(sd(`Δ47I.CDES....`), 3),
            `N Δ47 ICDES (‰)` = length(`Δ47I.CDES....`),
            na.action(na.omit)
            
  )

D47ICDES_shap <- D47ICDES %>%
  group_by(`Mass.Spectrometer`, Standard) %>%
  filter(n() >=3) %>%
  shapiro_test(`Δ47I.CDES....`)

D47ICDES <- D47ICDES %>%
  group_by(`Mass.Spectrometer`, Standard) %>%
  mutate(D47ICDESse = pracma::std_err(`Δ47I.CDES....`), na.action(na.omit))

# Make a nice summary table
D47ICDES_summary_final <- merge(D47ICDES_summary, D47ICDES_shap, all = TRUE)

write.csv(D47ICDES_summary_final, "D47ICDES_summary_final2.csv", row.names = FALSE)

# Shorten the names
D47ICDES$MassSpec2 <- recode(D47ICDES$`Mass.Spectrometer`, "Config 1a" = "C1a", "Config 1b" = "C1b", "Config 1c" = "C1c", "Config 1d" = "C1d", "Config 2" = "C2",
                              "Config 3" = "C3") 
names(D47ICDES)

# Set up lme
lme_D47ICDES <- lme(`Δ47I.CDES....` ~ `Mass.Spectrometer`+Standard, 
                    data=D47ICDES[D47ICDES$Mass.Spectrometer != "Multiple" & !is.na(D47ICDES$D47ICDESse),], random = ~1|D47ICDESse , method = "REML", na.action = na.omit)
summary(lme_D47ICDES)

# Get pairwise contrasts

modpairwise1 <- emmeans(lme_D47ICDES, pairwise ~ `Mass.Spectrometer`)
modpairwise2 <- emmeans(lme_D47ICDES, pairwise ~ `Mass.Spectrometer`|Standard)

plot(modpairwise1,comparisons = TRUE, adjust = "mvt", 
     horizontal = FALSE, colors = "darkgreen") + theme_bw() +
  xlab("Estimated marginal means") + ylab("")

ggsave("D47ICDES_emmeansplot.tiff", dpi = 800, compression = "lzw", scale = 1.25)

write.csv(modpairwise1, file = "D47ICDES_emmeans.csv", row.names = FALSE)

# Power analysis, Configs 1b and 3, D47 ICDES

# ETH-1
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"| D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-1"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-1"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-1"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 27, length = 907

# ETH-2
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"| D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-2"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-2"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-2"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 26, length = 894

# ETH-3
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"| D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-3"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-3"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-3"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 8, length = 833

# ETH-4
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"| D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-4"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-4"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "ETH-4"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 8, length = 831

# TV03
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"| D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "TV03"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "TV03"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "TV03"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 10, length = 703

# Veinstrom
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"| D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "Veinstrom"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "Veinstrom"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "Veinstrom"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 7, length = 927

# Carmel Chalk
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"| D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "Carmel Chalk"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "Carmel Chalk"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "Carmel Chalk"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 7, length = 905

# Carrara Marble
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"| D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "Carrara Marble"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "Carrara Marble"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 1b"|D47ICDES$`Mass.Spectrometer` == "Config 3" & D47ICDES$Standard == "Carrara Marble"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 13, length = 778

# Config 2, D47 ICDES

# ETH-1
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-1"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-1"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-1"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 3, length = 402

# ETH-2
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-2"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-2"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-2"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 3, length = 386

# ETH-3
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-3"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-3"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-3"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# ETH-4
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-4"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-4"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "ETH-4"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# TV03
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "TV03"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "TV03"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "TV03"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# Veinstrom
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "Veinstrom"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "Veinstrom"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "Veinstrom"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# Carmel Chalk
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "Carmel Chalk"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "Carmel Chalk"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "Carmel Chalk"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# Carrara Marble
pwr.t.test(d = (mean(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "Carrara Marble"])-sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "Carrara Marble"]))/sd(D47ICDES$`Δ47I.CDES....`[D47ICDES$`Mass.Spectrometer` == "Config 2" & D47ICDES$Standard == "Carrara Marble"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 3, length = 44

# Config 1a

# ETH-1
pwr.t.test(d = (mean(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-1"])-sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-1"]))/sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-1"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# n = 3, length = 36

# ETH-2
pwr.t.test(d = (mean(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-2"])-sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-2"]))/sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-2"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# ETH-3
pwr.t.test(d = (mean(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-3"])-sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-3"]))/sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-3"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# ETH-4
pwr.t.test(d = (mean(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-4"])-sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-4"]))/sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "ETH-4"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# TV03
pwr.t.test(d = (mean(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "TV03"])-sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "TV03"]))/sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "TV03"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# Veinstrom
pwr.t.test(d = (mean(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "Veinstrom"])-sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "Veinstrom"]))/sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "Veinstrom"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# Carmel Chalk
pwr.t.test(d = (mean(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "Carmel Chalk"])-sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "Carmel Chalk"]))/sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "Carmel Chalk"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign

# Carrara Marble
pwr.t.test(d = (mean(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "Carrara Marble"])-sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "Carrara Marble"]))/sd(D47CDES90$`Δ47CDES.90....`[D47CDES90$`Mass.Spectrometer` == "Config 1a" & D47CDES90$Standard == "Carrara Marble"]), sig.level = 0.05, power = 0.95, type = "two.sample")
# Error in uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09)) : 
# f() values at end points not of opposite sign