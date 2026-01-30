## libraries
library(here)
library(data.table)
library(ggplot2)
library(viridis)
library(glue)
library(ggthemes)

## utilities
rdb <- "https://raw.githubusercontent.com/petedodd/"
rdd <- "adotb/refs/heads/main/rawdata/" #ado TB repo
there <- function(x) glue("{rdb}{rdd}{x}")

## key to WHO regions
fn <- here("data/whokey.Rdata")
if (!file.exists(fn)) {
  load(url(there("whokey.Rdata"))) # from adotb repo
  save(whokey, file = fn)
} else {
  load(file = fn)
}

## --- mixing data
## contact data
load(url(there("synthetic_contacts_2021.csv.Rdata"))) #from adotb repo

## check
tmp <- synthetic_contacts_2021[
  iso3c == "IND" &
    location_contact == "all" &
    setting == "overall"
]

dcast(tmp, age_cotactee ~ age_contactor,
  value.var = "mean_number_of_contacts"
)


## --- TB estimates in 2024 for 2023
## read in WHO age-specific incidence
fn <- here("data/E.Rdata")
if (!file.exists(fn)) { # get if not there
  E <- fread(here("data/TB_burden_age_sex_2024-10-30.csv")) # from adotb repo
  E[, unique(age_group)]
  exa <- c("all", "0-14", "15plus", "18plus") # exlude age groups
  E <- E[
    !age_group %in% exa & risk_factor == "all",
    .(iso3, sex, age_group, TB = best, TB.sd = (hi - lo) / 3.92)
  ] # choose right age groups
  save(E, file = fn)
} else { # load if there
  load(fn)
}

## --- WPP24 demography for 2023
load(here("data/N80_2023.Rdata")) #already licked into shape
totpop <- N80[, sum(PopTotal)]

## saving out as used in step 2
if (!file.exists(here("data/totpop.Rdata"))) {
  totpop <- N80[, sum(PopTotal)]
  save(totpop, file = here("data/totpop.Rdata"))
}


## age key
akey <- data.table(AgeGrp = unique(N80$AgeGrp))
akey[, age_group := AgeGrp]
akey[age_group %in% c("65-69", "70-74", "75-79", "80+"), age_group := "65plus"]
akey[age_group %in% c("5-9", "10-14"), age_group := "5-14"]
akey[age_group %in% c("15-19", "20-24"), age_group := "15-24"]
akey[age_group %in% c("25-29", "30-34"), age_group := "25-34"]
akey[age_group %in% c("35-39", "40-44"), age_group := "35-44"]
akey[age_group %in% c("45-49", "50-54"), age_group := "45-54"]
akey[age_group %in% c("55-59", "60-64"), age_group := "55-64"]
akey[, acat := gsub("plus", "+", age_group)]
akey[, cage := AgeGrp]
akey[AgeGrp %in% c("75-79", "80+"), cage := "75+"] # for contacts
agz <- akey[, unique(acat)]

if (!file.exists(here("data/agz.Rdata"))) {
  save(agz, file = here("data/agz.Rdata"))
}

## ## summary(CD)

## CD <- synthetic_contacts_2021[
##   location_contact == "all",
##   .(
##     iso3 = iso3c, age_contactor,
##     age_contactee = age_cotactee,
##     ctx = mean_number_of_contacts
##   )
## ]


## ## aggregating contacts
## CD[, AO := gsub(" to ", "-", age_contactor)]
## CD[, AI := gsub(" to ", "-", age_contactee)]
## CD <- merge(CD, unique(akey[, .(cage, acato = acat)]),
##   by.x = "AO", by.y = "cage",
##   all.x = TRUE, all.y = FALSE
## )
## CD <- merge(CD, unique(akey[, .(cage, acati = acat)]),
##   by.x = "AI", by.y = "cage",
##   all.x = TRUE, all.y = FALSE
##   )

## CD <- CD[, .(ctx = sum(ctx)), by = .(iso3, acato, acati)]
## CD <- merge(CD, whokey, by = "iso3") # regions


## ## inspect
## CDR <- CD[, .(contacts = mean(ctx)),
##   by = .(g_whoregion, acato, acati)
## ] # regional ave
## CDR[, acato := factor(acato, levels = agz, ordered = TRUE)]
## CDR[, acati := factor(acati, levels = agz, ordered = TRUE)]


## GP <- ggplot(data = CDR, aes(x = acato, y = acati, fill = contacts)) +
##   geom_tile() +
##   scale_fill_viridis() +
##   theme(legend.position = "bottom") +
##   facet_wrap(~g_whoregion) +
##   xlab("Age of contactor") +
##   ylab("Age of contactee") +
##   ggtitle("Step 2: Regional average contact patterns")
## GP

## ggsave(GP, file = here("plots/ARI_step2_contacts.png"), w = 12, h = 8)

## include sex also
NS <- merge(N80, akey)
NS <- NS[, .(
  pop.total = sum(PopTotal),
  pop.female = sum(PopFemale),
  pop.male = sum(PopMale)
),
by = .(iso3, age_group)
]

## join demography + TB
EW <- dcast(E, iso3 + age_group ~ sex, value.var = c("TB", "TB.sd"))
EW <- merge(EW, NS, by = c("iso3", "age_group"))
EW[, c(
  "pcTB.total", "pcTB.male", "pcTB.female",
  "pcTB.total.sd", "pcTB.male.sd", "pcTB.female.sd"
) := .(
  (TB_m + TB_f) / pop.total,
  TB_m / pop.male, TB_f / pop.female,
  sqrt(TB.sd_m^2 + TB.sd_f^2) / pop.total,
  TB.sd_m / pop.male, TB.sd_f / pop.female
)] # per capita TB
EW <- merge(EW, unique(akey[, .(age_group, acat)]), by = "age_group")
EW$acat <- factor(EW$acat, levels = agz, ordered = TRUE)
refs <- EW[acat == "15-24", .(iso3, refpcTB = pcTB.female)] # ref cat
EW <- merge(EW, refs, by = c("iso3"))

## relative per capita TB
EW[, c("relpcTB.male", "relpcTB.female") :=
  .(pcTB.male / refpcTB, pcTB.female / refpcTB)]
EW[, c("relpcTB.male.sd", "relpcTB.female.sd") :=
  .(pcTB.male.sd / refpcTB, pcTB.female.sd / refpcTB)]
EW[!is.finite(relpcTB.male), relpcTB.male := NA]
EW[!is.finite(relpcTB.female), relpcTB.female := NA]

## look:
EW[iso3 == "VNM"]

EW <- merge(EW, whokey, by = "iso3") # merge region

## make regional version
ER <- EW[, .(
  relpcTB.male = mean(relpcTB.male, na.rm = TRUE),
  relpcTB.female = mean(relpcTB.female, na.rm = TRUE),
  relpcTB.male.sdpc = mean(relpcTB.male.sd / relpcTB.male, na.rm = TRUE),
  relpcTB.female.sdpc = mean(relpcTB.female.sd / relpcTB.female, na.rm = TRUE)
),
by = .(g_whoregion, acat)
]
ER[, iso3 := g_whoregion]

## reshape and plot this data
ERM <- melt(ER, id = c("g_whoregion", "iso3", "acat"))
ERM[, type := ifelse(grepl("sdpc", variable), "sdpc", "mid")]
ERM[, variable := gsub("\\.sdpc", "", variable)]
ERM[, c("variable", "sex") := tstrsplit(variable, split = "\\.")]
ERM <- dcast(ERM,
  g_whoregion + iso3 + acat + variable + sex ~ type,
  value.var = "value"
)

GP <- ggplot(
  ERM[!acat %in% c("0-4", "5-14")],
  aes(acat, mid, group = paste(iso3, sex), col = sex)
) +
  geom_line() +
  facet_wrap(~g_whoregion) +
  geom_hline(yintercept = 1, col = 2, lty = 2) +
  scale_y_sqrt() +
  ylab("Square-root of relative per capita TB incidence") +
  xlab("Age category") +
  ggtitle("Step 1: Relative per capita TB incidence (WHO estimates)") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
GP

ggsave(GP, file = here("output/step1_ARI_percapTB.png"), w = 12, h = 7)


## --- contacts
## restrict & rename
CD <- synthetic_contacts_2021[
  location_contact == "all" & setting == "overall",
  .(
    iso3 = iso3c, age_contactor,
    age_contactee = age_cotactee,
    ctx = mean_number_of_contacts
  )
]

## -- aggregating contacts
CD[, AO := gsub(" to ", "-", age_contactor)]
CD[, AI := gsub(" to ", "-", age_contactee)]
CD <- merge(CD, unique(akey[, .(cage, acato = acat)]),
  by.x = "AO", by.y = "cage",
  all.x = TRUE, all.y = FALSE
)
CD <- merge(CD, unique(akey[, .(cage, acati = acat)]),
  by.x = "AI", by.y = "cage",
  all.x = TRUE, all.y = FALSE
)

## aggregate 'to' contacts (over acati):
CD <- CD[, .(ctx = sum(ctx)), by = .(iso3, acati, AO)]

## compare age groups
setdiff(N80[, unique(AgeGrp)], CD[, unique(AO)])
N75 <- N80[, AO := ifelse(
  AgeGrp %in% c("75-79", "80+"),
  "75+",
  as.character(AgeGrp)
)]
unique(N75[, .(AgeGrp, AO)]) # check
## aggregate last cat:
N75 <- N75[, .(PopTotal = sum(PopTotal)), by = .(iso3, AO)]

## merge pops
setdiff(N75[, unique(AO)], CD[, unique(AO)])
setdiff(CD[, unique(AO)], N75[, unique(AO)])
CD <- merge(CD, N75, by = c("iso3", "AO"))

dcast(CD[iso3 == "IND"],
      acati ~ AO,
  value.var = "ctx"
  )                                        #cf same as sheet
unique(CD[iso3 == "IND", .(AO, PopTotal)]) #cf pops in validation sheet

## acato back in
CD <- merge(CD,
  unique(akey[, .(acato = acat, AO = cage)]),
  by = "AO"
)

## weighted mean over contactor
CD <- CD[, .(
  ctx = weighted.mean(ctx, PopTotal)
),
by = .(iso3, acato, acati)
]

## order factors
CD[, acato := factor(acato, levels = agz, ordered = TRUE)]
CD[, acati := factor(acati, levels = agz, ordered = TRUE)]

dcast(CD[iso3 == "IND"],
  acati ~ acato,
  value.var = "ctx"
) # cf same as sheet


## merging in aggregate population data
NS[, acato := age_group]
NS[acato == "65plus", acato := "65+"]
CD <- merge(CD, NS[, .(iso3, acato, pop.total)], by = c("iso3", "acato"))
CD <- merge(CD, whokey, by = "iso3") # regions

summary(CD)

## === inspect

## individual countries
CD[iso3 == "IND"][, sum(ctx), by = acato] #cf very close to validation sheet

## patterns by age
CDR <- CD[, .(contacts = mean(ctx)),
  by = .(g_whoregion, acato, acati)
] # regional average
CDR[, acato := factor(acato, levels = agz, ordered = TRUE)]
CDR[, acati := factor(acati, levels = agz, ordered = TRUE)]

## plot WAIFW matrices
GP <- ggplot(data = CDR, aes(x = acato, y = acati, fill = contacts)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.position = "bottom") +
  facet_wrap(~g_whoregion) +
  xlab("Age of contactor") +
  ylab("Age of contactee") +
  ggtitle("Step 2: Regional average WAIFW matrices")
GP

ggsave(GP, file = here("output/step2_ARI_WAIFW.png"), w = 12, h = 8)


## create and plot regional averages
cdr <- CD[, .(
  contacts = weighted.mean(ctx, w = pop.total),
  pop.total = sum(pop.total)
),
by = .(g_whoregion, acato, acati)
] # regional ave
cdr[, acato := factor(acato, levels = agz, ordered = TRUE)]
cdr[, acati := factor(acati, levels = agz, ordered = TRUE)]
cdr[, totar := sum(contacts),
  by = .(acato, g_whoregion)
] # normalize by contactor contacts
cdrr <- cdr[, .(contacts = sum(contacts * pop.total)),
  by = .(g_whoregion, acati)
]
cdrr[, totar := sum(contacts), by = .(g_whoregion)] # normalize
cdrr[, acato := NA]

GP <- ggplot(data = cdr, aes(
  x = acati, y = contacts / totar,
  col = acato, group = paste(g_whoregion, acato)
)) +
  geom_line() +
  geom_line(data = cdrr, col = 2, lty = 2, lwd = 2) +
  theme_linedraw() +
  scale_color_colorblind(name = "Age contactor") +
  theme(legend.position = "bottom") +
  facet_wrap(~g_whoregion) +
  xlab("Age of contactee") +
  ylab("Proportion of contacts made to each
contactee age group by contactor group") +
  ggtitle("Step 2: Regional average contact patterns") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
GP

ggsave(GP, file = here("output/step2_ARI_line_contacts.png"), w = 12, h = 8)

fwrite(cdrr[, .(prop = contacts / totar, acati, g_whoregion)],
  file = here("output/step2_contactprops.csv")
)


## including pops
cdr0 <- copy(cdr)
cdr0[acati %in% c("0-4", "5-14"), contacts := 0]
cdr0[, acato := factor(acato, levels = agz, ordered = TRUE)]
cdr0[, acati := factor(acati, levels = agz, ordered = TRUE)]
cdr0[, totar := sum(contacts),
  by = .(acato, g_whoregion)
] # normalize by contactor contacts
cdr0r <- cdr0[, .(contacts = sum(contacts * pop.total)),
  by = .(g_whoregion, acati)
]
cdr0r[, totar := sum(contacts), by = .(g_whoregion)] # normalize
cdr0r[, acato := NA]

GP <- ggplot(
  data = cdr0,
  aes(
    x = acati, y = contacts / totar,
    col = acato, group = paste(g_whoregion, acato)
  )
) +
  geom_line() +
  geom_line(data = cdr0r, col = 2, lty = 2, lwd = 2) +
  theme_linedraw() +
  scale_color_colorblind(name = "Age contactor") +
  theme(legend.position = "bottom") +
  facet_wrap(~g_whoregion) +
  xlab("Age of contactee") +
  ylab("Proportion of contacts made to each contactee
age group by contactor group") +
  ggtitle("Step 2: Regional average 'effective' contact patterns:
zeroing child contactees") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
GP

ggsave(GP,
  file = here("output/step2_ARI_line_contacts_inf.png"),
  w = 12, h = 8
) #


fwrite(cdr0r[, .(prop = contacts / totar, acati, g_whoregion)],
  file = here("output/step2_contactprops0.csv")
)


EW[acat == "0-4"]

## merge in contact data
EC <- merge(
  EW[, .(
    iso3, g_whoregion, acat,
    pop.female, pop.male,
    pcTB.female, pcTB.male,
    pcTB.female.sd, pcTB.male.sd
  )],
  CD,
  by.x = c("iso3", "g_whoregion", "acat"),
  by.y = c("iso3", "g_whoregion", "acato"), # contactor
  all.x = TRUE, all.y = FALSE
)
EC <- EC[!is.na(ctx)]
EC[, sum(pcTB.female)]
EC[, sum(pcTB.male)]

length(unique(EC$iso3)) # 177 countries


## === assortativity
## LN(-0.57,0.085)
ass <- 0.6 ## exp(-0.570 + 0.5*0.085^2)

## how does this data translate to: ass[s,A,A']?
## assuming:
## C[s,s',a,a'] = [ ass[s,A,A'](s==s') + (1-ass[s,A,A'])(s!=s') ] K[a,a']
## and:
## N[f,a] = N[m,a], aggregating to calculate fraction of contacts a\in A etc:
## F[s,s',A,A'] = ass[s,A,A'](s==s') + (1-ass[s,A,A'])(s!=s')
## in other words: aggregates consistently if N[f,a] = N[m,a]
ASS <- SASS <- array(0.5,
                     dim = c(2, 2, 2),
                     dimnames = list(fr = c("K", "A"),
                                     to = c("K", "A"),
                                     sex = c("M", "F")))
## data from Katherine:
## Boys w/ children: 	0.61 [0.59; 0.63]	0.61 [0.59; 0.63]
ASS["K", "K", "M"] <- 0.61
SASS["K", "K", "M"] <- abs(0.59 - 0.63) / 3.92
## Boys w/ adults: 	0.42 [0.40; 0.44]	0.41 [0.39; 0.43]
ASS["K", "A", "M"] <- 0.41
SASS["K", "A", "M"] <- abs(0.39 - 0.43) / 3.92
## Girls w/ children: 	0.60 [0.57; 0.62]	0.60 [0.57; 0.63]
ASS["K", "K", "F"] <- 0.60
SASS["K", "K", "F"] <- abs(0.57 - 0.63) / 3.92
## Girls w/ adults: 	0.62 [0.60; 0.64]	0.63 [0.61; 0.66]
ASS["K", "A", "F"] <- 0.63
SASS["K", "A", "F"] <- abs(0.61 - 0.66) / 3.92
## Men w/ children: 	0.53 [0.51; 0.55]	0.52 [0.51; 0.54]
ASS["A", "K", "M"] <- 0.52
SASS["A", "K", "M"] <- abs(0.5 - 0.54) / 3.92
## Men w/ adults: 	0.56 [0.54; 0.58]	0.54 [0.51; 0.57]
ASS["A", "A", "M"] <- 0.54
SASS["A", "A", "M"] <- abs(0.51 - 0.57) / 3.92
## Women w/ children: 0.52 [0.50; 0.53)	0.52 [0.50; 0.54]
ASS["A", "K", "F"] <- 0.52
SASS["A", "K", "F"] <- abs(0.50 - 0.54) / 3.92
## Women w/ adults: 	0.60 [0.58; 0.62]	0.61 [0.58; 0.64]
ASS["A", "A", "F"] <- 0.61
SASS["A", "A", "F"] <- abs(0.58 - 0.64) / 3.92

## check
ASS

## merge in ASS
EC[, fr := ifelse(acat %in% c("0-4", "5-14"), "K", "A")]
EC[, to := ifelse(acati %in% c("0-4", "5-14"), "K", "A")]
DASS <- as.data.table(ASS)
DASW <- dcast(DASS, fr + to ~ sex, value.var = "value")
SASW <- dcast(as.data.table(SASS), fr + to ~ sex, value.var = "value")
SASW <- SASW[, .(fr, to, F.sd = F, M.sd = M)]
EC <- merge(EC, DASW, by = c("fr", "to"))
EC <- merge(EC, SASW, by = c("fr", "to"))


## under 15 not infectious
EC[acat %in% c("0-4", "5-14"),
   c("pcTB.female", "pcTB.male", "pcTB.female.sd", "pcTB.male.sd") := 0]

## look:
EC[iso3 == "ZAF", .(acat, acati, ctx)]
CD[iso3 == "ZAF"]


## country weightings (TB)
ET <- E[, .(TB = sum(TB)), by = iso3]
EC <- merge(EC, ET, by = "iso3")


##   NOTE  to  <- from
EC[, c(
  "ari.male.male",
  "ari.female.male",
  "ari.male.female", # NOTE the below are read right to left
  "ari.female.female"
) := .(
  ctx * pcTB.male * M,
  ctx * pcTB.male * (1 - F), # NB assortativity of contactee
  ctx * pcTB.female * (1 - M),
  ctx * pcTB.female * F
)]

## counterfactual
EC[, c(
  "ari.male.male0",
  "ari.female.male0",
  "ari.male.female0", # NOTE the below are read right to left
  "ari.female.female0"
) := .(
  ctx * pcTB.female * M,
  ctx * pcTB.female * (1 - F), # NB assortativity of contactee
  ctx * pcTB.female * (1 - M),
  ctx * pcTB.female * F
)]


## quick check:
## ass=0.6 for all means 1.1 as expected
EC[, sum(ari.male.male + ari.male.female)] /
  EC[, sum(ari.female.male + ari.female.female)] # to M/F
EC[, sum(ari.male.male + ari.female.male)] /
  EC[, sum(ari.female.female + ari.male.female)] # from M/F


## checks
TAB <- array(0.5,
  dim = c(2, 2),
  dimnames = list(fr = c("F", "M"), to = c("F", "M"))
)
TAB["F", "M"] <- EC[, sum(ari.male.female)]
TAB["M", "F"] <- EC[, sum(ari.female.male)]
TAB["F", "F"] <- EC[, sum(ari.female.female)]
TAB["M", "M"] <- EC[, sum(ari.male.male)]
TAB <- TAB / sum(TAB)

colSums(TAB)                      # TO
rowSums(TAB)                      # FROM
sum(TAB["M", ]) / sum(TAB["F", ]) # from M/F
sum(TAB[, "M"]) / sum(TAB[, "F"]) # to M/F

write.csv(TAB, file = "output/TAB.csv")


kds <- c("0-4", "5-14")
E[sex == "m", sum(TB)] / E[sex == "f", sum(TB)] # 1.6
E[sex == "m" & !age_group %in% kds, sum(TB)] /
  E[sex == "f" & !age_group %in% kds, sum(TB)] # 1.7

## rough corrobrating calc
foitoM <- 1.7 * 0.6 + 1 * 0.4
foitoF <- 1.7 * 0.4 + 1 * 0.6
foitoM / foitoF # 1.1
foifrM <- 1.7 * 0.6 + 1.7 * 0.4
foifrF <- 1 * 0.4 + 1 * 0.6
foifrM / foifrF # 1.7


## d(A*B)/(A*B) = sqrt((dA/A)^2+(dB/B)^2)
EC[, c(
  "ari.male.male.sd",
  "ari.female.male.sd",
  "ari.male.female.sd", # NOTE the below are read right to left
  "ari.female.female.sd"
) := .(
       pcTB.male * M * ctx *
       sqrt((M.sd / M)^2 + (pcTB.male.sd / pcTB.male)^2),
  pcTB.male * (1 - F) * ctx *
  sqrt((F.sd / (1 - F))^2 + (pcTB.male.sd / pcTB.male)^2),
  pcTB.female * (1 - M) * ctx *
  sqrt((M.sd / (1 - M))^2 + (pcTB.female.sd / pcTB.female)^2),
  pcTB.female * F * ctx *
  sqrt((F.sd / F)^2 + (pcTB.female.sd / pcTB.female)^2)
  )]
## set child SDs
EC[!is.finite(ari.female.female.sd), c(
  "ari.male.male.sd",
  "ari.female.male.sd",
  "ari.male.female.sd",
  "ari.female.female.sd"
) := 0]


## counterfactual
EC[, c(
  "ari.male.male0.sd",
  "ari.female.male0.sd",
  "ari.male.female0.sd", # NOTE the below are read right to left
  "ari.female.female0.sd"
) := .(
  pcTB.female * M * ctx *
    sqrt((M.sd / M)^2 + (pcTB.male.sd / pcTB.male)^2),
  pcTB.female * (1 - F) * ctx *
    sqrt((F.sd / (1 - F))^2 + (pcTB.male.sd / pcTB.male)^2),
  pcTB.female * (1 - M) * ctx *
    sqrt((M.sd / (1 - M))^2 + (pcTB.female.sd / pcTB.female)^2),
  pcTB.female * F * ctx *
    sqrt((F.sd / F)^2 + (pcTB.female.sd / pcTB.female)^2)
)]

## set child SDs
EC[!is.finite(ari.female.female0.sd), c(
  "ari.male.male0.sd",
  "ari.female.male0.sd",
  "ari.male.female0.sd",
  "ari.female.female0.sd"
) := 0]

## checks
EC[iso3 == "AFG" & acat == "0-4"] # acat is from
EC[iso3 == "AFG" & acati == "0-4"] # acat is to
EC[, sum(pcTB.male)] / EC[, sum(pcTB.female)]


save(EC, file = here("data/EC.Rdata"))
