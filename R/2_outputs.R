## libraries
library(here)
library(data.table)
library(ggplot2)
library(glue)
library(ggthemes)
library(scales)
library(ggpubr)

load(file = here("data/EC.Rdata"))


## restructure
ECM <- melt(
  EC[, .(
    iso3, g_whoregion, acat, acati,
    ari.male.male, ari.female.male,
    ari.male.female, ari.female.female,
    ari.male.male.sd, ari.female.male.sd,
    ari.male.female.sd, ari.female.female.sd
  )],
  id = c("iso3", "g_whoregion", "acat", "acati")
)
ECM[, type := ifelse(grepl("sd", variable), "sd", "mn")]
ECM[, variable := gsub("\\.sd", "", variable)]
ECM <- dcast(ECM, iso3 + g_whoregion + acat + acati + variable ~ type, value.var = "value")
ECM[, c("qty", "to", "from") := tstrsplit(variable, split = "\\.")]
ECM$acat <- factor(ECM$acat, levels = agz, ordered = TRUE)
ECM$acati <- factor(ECM$acati, levels = agz, ordered = TRUE)



## merge in population for infectees
pop <- EW[, .(iso3, acati = age_group, pop.female, pop.male)]
pop[acati == "65plus", acati := "65+"]
pop <- melt(pop, id = c("iso3", "acati"))
pop[, c("vr", "sex") := tstrsplit(variable, split = "\\.")]
ECM <- merge(ECM, pop[, .(iso3, to = sex, acati, popto = value)], by = c("iso3", "to", "acati"))
ECM[!is.finite(sd), sd := 0.1] # NOTE small nation noisy pattern safety


## NOTE weighted by population - assumes same beta everywhere
## infectees
## most aggregate version
ECG <- ECM[, .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
),
by = .(acati, from, to)
]
ECG[, c("tot", "tot.sd") := .(sum(ari), ssum(ari.sd))]
ECG[, c("fari", "fari.sd") := .(ari / tot, (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2))]
ECG$acati <- factor(ECG$acati, levels = agz, ordered = TRUE)
TXT <- ECG[, .(tot = sum(ari)), by = .(acati, to)]
TXT <- merge(ECG[, .(acati, from, to, ari, fari, fari.sd)], TXT, by = c("acati", "to")) #
TXT <- TXT[from == "male"]
TXT[acati == "65+" & to == "male", fari := fari - 1e-3]
TXT[, pcnt := paste0(format(round(100 * ari / tot)), "%")]


## check
ECG[to == "male", sum(ari)] / ECG[to == "female", sum(ari)]
ECG[from == "male", sum(ari)] / ECG[from == "female", sum(ari)]
CD[, mean(ctx), by = acato]
CD[, mean(ctx), by = acati]


## CI data
ECGT <- ECG[, .(fari = sum(fari), fari.sd = ssum(fari.sd)), by = .(acati, to)]
ECGT[, from := "male"]

cl <- colorblind_pal()(2)[2] # relevant colot


## NOTE most favoured
GPa <- ggplot(ECG, aes(acati, fari, fill = from)) +
  geom_bar(stat = "identity") +
  scale_fill_colorblind() +
  geom_errorbar(
    data = ECGT,
    aes(ymin = fari - 1.96 * fari.sd, ymax = fari + 1.96 * fari.sd), width = 0, col = 2
  ) +
  facet_wrap(~to) +
  theme_linedraw() +
  scale_y_continuous(label = percent) +
  xlab("Age of infectee") +
  ylab("Proportion of all exposure to each group") +
  geom_text(data = TXT, aes(acati, fari + 3e-3, label = pcnt), col = cl) +
  guides(color = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
GPa

ggsave(file = here("output/ARI_to.png"), w = 12, h = 6)


ggplot(ECG, aes(acati, fari, fill = from)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_colorblind() +
  facet_wrap(~to) +
  theme_linedraw() +
  scale_y_continuous(label = percent) +
  xlab("Age") +
  ylab("Proportion of all exposure") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here("output/ARI_to_stretch.png"), w = 10, h = 5)



## ==== conditional by sex
ECG[, c("tots", "tots.sd") := .(sum(ari), ssum(ari.sd)), by = to]
ECG[, c("faris", "faris.sd") := .(ari / tots, (ari / tots) * sqrt((ari.sd / ari)^2 + (tots.sd / tots)^2))]
TXTs <- ECG[, .(tots = sum(ari)), by = .(acati, to)]
TXTs <- merge(ECG[, .(acati, from, to, ari, faris, faris.sd)], TXTs, by = c("acati", "to"))
TXTs <- TXTs[from == "male"]
TXTs[, pcnt := paste0(format(round(100 * ari / tots)), "%")]

## CI data
ECGTs <- ECG[, .(faris = sum(faris), faris.sd = ssum(faris.sd)), by = .(acati, to)]
ECGTs[, from := "male"]



## NOTE most favoured
ggplot(ECG, aes(acati, faris, fill = from)) +
  geom_bar(stat = "identity") +
  scale_fill_colorblind() +
  geom_errorbar(
    data = ECGTs,
    aes(ymin = faris - 1.96 * faris.sd, ymax = faris + 1.96 * faris.sd), width = 0, col = 2
  ) +
  facet_wrap(~to) +
  theme_linedraw() +
  scale_y_continuous(label = percent) +
  xlab("Age") +
  ylab("Proportion of all exposure") +
  geom_text(data = TXTs, aes(acati, faris + 3e-3, label = pcnt), col = cl) +
  guides(color = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here("output/ARI_to_conditional.png"), w = 10, h = 5)



## ---- counterfactual

## restructure
ECM0 <- melt(
  EC[, .(
    iso3, g_whoregion, acat, acati,
    ari.male.male0, ari.female.male0,
    ari.male.female0, ari.female.female0,
    ari.male.male0.sd, ari.female.male0.sd,
    ari.male.female0.sd, ari.female.female0.sd
  )],
  id = c("iso3", "g_whoregion", "acat", "acati")
)

ECM0[, type := ifelse(grepl("sd", variable), "sd", "mn")]
ECM0[, variable := gsub("\\.sd", "", variable)]
ECM0[, variable := gsub("0", "", variable)]
ECM0 <- dcast(ECM0, iso3 + g_whoregion + acat + acati + variable ~ type, value.var = "value")
ECM0[, c("qty", "to", "from") := tstrsplit(variable, split = "\\.")]
ECM0$acat <- factor(ECM0$acat, levels = agz, ordered = TRUE)
ECM0$acati <- factor(ECM0$acati, levels = agz, ordered = TRUE)
## merge in population for infectees
ECM0 <- merge(ECM0, pop[, .(iso3, to = sex, acati, popto = value)], by = c("iso3", "to", "acati"))
ECM0[!is.finite(sd), sd := 0.1] # NOTE small nation noisy pattern safety


## NOTE weighted by population - assumes same beta everywhere
## infectees
## most aggregate
ECG1 <- ECM[, .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
)]
ECG0 <- ECM0[, .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
)]

## ---------------- stats
ECGL <- ECM[, .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
),
by = .(acati, acat, from, to)
]

ECGL[, c("tot", "tot.sd") := .(sum(ari), ssum(ari.sd))]
ECGL[, c("fari", "fari.sd") := .(ari / tot, (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2))]
ECGL[!is.finite(fari.sd), fari.sd := 0] # kids manually 0

## utility function to format X / Y as percent
XoYpc <- function(num, den) {
  mid <- num$fari / den$fari
  mid.sd <- mid * sqrt((num$fari.sd / num$fari)^2 + (den$fari.sd / den$fari)^2)
  hi <- min(1, mid + 1.96 * mid.sd)
  lo <- max(0, mid - 1.96 * mid.sd)
  paste0(format(round(100 * mid)), "% (", format(round(100 * lo)), "% to ", format(round(100 * hi)), "%)")
}


## NOTE fraction of everyone that are male adults
NS[!age_group %in% c("0-4", "5-14"), sum(pop.male)] / NS[, sum(pop.total)] # 37 % of population


## out text
outtxt <- list()

## to men
den <- ECGL[to == "male" & !acati %in% c("0-4", "5-14"), .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  from == "male" & !acat %in% c("0-4", "5-14") & # from men
    to == "male" & !acati %in% c("0-4", "5-14"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[1]] <- data.table(qty = "% infections to men, from men", value = XoYpc(num, den))

## to women
den <- ECGL[to == "female" & !acati %in% c("0-4", "5-14"), .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  from == "male" & !acat %in% c("0-4", "5-14") & # from men
    to == "female" & !acati %in% c("0-4", "5-14"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[2]] <- data.table(qty = "% infections to women, from men", value = XoYpc(num, den))

## to kids
den <- ECGL[acati %in% c("0-4", "5-14"), .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  from == "male" & !acat %in% c("0-4", "5-14") & # from men
    acati %in% c("0-4", "5-14"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[3]] <- data.table(qty = "% infections to children, from men", value = XoYpc(num, den))

## to anyone
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  from == "male" & !acat %in% c("0-4", "5-14"), # from men
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[4]] <- data.table(qty = "% infections, from men", value = XoYpc(num, den))

## to anyone, from 25-44
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  from == "male" & acat %in% c("25-34", "35-44"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[5]] <- data.table(qty = "% infections, from M25-44", value = XoYpc(num, den))

## to men
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  to == "male",
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[6]] <- data.table(qty = "% infections to men ", value = XoYpc(num, den))

## population data: fraction male adults
val <- NS[!age_group %in% c("0-4", "5-14"), sum(pop.male)] / NS[, sum(pop.total)]
outtxt[[7]] <- data.table(
  qty = "% of population men >=15",
  value = paste0(format(round(100 * val)))
)

## population data: fraction male adults 15-44
val <- NS[age_group %in% c("25-34", "35-44"), sum(pop.male)] / NS[, sum(pop.total)]
outtxt[[8]] <- data.table(
  qty = "% of population M25-44",
  value = paste0(format(round(100 * val)))
)

## counterfactual reduction
den <- ECG1[, .(fari = ari, fari.sd = ari.sd)]
num <- ECG0[, .(fari = ari, fari.sd = ari.sd)]
outtxt[[9]] <- data.table(qty = "% exposure if F were M ", value = XoYpc(num, den))

## join
outtxt <- rbindlist(outtxt)
outtxt

fwrite(outtxt, file = here("output/outtxt.csv"))


## • regions (potentially in text)
## • children different?
## • patterns across countries: incidence?
## • some sort of amplifier estimate: eg (proportional change in transmission) / (proportional reduction in male prevalence, eg to F)?
pcm <- 0.64
1 - 2 * (1 - pcm) # 28% lower


## going by acato ------infector
ECM[acat == "65+"]
ECM[, mean(mn), by = acat] # from mean per cap contrib 20% higher in 65+
ECM[, sum(popto), by = acat] # const
ECM[, sum(popto), by = acati] # const
ECM[, unique(acati)]


tmp <- ECM[iso3 == "ZAF" & to == "female" & acati == "25-34"]
tmp <- ECM[iso3 == "VNM" & to == "female" & acati == "25-34"]
tmp <- ECM[iso3 == "VNM" & to == "female", .(mn = sum(mn * popto)), by = .(acat, from)]
tmp[, tot := sum(mn)]
ggplot(tmp, aes(acat, mn / tot, col = from, group = from)) +
  geom_line()


ECGi <- ECM[, .(ari = sum(mn * popto), ari.sd = ssum(sd * popto)),
  by = .(acat, from, to)
] # total exposure happening, by age of source
ECGi[, c("tots", "tots.sd") := .(sum(ari), ssum(ari.sd))]
ECGi[, c("faris", "faris.sd") := .(ari / tots, (ari / tots) * sqrt((ari.sd / ari)^2 + (tots.sd / tots)^2))]
ECGi$acat <- factor(ECGi$acat, levels = agz, ordered = TRUE)


## CI data
ECGTsi <- ECGi[, .(faris = sum(faris), faris.sd = ssum(faris.sd)), by = .(acat, to)]
ECGTsi[, sum(faris)]
ECGTsi[, from := "male"]


## NOTE most favoured 2
GPb <- ggplot(ECGi, aes(acat, faris, fill = from)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    data = ECGTsi,
    aes(ymin = faris - 1.96 * faris.sd, ymax = faris + 1.96 * faris.sd), width = 0, col = 2
  ) +
  scale_fill_colorblind() +
  facet_wrap(~to) +
  theme_linedraw() +
  scale_y_continuous(label = percent) +
  xlab("Age of infector") +
  ylab("Proportion of all exposure from each group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
GPb

ggsave(file = here("output/ARI_from.png"), w = 10, h = 5)


## ----------------
## combined plot:

ggarrange(GPa, GPb, ncol = 1, labels = c("A", "B"), common.legend = TRUE)
ggsave(file = here("output/ARI_BOTH.png"), w = 10, h = 10)


## for this plot:
## columns are representing sexes as exposure targets
## colors are representing sexes as sources
## ages are targets for A and sources for B
## A and B sum to 1


## -----------------
ECGir <- ECM[, .(ari = sum(mn * popto)),
  by = .(acat, from, to, g_whoregion)
] # total exposure happening, by age of source
ECGir[, ari := ari / sum(ari), by = g_whoregion]
ECGir$acat <- factor(ECGir$acat, levels = agz, ordered = TRUE)



ggplot(ECGir, aes(acat, ari, fill = from)) +
  geom_bar(stat = "identity") +
  scale_fill_colorblind() +
  facet_grid(to ~ g_whoregion) +
  theme_linedraw() +
  scale_y_continuous(label = percent) +
  xlab("Age of infector") +
  ylab("Proportion of all infection") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here("output/ARI_from_reg.png"), w = 20, h = 10)




## exposure weighting by region
ECGgr <- ECM[, .(ari = sum(mn * popto)),
  by = .(g_whoregion)
] # total exposure happening, by age of source
ECGgr[, ari := ari / sum(ari)]


ggplot(ECGgr, aes(g_whoregion, ari)) +
  geom_bar(stat = "identity") +
  theme_linedraw() +
  scale_y_continuous(label = percent) +
  xlab("WHO region") +
  ylab("Proportion of all exposure") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## TODO what is this?

## ggsave(file=here('../plots/ARI_inreg.png'),w=6,h=5)

CDS <- CD[!acati %in% kds, .(ctx = sum(ctx)), by = .(acato, g_whoregion)]
CDS[, tot := sum(ctx), by = g_whoregion]
CDS$acato <- factor(CDS$acato, levels = agz, ordered = TRUE)

ggplot(CDS, aes(acato, ctx / tot)) +
  geom_bar(stat = "identity") +
  theme_linedraw() +
  facet_wrap(~g_whoregion) +
  scale_y_continuous(label = percent) +
  xlab("Contactor age") +
  ylab("Proportion of all contacts to adults from each age") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(file = here("output/Contacts_reg.png"), w = 10, h = 7) # TODO BUG check


## pc TB by region?
EG <- merge(E, whokey, by = "iso3")
EG <- EG[, .(TB = sum(TB)), by = .(sex = toupper(sex), acat = age_group, g_whoregion)]
EG[, sex := ifelse(sex == "F", "female", "male")]
EG[acat == "65plus", acat := "65+"]
popr <- merge(pop, whokey, by = "iso3")
popr <- popr[, .(n = sum(value)), by = .(g_whoregion, sex, acat = acati)]

EG <- merge(EG, popr, by = c("g_whoregion", "sex", "acat"))
EG[, pcTB := 1e2 * TB / n] # pop in 1000s
EG$acat <- factor(EG$acat, levels = agz, ordered = TRUE)



ggplot(EG, aes(acat, pcTB)) +
  geom_bar(stat = "identity") +
  theme_linedraw() +
  facet_wrap(~g_whoregion, scale = "free") +
  scale_y_continuous() +
  xlab("Age") +
  ylab("TB incidence per 100,000") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here("output/pcTB_reg.png"), w = 10, h = 7)

