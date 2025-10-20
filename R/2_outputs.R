## libraries
library(here)
library(data.table)
library(ggplot2)
library(ggthemes)
library(scales)
library(ggpubr)

## utility
ssum <- function(x) sqrt(sum(x^2))
cl <- colorblind_pal()(2)[2] # relevant color (gold)
clz <- colorblind_pal()(8) #more colors

## data from file 1
load(file = here("data/totpop.Rdata"))
load(file = here("data/EC.Rdata"))
load(file = here("data/E.Rdata"))
load(file = here("data/agz.Rdata"))
load(file = here("data/whokey.Rdata"))
NS <- unique(EC[, .(iso3, acat, pop.female, pop.male)]) # convenient store

whoz <- c("AFR", "AMR", "EMR", "EUR", "SEA", "WPR")
whozt <- c(
  "Africa", "The Americas",
  "Eastern Mediterranean", "Europe", "South-East Asia",
  "Western Pacific"
)
regkey <- data.table(g_whoregion = whoz, region = whozt)


## NOTE
## output/Contacts_reg.png currently commented out:
## requires CD data saving from 1

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
ECM <- dcast(ECM,
  iso3 + g_whoregion + acat + acati + variable ~ type,
  value.var = "value"
)
ECM[, c("qty", "to", "from") := tstrsplit(variable, split = "\\.")]
ECM$acat <- factor(ECM$acat, levels = agz, ordered = TRUE)
ECM$acati <- factor(ECM$acati, levels = agz, ordered = TRUE)



## merge in population for infectees
pop <- NS[, .(iso3, acati = acat, pop.female, pop.male)]
pop <- melt(pop, id = c("iso3", "acati"))
pop[, c("vr", "sex") := tstrsplit(variable, split = "\\.")]
ECM <- merge(ECM,
  pop[, .(iso3, to = sex, acati, popto = value)],
  by = c("iso3", "to", "acati")
)
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
ECG[, c("fari", "fari.sd") :=
        .(ari / tot, (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2))]
ECG$acati <- factor(ECG$acati, levels = agz, ordered = TRUE)
TXT <- ECG[, .(tot = sum(ari)), by = .(acati, to)]
TXT <- merge(ECG[, .(acati, from, to, ari, fari, fari.sd)],
  TXT,
  by = c("acati", "to")
)
TXT <- TXT[from == "male"]
TXT[acati == "65+" & to == "male", fari := fari - 1e-3]
TXT[, pcnt := paste0(format(round(100 * ari / tot)), "%")]


## check
ECG[to == "male", sum(ari)] / ECG[to == "female", sum(ari)]
ECG[from == "male", sum(ari)] / ECG[from == "female", sum(ari)]
## CD[, mean(ctx), by = acato]
## CD[, mean(ctx), by = acati]


## CI data
ECGT <- ECG[, .(fari = sum(fari), fari.sd = ssum(fari.sd)), by = .(acati, to)]
ECGT[, from := "male"]


## NOTE most favoured
GPa <- ggplot(ECG, aes(acati, fari, fill = from)) +
  geom_bar(stat = "identity") +
  scale_fill_colorblind() +
  geom_errorbar(
    data = ECGT,
    aes(ymin = fari - 1.96 * fari.sd, ymax = fari + 1.96 * fari.sd),
    width = 0, col = 2
  ) +
  ggh4x::facet_wrap2(~to, axes = "all") +
  theme_classic() +
  ggpubr::grids() +
  scale_y_continuous(label = percent) +
  xlab("Age group of infectee (years)") +
  ylab("Proportion of all exposure to each group") +
  geom_text(data = TXT, aes(acati, fari + 3e-3, label = pcnt), col = cl) +
  guides(color = "none") +
  ggtitle("exposure") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(2, "lines"), # or 3
    strip.text = element_text(face = "italic", size = 12),
    strip.background = element_blank(),
    strip.placement = "outside"
  )
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
ECG[, c("faris", "faris.sd") :=
  .(ari / tots, (ari / tots) * sqrt((ari.sd / ari)^2 + (tots.sd / tots)^2))]
TXTs <- ECG[, .(tots = sum(ari)), by = .(acati, to)]
TXTs <- merge(ECG[, .(acati, from, to, ari, faris, faris.sd)],
  TXTs,
  by = c("acati", "to")
)
TXTs <- TXTs[from == "male"]
TXTs[, pcnt := paste0(format(round(100 * ari / tots)), "%")]

## CI data
ECGTs <- ECG[, .(faris = sum(faris), faris.sd = ssum(faris.sd)),
  by = .(acati, to)
]
ECGTs[, from := "male"]



## NOTE most favoured
ggplot(ECG, aes(acati, faris, fill = from)) +
  geom_bar(stat = "identity") +
  scale_fill_colorblind() +
  geom_errorbar(
    data = ECGTs,
    aes(ymin = faris - 1.96 * faris.sd, ymax = faris + 1.96 * faris.sd),
    width = 0, col = 2
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
ECM0 <- dcast(ECM0, iso3 + g_whoregion + acat + acati + variable ~ type,
  value.var = "value"
)
ECM0[, c("qty", "to", "from") := tstrsplit(variable, split = "\\.")]
ECM0$acat <- factor(ECM0$acat, levels = agz, ordered = TRUE)
ECM0$acati <- factor(ECM0$acati, levels = agz, ordered = TRUE)
## merge in population for infectees
ECM0 <- merge(ECM0, pop[, .(iso3, to = sex, acati, popto = value)],
  by = c("iso3", "to", "acati")
)
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
ECGL[, c("fari", "fari.sd") :=
  .(ari / tot, (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2))]
ECGL[!is.finite(fari.sd), fari.sd := 0] # kids manually 0

## utility function to format X / Y as percent
XoYpc <- function(num, den) {
  mid <- num$fari / den$fari
  mid.sd <- mid * sqrt((num$fari.sd / num$fari)^2 + (den$fari.sd / den$fari)^2)
  hi <- min(1, mid + 1.96 * mid.sd)
  lo <- max(0, mid - 1.96 * mid.sd)
  paste0(
    format(round(100 * mid)),
    "% (95%UI: ", format(round(100 * lo)),
    "% to ", format(round(100 * hi)), "%)"
  )
}


## NOTE fraction of everyone that are male adults
NS[!acat %in% c("0-4", "5-14"), sum(pop.male)] /
  NS[, sum(pop.male + pop.female)] # 37 % of population

## out text
ok <- 1
outtxt <- list()

## to men
den <- ECGL[
  to == "male" & !acati %in% c("0-4", "5-14"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
num <- ECGL[
  from == "male" & !acat %in% c("0-4", "5-14") & # from men
    to == "male" & !acati %in% c("0-4", "5-14"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% transmission to men, from men",
  value = XoYpc(num, den)
)
ok <- ok + 1

## to women
den <- ECGL[
  to == "female" & !acati %in% c("0-4", "5-14"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
num <- ECGL[
  from == "male" & !acat %in% c("0-4", "5-14") & # from men
    to == "female" & !acati %in% c("0-4", "5-14"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% transmission to women, from men",
  value = XoYpc(num, den)
)
ok <- ok + 1

## to kids
den <- ECGL[
  acati %in% c("0-4", "5-14"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
num <- ECGL[
  from == "male" & !acat %in% c("0-4", "5-14") & # from men
    acati %in% c("0-4", "5-14"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% transmission to children, from men",
  value = XoYpc(num, den)
)
ok <- ok + 1

## to anyone
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  from == "male" & !acat %in% c("0-4", "5-14"), # from men
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% transmission from men",
  value = XoYpc(num, den)
)
ok <- ok + 1

## to anyone, from M25-44
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  from == "male" & acat %in% c("25-34", "35-44"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% transmission from M25-44",
  value = XoYpc(num, den)
)
ok <- ok + 1

## to anyone, from M55+
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  from == "male" & acat %in% c("55-64", "65+"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% transmission from M55+",
  value = XoYpc(num, den)
)
ok <- ok + 1

## to anyone, from A15-24
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  acat %in% c("15-24"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% transmission from A15-24",
  value = XoYpc(num, den)
)
ok <- ok + 1

## to anyone, from A25-44
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  acat %in% c("25-34", "35-44"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% transmission from A25-44",
  value = XoYpc(num, den)
)
ok <- ok + 1

## to anyone, from A55+
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  acat %in% c("55-64", "65+"),
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% transmission from A55+",
  value = XoYpc(num, den)
)
ok <- ok + 1


## to men
den <- ECGL[, .(fari = sum(fari), fari.sd = ssum(fari.sd))]
num <- ECGL[
  to == "male",
  .(fari = sum(fari), fari.sd = ssum(fari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% exposure to men ",
  value = XoYpc(num, den)
)
ok <- ok + 1

## population data: fraction male adults
val <- NS[
  !acat %in% c("0-4", "5-14"),
  sum(pop.male)
] / NS[, sum(pop.male + pop.female)]
outtxt[[ok]] <- data.table(
  qty = "% of population men >=15",
  value = paste0(format(round(100 * val)))
)
ok <- ok + 1

## population data: fraction male adults 15-44
val <- NS[
  acat %in% c("25-34", "35-44"),
  sum(pop.male)
] / NS[, sum(pop.male + pop.female)]
outtxt[[ok]] <- data.table(
  qty = "% of population M25-44",
  value = paste0(format(round(100 * val)))
)
ok <- ok + 1

## counterfactual reduction
den <- ECG1[, .(fari = ari, fari.sd = ari.sd)]
num <- ECG0[, .(fari = ari, fari.sd = ari.sd)]
outtxt[[ok]] <- data.table(
  qty = "% exposure if M were F ",
  value = XoYpc(num, den)
)
ok <- ok + 1
outtxt[[ok]] <- data.table(
  qty = "% reduction if M were F ",
  value = XoYpc(den-num, den)
)
ok <- ok + 1

## number of countries
outtxt[[ok]] <- data.table(
  qty = "No. countries",
  value = length(unique(ECM$iso3))
)
ok <- ok + 1

## % population
outtxt[[ok]] <- data.table(
  qty = "% global population",
  value = round(
    1e2 * unique(ECM[, .(iso3, to, acati, popto)])[, sum(popto)] / totpop,
    digits = 1
  )
)
ok <- ok + 1

## % TB
outtxt[[ok]] <- data.table(
  qty = "% global TB",
  value = round(
    1e2 * E[, sum(TB)] / 10844410,
    digits = 1
  )
)
ok <- ok + 1

## exposure to A15-24
den <- ECGL[, .(fari = sum(ari), fari.sd = ssum(ari.sd))]
num <- ECGL[
  acati %in% c("15-24"),
  .(fari = sum(ari), fari.sd = ssum(ari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% exposure in A15-24",
  value = XoYpc(num, den)
)
ok <- ok + 1

## exposure A55+
den <- ECGL[, .(fari = sum(ari), fari.sd = ssum(ari.sd))]
num <- ECGL[
  acati %in% c("55-64", "65+"),
  .(fari = sum(ari), fari.sd = ssum(ari.sd))
]
outtxt[[ok]] <- data.table(
  qty = "% exposure in A55+",
  value = XoYpc(num, den)
)
ok <- ok + 1


## join
outtxt <- rbindlist(outtxt)
outtxt

fwrite(outtxt, file = here("output/outtxt.csv"))


## going by acato ------infector
ECM[acat == "65+"]
ECM[, mean(mn), by = acat] # from mean per cap contrib 20% higher in 65+
ECM[, sum(popto), by = acat] # const
ECM[, sum(popto), by = acati] # const
ECM[, unique(acati)]


tmp <- ECM[iso3 == "ZAF" & to == "female" & acati == "25-34"]
tmp <- ECM[iso3 == "VNM" & to == "female" & acati == "25-34"]
tmp <- ECM[iso3 == "VNM" & to == "female",
  .(mn = sum(mn * popto)),
  by = .(acat, from)
]
tmp[, tot := sum(mn)]
ggplot(tmp, aes(acat, mn / tot, col = from, group = from)) +
  geom_line()


ECGi <- ECM[, .(ari = sum(mn * popto), ari.sd = ssum(sd * popto)),
  by = .(acat, from, to)
] # total exposure happening, by age of source
ECGi[, c("tots", "tots.sd") := .(sum(ari), ssum(ari.sd))]
ECGi[, c("faris", "faris.sd") :=
  .(
    ari / tots,
    (ari / tots) * sqrt((ari.sd / ari)^2 + (tots.sd / tots)^2)
  )]
ECGi$acat <- factor(ECGi$acat, levels = agz, ordered = TRUE)


## CI data
ECGTsi <- ECGi[,
  .(faris = sum(faris), faris.sd = ssum(faris.sd)),
  by = .(acat, to)
]
ECGTsi[, sum(faris)]
ECGTsi[, from := "male"]

## make label text for plot
TXTi <- ECGi[, .(tot = sum(ari)), by = .(acat, to)]
TXTi <- merge(ECGi[, .(acat, from, to, ari, faris, faris.sd)],
  TXTi,
  by = c("acat", "to")
  )
TXTi <- TXTi[from == "male"]
TXTi[acat == "65+" & to == "male", faris := faris - 1e-3]
TXTi[, pcnt := paste0(format(round(100 * ari / tot)), "%")]
TXTi[!is.finite(faris.sd), pcnt := NA_character_]


## NOTE most favoured 2
GPb <- ggplot(ECGi, aes(acat, faris, fill = from)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    data = ECGTsi,
    aes(ymin = faris - 1.96 * faris.sd, ymax = faris + 1.96 * faris.sd),
    width = 0, col = 2
  ) +
  scale_fill_colorblind() +
  ggh4x::facet_wrap2(~to, axes = "all") +
  theme_classic() +
  ggpubr::grids() +
  scale_y_continuous(label = percent) +
  geom_text(data = TXTi, aes(acat, faris + 3e-3, label = pcnt), col = cl) +
  xlab("Age group of infector (years)") +
  ylab("Proportion of all transmission from each group") +
  ggtitle("transmission") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(2, "lines"), # or 3
    strip.text = element_text(face = "italic", size = 12),
    strip.background = element_blank(),
    strip.placement = "outside"
  )
GPb

ggsave(file = here("output/ARI_from.png"), w = 10, h = 5)


## ----------------
## combined plot:

ggarrange(GPa, GPb,
  ncol = 1, labels = c("A", "B"),
  common.legend = TRUE
)

ggsave(file = here("output/ARI_BOTH.png"), w = 10, h = 10)

fn <- here("output/figs")
if (!file.exists(fn)) dir.create(fn)

ggsave(file = here("output/figs/fig1.pdf"), w = 10, h = 10)


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


## --- regional aggregates over sex
ECGA <- ECM[, .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
),
by = .(acat, from, g_whoregion)
] # total exposure happening, by age of source
ECGA[, c("tot", "tot.sd") := .(sum(ari), ssum(ari.sd)), by = g_whoregion]
ECGA[, c("ari", "ari.sd") :=
  .(
    ari / tot,
    (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2)
  )]
ECGA <- merge(ECGA, regkey, by = "g_whoregion")
ECGA$region <- factor(ECGA$region, levels = whozt, ordered = TRUE)
ECGA$acat <- factor(ECGA$acat, levels = agz, ordered = TRUE)

## aggregated over from also
BRT1 <- ECM[, .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
),
by = .(acat, g_whoregion)
]
BRT1[, c("tot", "tot.sd") := .(sum(ari), ssum(ari.sd)), by = g_whoregion]
BRT1[, c("ari", "ari.sd") :=
  .(
    ari / tot,
    (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2)
  )]
BRT1 <- merge(BRT1, regkey, by = "g_whoregion")
BRT1$region <- factor(BRT1$region, levels = whozt, ordered = TRUE)
BRT1$acat <- factor(BRT1$acat, levels = agz, ordered = TRUE)


## plot
ggplot(ECGA, aes(acat, ari, fill = from)) +
  geom_bar(stat = "identity") +
  scale_fill_colorblind() +
  ggh4x::facet_wrap2(~region, axes = "all") +
  theme_classic() +
  ggpubr::grids() +
  scale_y_continuous(label = percent) +
  xlab("Age of infector") +
  ylab("Proportion of all infection") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(2, "lines"), # or 3
    strip.text = element_text(face = "italic"), # , size = 12
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "top"
  )

ggsave(file = here("output/ARIA_from_reg.png"), w = 7, h = 5)


## most aggregate version
ECGR <- ECM[, .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
),
by = .(g_whoregion, acati, from)
]
ECGR[, c("tot", "tot.sd") := .(sum(ari), ssum(ari.sd)), by = g_whoregion]
ECGR[, c("fari", "fari.sd") :=
  .(
    ari / tot,
    (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2)
  )]
ECGR <- merge(ECGR, regkey, by = "g_whoregion")
ECGR$region <- factor(ECGR$region, levels = whozt, ordered = TRUE)
ECGR$acati <- factor(ECGR$acati, levels = agz, ordered = TRUE)

## aggregating also over from
BRT2 <- ECM[, .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
),
by = .(g_whoregion, acati)
]
BRT2[, c("tot", "tot.sd") := .(sum(ari), ssum(ari.sd)), by = g_whoregion]
BRT2[, c("fari", "fari.sd") :=
  .(
    ari / tot,
    (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2)
  )]
BRT2 <- merge(BRT2, regkey, by = "g_whoregion")
BRT2$region <- factor(BRT2$region, levels = whozt, ordered = TRUE)
BRT2$acati <- factor(BRT2$acati, levels = agz, ordered = TRUE)

## plot
ggplot(ECGR, aes(acati, fari, fill = from)) +
  geom_bar(stat = "identity") +
  scale_fill_colorblind() +
  facet_wrap(~region, scales = "free") +
  theme_classic() +
  ggpubr::grids() +
  scale_y_continuous(label = percent) +
  xlab("Age of infectee") +
  ylab("Proportion of all exposure to each group") +
  ## geom_text(data = TXT, aes(acati, fari + 3e-3, label = pcnt), col = cl) +
  guides(color = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(2, "lines"), # or 3
    strip.text = element_text(face = "italic"), #, size = 12
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "top"
  )

ggsave(file = here("output/ARIA_to_reg.png"), w = 7, h = 5)

## both
BRT <- rbind(
  BRT1[, .(region, acat,
    quantity = "transmission", from = "all",
    value = ari, v.sd = ari.sd
  )],
  BRT2[, .(region,
    acat = acati,
    quantity = "exposure", from = "all",
    value = fari, v.sd = fari.sd
  )]
)
BRT[!is.finite(v.sd), c("value", "v.sd") := 0.0]

cvz <- c(6, 7)
ggplot(
  BRT,
  aes(acat, value,
    col = quantity,
    group = paste(region, quantity)
  )
) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = value - 1.96 * v.sd,
      ymax = value + 1.96 * v.sd,
      fill = quantity
    ),
    alpha = 0.3, col = NA
  ) +
  scale_color_manual(values = clz[cvz]) +
  scale_fill_manual(values = clz[cvz]) +
  facet_wrap(~region, scales = "free") +
  theme_classic() +
  ggpubr::grids() +
  scale_y_continuous(label = percent) +
  xlab("Age group (years)") +
  ylab("Proportion of all exposure to or transmission from each group") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(2, "lines"), # or 3
    strip.text = element_text(face = "italic"), # , size = 12
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "top",
    legend.title = element_blank()
  )


ggsave(file = here("output/ARIB_to_reg2.png"), w = 10, h = 7)

ggsave(file = here("output/figs/fig2.pdf"), w = 10, h = 7)




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

ggsave(file = here("output/ARI_inreg.png"), w = 6, h = 5)


## pc TB by region?
EG <- merge(E, whokey, by = "iso3")
EG <- EG[, .(TB = sum(TB)),
  by = .(sex = toupper(sex), acat = age_group, g_whoregion)
]
EG[, sex := ifelse(sex == "F", "female", "male")]
EG[acat == "65plus", acat := "65+"]

popr <- melt(NS, id = c("iso3", "acat"))
popr[, sex := ifelse(grepl("f", variable), "female", "male")]
popr <- merge(popr, whokey, by = "iso3")
popr <- popr[, .(n = sum(value)), by = .(g_whoregion, sex, acat)]

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


## ---- combined figure

BECR <- merge(ECGA[, .(region, from, acat, ari)],
  ECGR[, .(region, from, acat = acati, fari)],
  by = c("region", "from", "acat")
)
BRM <- melt(BECR, id = c("region", "from", "acat"))
BRM <- rbind(BRM[from == "male"], BRT, fill = TRUE)

BRM[, quantity := ifelse(variable == "fari", "exposure", "transmission")]



## male only data
## to:
BRM1 <- ECM[, .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
),
by = .(acat, g_whoregion, sex = from)
]


BRM1[, c("tot", "tot.sd") := .(sum(ari), ssum(ari.sd)), by = g_whoregion]
BRM1[, c("ari", "ari.sd") :=
  .(
    ari / tot,
    (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2)
  )]
BRM1 <- merge(BRM1, regkey, by = "g_whoregion")
BRM1$region <- factor(BRM1$region, levels = whozt, ordered = TRUE)
BRM1$acat <- factor(BRM1$acat, levels = agz, ordered = TRUE)
## from:
BRM2 <- ECM[from == "male", .(
  ari = sum(mn * popto),
  ari.sd = ssum(sd * popto)
),
by = .(g_whoregion, acati, sex = to)
]

BRM2[, c("tot", "tot.sd") := .(sum(ari), ssum(ari.sd)), by = g_whoregion]
BRM2[, c("fari", "fari.sd") :=
  .(
    ari / tot,
    (ari / tot) * sqrt((ari.sd / ari)^2 + (tot.sd / tot)^2)
  )]
BRM2 <- merge(BRM2, regkey, by = "g_whoregion")
BRM2$region <- factor(BRM2$region, levels = whozt, ordered = TRUE)
BRM2$acati <- factor(BRM2$acati, levels = agz, ordered = TRUE)

## both:
BRM <- rbind(
  BRM1[sex=='male', .(region, acat,
    quantity = "transmission", from = "male",
    value = ari, v.sd = ari.sd
  )],
  BRM2[sex=='male', .(region,
    acat = acati,
    quantity = "exposure", from = "male",
    value = fari, v.sd = fari.sd
  )]
)
BRM[!is.finite(v.sd), c("value", "v.sd") := 0.0]
## males and all
BRB <- rbind(BRT, BRM)

TXTb <- dcast(BRB, region + acat + quantity ~ from, value.var = "value")
TXTb[, txt := paste0(round(1e2 * male / all), "%")]
TXTb[all == 0, txt := NA_character_]
TXTb[, from := "male"]
TXTb[, v.sd := 0.0]
TXTb[, value := male]
TXTb[, c("all", "male") := NULL]

## plot:
symbs <- c("\u25CF", "\u2642")
GP <- ggplot(
  BRB,
  aes(acat, value,
    col = quantity, lty = from, shape = from,
    group = paste(region, quantity, from)
  )
) +
  geom_line() +
  geom_point(
    size = I(3),
    stroke = 3
  ) +
  scale_shape_manual(values = symbs) +
  geom_ribbon(
    aes(
      ymin = value - 1.96 * v.sd,
      ymax = value + 1.96 * v.sd,
      fill = quantity, alpha = from
    ),
    col = NA
  ) +
  geom_errorbar(
    data = BRB[from == "male"],
    aes(
      ymin = value - 1.96 * v.sd,
      ymax = value + 1.96 * v.sd,
      col = quantity
    ), width = 0.25
  ) +
  ## ggrepel::geom_text_repel(
  ##   data = TXTb,
  ##   aes(acat,
  ##     value,
  ##     label = txt
  ##   ),
  ##   point.padding = 0.2,
  ##   ## min.segment.length = 0,
  ##   box.padding = 0.3,
  ##   size = 2.5
  ## ) +
  scale_alpha_manual(values = c(0.3, 0)) +
  scale_color_manual(values = clz[cvz]) +
  scale_fill_manual(values = clz[cvz]) +
  facet_wrap(~region, scales = "free") +
  theme_classic() +
  ggpubr::grids() +
  scale_y_continuous(label = percent) +
  xlab("Age") +
  ylab("Proportion of all exposure to or transmission from each group") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(2, "lines"), # or 3
    strip.text = element_text(face = "italic"), # , size = 12
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "top"
  )
ggsave(GP, file = here("output/ARIB_to_reg.png"), w = 10, h = 7)
