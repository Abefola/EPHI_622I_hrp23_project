
# Abebe Fola
# Feb 2023


library(boot)
library(confintr)
set.seed(123)

# Calculate 95% CIs for weighted proportions using bootstrapping

 #### BY DISTRICT ####

#Metema
prop_res <- c(0.83,0.83,0.71)
pf_N <- c(24,30,178)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

metema <- ci_proportion(x= x_N, n = N, R = 2000, type = "bootstrap")
metema$estimate
metema$interval

# Quara
prop_res <- c(1,0.87,0.78)
pf_N <- c(42, 39, 435)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

quara <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
quara$estimate
quara$interval

#Tegede
prop_res <- c(0.67,1,0.83)
pf_N <- c(33,65,196)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

tegede <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
tegede$estimate
tegede$interval

#W Armachiho
prop_res <- c(0.67,0.80,0.73)
pf_N <- c(32,75,187)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

WA <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
WA$estimate
WA$interval

#Itang
prop_res <- c(1,1,0.857)
pf_N <- c(5,2,144)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

itang <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
itang$estimate
itang$interval

#Kule
prop_res <- c(0.80,0.67,0.87)
pf_N <- c(7,4,460)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

kule <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
kule$estimate
kule$interval

#Ahferom
prop_res <- c(0.67,1,0.80)
pf_N <- c(16,18,83)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

ahf <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
ahf$estimate
ahf$interval

#Atseged Tsimbila
prop_res <- c(0.75,0.95,0.89)
pf_N <- c(26,48,86)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

AT <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
AT$estimate
AT$interval

#Gulomekeda
prop_res <- c(1,0.67,0.91)
pf_N <- c(2,3,14)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

gulo <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
gulo$estimate
gulo$interval

#K/Humera
prop_res <- c(1,0.93,0.77)
pf_N <- c(21,39,112)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

KH <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
KH$estimate
KH$interval

#L. Adiabo
prop_res <- c(0,1,0.89)
pf_N <- c(11,22,111)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

LA <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
LA$estimate
LA$interval

#T/Adaibo
prop_res <- c(0.5,0.83,0.92)
pf_N <- c(6,10,51)
N = sum(pf_N)
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

TA <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
TA$estimate
TA$interval

#### BY REGION ####

# AMHARA
prop_res <- c(metema$estimate, quara$estimate, tegede$estimate, WA$estimate)
pf_N <- c(232,516,294,294)
N = sum(pf_N)
amhara_n = N
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

am <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
am$estimate
am$interval

# GAMBELA
prop_res <- c(itang$estimate, kule$estimate)
pf_N <- c(151, 471)
N = sum(pf_N)
gambela_n = N
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

ga <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
ga$estimate
ga$interval

# TIGRAY
prop_res <- c(ahf$estimate, AT$estimate, gulo$estimate, KH$estimate, 
              LA$estimate, TA$estimate)
pf_N <- c(117, 160, 19, 172, 144, 67)
N = sum(pf_N)
tigray_n = N
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

tg <- ci_proportion(x = x_N, n = N, R = 2000, type = "bootstrap")
tg$estimate
tg$interval

 #### OVERALL ####
prop_res <- c(am$estimate, ga$estimate, tg$estimate)
pf_N <- c(amhara_n, gambela_n, tigray_n)
N = sum(pf_N)
overall_n = N
x <- weighted.mean(prop_res, pf_N)
x_N <- x * N

ci_proportion(x = x_N, n = N, R = 3000, type = "bootstrap")
