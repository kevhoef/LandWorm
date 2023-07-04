##### MODELES AVEC EC #####

library(lme4)
library(car)
library("corrplot")
library(lmerTest)
var_quanti_cor <- SBT_error[,c("Nbr_EW", "ADJV","Div_LIEC")]
Matrice_cor <-cor(na.omit(var_quanti_cor), method = c("pearson"))
corrplot(Matrice_cor, type="upper", order="hclust", tl.col="black", tl.srt=45)
corrplot(Matrice_cor, method="number", type="upper")


SBT_error_MR= SBT_error_MR%>% 
  mutate(MR_rounded = MR*100) %>%
  mutate(MR_rounded = round(MR_rounded, 0))


any(is.na(SBT_error_train$Nbr_EW))==F # if true = OK
any(is.na(SBT_error_train$ID_Site))==F # if true = OK
any(is.na(SBT_error_train$Div_LIEC))==F # if true = OK
any(is.na(SBT_error_train$ADJV))==F # if true = OK
any(is.na(SBT_error_train$nuages))==F # if true = OK
any(is.na(SBT_error_train$Train_count))==F # if true = OK
any(is.na(SBT_error_train$observateur))==F # if true = OK




# Analyses satistiques
## Modèle à 4 facteurs fixes + 1 aléatoire, SANS interaction
library(gam)
SBT_error_MR$ID_Site <- as.ordered(SBT_error_MR$ID_Site)
model1<-gam(MR ~ EC* (Nbr_EW + Div_LIEC+ ADJV + nuages)+ s(ID_Site),data = SBT_error_MR)
summary(model1)

#MR
hist(SBT_error_MR$MR_rounded)

str(SBT_error_train)

model1<-glmer(MR_rounded ~ EC* (Nbr_EW + Div_LIEC+ ADJV + nuages) +(1|ID_Site),data = SBT_error_MR, family= "poisson"(link = "log"),nAGQ = 1)
model1<-stan_glmer(MR_rounded ~ EC* (Nbr_EW + Div_LIEC+ ADJV + nuages) +(1|ID_Site),data = SBT_error_MR, family= "poisson")
model1<-glmmTMB(MR_rounded ~ EC* (Nbr_EW + Div_LIEC+ ADJV + nuages) +(1|ID_Site),data = SBT_error_MR, family= "poisson")
model1<-glmmPQL(MR_rounded ~ EC* (Nbr_EW + Div_LIEC+ ADJV + nuages), random = ~1|ID_Site,data = SBT_error_MR, family= "quasipoisson")
Anova(model1, type = 3)
overdisp_fun(model1)
performance::check_overdispersion(model1)

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

str(SBT_error_UR)
summary(SBT_error_MR_EPI$Nbr_EW) # outlier ?
summary(SBT_error_MR_EPI$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_MR_EPI$ADJV)
summary(SBT_error_MR_EPI$nuages)
summary(model1)
model1
Anova(model1,type=3)

## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(MR_rounded~nuages,data=SBT_error_MR)
leveneTest(MR_rounded~nuages, data=SBT_error_MR)

bartlett.test(MR_rounded~ID_Site,data=SBT_error_MR)
leveneTest(log(MR_rounded+1)~ID_Site, data=SBT_error_MR)

bartlett.test(MR~ID_Site,data=SBT_error_MR)
leveneTest(MR~ID_Site, data=SBT_error_MR)

bartlett.test(sqrt(MR_rounded+1)~EC^2,data=SBT_error_MR)
leveneTest(log10(MR_rounded+1)~EC|nuages, data=SBT_error_MR)
#Independence
#car::durbinWatsonTest(model1) #  => error !!!!
#lmtest::dwtest(model1) # => error !!!!
#performance::check_autocorrelation(model1) # OK

vif(model1)
lmerTest::step(model1)

model1<-lmer(MR ~ EC * (Nbr_EW + Div_LIEC + ADJV + nuages) + (1 | ID_Site) ,data = SBT_error_MR)
#refaire tableau unique (une seule colonne MR + UR)
MuMIn::r.squaredGLMM(model1)


r_squared <- r.squaredGLMM(model1)
# Affichage du pourcentage de variance expliquée
variance_explained <- r_squared * 100
variance_explained

Anova(model1,type=3)
summary(model1)$adj.r.squared # 0.5583325


r_squared <- r.squaredGLMM(model1, type = c("marginal", "conditional"))

# Affichage des coefficients de détermination
r_squared

library(effects)
plot(Effect(c("Div_LIEC","EC"),model1))
plot(Effect(c("ADJV", "EC"),model1))
plot(Effect(c("EC", "nuages"),model1))
plot(Effect(c("EC", "Nbr_EW"),model1))

# Comparaisons multiples de moyennes 'pairwise test' (facteur digestat * sol)
emmeans::comp_tuckey <- emmeans(model1, pairwise ~ EC|ADJV)
multcomp::cld(comp_tuckey,comparison=TRUE, Letters=letters, reversed = T, type = "response")

means <- emmeans(model1, ~ Div_LIEC*EC)
contrasts <- contrast(means, interaction = "pairwise")
print(contrasts)
means <- aggregate(MR ~ EC + Nbr_EW, data = SBT_error_MR, FUN = mean)
print(means)


#### MODELE PAR EC ####


# MR EPIGES
## GAM à 4 facteurs fixes + 1 aléatoire, SANS interaction
library(gam)
SBT_error_MREPI$ID_Site <- as.ordered(SBT_error_MREPI$ID_Site)
model1<-gam(MR_EPI ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_MREPI)
step_model <- step.Gam(model1, scope = list(lower = model1, upper = ~ 1), direction = "both", k = 2, trace = FALSE)
stepcAIC(model1, data = SBT_error_MREPI)

summary(model1)

model <- gam(MR ~ EC * Nbr_EW + EC * Div_LIEC + EC * ADJV + EC * nuages + s(ID_Site), data = SBT_error_MR)

# LMM
hist(SBT_error_UR$UR)

model1<-lmer(log(MR_EPI+1) ~ Nbr_EW + Div_LIEC+ ADJV + nuages +(1|observateur),data = SBT_error_MREPI)

summary(SBT_error_UR_EPI$Nbr_EW) # outlier ?
summary(SBT_error_UR_EPI$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_UR_EPI$ADJV)
summary(SBT_error_UR_EPI$nuages)
summary(model1)
model1
Anova(model1,type=3)
## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(log(MR_EPI+1)~nuages,data=SBT_error_MREPI)
leveneTest(log(MR_EPI+1)~nuages, data=SBT_error_MREPI)

bartlett.test(log(MR_EPI+1)~ID_Site,data=SBT_error_MREPI)
leveneTest(log(MR_EPI+1)~ID_Site, data=SBT_error_MREPI)

#Independence
#car::durbinWatsonTest(model1) #  => error !!!!
#lmtest::dwtest(model1) # => error !!!!
#performance::check_autocorrelation(model1) # OK

vif(model1)
lmerTest::step(model1)

model1<-lmer(log(MR_EPI + 1) ~ Div_LIEC + ADJV + (1 | observateur),data = SBT_error_MREPI)
MuMIn::r.squaredGLMM(model1)


r_squared <- r.squaredGLMM(model1)
# Affichage du pourcentage de variance expliquée
variance_explained <- r_squared * 100
variance_explained

Anova(model1,type=2)



library(effects)
plot(Effect(c("Div_LIEC"),model1))
plot(Effect(c("ADJV"),model1))




#### MR EPA ####
## GAM à 4 facteurs fixes + 1 aléatoire, SANS interaction
library(gam)
SBT_error_MRANE_RH$ID_Site <- as.ordered(SBT_error_MRANE_RH$ID_Site)
model1<-gam(MR_ANE_RH ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_MRANE_RH)
step_model <- step.Gam(model1, scope = list(lower = model1, upper = ~ 1), direction = "both", k = 2, trace = FALSE)
stepcAIC(model1, data = SBT_error_MRANE_RH)

summary(model1)

model1<-gam(MR_ANE_RH ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_MRANE_RH)

# LMM
hist(SBT_error_UR$UR)

model1<-lmer(MR_ANE_RH ~ Nbr_EW + Div_LIEC+ ADJV + nuages +(1|observateur),data = SBT_error_MRANE_RH)

summary(SBT_error_UR_ANE_RH$Nbr_EW) # outlier ?
summary(SBT_error_UR_ANE_RH$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_UR_ANE_RH$ADJV)
summary(SBT_error_UR_ANE_RH$nuages)
summary(model1)
model1
Anova(model1,type=3)
## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(MR_ANE_RH~nuages,data=SBT_error_MRANE_RH)
leveneTest(MR_ANE_RH~nuages, data=SBT_error_MRANE_RH)

leveneTest(MR_ANE_RH~ID_Site, data=SBT_error_MRANE_RH)

#Independence
#car::durbinWatsonTest(model1) #  => error !!!!
#lmtest::dwtest(model1) # => error !!!!
#performance::check_autocorrelation(model1) # OK

vif(model1)
lmerTest::step(model1)

model1<-lmer(MR_ANE_RH ~ ADJV + (1 | ID_Site),data = SBT_error_MRANE_RH)
#refaire tableau unique (une seule colonne UR + UR)
MuMIn::r.squaredGLMM(model1)


r_squared <- r.squaredGLMM(model1)
# Affichage du pourcentage de variance expliquée
variance_explained <- r_squared * 100
variance_explained

Anova(model1,type=2)



library(effects)
plot(Effect(c("ADJV"),model1))



#### MR ANE_BH ####
## GAM à 4 facteurs fixes + 1 aléatoire, SANS interaction
library(gam)
SBT_error_MRANE_BH$ID_Site <- as.ordered(SBT_error_MRANE_BH$ID_Site)
model1<-gam(MR_ANE_BH ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_MRANE_BH)
step_model <- step.Gam(model1, scope = list(lower = model1, upper = ~ 1), direction = "both", k = 2, trace = FALSE)
stepcAIC(model1, data = SBT_error_MRANE_BH)

summary(model1)

model1<-gam(MR_ANE_BH ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_MRANE_BH)

# LMM
hist(SBT_error_UR$UR)

model1<-lmer(MR_ANE_BH ~ Nbr_EW + Div_LIEC+ ADJV + nuages +(1|observateur),data = SBT_error_MRANE_BH)

summary(SBT_error_UR_ANE_BH$Nbr_EW) # outlier ?
summary(SBT_error_UR_ANE_BH$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_UR_ANE_BH$ADJV)
summary(SBT_error_UR_ANE_BH$nuages)
summary(model1)
model1
Anova(model1,type=3)
## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(MR_ANE_BH~nuages,data=SBT_error_MRANE_BH)
leveneTest(MR_ANE_BH~nuages, data=SBT_error_MRANE_BH)

bartlett.test(MR_ANE_BH~ID_Site,data=SBT_error_MRANE_BH)
leveneTest(MR_ANE_BH~ID_Site, data=SBT_error_MRANE_BH)

#Independence
#car::durbinWatsonTest(model1) #  => error !!!!
#lmtest::dwtest(model1) # => error !!!!
#performance::check_autocorrelation(model1) # OK

vif(model1)
lmerTest::step(model1)

model1<-lmer(MR_ANE_BH ~ Div_LIEC + ADJV + (1 | observateur),data = SBT_error_MRANE_BH)
Anova(model1,type=2)
MuMIn::r.squaredGLMM(model1)


r_squared <- r.squaredGLMM(model1)
# Affichage du pourcentage de variance expliquée
variance_explained <- r_squared * 100
variance_explained

Anova(model1,type=2)



library(effects)
plot(Effect(c("ADJV"),model1))
plot(Effect(c("Div_LIEC"),model1))


#### MR END ####
## GAM à 4 facteurs fixes + 1 aléatoire, SANS interaction

SBT_error_MREND= SBT_error_MREND%>% 
  mutate(MR_rounded = MR_END*100) %>%
  mutate(MR_rounded = round(MR_rounded, 0))
library(gam)
SBT_error_MREND$ID_Site <- as.ordered(SBT_error_MREND$ID_Site)
model1<-gam(MR_END ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_MREND)
step_model <- step.Gam(model1, scope = list(lower = model1, upper = ~ 1), direction = "both", k = 2, trace = FALSE)
stepcAIC(model1, data = SBT_error_MREND)

summary(model1)

model1<-gam(MR_END ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_MREND)

# LMM
hist(SBT_error_UR$UR)

model1<-lmer(MR_END ~ Nbr_EW + Div_LIEC+ ADJV + nuages +(1|observateur),data = SBT_error_MREND)
model1<-glmmTMB(MR_rounded ~ EC* (Nbr_EW + Div_LIEC+ ADJV + nuages) +(1|ID_Site),data = SBT_error_MR, family= "poisson")

summary(SBT_error_UR_END$Nbr_EW) # outlier ?
summary(SBT_error_UR_END$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_UR_END$ADJV)
summary(SBT_error_UR_END$nuages)
summary(model1)
model1
Anova(model1,type=3)
## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(log(MR_END+1)~nuages,data=SBT_error_MREND)
leveneTest(log(MR_END+1)~nuages, data=SBT_error_MREND)

bartlett.test(log(MR_END+1)~ID_Site,data=SBT_error_MREND)
leveneTest(log(MR_END+1)~ID_Site, data=SBT_error_MREND)

#Independence
#car::durbinWatsonTest(model1) #  => error !!!!
#lmtest::dwtest(model1) # => error !!!!
#performance::check_autocorrelation(model1) # OK

vif(model1)
lmerTest::step(model1)

model1<-lmer(MR_END ~ Nbr_EW + ADJV + (1 | observateur),data = SBT_error_MREND)
Anova(model1,type=2)
MuMIn::r.squaredGLMM(model1)


r_squared <- r.squaredGLMM(model1)
# Affichage du pourcentage de variance expliquée
variance_explained <- r_squared * 100
variance_explained

Anova(model1,type=2)



library(effects)
plot(Effect(c("Nbr_EW"),model1))
plot(Effect(c("ADJV"),model1))


#### MODELE PAR EC ####


#### UR EPIGES ####
## GAM à 4 facteurs fixes + 1 aléatoire, SANS interaction
library(gam)
SBT_error_UREPI$ID_Site <- as.ordered(SBT_error_UREPI$ID_Site)
model1<-gam(UR_EPI ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_UREPI)
step_model <- step.Gam(model1, scope = list(lower = model1, upper = ~ 1), direction = "both", k = 2, trace = FALSE)
stepcAIC(model1, data = SBT_error_UREPI)

summary(model1)

model <- gam(UR ~ EC * Nbr_EW + EC * Div_LIEC + EC * ADJV + EC * nuages + s(ID_Site), data = SBT_error_UR)

# LMM
hist(SBT_error_UR$UR)

model1<-lmer(UR_EPI ~ Nbr_EW + Div_LIEC+ ADJV + nuages +(1|observateur),data = SBT_error_UREPI)

summary(SBT_error_UR_EPI$Nbr_EW) # outlier ?
summary(SBT_error_UR_EPI$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_UR_EPI$ADJV)
summary(SBT_error_UR_EPI$nuages)
summary(model1)
model1
Anova(model1,type=3)
## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(UR_EPI~nuages,data=SBT_error_UREPI)
leveneTest(UR_EPI~nuages, data=SBT_error_UREPI)

bartlett.test(UR_EPI~ID_Site,data=SBT_error_UREPI)
leveneTest(UR_EPI~ID_Site, data=SBT_error_UREPI)

#Independence
#car::durbinWatsonTest(model1) #  => error !!!!
#lmtest::dwtest(model1) # => error !!!!
#performance::check_autocorrelation(model1) # OK

vif(model1)
lmerTest::step(model1)

model1<-lmer(log(UR_EPI + 1) ~ (1 | observateur),data = SBT_error_UREPI)
#refaire tableau unique (une seule colonne UR + UR)
MuMIn::r.squaredGLMM(model1)


r_squared <- r.squaredGLMM(model1)
# Affichage du pourcentage de variance expliquée
variance_explained <- r_squared * 100
variance_explained

Anova(model1,type=2)



library(effects)
plot(Effect(c("Div_LIEC"),model1))
plot(Effect(c("ADJV"),model1))




#### UR EPA ####
## GAM à 4 facteurs fixes + 1 aléatoire, SANS interaction
library(gam)
SBT_error_URANE_RH$ID_Site <- as.ordered(SBT_error_URANE_RH$ID_Site)
model1<-gam(UR_ANE_RH ~ Nbr_EW+ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_URANE_RH)
step_model <- step.Gam(model1, scope = list(lower = model1, upper = ~ 1), direction = "both", k = 2, trace = FALSE)
stepcAIC(model1, data = SBT_error_URANE_RH)

summary(model1)

model1<-gam(UR_ANE_RH ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_URANE_RH)

# LMM
hist(SBT_error_UR$UR)

model1<-lmer(UR_ANE_RH ~ Nbr_EW + Div_LIEC+ ADJV + nuages +(1|ID_Site),data = SBT_error_URANE_RH)

summary(SBT_error_UR_ANE_RH$Nbr_EW) # outlier ?
summary(SBT_error_UR_ANE_RH$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_UR_ANE_RH$ADJV)
summary(SBT_error_UR_ANE_RH$nuages)
summary(model1)
model1
Anova(model1,type=3)
## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(UR_ANE_RH~nuages,data=SBT_error_URANE_RH)
leveneTest(UR_ANE_RH~nuages, data=SBT_error_URANE_RH)

bartlett.test(UR_ANE_RH~ID_Site,data=SBT_error_URANE_RH)
leveneTest(UR_ANE_RH~ID_Site, data=SBT_error_URANE_RH)

#Independence
#car::durbinWatsonTest(model1) #  => error !!!!
#lmtest::dwtest(model1) # => error !!!!
#performance::check_autocorrelation(model1) # OK

vif(model1)
lmerTest::step(model1)

model1<-lmer(UR_ANE_RH ~ Nbr_EW + ADJV + (1 | ID_Site),data = SBT_error_URANE_RH)
Anova(model1,type=2)
MuMIn::r.squaredGLMM(model1)


r_squared <- r.squaredGLMM(model1)
# Affichage du pourcentage de variance expliquée
variance_explained <- r_squared * 100
variance_explained

Anova(model1,type=2)



library(effects)
plot(Effect(c("ADJV"),model1))
plot(Effect(c("Nbr_EW"),model1))



#### UR ANE_BH ####
## GAM à 4 facteurs fixes + 1 aléatoire, SANS interaction
library(gam)
SBT_error_URANE_BH$ID_Site <- as.ordered(SBT_error_URANE_BH$ID_Site)
model1<-gam(UR_ANE_BH ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_URANE_BH)
step_model <- step.Gam(model1, scope = list(lower = model1, upper = ~ 1), direction = "both", k = 2, trace = FALSE)
stepcAIC(model1, data = SBT_error_URANE_BH)

summary(model1)

model1<-gam(UR_ANE_BH ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_URANE_BH)

# LMM
hist(SBT_error_UR$UR)

model1<-lmer(UR_ANE_BH ~ Nbr_EW + Div_LIEC+ ADJV + nuages +(1|ID_Site),data = SBT_error_URANE_BH)

summary(SBT_error_UR_ANE_BH$Nbr_EW) # outlier ?
summary(SBT_error_UR_ANE_BH$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_UR_ANE_BH$ADJV)
summary(SBT_error_UR_ANE_BH$nuages)
summary(model1)
model1
Anova(model1,type=3)
## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(UR_ANE_BH~nuages,data=SBT_error_URANE_BH)
leveneTest(UR_ANE_BH~nuages, data=SBT_error_URANE_BH)

bartlett.test(UR_ANE_BH~ID_Site,data=SBT_error_URANE_BH)
leveneTest(UR_ANE_BH~ID_Site, data=SBT_error_URANE_BH)

#Independence
#car::durbinWatsonTest(model1) #  => error !!!!
#lmtest::dwtest(model1) # => error !!!!
#performance::check_autocorrelation(model1) # OK

vif(model1)
lmerTest::step(model1)

model1<-lmer(UR_ANE_BH ~ ADJV + (1 | ID_Site),data = SBT_error_URANE_BH)
Anova(model1,type=2)
MuMIn::r.squaredGLMM(model1)


r_squared <- r.squaredGLMM(model1)
# Affichage du pourcentage de variance expliquée
variance_explained <- r_squared * 100
variance_explained

Anova(model1,type=2)



library(effects)
plot(Effect(c("ADJV"),model1))
plot(Effect(c("Div_LIEC"),model1))


#### UR END ####
## GAM à 4 facteurs fixes + 1 aléatoire, SANS interaction
SBT_error_UREND= SBT_error_UREND%>% 
  mutate(UR_rounded = UR_END*100) %>%
  mutate(UR_rounded = round(UR_rounded, 0))

library(gam)
SBT_error_UREND$ID_Site <- as.ordered(SBT_error_UREND$ID_Site)
model1<-gam(UR_rounded ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_UREND)
step_model <- step.Gam(model1, scope = list(lower = model1, upper = ~ 1), direction = "both", k = 2, trace = FALSE)
stepcAIC(model1, data = SBT_error_UREND)

summary(model1)

model1<-gam(UR_END ~ Div_LIEC+ ADJV + nuages+ s(ID_Site),data = SBT_error_UREND)

# LMM
hist(SBT_error_UR$UR)

model1<-lmer(UR_rounded ~ Nbr_EW + Div_LIEC+ ADJV + nuages +(1|ID_Site),data = SBT_error_UREND)
model1<-glmmTMB(UR_rounded ~ EC* (Nbr_EW + Div_LIEC+ ADJV + nuages) +(1|ID_Site),data = SBT_error_UR, family= "poisson")

summary(SBT_error_UR_END$Nbr_EW) # outlier ?
summary(SBT_error_UR_END$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_UR_END$ADJV)
summary(SBT_error_UR_END$nuages)
summary(model1)
model1
Anova(model1,type=3)
## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(log(UR_rounded+1)~nuages,data=SBT_error_UREND)
leveneTest(log(UR_rounded+1)~nuages, data=SBT_error_UREND)

bartlett.test(UR_END~ID_Site,data=SBT_error_UREND)
leveneTest(UR_END~ID_Site, data=SBT_error_UREND)

#Independence
#car::durbinWatsonTest(model1) #  => error !!!!
#lmtest::dwtest(model1) # => error !!!!
#performance::check_autocorrelation(model1) # OK

vif(model1)
lmerTest::step(model1)

model1<-lmer(UR_END ~ Div_LIEC + nuages + (1 | ID_Site),data = SBT_error_UREND)
Anova(model1,type=2)
MuMIn::r.squaredGLMM(model1)


r_squared <- r.squaredGLMM(model1)
# Affichage du pourcentage de variance expliquée
variance_explained <- r_squared * 100
variance_explained

Anova(model1,type=2)



library(effects)
plot(Effect(c("Div_LIEC"),model1))
plot(Effect(c("nuages"),model1))
comp_tuckey <- emmeans(model1, pairwise ~ nuages, adjust="none")
multcomp::cld(comp_tuckey,comparison=TRUE, Letters=letters, reversed = T, type = "response")
tapply(perte_pds,digestat,mean)
sort(tapply(perte_pds,digestat,mean))
hsd.out<-HSD.test(model1,"nuages")
hsd.out$groups

summary(glht(model1, linfct = mcp(nuages = "Tukey")), test = adjusted("bonferroni"))

difflsmeans: difflsmeans(model1,test.effs=NULL,ddf="Satterthwaite")
?adjust




test =lsmeans(model1, ~ nuages)
pairs(test, by="nuages")


#MR_EPA
SBT_error_MR_EPA = drop_na(SBT_error, "ID_Site","ADJV", "MR_EPA", "nuages", "ADJV", "Div_LIEC", "Nbr_EW")
hist(SBT_error_MR_EPA$MR_EPA)

comparisons <- glht(model1, linfct = mcp(EC = "Tukey"))
comparison_results <- summary(comparisons)

# Afficher les résultats des comparaisons
print(comparison_results)


model1<-lmer(MR_EPA~ Nbr_EW + Div_LIEC + ADJV + nuages + (1|ID_Site),data = SBT_error_MR_EPA)
summary(SBT_error_MR_EPA$Nbr_EW) # outlier ?
summary(SBT_error_MR_EPA$Div_LIEC) # LIEC 3 et 4 sur représenté ?
summary(SBT_error_MR_EPA$ADJV)
summary(SBT_error_MR_EPA$nuages)

Anova(model1,type=2)

## Vérification des conditions 
#Normalité
shapiro.test(residuals(model1))

#Homoscédasticité
bartlett.test(MR_EPI~nuages,data=SBT_error_MR_EPA)

#Independence
car::durbinWatsonTest(model1) #  => error !!!!
lmtest::dwtest(model1) # => error !!!!
performance::check_autocorrelation(model1) # PAS OK

vif(model1)
lmerTest::step(model1)

model2<-lmer(MR_EPA ~ ADJV + (1 | ID_Site),data = SBT_error_MR_EPI)
Anova(model2,type=2)


plot(Effect(c("ADJV"),model2))

# Comparaisons multiples de moyennes 'pairwise test' (facteur digestat * sol)
emmeans::comp_tuckey <- emmeans(model1, pairwise ~ nuages)
multcomp::cld(comp_tuckey,comparison=TRUE, Letters=letters, reversed = T, type = "response")

#################################################
#################################################
#################################################
#################################################

devtools::install_github("erblast/easyalluvial")


# Préparer les données des liens (links)
links <- data.frame(source = CS_error$FIEC,
                    target = CS_error$LIEC,
                    value = CS_error$Nbr_EW)

# Préparer les données des nœuds (nodes)
nodes <- data.frame(ID = unique(c(CS_error$FIEC, CS_error$LIEC)))

# Créer le diagramme Sankey
sankey <- ggplot() +
  geom_parallel_sets(data = links,
                     aes(axis = (source),
                         axis2 = (target),
                         weight = value),
                     alpha = 0.5,
                     edge_width = 0.5,
                     node_width = 0.5,
                     node_colour = "lightblue",
                     node_linetype = "solid",
                     node_labels = TRUE) +
  theme_minimal()
# Afficher le diagramme
print(sankey)


# Préparer les données des liens (links)
links <- CS_error %>%
  group_by(FIEC, LIEC) %>%
  summarise(Nbr_EW = sum(Nbr_EW)) %>%
  rename(source = FIEC, target = LIEC, value = Nbr_EW)

# Préparer les données des nœuds (nodes)
nodes <- data.frame(ID = unique(c(CS_error$FIEC, CS_error$LIEC)))
links <- as.data.frame(links)

# Create the Sankey diagram with networkD3
sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target", Value = "value")

# Display the diagram
sankey
# Créer le diagramme Sankey avec networkD3
sankey <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target", Value = "value")

# Afficher le diagramme
sankey


alluvial_wide(CS_error , bins = 5 , bin_labels = c('FIEC','LIEC'),fill_by = 'all_flows')
data_wide = data_wide %>%
  mutate_at( vars(categoricals), as.factor ) %>%
  mutate( car_id = row_number() )

library(networkD3)
library(htmlwidgets)
# Supposons que vos données sont stockées dans un objet appelé "donnees"
# Créez une liste unique de tous les nœuds (sources, cibles)
nodes <- data.frame(ID = unique(c(CS_error$FIEC, CS_error$LIEC)))

# Créez un tableau de liens (links) avec les colonnes "source", "target" et "value"
links <- data.frame(source = CS_error$FIEC,
                    target = CS_error$LIEC,
                    value = CS_error$Nbr_EW)
view(CS_error)
sankey <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "source", Target = "target",
                        Value = "value", NodeID = "ID")
# Afficher le diagramme
sankey


sankey <- sankeyNetwork(
  Links = CS_error,
  Source = "FIEC",
  Target = "LIEC",
  Value = "Nbr_EW",
  NodeID = unique(c(CS_error$FIEC, CS_error$LIEC))
)

Links <- as.data.frame(Links)

x11()


# Convert the Links object to a data frame
links <- data.frame(source = CS_error$FIEC,
                    target = CS_error$LIEC,
                    value = CS_error$Nbr_EW)




# Supposons que vous ayez 5 facteurs pour chaque colonne
factors_CE <- unique(CS_error$FIEC)
factors_LIEC <- unique(CS_error$LIEC)


sankey <- sankeyNetwork(Links = links,
                        Nodes = data.frame(ID = unique(c(factors_CE, factors_LIEC))),
                        Source = "source", Target = "target",
                        Value = "value", NodeID = "ID",
                        sinksRight = FALSE, height = 600, width = 800,
                        nodeWidth = 30, nodePadding = 10,
                        iterations = 0)
grDevices
plot(sankey)
?Devices
windows()
view(CS_error)
CS_error$Nbr_EW_OK <- round(CS_error$Nbr_EW, digits = 0)
round

x11()
length(unique(c(links$source, links$target)))
length(nodes$area)


summary(CS_error)

links <- data.frame(source = CS_error$FIEC,
                    target = CS_error$LIEC,
                    value = CS_error$Nbr_EW)

# Préparez les données des nœuds (nodes)
nodes <- data.frame(ID = unique(c(links$source, links$target)))
nodes$area <- rep("EPI", "EPA", "ANS", "END", "LIEC_X", nrow(nodes))  # Remplacez "Area A" par les valeurs appropriées pour chaque nœud
nodes$area <- c("EPI", "EPA", "ANS", "END", "LIEC_X")

sankey <- sankeyNetwork(Links = links,
                        Nodes = nodes,
                        Source = "source",
                        Target = "target",
                        Value = "value",
                        NodeID = "ID",
                        sinksRight = FALSE,
                        height = 600,
                        width = 800,
                        nodeWidth = 30,
                        nodePadding = 10,
                        iterations = 0)

# Create a data frame for the Nodes
NodesData <- data.frame(ID = unique(c(CS_error$FIEC, CS_error$LIEC)))

# Create the Sankey diagram
sankey <- sankeyNetwork(
  Links = Links,
  Nodes = NodesData,
  Source = "FIEC",
  Target = "LIEC",
  Value = "Nbr_EW",
  NodeID = "ID"
)

# Display the diagram
sankey



links <- data.frame(source = CS_error$FIEC,
                    target = CS_error$LIEC,
                    value = CS_error$Nbr_EW)

# Préparer les données des nœuds (nodes)
nodes <- data.frame(ID = unique(c(CS_error$FIEC, CS_error$LIEC)))

# Créer le diagramme Sankey
sankey <- sankeyNetwork(Links = links,
                        Nodes = nodes,
                        Source = "source",
                        Target = "target",
                        Value = "value",
                        NodeID = "ID",
                        sinksRight = FALSE,
                        height = 600,
                        width = 800,
                        nodeWidth = 30,
                        nodePadding = 10,
                        iterations = 0)


# Préparer les données des liens (links)
links <- data.frame(source = CS_error$FIEC,
                    target = CS_error$LIEC,
                    value = CS_error$Nbr_EW)

# Préparer les données des nœuds (nodes)
nodes <- data.frame(ID = unique(c(CS_error$FIEC, CS_error$LIEC)))

# Créer le diagramme Sankey
sankey <- SankeyDiagram(links = links, 
                        nodes = nodes, 
                        width = 800, 
                        height = 600)

# Afficher le diagramme
sankey



Links <- as.data.frame(CS_error)

# Create a data frame for the Nodes
NodesData <- data.frame(ID = unique(c(CS_error$FIEC, CS_error$LIEC)))

# Create the Sankey diagram
sankey <- sankeyNetwork(
  Links = Links,
  Nodes = NodesData,
  Source = "FIEC",
  Target = "LIEC",
  Value = "Nbr_EW",
  NodeID = "ID"
)

sankey <- sankeyNetwork(
  Links = CS_error,
  Source = "FIEC",
  Target = "LIEC",
  Value = "Nbr_EW",
  NodeID = unique(c(CS_error$FIEC, CS_error$LIEC))
)
??sankeyNetwork
sankeyNetwork(Links, Nodes, Source, Target, Value, NodeID, NodeGroup = NodeID,
              LinkGroup = NULL, units = "",
              colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"), fontSize = 7,
              fontFamily = NULL, nodeWidth = 15, nodePadding = 10, margin = NULL,
              height = NULL, width = NULL, iterations = 32, sinksRight = TRUE)

###################################################
##################################################
###################################################
###################################################














#TEST###

SBT_fusion_error_comptage_obs=right_join(sbt_eni_comptage_VDT, sbt_eni_VDT_aggreg, by ="id_obs")
anti_join(sbt_eni_comptage_VDT, sbt_eni_VDT_aggreg, by = "id_obs")
write.table(SBT_fusion_error_comptage_obs, "D:/Home/khoeffner/Downloads/sbt_test.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")

###### S?l?ction des variables d'int?r?ts ######  
SBT_error =select(SBT_fusion_error_comptage_obs,ID_Site,id_obs,Annee, MR_EPI, MR_ANE_RH, MR_ANE_BH, MR_END, UR_EPI, UR_ANE_RH, UR_ANE_BH, UR_END,Nbr_EW,AB_field,AD,SA,JV,ADJV, Div_LIEC, ID_Site,nuages,pluie,vent,temperature_air)
SBT_error_MR$ID_Site <- as.factor(SBT_error_MR$ID_Site)
SBT_error_MR$nuages <- as.factor(SBT_error_MR$nuages)
SBT_error_MR$ID_Site <- as.factor(SBT_error_MR$ID_Site)


write.table(SBT_error, "D:/Home/khoeffner/Downloads/SBT_error.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")

unique(SBT_error$pluie)
unique(SBT_error$ID_Site)
unique(SBT_error$ID_Site)
str(SBT_error_MR)


SBT_error = drop_na(SBT_error, "ID_Site","ADJV", "MR", "pluie", "ADJV", "Div_LIEC", "Nbr_EW")
SBT_error$ID_Site <- as.factor(SBT_error$ID_Site)
SBT_error$nuages <- as.factor(SBT_error$nuages)


unique(SBT_error_final$pluie)
unique(SBT_error_final$Nbr_EW)
unique(SBT_error_final$ADJV)
unique(SBT_error_final$MR_ANS)




sum(SBT_error_final$Nbr_EW)
sum(SBT_error_final$AB_field)

sum(SBT_error_final$Nbr_EW)
sum(SBT_error_final$AB_field)
SBT_error = filter(SBT_error,!is.na(AB_field))




test=table(SBT_error$ID_Site)
table(SBT_error_final$ID_Site,SBT_error_final$Annee)
write.table(test, "D:/Home/khoeffner/Downloads/test.csv", sep=",", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")


SBT_error_test = filter(SBT_error_final,!is.na(AB_field))
SBT_error_test = filter(SBT_error_test,!is.na(Nbr_EW))

taux_error = (sum(SBT_error_test$Nbr_EW) * 100) /sum(SBT_error_test$AB_field)

SBT_ACM =select(SBT_fusion_error_comptage,nuages,pluie,vent,temperature_air)
SBT_ACM <- within(SBT_ACM, {ID_Site <- NULL })

write_xlsx(SBT_fusion_error_comptage,"D:/Home/khoeffner/Downloads/sbt_database.xlsx")

res.mca = MFA(SBT_ACM, group=c(3,3,3), ncp = 5, graph = TRUE)

data_list <- list(airlines, nuages,pluie,vent,temperature_air)
