##########################################################
##########################################################
###########                                 ##############
#########     Starts here: CV 31 años big farms  #########
###########                                 ##############
##########################################################
##########################################################
#setwd("Z:\\Adrian\\2020\\Chp3_data")
setwd("C:\\Users\\Adrian\\Google Drive\\doutorado\\2020\\dataChp3")
spatial<-read.csv("CVz31anos.csv", sep=","); nrow(spatial)
spatialC<-read.csv("CVz31anos_sept2021.csv", sep=","); nrow(spatialC)
setwd("C:\\Users\\Adrian\\Google Drive\\doutorado\\2021\\Chp3_papers")
spatialB<-read.csv("CVz31anos_antigo.csv", sep=","); nrow(spatialB)

################################################
##3What the fu*** is going on:###########
nrow(spatialC)-nrow(spatialB)
subsetdiffer = spatialC[spatialC$codigo_ibg %in% spatialB$codigo_ibg,] 
nrow(subsetdiffer)
subsetdiffer$cc<-paste(subsetdiffer$crop, subsetdiffer$codigo_ibg)
spatialB$cc<-paste(spatialB$crop, spatialB$codigo_ibg)
spatial$cc<-paste(spatial$crop, spatial$codigo_ibg)

sb = subsetdiffer[!subsetdiffer$cc=="cafe 2910909",] ###Only 

plot(spatialB$CV~sb$CV)
#################################################
setwd("C:\\Users\\Adrian\\Google Drive\\doutorado\\2020\\dataChp3")
#spatial<-read.csv("CVz31anos.csv", sep=",")     ##Large crop area (>1000 ha) per crop
#spatial<-read.csv("CVz31anos_sept2021.csv", sep=","); nrow(spatial)
spatial<-read.csv("CVz31anos_052022.csv", sep=","); nrow(spatial)
head(spatial); hist(spatial$CV); summary(spatial$CV)

#spatial<-read.csv("CVz31anos_small.csv", sep=",")    ##Small crop area (>100 ha) per crop
#spatial<-read.csv("CVz31anos_Pol.csv", sep=",")  ##Only dependent crops, area >15 ha per crop

summary(spatial$P_dependency)
spatial<-spatial[complete.cases(spatial),]
nrow(spatial)
spatial<-spatial[spatial$CV>0,]  ##66 com valores iguais a 0
#spatial_spd<-spatial[spatial$P_dependency>=0.05,]

setwd("C:/Users/Adrian/Google Drive/doutorado/2020/dataChp3")
mesor<-read.csv("mesoRegiones.csv", sep=",");head(mesor)
levels(as.factor(mesor$NM_MESO)) ##139 meso regioes
mesor$meso<-mesor$NM_MESO; mesor$uf<-mesor$SIGLA_UF
spatial<- merge(spatial, mesor[,c(1,6,7)],by="codigo_ibg", all.x=T )
nrow(spatial) ##3973 for new Dataset
head(spatial)


### Agregando las coordenadas x , y del centor de cada municipio
setwd("C:/Users/Adrian/Google Drive/doutorado/2019/municipios")
coord<-read.csv("CoordenadasBrasil.csv", sep=";")
head(coord)
spatial<-merge(spatial, coord, by="codigo_ibg", all.X=T); head(spatial)


setwd("G:/Meu Drive/doutorado/2019/municipios")
area_muncp<-read.csv("area_muncp.csv", sep=";")
head(area_muncp)
spatial<-merge(spatial, area_muncp[,c(4,6)], by="codigo_ibg", all.X=T)
nrow(spatial)#3932
##Incorporating crop cover
spatial$crop_cover<-100*(spatial$area/spatial$area_muncp)
summary(spatial$crop_cover)
oq<-spatial[spatial$crop_cover>80,]
oq$oq<-(oq$area/oq$area_muncp)*100; oq

###Lets calculate Forest area and used to calcula new edge density (fragmentation)
##Forest area:
spatial$forest_area<-spatial$area_muncp*spatial$FC/100
#Fragmentation = frag
spatial$frag<-spatial$edgeT/spatial$forest_area*100 ##Multipliquei por 100 para ficar em km/km2
summary(spatial$frag)
##Fragmentation metric matheron:
spatial$mathe<-spatial$edgeT*1000/(sqrt(spatial$forest_area)*sqrt(spatial$area_muncp))
summary(spatial$mathe)

#####Identifying each coffee species as different crops:
canefora<-spatial[spatial$crop=="cafe"&spatial$contribution=="essential",]; nrow(canefora) #47
arabica<-spatial[spatial$crop=="cafe"&spatial$contribution=="modest",]; nrow(arabica) #274
##Substituting café for each coffee species
spatial[spatial$crop=="cafe"&spatial$contribution=="essential",3]<-"C.canephora"
spatial[spatial$crop=="cafe"&spatial$contribution=="modest",3]<-"C.arabica"


####Tranforming bean into the modest pollianotr dependency category:
spatial[spatial$crop=="feijao",4]<-0.25
spatial[spatial$crop=="feijao",5]<-"modest"

##Checking that the change was made
require(doBy)
summaryBy(FC+P_dependency~crop+contribution, FUN=mean, data=spatial)
cultivos<-summaryBy(FC+P_dependency~crop+contribution, FUN=mean, data=spatial)
# Distr de la var de respuesta CV (manual) y CV2 (summaryBy)
##############################################################
require(fitdistrplus)
spatial<-spatial[complete.cases(spatial),]; spatial<-subset(spatial, spatial$CV>0) ##hay 69 0, sin variación en el tiempo
spatial2<-spatial[spatial$contribution!="essential",];nrow(spatial); nrow(spatial2)
respuesta<-1/(spatial$CV)
respuesta<-spatial$z
respuesta<-spatial$CV
#respuesta<-mandioca$CV

fit.normal=fitdist(respuesta,"norm")
fit.lognormal=fitdist(respuesta,"lnorm")
fit.gamma=fitdist(respuesta,"gamma")
par(mfrow=c(1,2))
cdfcomp(list(fit.normal,fit.lognormal,fit.gamma),horizontals=F, addlegend=T,legendtext=c("Normal","log", "Gamma"))
qqcomp(list(fit.normal,fit.lognormal,fit.gamma),addlegend=T,legendtext=c("Normal","log", "Gamma"))
gofstat(list(fit.normal, fit.lognormal,fit.gamma))$aic
###mejor ajuste: gamma para CV
par(mfrow=c(1,1))
#############Best distribution Gamma despite eliminating non normal distributions###############################


##avaliemos la correlación de las variables restantes
##Only edgeD and FC correlated and prep -> Pr_seasonality
##############################################################
names(spatial)
var<-spatial[,c(4,6,7,9:18,29,34:37)]
names(var)
rav<-cor(var, method= "spearman")  ###corelación tipo pearson

round(rav,2)
library(Hmisc)
rav2<-rcorr(as.matrix(var, type="spearman"))
rav2

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

ad<-flattenCorrMatrix(round(rav2$r, 2), round(rav2$P, 2))
str(ad)
forest<-rbind(ad[ad$row=="FC",], ad[ad$column=="FC",]); forest
ad2<-ad[ad$cor>=0.5|ad$cor<=-0.5,]; ad2
#######################################################3

#####Transforming he predictive variables, by scaling########
names(spatial)

spatial$Pollination.s<-scale(spatial$P_dependency, scale=T)
spatial$edgeD.s<-scale(spatial$edgeD, scale=T) ##Correlacionada con FC?
spatial$edgeT.s<-scale(spatial$edgeT, scale=T)
spatial$patch.s<-scale(spatial$patchD, scale=T)
spatial$FC.s<-scale(spatial$FC, scale=T)
spatial$FC.log<-log(spatial$FC)
spatial$area.s<-scale(spatial$area, scale=T)
spatial$crop_cover.log<-log(spatial$crop_cover)
spatial$frag.log<-log(spatial$frag)

#climaticas
spatial$temp.s<-scale(spatial$temp, scale=T)
spatial$t_range.s<-scale(spatial$t_range, scale=T)
spatial$prec.s<-scale(spatial$prec, scale=T)
spatial$Pp_seasonality.s<-scale(spatial$Pp_seasonality, scale=T)
spatial$prepCV.s<-scale(spatial$prepCV, scale=T)
spatial$soil.s<-scale(spatial$soilQ, scale=T)
spatial$div_area.s<-scale(spatial$div_area.mean, scale=T) ###Ojo variable del paisaje


#Selección de modelos
require(MuMIn)
require(lme4)
require(nlme)
require(visreg)
require(car)
require(ggplot2)
require(ggeffects)

cont=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000))  ##### In order to increase the number of interactions#######


####Considering this full models after the second meeting with lucas###
#Suggestions:
#### 1) Substitute Patch density or Edge density for Total Forest edge
#### 2) Eliminate Pollinator dependency variable
#### 3) Allow for slope to vary for each crop regarding forest cover
CVlocas<-lme(CV~poly(log(FC),2)*(div_area.s+soil.s+log(crop_cover)+P_dependency)+prepCV+offset(area_muncp),
             random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
vif(CVlocas)
summary(CVlocas)

CVlocas2<-lme(CV~poly(log(edgeD),2)*(div_area.s+soil.s+log(crop_cover))+prepCV+offset(area_muncp),
             random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
vif(CVlocas2)
summary(CVlocas2)

CVlocas3<-lme(CV~poly(log(edgeD),2)*(soil.s+log(crop_cover)+div_area.s+P_dependency)+prepCV,
              random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
vif(CVlocas3)
summary(CVlocas3)
model.sel(CVlocas2, CVlocas3, CVlocas, rank="AIC")
##LEts test including municipality in the random variables
CVlocas3b<-lme(CV~poly(log(edgeD),2)*(soil.s+log(FC)+div_area.s)+log(crop_cover)+prepCV,
              random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
#Does not work
##Given that div_area is among the mos relevant varibles, i'll test the interactions with it:
CVdiv_Area<-lme(CV~div_area.s*(poly(log(edgeD),2)+soil.s+(crop_cover)+P_dependency)+prepCV,
               random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)

vif(CVdiv_Area) #how are we doing? 
model.sel(CVdiv_Area, CVlocas3, rank="AIC")


##Lets compare among forest features:
fc<-lme(CV~log(FC), random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
edgeD<-lme(CV~log(edgeD), random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
patch<-lme(CV~patch.s,random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
S_fc<-lme(CV~poly(log(FC), 2),
          random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
S_edgeD<-lme(CV~poly(log(edgeD), 2),
             random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
S_patch<-lme(CV~poly(patch.s, 2),
             random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)

model.sel(fc, edgeD, patch, S_fc, S_edgeD, S_patch, rank="AIC")

#Given that Patch density forest feature is a better predictor than the others, 
#we proceed to run a full model with forest interacting with demand landscape features:
CVmodel<-lme(CV~poly(patch.s,2)*(log(crop_cover)+div_area.s+P_dependency)+soil.s+prepCV,
             random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
vif(CVmodel)

#Allowing FC slop to vary for each crop not only the intercept
CVlucas<-lme(CV~poly(log(edgeD),2)*(log(crop_cover)+div_area.s+P_dependency+soil.s)+prepCV,
             random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
vif(CVlucas)
model.sel(CVlucas, CVlocas,CVlocas3, rank="AIC")
#Including FC to vary per crop makes the model much better


## normal plot of standardized residuals by gender
qqnorm(CVlucas, ~ resid(., type = "p") | crop, abline = c(0, 1))
qqnorm(CVlucas, ~ resid(., type = "p") , abline = c(0, 1))
qqnorm(CVlocas, ~ resid(., type = "p") , abline = c(0, 1))

## normal plots of random effects
qqnorm(CVlucas, ~ranef(.))

spatial$FC.log<-log(spatial$FC)
CVlucas_log<-lme(CV~poly(FC.log,2)*log(edgeT)+log(crop_cover)+prepCV,
             random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)



### Model selection accounting for multicollinearity ##21h56 ->  até as  =  min to run the model 
options(na.action = na.fail)
model <- CVmodel## Define model to be used
Allmodels <- dredge(model, rank = "AIC")#, m.lim=c(0, 3), extra= c(max.r)) ###Run dredge specifying the number of predictor variables and including the max.r function
#NCM <- get.models(Allmodels, subset = max.r<=0.6) ##Retrieve non-collinear models (max.r <=0.6)
model.sel(Allmodels, rank = "AIC") ##Final model selection table
impNCM<-importance(Allmodels); head(impNCM, 16)

par(mar = c(9,3,3,2))
layout(matrix(1:1,1,1))
barplot(impNCM[1:15], main= "Relative importance of the variables", las=2, 
        col=c("#1BD465","#E6800A", "#808080", "#00FFFF","#FFFF00",  "#808080","#E6800A",
              "#00FFFF","#E6800A"), cex.names=1)

#####After eliminating edgeT from the full model, this was the one selected:
selected<-lme(CV~poly(patch.s, 2)*(P_dependency+div_area.s)+prepCV,
              random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
vif(selected)
selected2<-lme(CV~patch.s*(P_dependency+div_area.s)+prepCV,
              random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
sselected<-lme(CV~poly(patch.s,2)*(P_dependency+div_area.s)+poly(log(FC),2)+prepCV,
               random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
vif(sselected)
model.sel(selected, selected2,sselected, rank="AIC")


##WHat if we allow the slopes to change for every crop:
#1st) using FC in the random slopes
slopes_F<-lme(CV~poly(patch.s,2)*div_area.s+prepCV,
           random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
#2nd) USing patch density in the random slopes and eliminating P_dependency
slopes_p<-lme(CV~poly(patch.s,2)*div_area.s+prepCV,
            random=list(~patch.s|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
#3rd) Using both patch density for slopes and maintaining P_dependency
slopes_Fp<-lme(CV~poly(patch.s,2)*(P_dependency+div_area.s)+prepCV,
              random=list(~patch.s|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)

model.sel(slopes_p,slopes_F,slopes_Fp,selected, rank="AIC") ##Using the variability for each crop is good,
##############################################################but no need for P_dependency
vif(slopes_Fp)


##Lets compare the effect of incorporating the variations in slopes in relationship to FC:
Final_wo_slopes<-lme(CV~poly(log(FC),2)*(log(edgeT))+log(crop_cover)+prepCV,
                     random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
slopes_Final<-lme(CV~poly(log(FC),2)*(log(edgeT))+log(crop_cover)+prepCV,
              random=list(~FC|crop,~FC|meso), correlation=corRatio(form=~x+y),  data=spatial)
vif(slopes_Final)
slopes_null<-lme(CV~1,
                  random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
wo_slopes_null<-lme(CV~1,
                 random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)

model.sel(Final_wo_slopes,slopes_Final, slopes_null ,wo_slopes_null,  rank="AIC")

r.squaredLR(slopes_Final)
r.squaredLR(Final_wo_slopes)
r.squaredLR(slopes_p)
r.squaredLR(selected)

###Lets plot the residuals
plot(selected, resid(., type = "p") ~ fitted(.)| crop, abline = 0)
plot(selected, resid(., type = "p") ~ fitted(.), abline = 0)
plot(selected, resid(., type = "p") ~ FC, abline = 0)
plot(selected, crop ~ resid(.))

qqnorm(selected, ~ resid(., type = "p") , abline = c(0, 1))
qqnorm(selected, ~ranef(.))

library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)

plot_model(selected, vline.color = "grey", show.values = TRUE, value.offset = .3)+theme_classic()


##What if we plot log(FC) instead of FC:
bosque<-ggpredict(CVlucas_log, terms = c("FC.log [all]","edgeT [20]", "crop"))
plot(bosque, add.data = T, ci=F)+ylab("Yield instability (CV)")+xlab("Forest cover (%)")+
  theme_ggeffects(base_size = 16)+theme_classic()

plot(patch)+theme_classic()

crop<-ggpredict(slopes_p, terms = c("patch.s [all]","div_area.s [0]", "crop"), type="re")
plot(crop, add.data = T, ci=F)+ylab("Yield instability (CV)")+xlab("Forest patch density")+
  theme_ggeffects(base_size = 16)+theme_classic()

plot(crop)+theme_classic()


#########################################
###   Plotting the fixed effects      ###
#########################################
##forest patch density
bosq<-ggpredict(selected, terms = c("patch.s [all]"))
plot(bosq, add.data=T, ci=T)+ylab("Yield instability (CV)")+xlab("Forest patch density (# of patches/ha)")+
  theme_ggeffects(base_size = 16)+theme_classic()+geom_line(col="dark green", size=1.3)

summary(spatial$patchD)

##Crop diversity (Shannon index):
cdiv<-ggpredict(selected, terms = c("div_area.s [all]"))
plot(cdiv, add.data=T, ci=F)+ylab("Yield instability (CV)")+xlab("Crop diversity (Shannon index)")+
  theme_ggeffects(base_size = 16)+theme_classic()+geom_line(col="light green", size=1.3)

##The interaction between patch density and crop diversity
mm<-ggpredict(selected, terms = c("patch.s [all]","div_area.s [-1.5,0,1.5]"))
plot(mm, add.data=T, line.size = 1.2, colors = "hero")+ylab("Yield instability (CV)")+xlab("Forest patch density (# of patches/ha)")+
  theme_ggeffects(base_size = 16)+theme_classic()
plot(mm,  line.size = 1.2, colors = "hero")+ylab("Yield instability (CV)")+xlab("Forest patch density (# of patches/ha)")+
  theme_ggeffects(base_size = 16)+theme_classic()


###Pollinator dependency interaction with patch density:
crops<-ggpredict(selected, terms = c("patch.s [all]","P_dependency [0,0.25,0.95]"))
plot(crops, add.data=T,line.size = 1.2, ci=T)+ylab("Yield instability (CV)")+xlab("Forest patch density (# patches/ha)")+
  theme_ggeffects(base_size = 24)+theme_classic() 
plot(crops, line.size = 1.2, ci=T)+ylab("Yield instability (CV)")+xlab("Forest patch density (# patches/ha)")+
  theme_ggeffects(base_size = 24)+theme_classic() 

edge<-ggpredict(slopes_Final, terms = c("edgeT [all]","FC [2, 8, 20]"))
plot(edge, add.data=T)+ylab("Yield instability (CV)")+xlab("Forest cover (%)")+theme_ggeffects(base_size = 16)+theme_classic() 
plot(edge)+ylab("Yield instability (CV)")+xlab("Forest cover (%)")+theme_ggeffects(base_size = 16)+theme_classic() 

#####What happens if we transform P_dependency into categorical
int2_contribution<-lme(CV~poly(log(FC),2)*(log(edgeT)+contribution)+log(crop_cover)+div_area.s+prepCV,
            random=list(~1|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
model.sel(test4, int1e, int2ep,int3ep,int2_contribution, rank="AIC")
pollination<-ggpredict(int2_contribution, terms = c("FC [all]","contribution"))
plot(pollination, add.data = T)+ylab("Yield instability (CV)")+xlab("Forest cover (%)")+
  theme_ggeffects(base_size = 16)+theme_classic() 

###plotting as continus variable
pollination<-ggpredict(int2ep, terms = c("FC [all]","P_dependency [0, 0.25, 0.95]"))
plot(pollination, add.data = T)+ylab("Yield instability (CV)")+xlab("Forest cover (%)")+
  theme_ggeffects(base_size = 16)+theme_classic() 


##Covariavles plots:
#crop_cover
ff<-ggpredict(slopes_Final, terms = c("crop_cover [all]"))
plot(ff)+ylab("Yield instability (CV)")+xlab("Crop Cover (%)")+theme_ggeffects(base_size = 16)
plot(ff, add.data = T, ci=T)+ylab("Yield instability (CV)")+xlab("Crop Cover (%)")+theme_ggeffects(base_size = 16) +
  theme_classic()+geom_line(col="orange", size=1.3)
#Precipitation instability
pp<-ggpredict(slopes_Final, terms = c("prepCV"))
plot(pp)+ylab("Yield instability (CV)")+xlab("Annual precipitation instability")+
  theme_ggeffects(base_size = 16)+theme_classic()
plot(pp, add.data = T, ci=T)+ylab("Yield instability (CV)")+xlab("Annual precipitation instability")+
  theme_ggeffects(base_size = 16)+theme_classic()+geom_line(col="red", size=1.3)


#############################3
##Extra questions:
####1)How bad is a model without poly(FC,2)
slopes_Final<-lme(CV~poly(log(FC),2)*(log(edgeT))+log(crop_cover)+prepCV,
                  random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
slopes_sem_poly<-lme(CV~log(FC)*(log(edgeT))+log(crop_cover)+prepCV,
                  random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
model.sel(slopes_Final, slopes_sem_poly, rank="AIC") ###AIC difference of 12.68 points!

##Eventually try running Dredge without poly to see if tendencies are mantained?


####2)What if we replace FC with edge density, will the model be equally good? 
####                                           are hipothesis mantain?
slopes_edge<-lme(CV~poly(log(edgeD),2)*(log(edgeT))+log(crop_cover)+prepCV,
                     random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
slopes_edge_sempoly<-lme(CV~log(edgeD)*(log(edgeT))+log(crop_cover)+prepCV,
                 random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
slopes_edge_semedgeT<-lme(CV~log(edgeD)+log(crop_cover)+prepCV,
                         random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)
model.sel(slopes_Final, slopes_sem_poly, slopes_edge, slopes_edge_sempoly,slopes_edge_semedgeT, rank="AIC")
vif(slopes_edge)
####3) What if we substitute, the term which is in polynomic form (FC -> edgeT)
slopes_edgeT<-lme(CV~poly(log(edgeT),2)*(log(FC))+log(crop_cover)+prepCV,
                  random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)

slopes_aver<-lme(CV~poly(log(FC),2)*(log(edgeT))+log(edgeD)+log(crop_cover)+prepCV,
                  random=list(~FC|crop,~1|meso), correlation=corRatio(form=~x+y),  data=spatial)

model.sel(slopes_Final, slopes_edgeT,slopes_aver,slopes_edge, rank="AIC") 
##sloped with edgeD is always better, and including edgeD as additive effects improves the model



##How does edge Total, and edge Density related to Forest Cover
plot(spatial$edgeD~spatial$FC)
plot(spatial$edgeD~spatial$edgeT)

##The interaction between forest cover and forest edge
ee<-ggpredict(slopes_aver, terms = c("edgeD [all]"))
plot(ee, add.data=T, line.size = 1.2, colors = "hero")+ylab("Yield instability (CV)")+xlab("edge density (km/ha)")+
  theme_ggeffects(base_size = 16)+theme_classic()
plot(ee, add.data=T,ci=T)+ylab("Yield instability (CV)")+xlab("Forest cover (%)")+
  theme_ggeffects(base_size = 16)+theme_classic() +ylim(c(0, 0.4))


###Plotting the fixed effects
crops<-ggpredict(slopes_edge, terms = c("edgeD [all]","edgeT [3]", "crop"),type="re")
plot(crops, add.data=T, ci=T)+ylab("Yield instability (CV)")+xlab("Forest cover (%)")+
  theme_ggeffects(base_size = 16)+theme_classic() +geom_line(col="green")
plot(crops)+ylab("Yield instability (CV)")+xlab("Forest cover (%)")+theme_ggeffects(base_size = 16)+theme_classic() 


