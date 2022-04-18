#-------------------------------------------------------------------------------
#1) Import data
#-------------------------------------------------------------------------------
#median climate data
datsc=read.csv("Scalesia_Data_Median_Climate_26Feb2022.csv")
#mean climate data
datsc=read.csv("Scalesia_Data_Mean_Climate_26Feb2022.csv")
datsc=datsc[,-1]

#Some summary stats
mean(datsc$mean.Tomeg); sd(datsc$mean.Tomeg) #Tomega (aka curve breadth)
mean(datsc$mean.Popt); sd(datsc$mean.Popt) #Popt (Optimum rate of photosynthesis)
mean(datsc$mean.Topt); sd(datsc$mean.Topt) #Topt (Optimum Temperature of photosynthesis)
#-------------------------------------------------------------------------------
#2) Load libraries
#-------------------------------------------------------------------------------
pkgs=c("phytools", "ape", "phangorn", "corrplot", "nlme", "cowplot", "MASS", "caper", "geiger", "emmeans", "piecewiseSEM", "sjPlot", "ggplot2")
ldpcks = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
}
lapply(pkgs, ldpcks)

#-------------------------------------------------------------------------------
#3) Make phylogenetic tree
#-------------------------------------------------------------------------------
#Import phylogenetic tree from Fernández-Mazuecos, M., P. Vargas, R. A. McCauley, D. Monjas, A. Otero, J. A. Chaves, J. E. Guevara Andino, et al. 2020. The Radiation of Darwin’s Giant Daisies in the Galápagos Islands. Current Biology 30:4989-4998.e7.
#Please cite authors above
sctree1=read.tree("DEF_c85m20r_TreePL.pruned.tre") 
#Label tree tips
sctree1$tip.label=gsub("S_", "Scalesia ", sctree1$tip.label)

#Drop species with no data
sctree=drop.tip(sctree1, sctree1$tip.label[c(2,3,6,11)])

#Make tree ultrmetric
sctree<-force.ultrametric(sctree)

#order datsc to match order of tree
datsc=datsc[match(sctree$tip.label, datsc$species),]

#Plot a tree with Omega values
ggPal=colorRampPalette(c("dark green","green"))
as.matrix(datsc$mean.Popt)
traits=as.matrix(datsc$mean.Popt)
row.names(traits)=datsc$species
kolors=ggPal(21)[as.numeric(cut((traits[,1]),breaks = 21))]
phylosig(sctree, datsc$mean.Tomeg, method="K", nsim=1000, test=T)
plotTree.wBars(sctree,(traits[,1]/5),fsize=1, tip.labels=T,  type="phylogram",  
               col=kolors, method="plotTree", lwd=2 , mar=c(5,6,5,1), args.boxplot=list((axes=T)))
axisPhylo(1, backward = T)
mtext("Ma",1,2, font = 2)
#dev.off()

#-------------------------------------------------------------------------------
#4) Correlation plot phylogenetically corrected data
#-------------------------------------------------------------------------------
#Create data frame with data of interest
datsp2=data.frame(species=datsc$species, 
                  mean.Topt=datsc$mean.Topt, #mean optimum tmeperaure for photosynthesis
                  mean.Popt=datsc$mean.Popt, #mean optimum rate of photosynthesis
                  mean.Tomeg=datsc$mean.Tomeg, #mean breadth of temperature response curve
                  mn.elw=datsc$mn.elw, #mean effective leaf width (cm)
                  area=datsc$area, #leaf area
                  bio3=datsc$bio3, 
                  bio5=datsc$bio5,
                  bio12=datsc$bio12, 
                  bio2=datsc$bio2,
                  bio4=datsc$bio4) 

#Re-order data according to tips of phylo. tree
datsp2=datsp2[match(sctree$tip.label, datsp2$species),] 

#Create matrix of data used to test for phylogenetically-correceted trait correlations
datmat2=matrix(as.numeric(unlist(datsp2[,-1])),nrow=nrow(datsp2[,-1]))
colnames(datmat2) <- colnames(datmat2[,-1])
rownames(datmat2) <- datsp2$species

#Make sure order of traits matches phylo. tree tips one more time, (all should be true)
datsp2$species==sctree$tip.label

#Create phylogenetic variance co-variance matrix
obj2=phyl.vcv(datmat2, vcv(sctree),1)
rmat2=cov2cor(obj2$R)
tmat2<-rmat2*sqrt((Ntip(sctree)-2)/(1-rmat2^2))
Pmat2<-2*pt(abs(tmat2),df=21-2,lower.tail=F)

#Name columns of matrix
var_names=c("=T[opt]", "=P[opt]", ":Omega",  "=ELW", "=Area", "=Bio3 (ISO)", "=Bio5 (MTWM)",  "=Bio12 (MAP)", "=Bio2", "=Bio4")
colnames(rmat2) = var_names
rownames(rmat2) = var_names
rownames(Pmat2) = var_names
colnames(Pmat2) = var_names

#pdf("your file path/FIG3.Part.Corrplot.pdf", width = 6, height =6 )
corrplot(rmat2, type = "upper", p.mat = Pmat2, sig.level = 0.05, tl.col="black")
#dev.off()
print(rmat2)
print(Pmat2)

#-------------------------------------------------------------------------------
#5) Correlations among traits plot(not phylogenetically corrected)
#-------------------------------------------------------------------------------
#Similar steps as above
datsp3=data.frame(species=datsc$species, 
                  mean.Topt=datsc$mean.Topt, 
                  mean.Popt=datsc$mean.Popt, 
                  mean.Tomeg=datsc$mean.Tomeg, 
                  mn.elw=datsc$mn.elw,
                  area=datsc$area, 
                  bio3=datsc$bio3, 
                  bio5=datsc$bio5, 
                  bio12=datsc$bio12, 
                  bio2=datsc$bio2,
                  bio4=datsc$bio4) 
datmat3=matrix(as.numeric(unlist(datsp3[,-1])),nrow=nrow(datsp3[,-1]))
colnames(datmat3) <- colnames(datsp3[,-1])
rownames(datmat3) <- datsp3$species

M <- cor(datmat3)
pmat1 <- cor.mtest(datmat3, conf.level = .95)
pmat1=pmat1$p

var_names2=c("=T[opt]", "=P[opt]", ":Omega",  "=ELW", "=Area",  "=Bio3 (ISO)", "=Bio5 (MTWM)",  "=Bio12 (MAP)", "=Bio2", "=Bio4")
colnames(M ) = var_names2
rownames(M ) = var_names2
rownames(pmat1) = var_names2
colnames(pmat1) = var_names2

corrplot(M, p.mat=pmat1, type = "upper", sig.level = 0.05, tl.col="black")

#-------------------------------------------------------------------------------
#6) Test of hypothesis 1 & figure 4
#-------------------------------------------------------------------------------
#Collapse and Resolve Multichotomies
sctree2=multi2di(sctree)


#Phylogenetic generalized least squares for making the figure 4
#Divide bio5 by 10 for correct units
datsc$bio5=datsc$bio5/10

pmodh1=gls(mn.elw~bio5, correlation=corBrownian(1, sctree2), data=datsc, method = "ML")
tab_model(pmodh1)

simpleploth1a=plot_model(pmodh1, type="pred", terms = c("bio5"), title="", geom.alpha=0.5, 
                         axis.title = c( expression("Bio5"~degree~"C"), expression("Effective leaf width (cm)")), 
                         axis.textsize.x=1, axis.textsize.y=1, show.data=T)+   theme_classic()+
  geom_point(aes(y=mn.elw, x=bio5), alpha=0.5, data=datsc)+
  coord_cartesian(ylim = c(0, 11), xlim = c(26.5, 31)) +
  annotate('text', y=0.5, x=26.5, label="Pearson's r = -0.62", hjust=0)+
  annotate('text', y=0, x=26.5, label="p = 0.03", hjust=0)+
  geom_errorbar(aes(ymin=mn.elw.lci, ymax=mn.elw.uci, x=bio5),  alpha=0.25, data=datsc, inherit.aes = FALSE) +
  geom_errorbar(aes(xmin=bio5.lci/10, xmax=bio5.uci/10, y=mn.elw),  alpha=0.25, data=datsc, inherit.aes = FALSE)+
  annotate('text', x=27, y=11, label="A", fontface="bold", hjust=2) 

#Results of phylogenetically correct correlation (see section 4 above) effective leaf width, bio5
print(rmat2[4,7])
print(Pmat2[4,7])

#Results of correlation not phylogenetically corrected (see section 4 above) effective leaf width, bio5
print(M[4,7])
print(pmat1[4,7])


pmodh1b=gls(mn.elw~bio12, correlation=corBrownian(1, sctree2), data=datsc, method = "ML")
tab_model(pmodh1b)
simpleploth1b=plot_model(pmodh1b, type="pred", terms = c("bio12"), title="", geom.alpha=0.5, 
                         axis.title = c( expression("Bio12 (mm)"), expression("Effective leaf width (cm)")), 
                         axis.textsize.x=1, axis.textsize.y=1, show.data=T)+   theme_classic()+
  geom_point(aes(y=mn.elw, x=bio12), alpha=0.5, data=datsc)+
  coord_cartesian(ylim = c(0, 11), xlim = c(250, 1200)) +
  annotate('text', y=0.5, x=750, label="Pearson's r = 0.47", hjust=0)+
  annotate('text', y=0, x=750, label="p = 0.13", hjust=0)+
  geom_errorbar(aes(ymin=mn.elw.lci, ymax=mn.elw.uci, x=bio12),  alpha=0.25, data=datsc, inherit.aes = FALSE) +
  geom_errorbar(aes(xmin=bio12.lci, xmax=bio12.uci, y=mn.elw),  alpha=0.25, data=datsc, inherit.aes = FALSE)+
  annotate('text', x=350, y=11, label="B", fontface="bold", hjust=2) 

#Results of phylogenetically correct correlation (see section 4 above) effective leaf width, bio12
print(rmat2[4,8])
print(Pmat2[4,8])
#Results of correlation not phylogenetically corrected (see section 4 above) effective leaf width, bio12
print(M[4,8])
print(pmat1[4,8])

#pdf("Hypothesis_1R.pdf", width = 8.5, height =4.25 )
cowplot::plot_grid(simpleploth1a, simpleploth1b) 
#dev.off()


#-------------------------------------------------------------------------------
#7) Test of hypothesis 2 & figure 5
#-------------------------------------------------------------------------------
#Dived by 10 for correct units
datsc$bio3=datsc$bio3/10

#Perform phylogenetic least squares
pmod=gls(mean.Tomeg~bio3, correlation=corBrownian(1, sctree2), data=datsc, method = "ML")
tab_model(pmod, show.aic = T)
simpleplot=plot_model(pmod, type="pred", 
                      terms = c("bio3"), 
                      title="", geom.alpha=0.5, 
                      axis.title = c( expression(paste("Isothermality")), expression(Omega~ "("~degree~"C)")), 
                      axis.textsize.x=1, axis.textsize.y=1, show.data=T)+  theme_classic()+
  coord_cartesian(ylim = c(5, 24), xlim = c(5.2, 6.4))+
  geom_point(aes(x=bio3, y=mean.Tomeg), alpha=0.5, data=datsc)+
  geom_errorbar(aes(xmin=bio3.lci/10, xmax=bio3.uci/10, y=mean.Tomeg),  alpha=0.25, data=datsc, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin=mean.Tomeg.lci, ymax=mean.Tomeg.uci, x=bio3),  alpha=0.25, data=datsc, inherit.aes = FALSE) 

#Results of phylogenetically correct correlation (see section 4 above) omega and bio3
print(rmat2[3,6])
print(Pmat2[3,6])

#Results of correlation not phylogenetically corrected (see section 4 above) omega and bio3
print(M[3,6])
print(pmat1[3,6])

#pdf("Isothermality.v.Omega.pdf", width = 4.25, height =4.25 )
simpleplot
#dev.off()


#-------------------------------------------------------------------------------
#8) Test of hypothesis 3 & figure 6
#-------------------------------------------------------------------------------
#Phylogenetic generalized least squares
pmod2=gls(mean.Tomeg~mean.Popt, correlation=corBrownian(1, sctree2), data=datsc, method = "ML")
tab_model(pmod2, show.aic = T)

simpleplot2=plot_model(pmod2, type="pred", terms = c("mean.Popt"), title="", geom.alpha=0.5, 
                       axis.title = c( expression(P[opt]~"("~mu~"mol m"^-2~"s"^-1~")"), expression(Omega~ "("~degree~"C)")), 
                       axis.textsize.x=1, axis.textsize.y=1, show.data=T)+   theme_classic()+
  geom_point(aes(x=mean.Popt, y=mean.Tomeg), alpha=0.5, data=datsc)+
  coord_cartesian(xlim = c(8, 28), ylim = c(3,25)) +
  annotate('text', x=25, y=23, label="Pearson's r = -0.82", hjust=1)+
  annotate('text', x=25, y=21, label="p < 0.01", hjust=1)+
  geom_errorbar(aes(xmin=mean.Popt.lci, xmax=mean.Popt.uci, y=mean.Tomeg),  alpha=0.25, data=datsc, inherit.aes = FALSE) +
  geom_errorbar(aes(ymin=mean.Tomeg.lci, ymax=mean.Tomeg.uci, x=mean.Popt),  alpha=0.25, data=datsc, inherit.aes = FALSE) 

#Results of phylogenetically correct correlation (see section 4 above) Popt and bio3
print(rmat2[2,6])
print(Pmat2[2,6])

#Results of correlation not phylogenetically corrected (see section 4 above) Popt and bio3
print(M[2,6])
print(pmat1[2,6])

#pdf("Popt.v.Omega.pdf", width = 4.25, height =4.25 )
simpleplot2
#dev.off()


#-------------------------------------------------------------------------------
#9) Figure 7 Climate interaction with leaf size 
#-------------------------------------------------------------------------------
#datsc$bio3=datsc$bio3/1000
pmod3a=gls(mean.Popt~bio3, correlation=corBrownian(1, sctree2), data=datsc, method = "ML")
tab_model(pmod3a, show.aic = T)
pmod3=gls(mean.Popt~bio3*mn.elw, correlation=corBrownian(1, sctree2), data=datsc, method = "ML")
summary(pmod3)
tab_model(pmod3, show.aic = T)
intervals(pmod3)
summary(datsc$mn.elw)
summary(datsc$bio3)

#The relationship between popt and bio3 may change based on the leaf morphology of each species.
#Predict you the relationship between popt and bio3 change for species woith different leaf sizes.
effa <- max(datsc$mn.elw) 
eff <- mean(datsc$mn.elw)
effb <- min(datsc$mn.elw) 
(effar <- round(effa,0.5))
(effr <- round(eff,0.5))
(effbr <- round(effb,0.5))
mylist <- list(mn.elw=c(effbr,effr,effar)) 

emtrends(pmod3, ~mn.elw, var="bio3", at=mylist)
(mylist <- list(bio3=seq(5.2,6.2,by=0.05),mn.elw=c(effbr,effar)))
emmip(pmod3,mn.elw~bio3,at=mylist, CIs=TRUE)
emtrends(pmod3, pairwise ~mn.elw, var="bio3",at=mylist, adjust="Tukey")
summary(pmod3$coefficients)

(mylist <- list(bio3=5.2, mn.elw=c(effbr,effar)))
emmeans(pmod3, pairwise ~ bio3*mn.elw, at=mylist)

(mylist <- list(bio3=6.2, mn.elw=c(effbr,effar)))
emmeans(pmod3, pairwise ~ bio3*mn.elw, at=mylist)

tab_model(pmod3)
#Plot data:
(mylist <- list(bio3=seq(5.1,6.4,by=0.05),mn.elw=c(effbr,effar)))
leafbio3dat <- emmip(pmod3, mn.elw~bio3,at=mylist, CIs=TRUE, plotit=FALSE)

leafbio3dat$fmnelw <- factor(leafbio3dat$mn.elw)
levels(leafbio3dat) <- c("small","large")

#Plot results of models: (Figure 6)
library(ggplot2)
p1=ggplot(data=leafbio3dat, aes(x=bio3, y=yvar)) +
  geom_line(data=leafbio3dat, aes(x=bio3, y=yvar, linetype=fmnelw), alpha=0.75) + 
  geom_ribbon(aes(ymax=UCL, ymin=LCL, fill=fmnelw), show.legend = FALSE, alpha=0.4)+ 
  scale_fill_manual(values=c("dark gray", "dark gray"))+
  scale_color_manual(values=c("black", "black"))+
  theme_classic()+
  coord_cartesian(ylim = c(0, 40), xlim = c(5.2,6.3)) +
  scale_linetype_manual(values=c("solid", "dotdash"))+
  labs(linetype = "leaf width (cm)")+
  xlab("Isothermality") + labs(y=expression(P[opt]~"("~mu~"mol m"^-2~"s"^-1~")"))+
  theme(legend.position = c(.85, 0.9))


rbPal=colorRampPalette(c('green','black'))
datsc$Col=rbPal(11)[as.numeric(cut(datsc$mean.Tomeg, breaks = 11))] #wue_to, #mean.gwsopt
p2= p1+geom_point(data=datsc, mapping=aes(bio3, mean.Popt, color=bio5), alpha=0.75, size=datsc$mn.elw) +
  scale_color_gradient(low="blue", high="red", name = expression("BIO5 ("~degree~"C)") )+ 
  guides(color = guide_colourbar(order = 1))+
  theme(legend.position = c(.7, 0.80), legend.box = 'horizontal')+
  annotate('text', x=5.25, y=40, label="B", fontface="bold", hjust=2)

#Bio3
(mylist <- list(bio3=seq(5.1,6.4,by=0.05),mn.elw=c(eff)))
leafbio3dat <- emmip(pmod3, mn.elw~bio3,at=mylist, CIs=TRUE, plotit=FALSE)
leafbio3dat$fmnelw <- factor(leafbio3dat$mn.elw)
library(ggplot2)
p3=ggplot(data=leafbio3dat, aes(x=bio3, y=yvar)) +
  geom_line(data=leafbio3dat, aes(x=bio3, y=yvar), alpha=0.5) + 
  geom_ribbon(aes(ymax=UCL, ymin=LCL), show.legend = FALSE, alpha=0.4)+ 
  scale_fill_manual(values=c("dark gray" ))+
  scale_color_manual(values=c("black"))+
  theme_classic()+
  coord_cartesian(ylim = c(0, 40), xlim = c(5.2,6.3)) +
  scale_linetype_manual(values=c("solid"))+
  labs(linetype = "leaf width (cm)")+
  xlab("Isothermality") + labs(y=expression(P[opt]~"("~mu~"mol m"^-2~"s"^-1~")"))+
  theme(legend.position = c(.85, 0.9))

rbPal=colorRampPalette(c('blue','red'))
datsc$Col=rbPal(11)[as.numeric(cut(datsc$bio5,breaks = 11))] #wue_to, #mean.gwsopt
p4=p3+geom_errorbar(aes(ymin=mean.Popt.lci, ymax=mean.Popt.uci, x=bio5),  alpha=0.25, data=datsc, inherit.aes = FALSE) +
  geom_errorbar(aes(xmin=bio3.lci, xmax=bio3.uci, y=mean.Popt),  alpha=0.25, data=datsc, inherit.aes = FALSE)+
  geom_point(data=datsc, mapping=aes(bio3, mean.Popt), alpha=0.5, size=datsc$mn.elw) +
  geom_point(aes(x=6, y=38), size=0.72, color="black", alpha=0.75)+
  annotate('text', x=6.1, y=38, label="S. helleri", fontface="italic", hjust=2) +
  annotate('text', x=6.1, y=38, label="0.72 cm", hjust=0) +
  geom_point(aes(x=6, y=36), size=8.6, color="black", alpha=0.75)+
  annotate('text', x=6.1, y=36, label="S. cordata", fontface="italic", hjust=2)+
  annotate('text', x=6.1, y=36, label="8.62 cm", hjust=0)+
  annotate('text', x=6, y=40, label="Leaf width:", fontface="bold", hjust=0.5)+
  annotate('text', x=5.25, y=40, label="A", fontface="bold", hjust=2)




#pdf("Fig6.Isothermality.v.Popt.pdf", width = 8.75, height =4.25 )
cowplot::plot_grid(p4, p2) 
#dev.off()


