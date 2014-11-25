library(lavaan)
library(ggm)

data<-read.csv("../R_scripts/Harman.csv")

data<-data[-18:-24,]
attach(data)
data<-data.frame(Visual.Perception,Cubes,Paper.Form.Board,Flags,General.Information,Paragraph.Comprehension,Sentence.Completion,
                 Word.Classification,Word.Meaning,Addition,Code,Counting.Dots,Straight.Curved.Capitals,Word.Recognition,
                 Number.Recognition,Figure.Recognition,Object.Number)
detach(data)
rownames(data)<-colnames(data)


model1<-'
Spatial.Relations =~ Visual.Perception + Cubes + Paper.Form.Board + Flags
Verbal            =~ General.Information + Paragraph.Comprehension + Sentence.Completion + Word.Classification + Word.Meaning
Perceptual.Speed  =~ Addition + Code + Counting.Dots + Straight.Curved.Capitals
Recognition       =~ Word.Recognition + Number.Recognition + Figure.Recognition + Object.Number
'
fit1<-cfa(model1, std.lv=T, sample.cov = as.matrix(data),sample.nobs = 145)
summary(fit1)

semPaths(fit1, "std", 
         style="lisrel", #残差間相関表示分かるようにlisrelスタイル
         fade=F,  #色が薄くなるのをやめる
         gray=T,  #モノクロ指定
         sizeLat =4,
         sizeMan =2,
         sizeInt =1)

amat<-UG(~Spatial.Relations*Verbal*Perceptual.Speed*Recognition)
fcor<-matrix(c(1.000, 0.561, 0.524, 0.544,
               0.561, 1.000, 0.458, 0.499,
               0.524, 0.458, 1.000, 0.527,
               0.544, 0.499, 0.527, 1.000), nrow=4)
rownames(fcor)<-colnames(fcor)<-rownames(amat)
parcormat<- (-1)*solve(fcor) /   (sqrt(diag(solve(fcor)) %*% t(diag(solve(fcor)))))
diag(parcormat)<- 1
parcormat

amat[2,3]<-amat[3,2]<-0
fitConGraph(amat,fcor,145)

model2<-'
Spatial.Relations =~ Visual.Perception + Cubes + Paper.Form.Board + Flags
Verbal            =~ General.Information + Paragraph.Comprehension + Sentence.Completion + Word.Classification + Word.Meaning
Perceptual.Speed  =~ Addition + Code + Counting.Dots + Straight.Curved.Capitals
Recognition       =~ Word.Recognition + Number.Recognition + Figure.Recognition + Object.Number
Verbal ~~ 0.3605668 * Perceptual.Speed
'
fit2<-cfa(model2, std.lv=T, sample.cov = as.matrix(data),sample.nobs = 145)
summary(fit2)
fitmeasures(fit2,fit.measures = "rmsea")
fitmeasures(fit2,fit.measures = "aic")

fcor<-matrix(c(1.000, 0.537, 0.485, 0.531,
               0.537, 1.000, 0.362, 0.472,
               0.485, 0.362, 1.000, 0.496,
               0.531, 0.472, 0.496, 1.000), nrow=4)
rownames(fcor)<-colnames(fcor)<-rownames(amat)
parcormat<- (-1)*solve(fcor) /   (sqrt(diag(solve(fcor)) %*% t(diag(solve(fcor)))))
diag(parcormat)<- 1
parcormat

amat[1,4]<-amat[4,1]<-0
fitConGraph(amat,fcor,145)

model3<-'
Spatial.Relations =~ Visual.Perception + Cubes + Paper.Form.Board + Flags
Verbal            =~ General.Information + Paragraph.Comprehension + Sentence.Completion + Word.Classification + Word.Meaning
Perceptual.Speed  =~ Addition + Code + Counting.Dots + Straight.Curved.Capitals
Recognition       =~ Word.Recognition + Number.Recognition + Figure.Recognition + Object.Number
Verbal ~~ 0.3632483 * Perceptual.Speed
Spatial.Relations ~~ 0.3625321 * Recognition
'
fit3<-cfa(model3, std.lv=T, sample.cov = as.matrix(data),sample.nobs = 145)
summary(fit3)
fitmeasures(fit3,fit.measures = "rmsea")
fitmeasures(fit3,fit.measures = "aic")



require(OpenMx)
data<-read.csv("../R_scripts/Harman.csv")

data<-data[-18:-24,]
attach(data)
data<-data.frame(Visual.Perception,Cubes,Paper.Form.Board,Flags,General.Information,Paragraph.Comprehension,Sentence.Completion,
                 Word.Classification,Word.Meaning,Addition,Code,Counting.Dots,Straight.Curved.Capitals,Word.Recognition,
                 Number.Recognition,Figure.Recognition,Object.Number)
detach(data)
rownames(data)<-colnames(data)
data<-as.matrix(data)
rownames(data)<-colnames(data)<-c('VisualPerception','Cubes','PaperFormBoard','Flags',
                                  'GeneralInformation','ParagraphComprehension','SentenceCompletion',
                                  'WordClassification','WordMeaning','Addition','Code','CountingDots',
                                  'StraightCurvedCapitals','WordRecognition',
                                  'NumberRecognition','FigureRecognition','ObjectNumber')

testDataCor<-mxData(observed=as.matrix(data),type="cor",numObs=145)

latents = c("SpatialRelations","Verbal","PerceptualSpeed","Recognition")
manifests = colnames(data) 
numSubjects = 145
cfa1<-mxModel("Model1",type="RAM",manifestVars=manifests, latentVars=latents,
              mxPath(from=manifests, arrows=2,free=T, values=1,labels=paste("error",1:17,sep="")),
              mxPath(from=latents,arrows=2,free=F,values=1,labels=c("varF1","varF2","varF3","varF4")),
              mxPath(from="SpatialRelations", to="Verbal", arrows=2, free=T,          values=1, labels="cov1"),
              mxPath(from="SpatialRelations", to="PerceptualSpeed", arrows=2, free=T, values=1, labels="cov2"),
              mxPath(from="SpatialRelations", to="Recognition", arrows=2, free=T,     values=1, labels="cov3"),
              mxPath(from="Verbal"          , to="PerceptualSpeed", arrows=2, free=T, values=1, labels="cov4"),
              mxPath(from="Verbal"          , to="Recognition", arrows=2, free=T,     values=1, labels="cov5"),
              mxPath(from="PerceptualSpeed" , to="Recognition", arrows=2, free=T,     values=1, labels="cov6"),
              mxPath(from="SpatialRelations", to= c("VisualPerception","Cubes","PaperFormBoard","Flags"), arrows=1,free=T,values=1,labels=paste("l",1:4,sep="")),
              mxPath(from="Verbal"          , to= c("GeneralInformation","ParagraphComprehension","SentenceCompletion","WordClassification","WordMeaning"), arrows=1,free=T,values=1,labels=paste("l",5:9,sep="")),
              mxPath(from="PerceptualSpeed" , to= c("Addition","Code","CountingDots","StraightCurvedCapitals"), arrows=1,free=T,values=1,labels=paste("l",10:13,sep="")),
              mxPath(from="Recognition"     , to= c("WordRecognition","NumberRecognition","FigureRecognition","ObjectNumber"), arrows=1,free=T,values=1,labels=paste("l",14:17,sep="")),
              mxData(observed=data,type="cor",numObs=numSubjects)
              )
cfa1<-mxRun(cfa1)
summary(cfa1)

cfa2<-mxModel("Model1",type="RAM",manifestVars=manifests, latentVars=latents,
              mxPath(from=manifests, arrows=2,free=T, values=1,labels=paste("error",1:17,sep="")),
              mxPath(from=latents,arrows=2,free=F,values=1,labels=c("varF1","varF2","varF3","varF4")),
             mxPath(from="SpatialRelations", to= c("VisualPerception","Cubes","PaperFormBoard","Flags"), arrows=1,free=T,values=1,labels=paste("l",1:4,sep="")),
              mxPath(from="Verbal"          , to= c("GeneralInformation","ParagraphComprehension","SentenceCompletion","WordClassification","WordMeaning"), arrows=1,free=T,values=1,labels=paste("l",5:9,sep="")),
              mxPath(from="PerceptualSpeed" , to= c("Addition","Code","CountingDots","StraightCurvedCapitals"), arrows=1,free=T,values=1,labels=paste("l",10:13,sep="")),
              mxPath(from="Recognition"     , to= c("WordRecognition","NumberRecognition","FigureRecognition","ObjectNumber"), arrows=1,free=T,values=1,labels=paste("l",14:17,sep="")),
             mxMatrix(type = "Symm", nrow = 4, ncol = 4,
                      free = c(F, T, T, T, T, F, T, T, T, T, F, T, T, T, T, F),
                      values = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                      name = "cov"),
             mxData(observed=data,type="cor",numObs=numSubjects)
)



cfa2<-mxRun(cfa2)
summary(cfa2)

semPaths(factorModel, "std", 
         style="lisrel", #残差間相関表示分かるようにlisrelスタイル
         #         layout="circle",
         fade=F,  #色が薄くなるのをやめる
         gray=T,  #モノクロ指定
         sizeLat =4,
         sizeMan =2,
         sizeInt =1)

fcor<-matrix(c(0,0,0,0,
               0.5613780,0,0,0,
               0.5239843,0.4584148,0,0,
               0.5435321,0.4994483,0.5273366,0),nrow=4)
fcor<-fcor+t(fcor)
diag(fcor)<-1
fcor
solve(fcor)


require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents     <- c("G")
factorModel1 <- mxModel("One Factor", type="RAM",
                       manifestVars = manifests,
                       latentVars = latents,
                       mxPath(from=latents, to=manifests),
                       mxPath(from=manifests, arrows=2),
                       mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
                       mxData(cov(demoOneFactor), type="cov", numObs=500))
summary(mxRun(factorModel1))

data(demoOneFactor)
factorModel2 <- mxModel("One Factor",
                       mxMatrix("Full", nrow=5, ncol=1, values=0.2, free=TRUE, name="A"),
                       mxMatrix("Symm", nrow=1, ncol=1, values=1, free=FALSE, name="L"),
                       mxMatrix("Diag", nrow=5, ncol=5, values=1, free=TRUE, name="U"),
                       mxAlgebra(A %*% L %*% t(A) + U, name="R"),
                       mxMLObjective("R", dimnames = names(demoOneFactor)),
                       mxData(cov(demoOneFactor), type="cov", numObs=500))
summary(mxRun(factorModel2))
