
library(nonlinearTseries)

#RQA Official trajectory matrix

windows()
rqa.analysis3=rqa(as.matrix(official_trajectory_matrix), radius=0.2 ,do.plot=TRUE)

#################################RQA-GROUPS##################
windows()
rqa.analysis1=rqa(clus1, radius=0.02 ,do.plot=TRUE)

windows()
rqa.analysis2=rqa(clus2, radius=0.02 ,do.plot=TRUE)

windows()
rqa.analysis3=rqa(clus3, radius=0.02 ,do.plot=TRUE)

windows()
rqa.analysis4=rqa(clus4, radius=0.02 ,do.plot=TRUE)

windows()
rqa.analysis5=rqa(clus5, radius=0.02 ,do.plot=TRUE)


#########################################

windows()
rqa.analysis11=rqa(clus1, radius=0.01 ,do.plot=TRUE)

windows()
rqa.analysis22=rqa(clus2, radius=0.01 ,do.plot=TRUE)

windows()
rqa.analysis33=rqa(clus3, radius=0.01 ,do.plot=TRUE)

windows()
rqa.analysis44=rqa(clus4, radius=0.01 ,do.plot=TRUE)

windows()
rqa.analysis55=rqa(clus5, radius=0.01 ,do.plot=TRUE)

