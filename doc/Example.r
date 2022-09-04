options(width=2000)
Stval<-0.85
ModelType=3;
cv_num=10; 

################ Loading the necessary codes  ##################
load("HeterogeneityCodes.RData")

#################################################################
################# Create the output dir #########################
##### All groups information will stored in this folder #########
#################################################################
outdir=getwd();
dir.create("MLRoutput",showWarnings=FALSE)
if(file.exists(paste(outdir,'groups1.txt',sep='/'))) file.remove(paste(outdir,'groups1.txt',sep='/'))
if(file.exists(paste(outdir,'groups2.txt',sep='/'))) file.remove(paste(outdir,'groups2.txt',sep='/'))
if(file.exists(paste(outdir,'groups3.txt',sep='/'))) file.remove(paste(outdir,'groups3.txt',sep='/'))

############### Simulate Disease #########################################
#### Dis1, Dis2 shared teh same mechanisms ###############################
#### Dis3-6 shared the same mechanisms ###################################
#### Disease model: 2-locus threshold interaction + 1 additive effect ####
##########################################################################


source("Disease_model_demo.r")
Dtpop=Disease_model_demo1()

NTrait=6
size1=size2=size3=size4=size5=size6=sizecont=300
size1=size5=150;
GROUPS=1:NTrait

RESULT=NULL;

NTotal=1000
for(rep in 1:NTotal)
{
	cat("REP=",rep,'\n')
	Dtnew=SampleSim(Dtpop,NTrait,c(size1,size2,size3,size4,size5,size6,sizecont))
	Dt=Dtnew$Dt
	Dt_test=Dtnew$Dt_test;
	samples=getsample(Dt,NTrait,Dt_test);
	aa=CROSS(samples$Training)

	DTSCASEss=DTSCASE=samples$Training$DTSCASE;
	DTSCASEsstest=DTSCASEtest=samples$Testing$DTSCASEtest;

	AUCcross1=aa$AUCcross1;
	AUCcross2=aa$AUCcross2;
	AUCcross3=aa$AUCcross3;

	t1=apply(AUCcross1,2,mean)
	groups1=which(t1==max(t1))[1];
	t1=apply(AUCcross2,2,mean)
	groups2=which(t1==max(t1))[1];
	t1=apply(AUCcross3,2,mean)
	groups3=which(t1==max(t1))[1];

	print(paste("Groups1=",groups1,"Groups2=",groups2,"Groups3=",groups3));

	total_group=NTrait;
	MODELsave=list()
	combine=list()
	GROUPStemp=GROUPS;
	bb=getalltestauc(groups1,groups2,groups3,samples$Training, samples$Testing,NTrait)
	RESULT=rbind(RESULT,bb$AUCs)
}

