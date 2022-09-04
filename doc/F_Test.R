library(fda)

### function for smoothing the data###
smoothing<-function(geno,pos)
{
  nN=nrow(geno)
  pos <- (pos-pos[1])/(pos[length(pos)]-pos[1])
  #pos2 = seq(min(pos), max(pos), length=2*length(pos))
  intrvllen<-0.0001
  pos3<-seq(0,1,intrvllen)
  shouldContinue=TRUE;
  
  if(shouldContinue==TRUE)
  {
    lambda<-rep(0,nN)
    for(j in 1:nN)
    {
      genospline<-try(smooth.spline(pos,geno[j,],nknots=length(pos)))
      if("try-error" %in% class(genospline)){
        lambda[j]<-NA}
      else {lambda[j]<-genospline$lambda}
    }
  }
  
  lambda_all<-mean(lambda,na.rm=TRUE)
  Knots = pos
  norder = 4
  nbasis=length(Knots) + norder - 2
  ybasis=create.bspline.basis(range(Knots), nbasis, norder, Knots)
  Lfdobj= 2
  yfdPar = fdPar(ybasis, Lfdobj, lambda=lambda_all)
  yfd = smooth.basis(pos,t(geno),yfdPar)$fd
  yfd<-center.fd(yfd)
  yij.f = eval.fd(pos3,yfd)
  z1=NULL;
  z1<-rbind(z1,apply(yij.f*intrvllen,MARGIN=2,sum))
  return(z1)
}



### function for calculating the U-Statistic
Mtest<-function(z,Tn,nN)
{
  ubar<-matrix(data=0,nrow=nN,ncol=ncol(Tn))
  A<-NULL;
  
  for(i1 in 1:nN)
  {
    u<-NULL;
    for(j1 in 1:nN) 
    {
      u<-rbind(u,Tn[i1,]-Tn[j1,]);      
    }
    ubar[i1,]<-apply(u,2,mean)
  }
  
  for(i1 in 1:ncol(Tn))
  {
    A<-cbind(A,(2*sqrt(nN)/(nN-1))*( t(z)*ubar[,i1]));
  }
  
  Astat<-apply(A,2,sum)  ############Thiis is sqrt(N)U
  
  ################# Covariance of Astat############
  #correction_fac<-(nN-1)/(nN-1.5)
  va<-as.numeric(var(t(z)))
  cov_Astat<-4*cov(Tn)*va;
  rankA<-rankMatrix(cov_Astat)[1];
  Astattr<-as.matrix(Astat,nrow=ntrait)
  chisqst<-(Astat%*%ginv( cov_Astat)%*%Astattr);
  return(pchisq( chisqst,rankA, ncp = 0, lower.tail = FALSE, log.p = FALSE));
  
}



###data###
## geno is seq data where rows correspond to individuals and columns to SNPs###
## no missing values allowed##
geno=read.table("geno.txt",sep=",") 

###y denotes the triats where rows correspond to individuals and columns to different traits###
## no missing values allowed##
y=read.table("trait.txt",sep=",")

###pos refers to positiond of the SNPs in the geno data###
## no missing values allowed##
pos=as.matrix(read.table("pos.txt",sep=","))
pos=as.vector(pos)
nN=nrow(geno)	
z=smoothing(geno,pos)
pval=Mtest(z,y,nN)
