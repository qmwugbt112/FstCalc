
# Illustration of method to get the Maximum likelihood estimate of Fst 
# given a list of local allele counts for each locus (alist) 
# & a list of global allele frequencies for each locus (plist)

########################################################################################
# Function to return log( P(α | p,Fst) ), the log multinomial dirichlet calculation of the likelihood of 
# 	allelic composition α, 
#	given Fst f, 
#	and allele frequencies p
# n.b. this log version deals well with small values of f and p
# Warning: p values must be non-zero
########################################################################################

lmd<-function(a=c(2,8),f=0.1,p=c(0.3,0.7)){	#n.b. default values illustrating locus with two alleles 2 alleles
	l<-1/f-1;n=sum(a);x<-l*p
	return(lgamma(l)+lgamma(n+1)-lgamma(n+l)	
			+sum(lgamma(x+a))	
			-sum(lgamma(a+1))
			-sum(lgamma(x))
			)	
	}	

# Vectorize the function so it can act over a vector different f values
lmdL<-Vectorize(lmd,'f')

########################################################################################



########################################################################################
# Function to return sum(log( P(α | p,Fst) )), for
#	al: a list of allelic counts for different loci (α)
#	ap: a list of allele frequencies (p)
#	f:	a vector of Fst values 	(f)
# Warning p-values should be non-zero
########################################################################################
combinedL<-function(al,pl,f){
	L<-rep(0,length(f))
	for (i in length(al)){
		L<-L+lmdL(al[[i]],f,pl[[i]])	
		}
	return(L)
	}

########################################################################################
#						Main Block										
########################################################################################

# Simple example data with very large deviations from expected frequencies
Locus1Counts<-c(1,4,30); Locus2Counts<-c(20,3,2,7);Locus3Counts<-c(30,2)
alist<-list(Locus1Counts,Locus2Counts,Locus3Counts)

Locus1p<-c(.3,.6,.1);Locus2p<-c(.2,.2,.2,.4);Locus3p<-c(.5,.5)
plist<-list(Locus1p,Locus2p,Locus3p)


# Obtaining the likelihood curve 
xvals<-1:999/1000

Lcurve<-combinedL(alist,f=xvals,plist)

# Get ML value
peakV<-max(Lcurve,na.rm=T)
MLv<-xvals[which(Lcurve==peakV)]

# Get 2 unit support limits
inside<-which(Lcurve>=(peakV-2))
supportLims<-c(xvals[inside][1],tail(xvals[inside],1))


plot(	xvals,
		Lcurve,
		type='l',
		xlab='Fst Value',
		ylab='Log likelihood',
		main=paste(	'Likelihood curve, ML value (',
					round(MLv,2),
					'*), & support limits (red)'
				)
	)


lines(	rep(supportLims[1],2),
		c(min(Lcurve),peakV-2),
		col='red')
lines(	rep(supportLims[2],2),
		c(min(Lcurve),peakV-2),
		col='red')		
lines(supportLims,rep(peakV-2,2),col='red')

points(	MLv,
		peakV,
		pch='*',
		cex=2,
		col='blue')


