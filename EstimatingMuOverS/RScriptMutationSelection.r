
#pdf("Sims_estimating_mu_over_s_DIFF_s_values.pdf")
setwd("~/Dropbox/SimulationsEstimatingSelCoeffSims/EstimatingMuOverS")
system("./Code_and_shellscript/make_HIV1site")    #compile the code 

#DataOverview<-data.frame("N"=Ne,"mu"=0,"cost"=0,"num_runs"=0,"datapointsperrun"=0,"equiPi"=0,"VarPi"=0,"expectedPi"=0,"t_half"=0)

Ne = 10000 #currently Ne cannot be changed in the sims 
nmax = 200 #(max num of patients, should be around 200)
NUMRUNS=1; numoutputs=nmax
theta = 10
mu = theta / Ne
numsites = 300
costs = sort(round(runif(numsites),2))
seed =1

for (site in 1:numsites){
print(paste("site",site))
cost = costs[site]
#make script 
x<-"#!/bin/bash"
x<-c(x,paste("mu=",mu,sep=""))
outputfrequency=min(c(2*Ne,ceiling(5/cost)))
#outputfrequency=2*Ne
x<-c(x,paste("output_every_Xgen=",outputfrequency,sep=""))
x<-c(x,paste("numgen_inN=",(numoutputs+2)*outputfrequency/Ne,sep=""))
x<-c(x,paste("start_output=",2*outputfrequency/Ne,sep=""))
x<-c(x,paste("cost=",cost,sep=""))
#x<-c(x,paste("for cost in", paste(costs,collapse = " ")))
#x<-c(x,"do")
x<-c(x,paste("for seed in",seed))
x<-c(x,"do",
     "echo \"", "$seed", "$mu", "$cost",
     "$output_every_Xgen", "$numgen_inN", "$start_output",
     paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site >../Data/Data_s_", site, "_T_", theta, "_cost_", cost,".txt",sep=""), 
     "done")
#x<-c(x,"done")
write(x,file="./Code_and_shellscript/tempscript.sh")
system("chmod 775 ./Code_and_shellscript/tempscript.sh")
#Run tempscript.sh
system("./Code_and_shellscript/tempscript.sh")
}

##########Read data
ListFiles<-list.files("../Data/")
ListCosts<-rep(0,length(ListFiles))
ListMeanFreqs<-rep(0,length(ListFiles))
for (i in 1:length(ListFiles)){
    ListCosts[i]<-as.numeric(substr(ListFiles[i], regexpr("cost",ListFiles[i])[1]+5, regexpr(".txt",ListFiles[i])[1]+-1))
    #read the file
    ListMeanFreqs[i]<-mean(read.csv(paste("../Data/",ListFiles[i],sep=""),sep = "\t")$freq)
}
plot(ListCosts,mu/ListMeanFreqs)
print(cor.test(ListCosts,mu/ListMeanFreqs))
#> plot(ListCosts,mu/ListMeanFreqs,ylim=c(0,1))
#> plot(ListCosts,mu/ListMeanFreqs,ylim=c(0,2))
#> plot(ListCosts,ListMeanFreqs,ylim=c(0,1))
#> plot(ListCosts,ListMeanFreqs,log="y")

###Hm, it looks like the very deleterious mutations have way too low frequencies. Fixed. 
###I need to check how i did the simulations SWAPPED ORDER OF MUTATION AND SELECTION. Now I sample after mutation, but before selection. 
###Also create a github repository. DONE

#####Looking at smaller number of patients. 

CorList<-c()
NumPatsList<-c(1,5,seq(10,200,by=20))   
for (Numpats in NumPatsList){
    ListCosts<-rep(0,length(ListFiles))
    ListMeanFreqs<-rep(0,length(ListFiles))
    Pats<-sample(200,Numpats)
    for (i in 1:length(ListFiles)){
     ListCosts[i]<-as.numeric(substr(ListFiles[i], regexpr("cost",ListFiles[i])[1]+5, regexpr(".txt",ListFiles[i])[1]+-1))
    #read the file
     ListMeanFreqs[i]<-mean(read.csv(paste("../Data/",ListFiles[i],sep=""),sep = "\t")$freq[Pats])
    }
plot(ListCosts,mu/ListMeanFreqs,main=paste("Numpats",Numpats),ylim=c(0,1))
Cor=round(cor.test(ListCosts,mu/ListMeanFreqs)$estimate,3)
text(0.6,0.2,paste("cor = ", Cor))
print(paste("numpats",Numpats,"cor",Cor))
CorList<-c(CorList,Cor)
}

plot(NumPatsList,CorList)

###Next: what if sequencing less precise! 


#####Old code from 2014, I don't think I need this. 
if (FALSE){
#READ FREQS FILE
read.csv(paste("../Data/Link",cost,"_",seed,".txt",sep=""),sep="\t",header=TRUE)->simdata
#in stead of the real frequencies, lets assume we have a sample from each patient of, say, 100, seqs. 
plot(c(0,0),col=0,xlim=c(-0.5,3.5),ylim=c(-1*mu/cost,8*mu/cost),xlab="Num Patients",xaxt="n",ylab="95% of observed average freq of allele",main=paste("Ne",Ne,", mu",mu,", Theta",2*Ne*mu ,", cost",cost))
axis(1, at=log10(samplesizes), labels=samplesizes)
abline(h=mu/cost,lty=2)

diffnumseqs <- c(50,100,200)
for (num_seqs_per_patient in diffnumseqs){
    for (i in 1:length(simdata$freq)){
        simdata$est_freq[i]<-rbinom(1,num_seqs_per_patient,simdata$freq[i])/num_seqs_per_patient}
    
    #system(paste("rm ","../Data/Link",cost,"_",seed,".txt",sep=""))
    samplesizes<-c(1,3,10,30,100,300,1000)
    for (num_patients in samplesizes){
        list_averages<-vector()
        for (i in 1:1000){ 
            list_averages<-c(list_averages,mean(sample(simdata$est_freq,num_patients)))}
        co=which(diffnumseqs==num_seqs_per_patient)
        X=(which(diffnumseqs==num_seqs_per_patient)-2)*0.1
        print(paste(num_seqs_per_patient,num_patients,log10(num_patients)-0.02+X))
        rect(log10(num_patients)-0.02+X,sort(list_averages)[25],log10(num_patients)+0.02 +X,sort(list_averages)[975],col=co)
        #	text(log10(num_patients)+X,sort(list_averages)[975]+0.15*mu/cost+(X+0.2)/200,paste(round(sort(list_averages)[25]/(mu/cost),2),"-",round(sort(list_averages)[975]/(mu/cost),2)),cex=0.6)
        
        text(log10(num_patients)+X,(-1+co/8)*mu/cost,paste(round(sort(list_averages)[25]/(mu/cost),2),"-",round(sort(list_averages)[975]/(mu/cost),2)),cex=0.5)	
    }
    text(log10(num_patients),6.5*mu/cost,"bars and",cex=0.8)
    text(log10(num_patients),6*mu/cost,"numbers indicate range",cex=0.8)
    text(log10(num_patients),5.5*mu/cost,"of 95% of estimates",cex=0.8)
    text(log10(num_patients),5.*mu/cost,"black:50, red: 100, ",cex=0.8)
    text(log10(num_patients),4.5*mu/cost,"green: 200 sequences/pat",cex=0.8)
    text(-0.3,1.1*mu/cost,"mu/s",cex=0.8)
    
}

dev.off()
}
#rbinom(n, size, prob)

