
plot.geno.vs.HW<-function(file,title=""){


    #read in the HW file from plink
        plink.hwe<-read.table(file,as.is=TRUE)

    names(plink.hwe)<-c("chr","SNP.id","which.inds","a1","a2","genotype","obs.het","exp.het","HWE.pval")

        counts<-sapply(plink.hwe$genotype,function(x){as.numeric(strsplit(x,"/")[[1]])})
        counts<-t(counts)
        tot.counts<-rowSums(counts)
        geno.freq<-counts/tot.counts
        allele.freq<-(geno.freq[,1]+.5*geno.freq[,2])

        these.minor<-sample(1:nrow(geno.freq),3000)
        these.major<-sample(1:nrow(geno.freq),3000)
        ss.allele<-c(allele.freq[these.minor],1-allele.freq[these.major])
        ss.geno<-rbind(geno.freq[these.minor,],geno.freq[these.major,c(3,2,1)])

        plot(ss.allele,ss.geno[,1],xlim=c(0,1),ylim=c(0,1),col=adjustcolor("red",0.1),
           xlab="allele frequency",ylab="genotype frequency",main=title)
        points(ss.allele,ss.geno[,3],xlim=c(0,1),ylim=c(0,1),col=adjustcolor("blue",0.1))
        points(ss.allele,ss.geno[,2],xlim=c(0,1),ylim=c(0,1),col=adjustcolor("green",0.1))
        smooth=1/5
        lines(lowess(ss.geno[,1]~ss.allele,f = smooth),col="black")
        lines(lowess(ss.geno[,3]~ss.allele,f = smooth),col="black")
        lines(lowess(ss.geno[,2]~ss.allele,f = smooth),col="black")

        x=1:1000/1000
        lines(x,x^2,lty=2)
        lines(x,2*x*(1-x),lty=2)
        lines(x,(1-x)^2,lty=2)
        legend(x=0.3,y=1,col=c("red","blue","green",rep("black",2)),
           legend=c("Homozygote AA","Homozygote aa","Heterozygote
        Aa","Mean","Hardy Weinberg Expectation"),pch=c(rep(1,3),rep(NA,2)),lty=c(rep(NA,3),1,2))
}

png(file="HGDP_HWE.png")
plot.geno.vs.HW(file="H938.hwe_reduced",title="HGDP")
dev.off()

file="H938.hwe_reduced"

# Read in the HWE table and compute counts and allele frequencies
hwe<-read.table(file,as.is=TRUE)
names(hwe)<-c("chr","SNP.id","which.inds","a1","a2","genotype","obs.het","exp.het","HWE.pval")
counts<-sapply(hwe$genotype,function(x){as.numeric(strsplit(x,"/")[[1]])})
counts<-t(counts)
tot.counts<-rowSums(counts)
allele.counts<-(2*counts[,1]+counts[,2])

# Flip allele counts so that we are sure we always have the minor
# allele frequency
# (Note: this uses a trick based on boolean math where true/false = 1/0).
counts.maf = allele.counts*(allele.counts<2*tot.counts-allele.counts)
+(2*tot.counts-allele.counts)*(allele.counts<2*tot.counts-allele.counts)

# Set the number of individuals by looking at the sites w/ the most
# observed data
n=max(tot.counts)

# Make the plot but filter on using only sites with fully observed
# data (i.e. totcounts==n)
hist(counts.maf[tot.counts==n],xlab="Minor allele count",
ylab="# of SNPs",main="Allele frequency spectra",breaks=n)

# Plot the expected minor allele frequency spectra for the standard
# neutral model (i.e. constant size population, all loci neutral)
# To do so we compute, Watterson's estimator of Theta
S=sum(tot.counts==n & counts.maf>0)
thetaW=S/sum(1/seq(1,2*n-1))
# Which determines the expected AFS
expectedAFS=(1/seq(1,n)+1/(n*2-seq(1,n))) * thetaW
# And then plot
lines(seq(1,n),expectedAFS,col=2)
# Note: This adds a red line displaying the expected AFS shape
# controlled to match the data w.r.t to Watterson's Theta (i.e. the total number of SNPs).

# Read in matrix of inferred ancestry coefficients for each individual.
Q=read.table("H938_Euro.LDprune.6.Q")
Qmat=as.matrix(Q)
barplot(t(Qmat),col=c("red","blue","gold","orange","purple","brown"),border=NA,space=0)

# To be able to label the graph we read in a
.clst file with population "cluster" labels for each indiv
clst=read.table("../data/Euro.clst.txt")
# And a fam file from the plink data
fam=read.table("../data/H938_Euro.LDprune.fam")

# Use the match function to link the family ids with the cluster id
clst_unord=clst$V3[match(fam$V2,clst$V2)]
# Re-order alphabetically
ordered_indices=order(clst_unord)
QmatO=Qmat[ordered_indices,]

# Compute where we will place the population labels in the barplot
n=length(ordered_indices)
clst_ord=clst_unord[ordered_indices]
breaks=c(0,which(clst_ord[1:(n-1)]!=clst_ord[2:n]),n)
nbrks=length(breaks)
midpts=(breaks[1:(nbrks-1)]+breaks[2:nbrks])/2

# Make the barplot
barplot(t(QmatO),col=c("red","blue","yellow","orange","purple","brown"),border=NA,space=0,inside=TRUE)
abline(v=breaks,lwd=2)
mtext(levels(clst_ord),side=1,at=midpts,las=2)

# Read in eigenvectors file
PCA=read.table("H938_Euro.LDprune.eigs")
names(PCA)=c("ID",paste("PC",(1:20),sep=""),"CaseControl")

# Note smartpca pushes the plink family and individual ids together so
# we need to extract out the ids afresh; note this code works just for
# this case
ids=substr(PCA$ID,start=6,stop=20)

# Read in clst table and fam file
clst=read.table("../../data/Euro.clst.txt")
# The list of countries as ordered in the fam file
clst_unord=clst$V3[match(ids,clst$V2)]

# Make a blank plot of the right size
plot(PCA$PC2,PCA$PC1,type="n",xlab="PC2",ylab="PC1")
# Add text labels at PC positions with abbreviations of the full
# labels.  The substr function will be used to make automatic abbreviations.
text(PCA$PC2,PCA$PC1,substr(clst_unord,1,2))

eval=read.table("H938_Euro.LDprune.eval")
plot(1:length(eval$V1),eval$V1/sum(eval$V1),xlab="PC",ylab="Eigenvalue")

# Read in location file
loc=read.table("H938_Euro.LDprune.loc")
x=loc$V7
y=loc$V8

# Read in clst table and fam file
clst=read.table("../data/Euro.clst.txt")


# The list of countries as ordered in the fam file
clst_unord=clst$V3[match(loc$V2,clst$V2)]


plot(x,y,type="n")
text(x,y,substr(clst_unord,1,2))

file1="H938_Italian.hwe_reduced"
file2="H938_Orcadian.hwe_reduced"

file1="H938_SardBasqItal.hwe_reduced"
file2="H938_OrcadianFrench.hwe_reduced"

# Extract pop1 MAF
pop1.hwe<-read.table(file1,as.is=TRUE)
names(pop1.hwe)<-c("chr","SNP.id","which.inds","a1","a2","genotype","obs.het","exp.het","HWE.pval")
pop1.counts<-sapply(pop1.hwe$genotype,function(x){as.numeric(strsplit(x,"/")[[1]])})
pop1.counts<-t(pop1.counts)
pop1.tot.counts<-rowSums(pop1.counts)
pop1.geno.freq<-pop1.counts/pop1.tot.counts
pop1.allele.freq<-(pop1.geno.freq[,1]+.5*pop1.geno.freq[,2])

# Extract pop2 MAF
pop2.hwe<-read.table(file2,as.is=TRUE)
names(pop2.hwe)<-c("chr","SNP.id","which.inds","a1","a2","genotype","obs.het","exp.het","HWE.pval")
pop2.counts<-sapply(pop2.hwe$genotype,function(x){as.numeric(strsplit(x,"/")[[1]])})
pop2.counts<-t(pop2.counts)
pop2.tot.counts<-rowSums(pop2.counts)
pop2.geno.freq<-pop2.counts/pop2.tot.counts
pop2.allele.freq<-(pop2.geno.freq[,1]+.5*pop2.geno.freq[,2])

# Note: we need to flip alleles because the minor allele in one pop
# isn't necessarily the minor allele in another.
pop2.allele.freq.flip=
pop2.allele.freq*(pop1.hwe$a1==pop2.hwe$a1)+(1-pop2.allele.freq)*(pop1.hwe$a1!=pop2.hwe$a1)

# Read in SNP positions table
pos=read.table("H938_Euro_chr2.pos")

# Calculate Fst w/ most basic version of F_{ST}: F_{ST} = Var(p) / (pbar*(1-pbar))
pbar=(pop1.allele.freq+pop2.allele.freq)/2
fst=(0.5*(pop1.allele.freq-pbar)^2+0.5*(pop2.allele.freq.flip-pbar)^2)/(pbar*(1-pbar))

# Plot of fst vs. position
plot(pos$V4,fst,col=adjustcolor("black",0.1),pch=16,xlab="Position")

# Identify the most extreme outlier
pos[which.max(fst),]

pop1.hwe[pop1.hwe$SNP.id=="rs10498067",]
pop2.hwe[pop2.hwe$SNP.id=="rs10498067",]

awk '$3=="Adygei"{print $1,$2,5}\
$3=="Basque"{print $1,$2,2}\
$3=="Italian"{print $1,$2,5}\
$3=="Sardinian"{print $1,$2,0}\
$3=="French"{print $1,$2,8}\
$3=="Orcadian"{print $1,$2,10}' \
../data/Euro.clst.txt > pheno.base.txt

pheno.base=read.table("pheno.base.txt")
# Scale the base phenotype to mean 0, sd 1
pheno.base.scale=scale(pheno.base$V3)
# Add some normally distributed noise (with as much variance as the base phenotype itself already has)
pheno.sim=rnorm(length(pheno.base$V3),mean=pheno.base$V3,sd=1)
# Output the phenotype to a file
write.table(cbind(pheno.base[,1:2],pheno.sim),"pheno.sim.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

qassoc=read.table("H938_Euro_sim.pheno.qassoc",header=TRUE)

# Make a Manhattan plot
# First set-up a plot of the right size
plot(1:length(qassoc$SNP),type="n",xlab="SNP index",ylab="-log10(p-value)",ylim=c(3,max(-log10(qassoc$P),na.rm=TRUE)+1))
# Next add the points (note: we only plot points with
# -log10(p-value)>3 to minimze the number of points plotted)
plot.these=which(-log10(qassoc$P)>3)
points(plot.these,-log10(qassoc$P[plot.these]),col=1+qassoc[plot.these,"CHR"]%%2,pch=16)
# Put in a line for a Bonferroni correction (0.05 / length(qassoc$SNP)
abline(h=-log10(0.05/length(qassoc$SNP)),lty=2,col="gray")

print(qassoc[head(order(qassoc$P),n=20),])

# Read in the p-values
qassoc=read.table("H938_Euro_sim.pheno.qassoc",header=TRUE)
# Produce expected p-values from the null (i.e. perfectly uniformly
# distributed).
nTests=length(qassoc$SNP)
Unif=seq(1/nTests,1-1/nTests,length=nTests)
# Sort the -log10 p-values (i.e. match on quantile)
logUnifOrder=order(-log10(Unif),decreasing=TRUE)
SNPorder=order(-log10(qassoc$P),decreasing=TRUE)
# Plot the p-values against against each other (Note: we do for only
# the top 150K SNPs to make the number of points plotted smaller)
qmax=max(-log10(qassoc$P),na.rm=TRUE)
plot(-log10(Unif[logUnifOrder][1:150e3]),-log10(qassoc$P[SNPorder][1:150e3]),pch=16,cex=0.5,xlab="-log(p)
Expected",ylab="-log(p) Observed",,ylim=c(0,qmax));
# put in a line for the expected relationship (y=x)
abline(0,1);
# Put in a line for a Bonferroni correction (0.05 / length(qassoc$SNP)
abline(h=-log10(0.05/length(qassoc$SNP)),lty=2,col="gray")
