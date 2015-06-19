
export PATH=$PATH:~/Escritorio/bin/

awk '$3=="Sardinian"||$3=="Basque"||$3=="Italian"||$3=="Adygei"||$3=="Orcadian"||$3=="French"\
     {print $0}' H938.clst.txt > Euro.clst.txt
plink --bfile H938 --keep Euro.clst.txt --make-bed --out H938_Euro

plink --bfile H938 --indep-pairwise 50 10 0.1
plink --bfile H938 --extract plink.prune.in --make-bed --out H938.LDprune

plink --bfile H938_Euro --indep-pairwise 50 10 0.1
plink --bfile H938_Euro --extract plink.prune.in --make-bed --out H938_Euro.LDprune

plink --bfile ../data/H938 --hardy --chr 2 --out H938

grep ALL H938.hwe > H938.hwe_reduced

../bin/admixture ../data/H938_Euro.LDprune.bed 6

cat > H938_Euro.LDprune.par
genotypename: ../data/H938_Euro.LDprune.bed
snpname: ../data/H938_Euro.LDprune.bim
indivname: ../data/H938_Euro.LDprune.PCA.fam
snpweightoutname: ./H938_Euro.LDprune.snpeigs
evecoutname: ./H938_Euro.LDprune.eigs
evaloutname: ./H938_Euro.LDprune.eval
phylipoutname: ./H938_Euro.LDprune.fst
numoutevec: 20
numoutlieriter: 0
outlieroutname: ./H938_Euro.LDprune.out
altnormstyle: NO
missingmode: NO
nsnpldregress: 0
noxdata: YES
nomalexhet: YES

# Deal with pesky smartpca ignore issue by creating new .fam file
awk '{print $1,$2,$3,$4,$5,1}' ../data/H938_Euro.LDprune.fam > ../data/H938_Euro.LDprune.PCA.fam

smartpca -p H938_Euro.LDprune.par

spa --bfile ../data/H938_Euro.LDprune --location-output
H938_Euro.LDprune.loc --model-output H938_Euro.LDprune.model

grep Italian ../data/Euro.clst.txt > Italian.clst.txt
grep Orcadian ../data/Euro.clst.txt > Orcadian.clst.txt

../bin/plink --bfile ../data/H938_Euro --keep Italian.clst.txt --hardy --chr2 --out H938_Italian
../bin/plink --bfile ../data/H938_Euro --keep Orcadian.clst.txt --hardy --chr 2 --out H938_Orcadian

grep ALL H938_Italian.hwe > H938_Italian.hwe_reduced
grep ALL H938_Orcadian.hwe > H938_Orcadian.hwe_reduced

# Build files for each grouping
awk '$3=="Sardinian"||$3=="Basque"||$3=="Italian"{print $0}' ../data/Euro.clst.txt > SardBasqItal.clst.txt
awk '$3=="Orcadian"||$3=="French"{print $0}' ../data/Euro.clst.txt > OrcadianFrench.clst.txt

# Run through plink
../bin/plink --bfile ../data/H938_Euro --keep SardBasqItal.clst.txt --hardy --chr 2 --out H938_SardBasqItal
../bin/plink --bfile ../data/H938_Euro --keep OrcadianFrench.clst.txt --hardy --chr 2 --out H938_OrcadianFrench

# Deal with the pesky --hardy output
grep ALL H938_SardBasqItal.hwe > H938_SardBasqItal.hwe_reduced
grep ALL H938_OrcadianFrench.hwe > H938_OrcadianFrench.hwe_reduced

# It's worth having a file of all the chr 2 SNP positions for plotting later
awk '$1==2{print $0}' ../data/H938_Euro.bim > H938_Euro_chr2.pos

../bin/plink --bfile ../data/H938_Euro --pheno .../data/pheno.sim.txt --assoc --out H938_Euro_sim.pheno
