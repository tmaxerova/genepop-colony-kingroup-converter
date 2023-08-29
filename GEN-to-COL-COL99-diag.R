f<-1 #leave this as it is
m<-0 #leave this as it is

##Parameters you must update:
path_to_input_genepop_file<-"C:\\Users\\Maxerová\\Desktop\\Misa\\experiment.txt"
output_folder_location<-"C:\\Users\\Maxerová\\Desktop\\Misa\\"
sex_vector<-c(f,f,m,m,f) #vector containing information about sex of the individual (from first to last sample), f for female, m for male
n_initial_empty_rows<-1
n_snps<-6
n_rubbish_rows<-1

##Parameters you can update (but there is no real reason for that):
output_txtfile_name<-"colony.txt"
output_csvfile_name<-"colony.csv"
males99_output_txtfile_name<-"colony_males99.txt"
males99_output_csvfile_name<-"colony_males99.csv"
output_sample_hetero_diagnostics_file<-"sample_hetero_proportion.csv"
output_snp_hetero_diagnostics_file<-"snp_hetero_proportion.csv"

##OUTPUT FILE SPECIFICATION.
complete_name<-paste0(output_folder_location,"ResultsFolder","_",substr(Sys.time(),1,10),"_",substr(Sys.time(),12,13),"-",substr(Sys.time(),15,16),"-",substr(Sys.time(),18,19))
dir.create(complete_name)
setwd(complete_name)

#GENEPOP TO COLONY CONVERTER.
nskippedrows<-n_initial_empty_rows+n_snps+n_rubbish_rows
origfile<-read.table(path_to_input_genepop_file,skip=nskippedrows)
snps_names<-t(read.table(path_to_input_genepop_file,skip=n_initial_empty_rows,nrows=n_snps))
snps_names_double<-rep(c(snps_names),each=2)
allele<-rep(c("_al1","_al2"),times=n_snps)
coln<-paste0(snps_names_double,allele)
rawdata<-origfile[,-1:-2]
rdnrow<-nrow(rawdata)
rdncol<-ncol(rawdata)
znrow<-nrow(rawdata)
zncol<-(ncol(rawdata)*2)
zmatrix<-matrix(0,znrow,zncol)
for(j in 1:rdncol)
{
  for(i in 1:rdnrow)
  {
    lok1<-substring(rawdata[i,j],1,3)
    zmatrix[i,((2*j)-1)]<-lok1
    lok2<-substring(rawdata[i,j],4,6)
    zmatrix[i,(2*j)]<-lok2
  }
}
rownames(zmatrix)<-origfile[,1]
colnames(zmatrix)<-coln
finalmatrix<-noquote(zmatrix)
write.table(finalmatrix, file=output_txtfile_name, quote=FALSE,row.names=TRUE,col.names=TRUE)
write.csv(finalmatrix, file=output_csvfile_name, quote=FALSE,row.names=TRUE)

#HETEROZYGOSITY DIAGNOSTICS GENERATOR (both for samples and SNPs).
hapdip_matrix<-matrix(0,rdnrow,rdncol)
for(b in 1:rdncol)
{
  for(a in 1:rdnrow)
  {
    if (substring(rawdata[a,b],1,3) != substring(rawdata[a,b],4,6))
    {
    hapdip_matrix[a,b]<-1
    }
    
  }
}
samp_diag<-numeric(rdnrow)
for (k in 1:rdnrow)
{
  samp_diag[k]<-(sum(hapdip_matrix[k,])/rdncol)
}
final_samp_diag<-rbind(origfile[,1],samp_diag)
snp_diag<-numeric(rdncol)
for (l in 1:rdncol)
{
  snp_diag[l]<-(sum(hapdip_matrix[,l])/rdnrow)
}
final_snp_diag<-rbind(snps_names,snp_diag)
write.table(final_samp_diag, file=output_sample_hetero_diagnostics_file, sep=",", col.names=FALSE)
write.table(final_snp_diag, file=output_snp_hetero_diagnostics_file, sep=",", col.names=FALSE)

#COLONY MALE REWRITER (in males, all second alleles are rewritten to -99).
malematrix<-finalmatrix
even_columns<-seq(from=2,to=zncol,by=2)
males_position<-which(sex_vector == 0)
if (length(sex_vector) != rdnrow){
  print("ERROR MESSAGE: Number of elements entered to the sex vector does not match number of samples in the input file.")
}else{
  for (s in 1:length(sex_vector)){
    if (sex_vector[s] == 0){
      for (t in even_columns)
      {
        for (u in males_position)
        {
         malematrix[u,t]<-(-99)
        }
      }
    }
  }
}
write.table(malematrix, file=males99_output_txtfile_name, quote=FALSE,row.names=TRUE,col.names=TRUE)
write.csv(malematrix, file=males99_output_csvfile_name, quote=FALSE,row.names=TRUE)