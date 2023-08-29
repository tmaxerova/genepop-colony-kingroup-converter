f<-1 #leave this as it is
m<-0 #leave this as it is

##Parameters you must update:
path_to_input_colony99_file<-"C:\\Users\\Maxerová\\Desktop\\Misa\\colony_males99.txt"
output_folder_location<-"C:\\Users\\Maxerová\\Desktop\\Misa\\"
sex_vector<-c(f,f,m,m,f) #vector containing information about sex of the individual (from first to last sample), f for female, m for male

##Parameters you can update (but there is no real reason for that):
output_kingroup_txtfile_name<-"kingroup.txt"
output_kingroup_maledot_txtfile_name<-"kingroup_maledot.txt"
output_sample_hetero_diagnostics_file<-"sample_hetero_proportion.csv"
output_snp_hetero_diagnostics_file<-"snp_hetero_proportion.csv"

##OUTPUT FILE SPECIFICATION.
complete_name<-paste0(output_folder_location,"ResultsFolder","_",substr(Sys.time(),1,10),"_",substr(Sys.time(),12,13),"-",substr(Sys.time(),15,16),"-",substr(Sys.time(),18,19))
dir.create(complete_name)
setwd(complete_name)

##COLONY99 TO COLONY CONVERTER.
colony99_data<-read.table(path_to_input_colony99_file)
datacol_number<-ncol(colony99_data)
datarow_number<-nrow(colony99_data)
position99<-matrix(0,datarow_number,datacol_number)
secondalleles_columns<-seq(from=2,to=datacol_number,by=2)
for (v in secondalleles_columns)
{
  for (w in 1:datarow_number)
  {
    if (colony99_data[w,v] == "-99")
    {
      colony99_data[w,v]<-colony99_data[w,v-1]
    }
  }
}
colony_data<-colony99_data

##COLONY TO KINGROUP CONVERTER.
colony_nrow<-nrow(colony_data)
colony_ncol<-ncol(colony_data)
kingnrow<-colony_nrow
kingncol<-(2*colony_ncol)+4-1
kingmatrix<-matrix(0,kingnrow,kingncol)
kingcolnames<-matrix(0,1,kingncol)
for(b in 1:colony_ncol)
{
  for(a in 1:colony_nrow)
  {
    ver1<-colony_data[a,b]
    kingmatrix[a,(3+2*b)]<-ver1
  }
}
comma_cols<-c(2,seq(from=4,to=kingncol,by=4))
slash_cols<-seq(from=6,to=kingncol,by=4)
rowheader<-rownames(colony_data)
for (g in 1:kingnrow)
{
  kingmatrix[g,1]<-substring(rowheader[g],1,5)
  kingmatrix[g,3]<-rowheader[g]
}
kingcolnames[1,1]<-"group"
kingcolnames[1,3]<-"sample"
for (d in comma_cols)
{
  for (c in 1:kingnrow)
  {
    kingmatrix[c,d]<-","
    kingcolnames[1,d]<-","
  }
}
for (ff in slash_cols)
{
  for (e in 1:kingnrow)
  {
    kingmatrix[e,ff]<-"/"
    kingcolnames[1,ff]<-"/"
  }
}
for (h in 1:length(colnames(colony_data)))
{
  kingcolnames[3+2*h]<-colnames(colony_data)[h]
}
colnames(kingmatrix)<-kingcolnames
write.table(kingmatrix, file=output_kingroup_txtfile_name, quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

#KINGROUP MALE REWRITER (in males, all second alleles are rewritten to ".").
malematrix_k<-kingmatrix
rewrite_columns<-seq(from=7,to=kingncol,by=4)
males_position<-which(sex_vector == 0)
if (length(sex_vector) != colony_nrow){
  print("ERROR MESSAGE: Number of elements entered to the sex vector does not match number of samples in the input file.")
}else{
  for (p in 1:length(sex_vector)){
    if (sex_vector[p] == 0){
      for (r in rewrite_columns)
      {
        for (q in males_position)
        {
          malematrix_k[q,r]<-"."
        }
      }
    }
  }
}
write.table(malematrix_k, file=output_kingroup_maledot_txtfile_name, quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

##HETEROZYGOSITY DIAGNOSTICS GENERATOR (both for samples and SNPs).
rdnrow<-(nrow(colony_data))
rdncol<-(ncol(colony_data))/2
hapdip_matrix<-matrix(0,rdnrow,rdncol)
for(b in seq(from=1,to=colony_ncol,by=2))
{
  for(a in 1:rdnrow)
  {
    if (colony_data[a,b] != colony_data[a,b+1])
    {
      hapdip_matrix[a,(b+1)/2]<-1
    }
    
  }
}
samp_diag<-numeric(rdnrow)
for (k in 1:rdnrow)
{
  samp_diag[k]<-(sum(hapdip_matrix[k,])/rdncol)
}
final_samp_diag<-rbind(rownames(colony_data),samp_diag)
snp_diag<-numeric(rdncol)
for (l in 1:rdncol)
{
  snp_diag[l]<-(sum(hapdip_matrix[,l])/rdnrow)
}
diag_names<-matrix(0,1,length(colnames(colony_data))/2)
for (x in 1:length(diag_names))
{
  diag_names[x]<-substring(colnames(colony_data)[2*x-1],1,5)
}
final_snp_diag<-rbind(diag_names,snp_diag)
write.table(final_samp_diag, file=output_sample_hetero_diagnostics_file, sep=",", col.names=FALSE)
write.table(final_snp_diag, file=output_snp_hetero_diagnostics_file, sep=",", col.names=FALSE)