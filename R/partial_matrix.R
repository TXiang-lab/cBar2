makeA_partial<-function(input_pedigree,exclude_missing_parents=TRUE,output_matrix_type="col_all",IND_rename=FALSE,
						matrix_log_det=FALSE,cpu_cores=1,col_three_threshold=0,trace_direction="backward",
						priority_rename_id=NULL,
			full_rank=TRUE #make sure A1 and A2 are full rank matrices 
		       ){

cat("Please make sure all animals in the pedigree have breed record!  \n") #所有个体均出现在第一列
Pedigree=data.frame(input_pedigree[,1:3],stringsAsFactors=F)
Pedigree[is.na(Pedigree)]=0

Breed=data.frame(input_pedigree[,c(1,4)],stringsAsFactors=F)
colnames(Breed)=c("Id","Breed")

if("col_three" %in% output_matrix_type&IND_rename==FALSE){
if(NA%in%as.numeric(input_pedigree[,1])){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
}

#针对 父母有一方为缺失值的情况
sire_missing=(Pedigree[,2]==0)&(Pedigree[,3]!=0)
dam_missing=(Pedigree[,2]!=0)&(Pedigree[,3]==0)
n1=sum(sire_missing)  #父亲缺失的数目
n2=sum(dam_missing)	  #母亲缺失的数目

sire_record=NULL
dam_record=NULL
if((n1+n2)>0){
cat("Found animals in pedigree only have sire or dam records,these missing records will automatically recoded by software !\n")

if(n1>0){

sire_record=-1*(1:n1)

Pedigree[sire_missing,2]=sire_record

#计算缺失父亲的品种记录，如果个体的品种记录为纯种，那么缺失的父亲的品种就为相应的纯种，反之则为另一个品种
breed_missing=data.frame(Id=sire_record,Breed=Breed[sire_missing,2]) #默认父亲与后代的品种一样
status_crossbreed=(Breed[sire_missing,2]==0)
breed_missing[status_crossbreed,2]=3-as.numeric(Breed[match(Pedigree[sire_missing,3][status_crossbreed],Breed[,1]),2])
Breed=rbind(Breed,breed_missing)
}

if(n2>0){
dam_record=-1*((n1+1):(n1+n2))
Pedigree[dam_missing,3]=dam_record
#计算缺失母亲的品种记录，如果个体的品种记录为纯种，那么缺失的母亲的品种就为相应的纯种，反之则为另一个品种
breed_missing=data.frame(Id=dam_record,Breed=Breed[dam_missing,2])
status_crossbreed=(Breed[dam_missing,2]==0)
breed_missing[status_crossbreed,2]=3-as.numeric(Breed[match(Pedigree[dam_missing,2][status_crossbreed],Breed[,1]),2])
Breed=rbind(Breed,breed_missing)
}

}

Pedigree=trace_pedigree(Pedigree,display_message=F,trace_direction=trace_direction,priority_rename_id=priority_rename_id)

error_id=do.call(c,Pedigree$error_id_set)
Pedigree=Pedigree$rename_ped

if(length(error_id)>0){stop("Please check provided input_pedigree carefully!")}

Pedigree$Breed=as.numeric(Breed[match(Pedigree[,1],Breed[,1]),2])
IND_name=Pedigree[,1]

num_ped=as.matrix(Pedigree[,c(3,4,5,7)])
num_ped=apply(num_ped,2,as.numeric)

if(exclude_missing_parents==TRUE){
record_pos=match(c(sire_record,dam_record),IND_name)
}else{
record_pos=as.numeric(NULL)
}

ped_result=makeA_partial_cpp(num_ped,IND_name,record_pos-1,full_rank)

IND_Breed1=ped_result[[3]]
IND_Breed2=ped_result[[4]]

Breed1=ped_result[[1]]
Breed2=ped_result[[2]]
rm(ped_result);gc();


#record the pedigree initial level to 1
if(length(record_pos)==0){
rename_ped=Pedigree[,c(1,3,4,5,6)]
}else{
rename_ped=Pedigree[-record_pos,c(1,3,4,5,6)]
}
rename_ped[,2:5]=rename_ped[,2:5]-length(record_pos)
rename_ped[rename_ped<0]=0
#

IND_Breed1_num=rename_ped[match(IND_Breed1,rename_ped[,1]),2]
IND_Breed2_num=rename_ped[match(IND_Breed2,rename_ped[,1]),2]

if(IND_rename==TRUE){
IND_Breed1=IND_Breed1_num
IND_Breed2=IND_Breed2_num
}else{
IND_Breed1_num=IND_Breed1
IND_Breed2_num=IND_Breed2
}

rownames(Breed1)=colnames(Breed1)=IND_Breed1
rownames(Breed2)=colnames(Breed2)=IND_Breed2


Breed1_three=NULL
Breed2_three=NULL

if("col_three" %in% output_matrix_type){
if(NA%in%as.numeric(IND_Breed1)|NA%in%as.numeric(IND_Breed2)){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
Breed1_three=matrix_col3(Breed1,IND_geno=as.numeric(IND_Breed1),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
Breed2_three=matrix_col3(Breed2,IND_geno=as.numeric(IND_Breed2),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
}

return(list(Breed1_A=Breed1,Breed1_A_three=Breed1_three,Breed2_A=Breed2,Breed2_A_three=Breed2_three,rename_ped=rename_ped))
}	   


makeAinv_partial<-function(input_pedigree,exclude_missing_parents=TRUE,output_matrix_type="col_all",IND_rename=FALSE,
						matrix_log_det=FALSE,cpu_cores=1,col_three_threshold=0,trace_direction="backward",
			   			priority_rename_id=NULL,
			  full_rank=TRUE #make sure A1 and A2 are full rank matrices 
			  ){

cat("Please make sure all animals in the pedigree have breed record!  \n") #所有个体均出现在第一列
Pedigree=data.frame(input_pedigree[,1:3],stringsAsFactors=F)
Pedigree[is.na(Pedigree)]=0

Breed=data.frame(input_pedigree[,c(1,4)],stringsAsFactors=F)
colnames(Breed)=c("Id","Breed")

if("col_three" %in% output_matrix_type&IND_rename==FALSE){
if(NA%in%as.numeric(input_pedigree[,1])){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
}

#针对 父母有一方为缺失值的情况
sire_missing=(Pedigree[,2]==0)&(Pedigree[,3]!=0)
dam_missing=(Pedigree[,2]!=0)&(Pedigree[,3]==0)
n1=sum(sire_missing)  #父亲缺失的数目
n2=sum(dam_missing)	  #母亲缺失的数目

sire_record=NULL
dam_record=NULL
if((n1+n2)>0){
cat("Found animals in pedigree only have sire or dam records,these missing records will automatically recoded by software !\n")

if(n1>0){

sire_record=-1*(1:n1)

Pedigree[sire_missing,2]=sire_record

#计算缺失父亲的品种记录，如果个体的品种记录为纯种，那么缺失的父亲的品种就为相应的纯种，反之则为另一个品种
breed_missing=data.frame(Id=sire_record,Breed=Breed[sire_missing,2]) #默认父亲与后代的品种一样
status_crossbreed=(Breed[sire_missing,2]==0)
breed_missing[status_crossbreed,2]=3-as.numeric(Breed[match(Pedigree[sire_missing,3][status_crossbreed],Breed[,1]),2])
Breed=rbind(Breed,breed_missing)
}

if(n2>0){
dam_record=-1*((n1+1):(n1+n2))
Pedigree[dam_missing,3]=dam_record
#计算缺失母亲的品种记录，如果个体的品种记录为纯种，那么缺失的母亲的品种就为相应的纯种，反之则为另一个品种
breed_missing=data.frame(Id=dam_record,Breed=Breed[dam_missing,2])
status_crossbreed=(Breed[dam_missing,2]==0)
breed_missing[status_crossbreed,2]=3-as.numeric(Breed[match(Pedigree[dam_missing,2][status_crossbreed],Breed[,1]),2])
Breed=rbind(Breed,breed_missing)
}

}

Pedigree=trace_pedigree(Pedigree,display_message=F,trace_direction=trace_direction,priority_rename_id=priority_rename_id)

error_id=do.call(c,Pedigree$error_id_set)
Pedigree=Pedigree$rename_ped

if(length(error_id)>0){stop("Please check provided input_pedigree carefully!")}

Pedigree$Breed=as.numeric(Breed[match(Pedigree[,1],Breed[,1]),2])
IND_name=Pedigree[,1]

num_ped=as.matrix(Pedigree[,c(3,4,5,7)])
num_ped=apply(num_ped,2,as.numeric)

if(exclude_missing_parents==TRUE){
record_pos=match(c(sire_record,dam_record),IND_name)
}else{
record_pos=as.numeric(NULL)
}

#return(list(num_ped,IND_name,record_pos))
ped_result=makeAinv_partial_cpp(num_ped,IND_name,record_pos-1,full_rank)

IND_Breed1=ped_result[[3]]
IND_Breed2=ped_result[[4]]

Breed1_Ainv=ped_result[[1]]
Breed2_Ainv=ped_result[[2]]
rm(ped_result);gc();


#record the pedigree initial level to 1
#record the pedigree initial level to 1
if(length(record_pos)==0){
rename_ped=Pedigree[,c(1,3,4,5,6)]
}else{
rename_ped=Pedigree[-record_pos,c(1,3,4,5,6)]
}
rename_ped[,2:5]=rename_ped[,2:5]-length(record_pos)
rename_ped[rename_ped<0]=0
#

IND_Breed1_num=rename_ped[match(IND_Breed1,rename_ped[,1]),2]
IND_Breed2_num=rename_ped[match(IND_Breed2,rename_ped[,1]),2]

if(IND_rename==TRUE){
IND_Breed1=IND_Breed1_num
IND_Breed2=IND_Breed2_num
}else{
IND_Breed1_num=IND_Breed1
IND_Breed2_num=IND_Breed2
}

rownames(Breed1_Ainv)=colnames(Breed1_Ainv)=IND_Breed1
rownames(Breed2_Ainv)=colnames(Breed2_Ainv)=IND_Breed2


Breed1_inv_three=NULL
Breed2_inv_three=NULL

if("col_three" %in% output_matrix_type){
if(NA%in%as.numeric(IND_Breed1)|NA%in%as.numeric(IND_Breed2)){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
Breed1_inv_three=matrix_col3(Breed1_Ainv,IND_geno=as.numeric(IND_Breed1),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
Breed2_inv_three=matrix_col3(Breed2_Ainv,IND_geno=as.numeric(IND_Breed2),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
}



return(list(Breed1_Ainv=Breed1_Ainv,Breed1_Ainv_three=Breed1_inv_three,Breed2_Ainv=Breed2_Ainv,Breed2_Ainv_three=Breed2_inv_three,rename_ped=rename_ped))
}


makeHAinv_partial<-function(input_pedigree,offspring_boa,
						 haplotype_hap,haplotype_sample,
						 output_matrix_type="col_all",IND_rename=FALSE,
						matrix_log_det=FALSE,cpu_cores=1,col_three_threshold=0,trace_direction="backward",
			                        priority_rename_id=NULL){


if("col_three" %in% output_matrix_type&IND_rename==FALSE){
if(NA%in%as.numeric(input_pedigree[,1])){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
}

A=makeA_partial(input_pedigree,IND_rename=IND_rename,trace_direction=trace_direction,priority_rename_id=priority_rename_id)
A_Sire=A[[1]]
A_Dam=A[[3]]
rename_ped=A$rename_ped
rm(A);gc();

#get the seperated haplotype data 
ind_breed1=input_pedigree[input_pedigree[,4]==1,1]
ind_breed2=input_pedigree[input_pedigree[,4]==2,1]
ind_cross=input_pedigree[input_pedigree[,4]==0,1]

ped=trace_pedigree(input_pedigree[,1:3],display_message=F,trace_direction=trace_direction,priority_rename_id=priority_rename_id)$ped[,1:3]
cross_ped=ped[ped[,1]%in%ind_cross,]
cross_ped=cross_ped[!duplicated(cross_ped[,1]),]

haplotype_sample=as.character(haplotype_sample[,1])				  
haplotype_sample=rep(haplotype_sample,each=2)

cross_hap=haplotype_hap[,haplotype_sample%in%ind_cross]
sire_hap=haplotype_hap[,haplotype_sample%in%ind_breed1]		
dam_hap=haplotype_hap[,haplotype_sample%in%ind_breed2]	

ind_cross=unique(haplotype_sample[haplotype_sample%in%ind_cross])
ind_breed1=unique(haplotype_sample[haplotype_sample%in%ind_breed1])
ind_breed2=unique(haplotype_sample[haplotype_sample%in%ind_breed2])


real_ind_cross=unique(colnames(offspring_boa))
cross_hap=cross_hap[,match_haplotype(real_ind_cross,ind_cross)]
ind_cross=real_ind_cross

ind_cross_num=rename_ped[match(unique(ind_cross),rename_ped[,1]),2]
ind_breed1_num=rename_ped[match(unique(ind_breed1),rename_ped[,1]),2]
ind_breed2_num=rename_ped[match(unique(ind_breed2),rename_ped[,1]),2]

#不进行重命名
if(IND_rename==FALSE){
ind_cross_num=ind_cross
ind_breed1_num=ind_breed1
ind_breed2_num=ind_breed2
}



G=G_matrix_partial_cpp(t(as.matrix(cross_hap)),t(as.matrix(offspring_boa)),
					   t(as.matrix(sire_hap)),t(as.matrix(dam_hap)))
G_Sire=G$G_Sire
G_Dam=G$G_Dam
rm(G);gc();

#for Sire
IND_pedigree=rownames(A_Sire)

IND_geno=c(ind_breed1_num,ind_cross_num)


IND_geno2=setdiff(IND_geno,IND_pedigree)    #有基因型但是无系谱
IND_A22=intersect(IND_pedigree,IND_geno)    #既有基因型又有系谱
IND_A11=setdiff(IND_pedigree,IND_A22)        #有系谱但是无基因型
IND_Additive=c(IND_A11,IND_A22)               #在H矩阵中的个体
pos_geno=match(IND_A22,IND_geno)-1            #既有基因型又有系谱的个体在 基因型数据中的位置
pos_A11=match(IND_A11,IND_pedigree)-1         #有系谱但是无基因型个体在  在系谱中的位置
pos_A22=match(IND_A22,IND_pedigree)-1         #既有基因型又有系谱的个体在  在系谱中的位置
pos_A=match(IND_Additive,IND_pedigree)-1      #H矩阵的个体  在系谱中的位置 
pos_H22=match(IND_A22,IND_Additive)-1         #既有基因型又有系谱的个体在  在H矩阵中的位置

real_geno=IND_geno[pos_geno+1] #真正在H矩阵中有基因型数据的个体号
n_pure=sum(real_geno%in%ind_breed1_num)
n_cross=sum(real_geno%in%ind_cross_num)

cat("For breed 1......\n")
H_Ainv_Sire=makeHA_partial_cpp(A_Sire, G_Sire, n_pure,n_cross, 
					 IND_geno,pos_A11,pos_A22,
					 pos_geno, pos_A,pos_H22,
					direct=FALSE,inverse=TRUE,omega=0.05)$Hinv


IND_Sire=IND_Additive
	
cat("\n")
cat("For breed 2......\n")
#for Dam
IND_pedigree=rownames(A_Dam)
IND_geno=c(ind_breed2_num,ind_cross_num)

IND_geno2=setdiff(IND_geno,IND_pedigree)    #有基因型但是无系谱
IND_A22=intersect(IND_pedigree,IND_geno)    #既有基因型又有系谱
IND_A11=setdiff(IND_pedigree,IND_A22)        #有系谱但是无基因型
IND_Additive=c(IND_A11,IND_A22)               #在H矩阵中的个体
pos_geno=match(IND_A22,IND_geno)-1            #既有基因型又有系谱的个体在 基因型数据中的位置
pos_A11=match(IND_A11,IND_pedigree)-1         #有系谱但是无基因型个体在  在系谱中的位置
pos_A22=match(IND_A22,IND_pedigree)-1         #既有基因型又有系谱的个体在  在系谱中的位置
pos_A=match(IND_Additive,IND_pedigree)-1      #H矩阵的个体  在系谱中的位置 
pos_H22=match(IND_A22,IND_Additive)-1         #既有基因型又有系谱的个体在  在H矩阵中的位置

real_geno=IND_geno[pos_geno+1] #真正在H矩阵中有基因型数据的个体号
n_pure=sum(real_geno%in%ind_breed2_num)
n_cross=sum(real_geno%in%ind_cross_num)


H_Ainv_Dam=makeHA_partial_cpp(A_Dam, G_Dam, n_pure,n_cross, 
					 IND_geno,pos_A11,pos_A22,
					 pos_geno, pos_A,pos_H22,
					direct=FALSE,inverse=TRUE,omega=0.05)$Hinv
					
IND_Dam=IND_Additive							
					
rownames(H_Ainv_Dam)=colnames(H_Ainv_Dam)=IND_Dam
rownames(H_Ainv_Sire)=colnames(H_Ainv_Sire)=IND_Sire
Breed1_inv_three=NULL
Breed2_inv_three=NULL
if("col_three" %in% output_matrix_type){
if(NA%in%as.numeric(IND_Sire)|NA%in%as.numeric(IND_Dam)){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
Breed1_inv_three=matrix_col3(H_Ainv_Sire,IND_geno=as.numeric(IND_Sire),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
Breed2_inv_three=matrix_col3(H_Ainv_Dam,IND_geno=as.numeric(IND_Dam),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
}


return(list(Breed1_H_Ainv=H_Ainv_Sire,Breed1_H_Ainv_three=Breed1_inv_three,Breed2_H_Ainv=H_Ainv_Dam,Breed2_H_Ainv_three=Breed2_inv_three,rename_ped=rename_ped))

}



makeHA_partial<-function(input_pedigree,offspring_boa,
						 haplotype_hap,haplotype_sample,
						 output_matrix_type="col_all",IND_rename=FALSE,
						matrix_log_det=FALSE,cpu_cores=1,col_three_threshold=0,trace_direction="backward",
			                        priority_rename_id=NULL){


if("col_three" %in% output_matrix_type&IND_rename==FALSE){
if(NA%in%as.numeric(input_pedigree[,1])){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
}

A=makeA_partial(input_pedigree,IND_rename=IND_rename,trace_direction=trace_direction,priority_rename_id=priority_rename_id)
A_Sire=A[[1]]
A_Dam=A[[3]]
rename_ped=A$rename_ped
rm(A);gc();

#get the seperated haplotype data 
ind_breed1=input_pedigree[input_pedigree[,4]==1,1]
ind_breed2=input_pedigree[input_pedigree[,4]==2,1]
ind_cross=input_pedigree[input_pedigree[,4]==0,1]

ped=trace_pedigree(input_pedigree[,1:3],display_message=F,trace_direction=trace_direction,priority_rename_id=priority_rename_id)$ped[,1:3]
cross_ped=ped[ped[,1]%in%ind_cross,]
cross_ped=cross_ped[!duplicated(cross_ped[,1]),]

haplotype_sample=as.character(haplotype_sample[,1])				  
haplotype_sample=rep(haplotype_sample,each=2)

cross_hap=haplotype_hap[,haplotype_sample%in%ind_cross]
sire_hap=haplotype_hap[,haplotype_sample%in%ind_breed1]		
dam_hap=haplotype_hap[,haplotype_sample%in%ind_breed2]	

ind_cross=unique(haplotype_sample[haplotype_sample%in%ind_cross])
ind_breed1=unique(haplotype_sample[haplotype_sample%in%ind_breed1])
ind_breed2=unique(haplotype_sample[haplotype_sample%in%ind_breed2])


real_ind_cross=unique(colnames(offspring_boa))
cross_hap=cross_hap[,match_haplotype(real_ind_cross,ind_cross)]
ind_cross=real_ind_cross

ind_cross_num=rename_ped[match(unique(ind_cross),rename_ped[,1]),2]
ind_breed1_num=rename_ped[match(unique(ind_breed1),rename_ped[,1]),2]
ind_breed2_num=rename_ped[match(unique(ind_breed2),rename_ped[,1]),2]

#不进行重命名
if(IND_rename==FALSE){
ind_cross_num=ind_cross
ind_breed1_num=ind_breed1
ind_breed2_num=ind_breed2
}



G=G_matrix_partial_cpp(t(as.matrix(cross_hap)),t(as.matrix(offspring_boa)),
					   t(as.matrix(sire_hap)),t(as.matrix(dam_hap)))
G_Sire=G$G_Sire
G_Dam=G$G_Dam
rm(G);gc();

#for Sire
IND_pedigree=rownames(A_Sire)

IND_geno=c(ind_breed1_num,ind_cross_num)


IND_geno2=setdiff(IND_geno,IND_pedigree)    #有基因型但是无系谱
IND_A22=intersect(IND_pedigree,IND_geno)    #既有基因型又有系谱
IND_A11=setdiff(IND_pedigree,IND_A22)        #有系谱但是无基因型
IND_Additive=c(IND_A11,IND_A22)               #在H矩阵中的个体
pos_geno=match(IND_A22,IND_geno)-1            #既有基因型又有系谱的个体在 基因型数据中的位置
pos_A11=match(IND_A11,IND_pedigree)-1         #有系谱但是无基因型个体在  在系谱中的位置
pos_A22=match(IND_A22,IND_pedigree)-1         #既有基因型又有系谱的个体在  在系谱中的位置
pos_A=match(IND_Additive,IND_pedigree)-1      #H矩阵的个体  在系谱中的位置 
pos_H22=match(IND_A22,IND_Additive)-1         #既有基因型又有系谱的个体在  在H矩阵中的位置

real_geno=IND_geno[pos_geno+1] #真正在H矩阵中有基因型数据的个体号
n_pure=sum(real_geno%in%ind_breed1_num)
n_cross=sum(real_geno%in%ind_cross_num)

cat("For breed 1......\n")
H_Ainv_Sire=makeHA_partial_cpp(A_Sire, G_Sire, n_pure,n_cross, 
					 IND_geno,pos_A11,pos_A22,
					 pos_geno, pos_A,pos_H22,
					direct=TRUE,inverse=FALSE,omega=0.05)$H


IND_Sire=IND_Additive
	
cat("\n")
cat("For breed 2......\n")
#for Dam
IND_pedigree=rownames(A_Dam)
IND_geno=c(ind_breed2_num,ind_cross_num)

IND_geno2=setdiff(IND_geno,IND_pedigree)    #有基因型但是无系谱
IND_A22=intersect(IND_pedigree,IND_geno)    #既有基因型又有系谱
IND_A11=setdiff(IND_pedigree,IND_A22)        #有系谱但是无基因型
IND_Additive=c(IND_A11,IND_A22)               #在H矩阵中的个体
pos_geno=match(IND_A22,IND_geno)-1            #既有基因型又有系谱的个体在 基因型数据中的位置
pos_A11=match(IND_A11,IND_pedigree)-1         #有系谱但是无基因型个体在  在系谱中的位置
pos_A22=match(IND_A22,IND_pedigree)-1         #既有基因型又有系谱的个体在  在系谱中的位置
pos_A=match(IND_Additive,IND_pedigree)-1      #H矩阵的个体  在系谱中的位置 
pos_H22=match(IND_A22,IND_Additive)-1         #既有基因型又有系谱的个体在  在H矩阵中的位置

real_geno=IND_geno[pos_geno+1] #真正在H矩阵中有基因型数据的个体号
n_pure=sum(real_geno%in%ind_breed2_num)
n_cross=sum(real_geno%in%ind_cross_num)


H_Ainv_Dam=makeHA_partial_cpp(A_Dam, G_Dam, n_pure,n_cross, 
					 IND_geno,pos_A11,pos_A22,
					 pos_geno, pos_A,pos_H22,
					direct=TRUE,inverse=FALSE,omega=0.05)$H
					
IND_Dam=IND_Additive							
					
rownames(H_Ainv_Dam)=colnames(H_Ainv_Dam)=IND_Dam
rownames(H_Ainv_Sire)=colnames(H_Ainv_Sire)=IND_Sire
Breed1_inv_three=NULL
Breed2_inv_three=NULL
if("col_three" %in% output_matrix_type){
if(NA%in%as.numeric(IND_Sire)|NA%in%as.numeric(IND_Dam)){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
Breed1_inv_three=matrix_col3(H_Ainv_Sire,IND_geno=as.numeric(IND_Sire),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
Breed2_inv_three=matrix_col3(H_Ainv_Dam,IND_geno=as.numeric(IND_Dam),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
}


return(list(Breed1_H_A=H_Ainv_Sire,Breed1_H_A_three=Breed1_inv_three,Breed2_H_A=H_Ainv_Dam,Breed2_H_A_three=Breed2_inv_three,rename_ped=rename_ped))

}





###########################################################G_A 
makeGA_partial<-function(input_pedigree,offspring_boa,
						 haplotype_hap,haplotype_sample,
						 output_matrix_type="col_all",IND_rename=FALSE,
						matrix_log_det=FALSE,cpu_cores=1,col_three_threshold=0,trace_direction="backward",
						priority_rename_id=NULL){


if("col_three" %in% output_matrix_type&IND_rename==FALSE){
if(NA%in%as.numeric(input_pedigree[,1])){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
}

rename_ped=trace_pedigree(input_pedigree[,1:3],display_message=F,trace_direction=trace_direction,priority_rename_id=priority_rename_id)$rename_ped[,c(1,3,4,5)]

#get the seperated haplotype data 
ind_breed1=input_pedigree[input_pedigree[,4]==1,1]
ind_breed2=input_pedigree[input_pedigree[,4]==2,1]
ind_cross=input_pedigree[input_pedigree[,4]==0,1]

ped=trace_pedigree(input_pedigree[,1:3],display_message=F,trace_direction=trace_direction,priority_rename_id=priority_rename_id)$ped[,1:3]
cross_ped=ped[ped[,1]%in%ind_cross,]
cross_ped=cross_ped[!duplicated(cross_ped[,1]),]

haplotype_sample=as.character(haplotype_sample[,1])				  
haplotype_sample=rep(haplotype_sample,each=2)

cross_hap=haplotype_hap[,haplotype_sample%in%ind_cross]
sire_hap=haplotype_hap[,haplotype_sample%in%ind_breed1]		
dam_hap=haplotype_hap[,haplotype_sample%in%ind_breed2]	

ind_cross=unique(haplotype_sample[haplotype_sample%in%ind_cross])
ind_breed1=unique(haplotype_sample[haplotype_sample%in%ind_breed1])
ind_breed2=unique(haplotype_sample[haplotype_sample%in%ind_breed2])


real_ind_cross=unique(colnames(offspring_boa))
cross_hap=cross_hap[,match_haplotype(real_ind_cross,ind_cross)]
ind_cross=real_ind_cross


ind_cross_num=rename_ped[match(unique(ind_cross),rename_ped[,1]),2]
ind_breed1_num=rename_ped[match(unique(ind_breed1),rename_ped[,1]),2]
ind_breed2_num=rename_ped[match(unique(ind_breed2),rename_ped[,1]),2]

#不进行重命名
if(IND_rename==FALSE){
ind_cross_num=ind_cross
ind_breed1_num=ind_breed1
ind_breed2_num=ind_breed2
}



G=G_matrix_partial_cpp(t(as.matrix(cross_hap)),t(as.matrix(offspring_boa)),
					   t(as.matrix(sire_hap)),t(as.matrix(dam_hap)))
G_Sire=G$G_Sire
G_Dam=G$G_Dam
rm(G);gc();

#
rownames(G_Dam)=colnames(G_Dam)=c(ind_breed2_num,ind_cross_num)
rownames(G_Sire)=colnames(G_Sire)=c(ind_breed1_num,ind_cross_num)

Breed1_inv_three=NULL
Breed2_inv_three=NULL

if("col_three" %in% output_matrix_type){
if(NA%in%as.numeric(c(ind_breed2_num,ind_cross_num))|NA%in%as.numeric(c(ind_breed1_num,ind_cross_num))){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
Breed1_inv_three=matrix_col3(G_Sire,IND_geno=as.numeric(c(ind_breed1_num,ind_cross_num)),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
Breed2_inv_three=matrix_col3(G_Dam,IND_geno=as.numeric(c(ind_breed2_num,ind_cross_num)),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
}


return(list(Breed1_G_A=G_Sire,Breed1_G_A_three=Breed1_inv_three,Breed2_G_A=G_Dam,Breed2_G_A_three=Breed2_inv_three,rename_ped=rename_ped))

}


###########################################################G_Ainv 
makeGAinv_partial<-function(input_pedigree,offspring_boa,
						 haplotype_hap,haplotype_sample,
						 output_matrix_type="col_all",IND_rename=FALSE,
						matrix_log_det=FALSE,cpu_cores=1,col_three_threshold=0,trace_direction="backward"){


if("col_three" %in% output_matrix_type&IND_rename==FALSE){
if(NA%in%as.numeric(input_pedigree[,1])){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
}

rename_ped=trace_pedigree(input_pedigree[,1:3],display_message=F,trace_direction=trace_direction,priority_rename_id=priority_rename_id)$rename_ped[,c(1,3,4,5)]

#get the seperated haplotype data 
ind_breed1=input_pedigree[input_pedigree[,4]==1,1]
ind_breed2=input_pedigree[input_pedigree[,4]==2,1]
ind_cross=input_pedigree[input_pedigree[,4]==0,1]

ped=trace_pedigree(input_pedigree[,1:3],display_message=F,trace_direction=trace_direction,priority_rename_id=priority_rename_id)$ped[,1:3]
cross_ped=ped[ped[,1]%in%ind_cross,]
cross_ped=cross_ped[!duplicated(cross_ped[,1]),]

haplotype_sample=as.character(haplotype_sample[,1])				  
haplotype_sample=rep(haplotype_sample,each=2)

cross_hap=haplotype_hap[,haplotype_sample%in%ind_cross]
sire_hap=haplotype_hap[,haplotype_sample%in%ind_breed1]		
dam_hap=haplotype_hap[,haplotype_sample%in%ind_breed2]	

ind_cross=unique(haplotype_sample[haplotype_sample%in%ind_cross])
ind_breed1=unique(haplotype_sample[haplotype_sample%in%ind_breed1])
ind_breed2=unique(haplotype_sample[haplotype_sample%in%ind_breed2])


real_ind_cross=unique(colnames(offspring_boa))
cross_hap=cross_hap[,match_haplotype(real_ind_cross,ind_cross)]
ind_cross=real_ind_cross


ind_cross_num=rename_ped[match(unique(ind_cross),rename_ped[,1]),2]
ind_breed1_num=rename_ped[match(unique(ind_breed1),rename_ped[,1]),2]
ind_breed2_num=rename_ped[match(unique(ind_breed2),rename_ped[,1]),2]

#不进行重命名
if(IND_rename==FALSE){
ind_cross_num=ind_cross
ind_breed1_num=ind_breed1
ind_breed2_num=ind_breed2
}



G=G_matrix_partial_cpp(t(as.matrix(cross_hap)),t(as.matrix(offspring_boa)),
					   t(as.matrix(sire_hap)),t(as.matrix(dam_hap)))
G_Sire=solve(G$G_Sire)
G_Dam=solve(G$G_Dam)
rm(G);gc();

#
rownames(G_Dam)=colnames(G_Dam)=c(ind_breed2_num,ind_cross_num)
rownames(G_Sire)=colnames(G_Sire)=c(ind_breed1_num,ind_cross_num)

Breed1_inv_three=NULL
Breed2_inv_three=NULL

if("col_three" %in% output_matrix_type){
if(NA%in%as.numeric(c(ind_breed2_num,ind_cross_num))|NA%in%as.numeric(c(ind_breed1_num,ind_cross_num))){stop("Provided pedigree is not numeric format, please set IND_rename=TRUE for outputing col_three format matrix!")}
Breed1_inv_three=matrix_col3(G_Sire,IND_geno=as.numeric(c(ind_breed1_num,ind_cross_num)),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
Breed2_inv_three=matrix_col3(G_Dam,IND_geno=as.numeric(c(ind_breed2_num,ind_cross_num)),det=matrix_log_det,cpu_cores=cpu_cores,threshold=col_three_threshold) 
}

return(list(Breed1_G_Ainv=G_Sire,Breed1_G_Ainv_three=Breed1_inv_three,Breed2_G_Ainv=G_Dam,Breed2_G_Ainv_three=Breed2_inv_three,rename_ped=rename_ped))

}


match_haplotype<-function(x,y){

pos=match(x,y)

pos_final=NULL
for(i in 1:length(pos)){
pos_final=c(pos_final,2*pos[i]-1,2*pos[i])
}
return(pos_final)
}


#A=makeA_partial(ped2)
#Ainv=makeAinv_partial(ped2)
#
#
#A1=A[[1]]
#A1_inv=Ainv[[1]]
#
#identical(round(A1,2),round(solve(A1_inv),2))




##########test_pedigree
#ped1=data.frame(
#			   V1=1:8,
#			   V2=c(0,0,1,0,0,4,2,3),
#			   V3=c(0,0,2,0,0,5,6,6),
#			   V4=c(1,1,1,2,2,2,0,0)
#			   )
#			   
#ped2=data.frame(
#			   V1=1:11,
#			   V2=c(0,0,1,1,0,0,6,6,7,4,1),
#			   V3=c(0,0,2,3,0,0,5,5,8,9,9),
#			   V4=c(1,1,0,0,2,2,0,0,0,0,0)
#			   )			   
#			   
#ped3=data.frame(
#				V1=1:11,
#				V2=c(0,0,0,0,1,3,3,5,7,9,5),
#				V3=c(0,0,0,0,2,2,4,6,6,8,8),
#				V4=c(1,1,2,2,1,0,2,0,0,0,0)
#				)
#			   
#		
#ped=as.matrix(data.frame(V1=1:8,
#			   V2=c(0,0,1,0,0,4,2,3),
#			   V3=c(0,0,2,0,0,5,6,6),
#			   V4=c(1,1,1,2,2,2,0,0))
#			   )
#
#input_pedigree=data.frame(V1=1:8,
#			   V2=c(0,0,0,0,0,4,2,3),
#			   V3=c(0,0,0,0,0,5,6,6),
#			   V4=c(1,1,1,2,2,2,0,0)
#			   )
#			   
#test_pedigree=data.frame(V1=1:11,
#						 V2=c(0,0,0,0,1,3,3,5,1,3,3),
#						 V3=c(0,0,0,0,2,2,4,7,7,2,2),
#						 V4=c(1,1,2,2,1,0,2,0,0,0,0))
#						 
#test_pedigree=data.frame(V1=1:11,
#						 V2=c(0,0,0,0,1,3,3,5,7,9,5),
#						 V3=c(0,0,0,0,2,2,4,6,6,8,8),
#						 V4=c(1,1,2,2,0,0,0,0,0,0,0))			   			   
#
#		
#
#a=makeA_partial_D(num_ped,IND_name,as.numeric(NULL)-1)
#solve(a$T1)
#A1=a$A1
#D1=a$D1
#T1=a$T1		
#identical(A1,T1%*%D1%*%t(T1))
#
#A2=a$A2
#D2=a$D2
#T2=a$T2		
#identical(A2,T2%*%D2%*%t(T2))
#
#
#
#dynmix::ldlt(A1)
#solve(dynmix::ldlt(A1)$L)
#DD=dynmix::ldlt(A1)$L%*%dynmix::ldlt(A1)$D
#num_ped[num_ped[,4]!=2,]






	   
