allele_tracing<-function(input_pedigree=NULL,
                         hap_win=50,
						 cpu_cores=1,					
						 haplotype_hap=NULL,
						 haplotype_map=NULL,
						 haplotype_sample=NULL, #can be obtained by blupADC directly
						 trace_direction="backward",
						 ){

ind_breed1=input_pedigree[input_pedigree[,4]==1,1]
ind_breed2=input_pedigree[input_pedigree[,4]==2,1]
ind_cross=input_pedigree[input_pedigree[,4]==0,1]


haplotype_sample=as.character(haplotype_sample[,1])				  
haplotype_sample=rep(haplotype_sample,each=2)

cross_hap=haplotype_hap[,haplotype_sample%in%ind_cross]
sire_hap=haplotype_hap[,haplotype_sample%in%ind_breed1]		
dam_hap=haplotype_hap[,haplotype_sample%in%ind_breed2]	

ind_cross=haplotype_sample[haplotype_sample%in%ind_cross]
ind_breed1=haplotype_sample[haplotype_sample%in%ind_breed1]
ind_breed2=haplotype_sample[haplotype_sample%in%ind_breed2]


ped=trace_pedigree(input_pedigree[,1:3],trace_direction=trace_direction)$ped[,1:3]
cross_ped=ped[ped[,1]%in%ind_cross,]
cross_ped=cross_ped[!duplicated(cross_ped[,1]),]

#已知parents 
Father_cross_ped=cross_ped[cross_ped[,2]%in%ind_breed1,]
Mother_cross_ped=cross_ped[!cross_ped[,2]%in%ind_breed1&cross_ped[,3]%in%ind_breed2,]

#未知parents
Unknow_cross_ped=cross_ped[!cross_ped[,2]%in%ind_breed1&!cross_ped[,3]%in%ind_breed2,]


#for animals without genotyped parents, to get pos_start and pos_end based on hap_win
block_pos=get_haplotype_block(haplotype_map,haplotype_window_nSNP=hap_win)$block$block_pos[,2:3]

#for animals with genotyped sires or dams, using the whole chromosome to get pos_end and pos_start
chr_name=unique(haplotype_map[,1])
chr_pos=matrix(NA,nrow=length(chr_name),ncol=2)
for(i in 1:length(chr_name)){
chr_pos[i,1]=match(chr_name[i],haplotype_map[,1])
chr_pos[i,2]=match(chr_name[i+1],haplotype_map[,1])-1
}
chr_pos[nrow(chr_pos),2]=nrow(haplotype_map)

offspring_status=allele_tracing_cpp(sire_hap,dam_hap,cross_hap,chr_pos-1,as.matrix(block_pos)-1,
			  as.matrix(Father_cross_ped),as.matrix(Mother_cross_ped), as.matrix(Unknow_cross_ped),
			  ind_breed1,ind_breed2,ind_cross,hap_win=hap_win,cpu_cores=cpu_cores)
			  
colnames(offspring_status)=c(rep(Father_cross_ped[,1],each=2),rep(Mother_cross_ped[,1],each=2),rep(Unknow_cross_ped[,1],each=2))
return(offspring_status)
}


#获取单倍型 block
get_haplotype_block<-function(data_map=NULL,
							  haplotype_window_nSNP=NULL,
							  haplotype_window_kb=NULL,
							  haplotype_window_block=NULL){

#同时提供三种参数中至少两种会报错
if(sum(c(!is.null(haplotype_window_nSNP),!is.null(haplotype_window_kb),!is.null(haplotype_window_block)))>=2){stop("Provided too much parateters: haplotype_window_nSNP, haplotype_window_kb, haplotype_window_block!")}
if(sum(c(!is.null(haplotype_window_nSNP),!is.null(haplotype_window_kb),!is.null(haplotype_window_block)))==0){stop("Provided no parateters: haplotype_window_nSNP, haplotype_window_kb, haplotype_window_block!")}

#获取SNP位置-考虑过了染色体
data_map=data.frame(data_map,stringsAsFactors=F)
data_map[,1]=as.numeric(data_map[,1])
data_map[,3]=as.numeric(data_map[,3])
chr_set=unique(data_map[,1])
chr_set_snp_n=0
for(i in 1:length(chr_set)){
tmp_data_map=data_map[data_map[,1]==chr_set[i],]
chr_set_snp_n=c(chr_set_snp_n,nrow(tmp_data_map))
}
chr_set_snp_n=cumsum(chr_set_snp_n)
chr_set_snp_n=chr_set_snp_n[-length(chr_set_snp_n)]

if(!is.null(haplotype_window_block)){

block_start=haplotype_window_block[,1]-1  #window_block为两列数据，均为位置信息
block_end=haplotype_window_block[,2]-1

}else if(!is.null(haplotype_window_nSNP)){ #根据SNP数目划分block

block_start=NULL
block_end=NULL
for(i in 1:length(chr_set)){
tmp_data_map=data_map[data_map[,1]==chr_set[i],]
block_start_tmp=seq(1,nrow(tmp_data_map),haplotype_window_nSNP)-1
block_end_tmp=seq(1,nrow(tmp_data_map),haplotype_window_nSNP)-2
block_end_tmp=c(block_end_tmp,nrow(tmp_data_map)-1)
block_end_tmp=block_end_tmp[-1]

#添加位置信息
block_start_tmp=block_start_tmp+chr_set_snp_n[i]
block_end_tmp=block_end_tmp+chr_set_snp_n[i]
block_start=c(block_start,block_start_tmp)
block_end=c(block_end,block_end_tmp)
}

}else if(!is.null(haplotype_window_kb)){  #根据物理距离划分block

block_start=NULL
block_end=NULL
for(i in 1:length(chr_set)){
tmp_data_map=data_map[data_map[,1]==chr_set[i],]
block=seq(min(tmp_data_map[,3]),max(tmp_data_map[,3]),haplotype_window_kb*1000)

block_1=block
block_2=block-1;
block_2=block_2[-1]; #去除最后一列
block_2=c(block_2,max(tmp_data_map[,3]))

#R-function too slow
#block_start_tmp=NULL
#block_end_tmp=NULL
#for(j in 1:length(block)){
#block_start_tmp=c(block_start_tmp,min(which(tmp_data_map[,3]>=block_1[j])))
#block_end_tmp=c(block_end_tmp,max(which(tmp_data_map[,3]<=block_2[j])))
#}

#Rcpp function 
block_result=define_block_window_kb_cpp(block_1,block_2,tmp_data_map[,3])
block_start_tmp=block_result[[1]]
block_end_tmp=block_result[[2]]

pos_status=block_start_tmp<=block_end_tmp
block_start_tmp=block_start_tmp[pos_status]
block_end_tmp=block_end_tmp[pos_status]
block_start_tmp=block_start_tmp-1
block_end_tmp=block_end_tmp-1
#添加位置信息
block_start_tmp=block_start_tmp+chr_set_snp_n[i]
block_end_tmp=block_end_tmp+chr_set_snp_n[i]
block_start=c(block_start,block_start_tmp)
block_end=c(block_end,block_end_tmp)
}
}
#进行单倍型构建分析
#写出block的信息
snp_map=as.character(data_map[,2])
block_pos=data.frame(block_num=paste0("block",1:length(block_start)),window_start_pos=block_start+1,window_end_pos=block_end+1)
block_snp=data.frame(block_num=paste0("block",1:length(block_start)),window_start_snp=snp_map[block_start+1],window_end_snp=snp_map[block_end+1])
return(list(block=list(block_pos=block_pos,block_snp=block_snp),
			block_rcpp=cbind(block_start,block_end) #c++程序的输入文件
			))
}		
