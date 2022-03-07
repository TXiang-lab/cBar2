#include "shared_function.h"

// [[Rcpp::export]]
arma::mat allele_tracing_cpp(arma::mat & sire_hap,arma::mat & dam_hap,arma::mat & cross_hap, arma::mat & chr_pos,arma::mat & block_pos,
						 CharacterMatrix Father_cross_ped,CharacterMatrix Mother_cross_ped, CharacterMatrix Unknow_cross_ped,
						 std::vector<std::string> sire_hap_ind,std::vector<std::string> dam_hap_ind,std::vector<std::string> cross_hap_ind,
						 int hap_win=50,int cpu_cores=1){

 omp_set_num_threads(cpu_cores);

						
int snp_size=cross_hap.n_rows,ind1_size=Father_cross_ped.nrow(),ind2_size=Mother_cross_ped.nrow(),ind3_size=Unknow_cross_ped.nrow();
std::string ind_offspring,ind_sire,ind_dam;
arma::mat hap_offspring,hap_sire,hap_dam;
//arma::mat Father_offspring_status(snp_size,ind1_size*2),Mother_offspring_status(snp_size,ind2_size*2),Unknow_offspring_status(snp_size,ind3_size*2);
arma::mat offspring_status(snp_size,(ind1_size+ind2_size+ind3_size)*2);
arma::vec offspring_L,offspring_R,sire_L,sire_R,dam_L,dam_R;
arma::uvec offsprint_L_sire_L,offsprint_L_sire_R,offsprint_R_sire_L,offsprint_R_sire_R;
arma::uvec offsprint_L_dam_L,offsprint_L_dam_R,offsprint_R_dam_L,offsprint_R_dam_R;

int i,j,m,n,index_offspring,index_sire,index_dam,index_max,pos_start,pos_end;
int status_L,status_R;
std::vector<double> v;
std::vector<string>::iterator iter_offspring,iter_sire,iter_dam;
std::vector<double>::iterator iter_maximum;

//situation1- Sire have been genotyped
Rcout<<"Start allele tracing for "<<Father_cross_ped.nrow()<<" animals with genotyped sires......"<<endl; 
Progress p1(ind1_size*chr_pos.n_rows,true);	
#pragma omp parallel for private(i,j,ind_offspring,ind_sire,iter_offspring,index_offspring,iter_sire,index_sire,hap_offspring,hap_sire,pos_start,pos_end,status_R,status_L,index_max,iter_maximum,v,offsprint_R_sire_R,offsprint_R_sire_L,offsprint_L_sire_R,offsprint_L_sire_L,sire_R,sire_L,offspring_R,offspring_L) 
 
for(i=0; i < ind1_size; i++){

	ind_offspring=Father_cross_ped(i,0);
    ind_sire=Father_cross_ped(i,1);

	iter_offspring = std::find(cross_hap_ind.begin(), cross_hap_ind.end(), ind_offspring);
    index_offspring = std::distance(cross_hap_ind.begin(), iter_offspring);

	iter_sire = std::find(sire_hap_ind.begin(), sire_hap_ind.end(), ind_sire);
    index_sire = std::distance(sire_hap_ind.begin(), iter_sire);

    hap_offspring=cross_hap.cols( index_offspring, index_offspring+1 );
	hap_sire=sire_hap.cols( index_sire, index_sire+1 );
	
    for(j=0; j < chr_pos.n_rows; j++){
	p1.increment();
	pos_start=chr_pos(j,0);
	pos_end=chr_pos(j,1);	

	offspring_L=hap_offspring(span(pos_start,pos_end),0).col(0);
	offspring_R=hap_offspring(span(pos_start,pos_end),1).col(0);
	sire_L=hap_sire(span(pos_start,pos_end),0).col(0);
	sire_R=hap_sire(span(pos_start,pos_end),1).col(0);


	offsprint_L_sire_L=find((offspring_L-sire_L)==0);
	offsprint_L_sire_R=find((offspring_L-sire_R)==0);
	offsprint_R_sire_L=find((offspring_R-sire_L)==0);
	offsprint_R_sire_R=find((offspring_R-sire_R)==0);

	v={offsprint_L_sire_L.size(),offsprint_L_sire_R.size(),offsprint_R_sire_L.size(),offsprint_R_sire_R.size()};
	iter_maximum = std::max_element(std::begin(v), std::end(v));  
	index_max=std::distance(std::begin(v), iter_maximum);

	if(index_max<=1){
		
		status_L=1;
		status_R=2;
		
	}else{
		
		status_L=2;
		status_R=1;		
	}
	
	offspring_status.submat(span(pos_start,pos_end),span(2*i,2*i)).fill(status_L);
	offspring_status.submat(span(pos_start,pos_end),span(2*i+1,2*i+1)).fill(status_R);
	}

}


//situation2- Mother have been genotyped
Rcout<<""<<endl;
Rcout<<"Start allele tracing for "<<Mother_cross_ped.nrow()<<" animals with genotyped dams......"<<endl; 
Progress p2(ind2_size*chr_pos.n_rows,true);	 
#pragma omp parallel for private(i,j,ind_offspring,ind_dam,iter_offspring,index_offspring,iter_dam,index_dam,hap_offspring,hap_dam,pos_start,pos_end,status_R,status_L,index_max,iter_maximum,v,offsprint_R_dam_R,offsprint_R_dam_L,offsprint_L_dam_R,offsprint_L_dam_L,dam_R,dam_L,offspring_R,offspring_L) 

for(i=0; i < ind2_size; i++){

	ind_offspring=Mother_cross_ped(i,0);
    ind_dam=Mother_cross_ped(i,2);

	iter_offspring = std::find(cross_hap_ind.begin(), cross_hap_ind.end(), ind_offspring);
    index_offspring = std::distance(cross_hap_ind.begin(), iter_offspring);

	iter_dam = std::find(dam_hap_ind.begin(), dam_hap_ind.end(), ind_dam);
    index_dam = std::distance(dam_hap_ind.begin(), iter_dam);

    hap_offspring=cross_hap.cols( index_offspring, index_offspring+1 );
	hap_dam=dam_hap.cols( index_dam, index_dam+1 );
	
    for(j=0; j < chr_pos.n_rows; j++){
	p2.increment();	
	pos_start=chr_pos(j,0);
	pos_end=chr_pos(j,1);	

	offspring_L=hap_offspring(span(pos_start,pos_end),0).col(0);
	offspring_R=hap_offspring(span(pos_start,pos_end),1).col(0);
	dam_L=hap_dam(span(pos_start,pos_end),0).col(0);
	dam_R=hap_dam(span(pos_start,pos_end),1).col(0);


	offsprint_L_dam_L=find((offspring_L-dam_L)==0);
	offsprint_L_dam_R=find((offspring_L-dam_R)==0);
	offsprint_R_dam_L=find((offspring_R-dam_L)==0);
	offsprint_R_dam_R=find((offspring_R-dam_R)==0);

	v={offsprint_L_dam_L.size(),offsprint_L_dam_R.size(),offsprint_R_dam_L.size(),offsprint_R_dam_R.size()};
	iter_maximum = std::max_element(std::begin(v), std::end(v));  
	index_max=std::distance(std::begin(v), iter_maximum);

	if(index_max<=1){
		
		status_L=2;
		status_R=1;
		
	}else{
		
		status_L=1;
		status_R=2;		
	}
	
	offspring_status.submat(span(pos_start,pos_end),span(2*i+2*ind1_size,2*i+2*ind1_size)).fill(status_L);
	offspring_status.submat(span(pos_start,pos_end),span(2*i+1+2*ind1_size,2*i+1+2*ind1_size)).fill(status_R);
	}

}
Rcout<<""<<endl;
//situation3- Non of Father and Mother have been genotyped 
Rcout<<"Start allele tracing for "<<Unknow_cross_ped.nrow()<<" animals without genotyped parents......"<<endl; 
arma::vec sire_both,dam_both,offsprint_L_sire(sire_hap.n_cols),offsprint_R_sire(sire_hap.n_cols),offsprint_L_dam(dam_hap.n_cols),offsprint_R_dam(dam_hap.n_cols);
arma::vec a(sire_hap.n_cols),b(sire_hap.n_cols),c(dam_hap.n_cols),d(dam_hap.n_cols);
arma::uvec offsprint_L_sire_both,offsprint_R_sire_both,offsprint_L_dam_both,offsprint_R_dam_both;
arma::uvec offsprint_L_sire_mean,offsprint_R_dam_mean,offsprint_L_dam_mean,offsprint_R_sire_mean;
arma::vec tmp1(snp_size),tmp2(snp_size);
//situation3- Non of Father and Mother have been genotyped
Progress p3(ind3_size*block_pos.n_rows,true);	
#pragma omp parallel for private(i,j,m,n,offsprint_L_sire,offsprint_R_sire,offsprint_L_dam,offsprint_R_dam,ind_offspring,iter_offspring,v,status_R,status_L,index_max,iter_maximum,offsprint_R_sire_mean,offsprint_L_dam_mean,offsprint_R_dam_mean,offsprint_L_sire_mean,offsprint_L_dam_both,offsprint_R_dam_both,offsprint_L_sire_both,offsprint_R_sire_both,sire_both,dam_both,offspring_L,offspring_R,pos_start,pos_end,hap_offspring,index_offspring) 
for(i=0; i < ind3_size; i++){

	ind_offspring=Unknow_cross_ped(i,0);

	iter_offspring = std::find(cross_hap_ind.begin(), cross_hap_ind.end(), ind_offspring);
    index_offspring = std::distance(cross_hap_ind.begin(), iter_offspring);

    hap_offspring=cross_hap.cols( index_offspring, index_offspring+1 );
	
    for(j=0; j < block_pos.n_rows; j++){
	p3.increment();	
	//Rcout<<j<<endl;
	offsprint_L_sire=a;
	offsprint_R_sire=b;


	offsprint_L_dam=c;
	offsprint_R_dam=d;	
	
	pos_start=block_pos(j,0);
	pos_end=block_pos(j,1);	

	offspring_L=hap_offspring(span(pos_start,pos_end),0).col(0);
	offspring_R=hap_offspring(span(pos_start,pos_end),1).col(0);
	
	
    for(m=0;m<sire_hap.n_cols;m++){
	sire_both=sire_hap(span(pos_start,pos_end),m).col(0);

	offsprint_L_sire_both=find((offspring_L-sire_both)==0);
	offsprint_R_sire_both=find((offspring_R-sire_both)==0);
		
    offsprint_L_sire[m]=offsprint_L_sire_both.size();
    offsprint_R_sire[m]=offsprint_R_sire_both.size();
    }

    for(n=0;n<dam_hap.n_cols;n++){
	dam_both=dam_hap(span(pos_start,pos_end),n).col(0);
	offsprint_L_dam_both=find((offspring_L-dam_both)==0);
	offsprint_R_dam_both=find((offspring_R-dam_both)==0);
		
    offsprint_L_dam[n]=offsprint_L_dam_both.size();
    offsprint_R_dam[n]=offsprint_R_dam_both.size();
    }	
	
	
	offsprint_L_sire_mean=find(offsprint_L_sire==hap_win);
	offsprint_R_dam_mean=find(offsprint_R_dam==hap_win);
	offsprint_L_dam_mean=find(offsprint_L_dam==hap_win);
	offsprint_R_sire_mean=find(offsprint_R_sire==hap_win);
	
	v={1.0 *offsprint_L_sire_mean.size()/offsprint_L_sire.size(),
	   1.0 *offsprint_R_dam_mean.size()/offsprint_R_dam.size(),
	   1.0 *offsprint_L_dam_mean.size()/offsprint_L_dam.size(),
	   1.0 *offsprint_R_sire_mean.size()/offsprint_R_sire.size()};
	   
	iter_maximum = std::max_element(std::begin(v), std::end(v));  
	index_max=std::distance(std::begin(v), iter_maximum);

	if(index_max<=1){
		
		status_L=1;
		status_R=2;
		
	}else{
		
		status_L=2;
		status_R=1;		
	}
	
	offspring_status.submat(span(pos_start,pos_end),span(2*i+(ind1_size+ind2_size)*2,2*i+(ind1_size+ind2_size)*2)).fill(status_L);
	offspring_status.submat(span(pos_start,pos_end),span(2*i+1+(ind1_size+ind2_size)*2,2*i+1+(ind1_size+ind2_size)*2)).fill(status_R);
	//tmp1.subvec(pos_start,pos_end).fill(status_L);
	//tmp2.subvec(pos_start,pos_end).fill(status_R);
	}
	//Unknow_offspring_status.col(2*i)=tmp1;
	//Unknow_offspring_status.col(2*i+1)=tmp2;
	

}

arma::uvec pos1,pos2;
//situation3-final result
for(i=0; i < ind3_size; i++){

for(j=0; j < chr_pos.n_rows; j++){

pos_start=chr_pos(j,0);
pos_end=chr_pos(j,1);

pos1=find(offspring_status(span(pos_start,pos_end),2*i+(ind1_size+ind2_size)*2).col(0)==1);
pos2=find(offspring_status(span(pos_start,pos_end),2*i+1+(ind1_size+ind2_size)*2).col(0)==1);


	if(pos1.size()>pos2.size()){
		
		status_L=1;
		status_R=2;
		
	}else{
		
		status_L=2;
		status_R=1;		
	}
	
	offspring_status.submat(span(pos_start,pos_end),span(2*i+(ind1_size+ind2_size)*2,2*i+(ind1_size+ind2_size)*2)).fill(status_L);
	offspring_status.submat(span(pos_start,pos_end),span(2*i+1+(ind1_size+ind2_size)*2,2*i+1+(ind1_size+ind2_size)*2)).fill(status_R);

}

}

	return offspring_status;
}
