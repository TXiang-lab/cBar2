#include "shared_function.h"

//计算杂交的A矩阵
// [[Rcpp::export]]
List makeA_partial_cpp(arma::mat Pedigree, 
					   std::vector<std::string> IND_name,  //重命名前的个体号，与Pedigree一一对应
					   IntegerVector record_pos,bool full_rank){   //前三列为个体号，父亲，母亲，第四列为品种(1,2,未知品种为0)


	int n_ind=Pedigree.n_rows;

	int Sire_id,Dam_id;
	double A1_svalue,A1_dvalue,A2_svalue,A2_dvalue,A12_svalue,A12_dvalue;
	arma::mat A1(n_ind,n_ind),A2(n_ind,n_ind); //品种的亲缘矩阵

    arma::mat Breed_matrix(n_ind,2);	
	A1.fill(0);
	A2.fill(0);
	Breed_matrix.fill(0);


	arma::vec Breed=Pedigree.col(3);
	arma::vec Animal=Pedigree.col(0);
	arma::vec Sire=Pedigree.col(1);
	arma::vec Dam=Pedigree.col(2);

	for(int i=0; i < n_ind; i++){
		Sire_id=Sire[i];
		Dam_id=Dam[i];
		
		if(Sire_id==0&&Dam_id==0){  //父母均未知,
		
			Breed_matrix(i,Breed[i]-1)=1; //父母未知的情况下,此时品种必须是已知的
			
			if(Breed[i]==1){			
				A1(i,i)=1;	
	
			}else if(Breed[i]==2){				
				A2(i,i)=1;					
			}	
		}else{                     //父母均已知
			
			Breed_matrix(i,0)=(Breed_matrix(Sire_id-1,0)+Breed_matrix(Dam_id-1,0))/2; //父母均已知,计算品种1的比例
			
			Breed_matrix(i,1)=(Breed_matrix(Sire_id-1,1)+Breed_matrix(Dam_id-1,1))/2; //父母均已知,计算品种1的比例
			
			A1(i,i)=Breed_matrix(i,0)+(A1(Sire_id-1,Dam_id-1))/2;	
			A2(i,i)=Breed_matrix(i,1)+(A2(Sire_id-1,Dam_id-1))/2;	

				for(int j=0;j<i;j++){
				
					A1_svalue=A1(j,Sire_id-1);	
					A2_svalue=A2(j,Sire_id-1);	

					A1_dvalue=A1(j,Dam_id-1);	
					A2_dvalue=A2(j,Dam_id-1);	
				
				A1(j,i)=A1(i,j)=(A1_svalue+A1_dvalue)/2;	
				A2(j,i)=A2(i,j)=(A2_svalue+A2_dvalue)/2;				
				}
			
		}
	}
	
	//删除 父母中有一个缺失的个体
	if(record_pos.size()>0){
		
		for(int t=0;t<record_pos.size();t++){
			Breed_matrix(record_pos[t],0)=0;
			Breed_matrix(record_pos[t],1)=0;
		}
	}
	
	 std::vector<std::string>  id1=IND_name;
	 std::vector<std::string>  id2=IND_name;


	arma::vec breed1=Breed_matrix.col(0);
	arma::vec breed2=Breed_matrix.col(1);
	
	if(full_rank==true){  // make sure A1 and A2 are full rank matrices	
	arma::uvec row_condition1=arma::find(breed1!=0);  //A1矩阵中为零的行
	arma::uvec row_condition2=arma::find(breed2!=0);  //A2矩阵中为零的行
	
	A1=A1.submat(row_condition1,row_condition1);
	A2=A2.submat(row_condition2,row_condition2);
	
	
	for(int k=breed1.size()-1;k>-1;k--){
		if(breed1[k]==0){
		id1.erase(id1.begin()+k);		
		}
		
		if(breed2[k]==0){
		id2.erase(id2.begin()+k);		
		}		
	}
	}
	
	return List::create(Named("Breed1_A") = A1,
						Named("Breed2_A") = A2,
						Named("id_Breed1_A")=id1,
						Named("id_Breed2_A")=id2,
						Named("Breed")=Breed_matrix); 
}


//计算杂交的A矩阵
// [[Rcpp::export]]
List makeAinv_partial_cpp(arma::mat Pedigree, 
					   std::vector<std::string> IND_name,  //重命名前的个体号，与Pedigree一一对应
					   IntegerVector record_pos,bool full_rank){   //前三列为个体号，父亲，母亲，第四列为品种(1,2,未知品种为0)


	int n_ind=Pedigree.n_rows;

	int Sire_id,Dam_id;
	double A1_svalue,A1_dvalue,A2_svalue,A2_dvalue,A12_svalue,A12_dvalue;
	arma::mat A1(n_ind,n_ind),A2(n_ind,n_ind); //品种的亲缘矩阵
	arma::mat D1_inv(n_ind,n_ind,fill::eye),D2_inv(n_ind,n_ind,fill::eye); 
	arma::mat T1_inv(n_ind,n_ind,fill::eye),T2_inv(n_ind,n_ind,fill::eye); 

    arma::mat Breed_matrix(n_ind,2);	
	A1.fill(0);
	A2.fill(0);
	Breed_matrix.fill(0);

	
	arma::vec Breed=Pedigree.col(3);
	arma::vec Animal=Pedigree.col(0);
	arma::vec Sire=Pedigree.col(1);
	arma::vec Dam=Pedigree.col(2);

	
	for(int i=0; i < n_ind; i++){
		
		Sire_id=Sire[i];
		Dam_id=Dam[i];
		
		if(Sire_id==0&&Dam_id==0){  //父母均未知,
		
			Breed_matrix(i,Breed[i]-1)=1; //父母未知的情况下,此时品种必须是已知的
			
			if(Breed[i]==1){			
				A1(i,i)=1;
	
			}else if(Breed[i]==2){				
				A2(i,i)=1;
			}	
		}else{                     //父母均已知		
			Breed_matrix(i,0)=(Breed_matrix(Sire_id-1,0)+Breed_matrix(Dam_id-1,0))/2; //父母均已知,计算品种1的比例
			Breed_matrix(i,1)=(Breed_matrix(Sire_id-1,1)+Breed_matrix(Dam_id-1,1))/2; //父母均已知,计算品种1的比例
			
			A1(i,i)=Breed_matrix(i,0)+(A1(Sire_id-1,Dam_id-1))/2;	
			A2(i,i)=Breed_matrix(i,1)+(A2(Sire_id-1,Dam_id-1))/2;	


			D1_inv(i,i)=Breed_matrix(i,0)-A1(Sire_id-1,Sire_id-1)/4-A1(Dam_id-1,Dam_id-1)/4;
			D2_inv(i,i)=Breed_matrix(i,1)-A2(Dam_id-1,Dam_id-1)/4-A2(Sire_id-1,Sire_id-1)/4;
			T1_inv(i,Sire_id-1)=-0.5;
			T2_inv(i,Dam_id-1)=-0.5;
			T1_inv(i,Dam_id-1)=-0.5;
			T2_inv(i,Sire_id-1)=-0.5;			


				for(int j=0;j<i;j++){
				
					A1_svalue=A1(j,Sire_id-1);	
					A2_svalue=A2(j,Sire_id-1);	

					A1_dvalue=A1(j,Dam_id-1);	
					A2_dvalue=A2(j,Dam_id-1);	
					
				A1(j,i)=A1(i,j)=(A1_svalue+A1_dvalue)/2;	
				A2(j,i)=A2(i,j)=(A2_svalue+A2_dvalue)/2;				
				}

		}
	}

	
	//删除 父母中有一个缺失的个体
	if(record_pos.size()>0){
		
		for(int t=0;t<record_pos.size();t++){
			Breed_matrix(record_pos[t],0)=0;
			Breed_matrix(record_pos[t],1)=0;
		}
	}

	 std::vector<std::string>  id1=IND_name;
	 std::vector<std::string>  id2=IND_name;

	arma::vec breed1=Breed_matrix.col(0);
	arma::vec breed2=Breed_matrix.col(1);


	arma::uvec row_condition1=arma::find(breed1!=0);  //A1矩阵中为零的行
	arma::uvec row_condition2=arma::find(breed2!=0);  //A2矩阵中为零的行
	

	D1_inv=D1_inv.submat(row_condition1,row_condition1);
	T1_inv=T1_inv.submat(row_condition1,row_condition1);
	A1=A1.submat(row_condition1,row_condition1);
	
	D2_inv=D2_inv.submat(row_condition2,row_condition2);
	T2_inv=T2_inv.submat(row_condition2,row_condition2);
	A2=A2.submat(row_condition2,row_condition2);
	
if(full_rank==true){		
	for(int k=breed1.size()-1;k>-1;k--){
		if(breed1[k]==0){
		id1.erase(id1.begin()+k);		
		}
		
		if(breed2[k]==0){
		id2.erase(id2.begin()+k);		
		}		
	}
}	
	D1_inv.diag()=1/D1_inv.diag();
	D2_inv.diag()=1/D2_inv.diag();	
	
//if(full_rank==false){
//
//	arma::uvec row_condition11=arma::find(breed1==0);  //A1矩阵中为零的行
//	arma::uvec row_condition22=arma::find(breed2==0);  //A2矩阵中为零的行
//	
//	D1_inv.submat(row_condition11,row_condition11)=arma::mat(row_condition11.size(),row_condition11.size(),fill::zeros);
//	T1_inv.submat(row_condition11,row_condition11)=arma::mat(row_condition11.size(),row_condition11.size(),fill::zeros);
//
//	D2_inv.submat(row_condition22,row_condition22)=arma::mat(row_condition22.size(),row_condition22.size(),fill::zeros);
//	T2_inv.submat(row_condition22,row_condition22)=arma::mat(row_condition22.size(),row_condition22.size(),fill::zeros);
//}


	arma::mat A1_inv(id1.size(),id1.size(),fill::zeros),A2_inv(id2.size(),id2.size(),fill::zeros);

	if(full_rank==true){  // make sure A1 and A2 are full rank matrices
	
	
	A1_inv=T1_inv.t()*D1_inv*T1_inv;
	A2_inv=T2_inv.t()*D2_inv*T2_inv;
	}else{
		
	A1_inv.submat(row_condition1,row_condition1)=T1_inv.t()*D1_inv*T1_inv;	
	A2_inv.submat(row_condition2,row_condition2)=T2_inv.t()*D2_inv*T2_inv;	
	}

	
	return List::create(Named("A1_inv") = A1_inv,
						Named("A2_inv") = A2_inv,
						Named("id_Breed1_A")=id1,
						Named("id_Breed2_A")=id2,
						Named("A1")=A1,
						Named("A2")=A2); 
}

//得到offsping-allele-line matrix for two-way cross
// [[Rcpp::export]]
List G_matrix_partial_cpp(arma::Mat<int> & cross_hap, arma::Mat<int> & cross_boa,
						  arma::Mat<int> & sire_hap,arma::Mat<int> & dam_hap){ //提供连个矩阵，一个是单倍型矩阵，一个是origin信息
	
	if((cross_hap.n_cols!=cross_boa.n_cols)||(cross_hap.n_rows!=cross_boa.n_rows)){	
		throw Rcpp::exception("Provided two matrix should have the same dimension!");	
	}
	
	int ind_cross_size=cross_hap.n_rows/2,snp_size=cross_hap.n_cols,i,j;
	int ind_sire_size=sire_hap.n_rows/2,ind_dam_size=dam_hap.n_rows/2;
	arma::Mat<double> M_cross_sire(ind_cross_size,snp_size,fill::zeros),M_cross_dam(ind_cross_size,snp_size,fill::zeros);
	
	for(i=0; i < ind_cross_size; i++){
		
		for(j=0; j < snp_size; j++){
			
			if(cross_boa(2*i,j)==1){  //来自父亲
				M_cross_sire(i,j)=cross_hap(2*i,j);
				M_cross_dam(i,j)=cross_hap(2*i+1,j);
			}else{
				M_cross_sire(i,j)=cross_hap(2*i+1,j);
				M_cross_dam(i,j)=cross_hap(2*i,j);	
			}
		
		}
	}


//pure sire and pure dam
	arma::Mat<double> M_sire(sire_hap.n_rows/2,snp_size,fill::zeros),M_dam(dam_hap.n_rows/2,snp_size,fill::zeros);
	for(i=0; i < sire_hap.n_rows/2; i++){	
		for(j=0; j < snp_size; j++){				
		M_sire(i,j)=sire_hap(2*i,j)+sire_hap(2*i+1,j);		
		}
	}	
	
	for(i=0; i < dam_hap.n_rows/2; i++){	
		for(j=0; j < snp_size; j++){				
		M_dam(i,j)=dam_hap(2*i,j)+dam_hap(2*i+1,j);		
		}
	}	
//arma::uvec sire_odd_pos=regspace<uvec>(0,2,sire_hap.n_rows-1);	
//arma::uvec sire_even_pos=regspace<uvec>(1,2,sire_hap.n_rows-1);	
//arma::uvec dam_odd_pos=regspace<uvec>(0,2,dam_hap.n_rows-1);	
//arma::uvec dam_even_pos=regspace<uvec>(1,2,dam_hap.n_rows-1);	
//
//arma::Mat<double> M_sire=conv_to<arma::Mat<double>>::from(sire_hap.rows(sire_odd_pos)+sire_hap.rows(sire_even_pos));
//arma::Mat<double> M_dam=conv_to<arma::Mat<double>>::from(dam_hap.rows(dam_odd_pos)+dam_hap.rows(dam_even_pos));	

arma::Mat<double> t_sire=M_sire;
arma::Mat<double> t_dam=M_dam;

	arma::Row<double> p_A(snp_size),p_B(snp_size);  // p_A ：allele frequency in breedA, p_B: allele frequency in breedB
	//construct partial genomic additive matrix
	for(j=0;j < snp_size;j++){
	p_A[j]=mean(M_sire.col(j))/2;
	p_B[j]=mean(M_dam.col(j))/2;
	M_cross_sire.col(j)=M_cross_sire.col(j)-p_A[j];
	M_cross_dam.col(j)=M_cross_dam.col(j)-p_B[j];
	M_sire.col(j)=M_sire.col(j)-mean(M_sire.col(j));
	M_dam.col(j)=M_dam.col(j)-mean(M_dam.col(j));
	}
	double s1=sum(2*p_A%(1-p_A)),s2=sum(2*p_B%(1-p_B));
	//double s1=sum(2*p_A%p_B),s2=sum(2*p_A%p_B);
	arma::Mat<double> G_c_line_s=M_cross_sire*M_cross_sire.t()/s1;
	arma::Mat<double> G_c_line_d=M_cross_dam*M_cross_dam.t()/s2;
	
	arma::Mat<double> G_s=M_sire*M_sire.t()/s1;
	arma::Mat<double> G_d=M_dam*M_dam.t()/s2;
	
	arma::Mat<double> G_s_c_line_s=M_sire*M_cross_sire.t()/s1;
	arma::Mat<double> G_d_c_line_d=M_dam*M_cross_dam.t()/s2;
	
	arma::Mat<double> G_Sire(ind_sire_size+ind_cross_size,ind_sire_size+ind_cross_size,fill::zeros);
	arma::Mat<double> G_Dam(ind_dam_size+ind_cross_size,ind_dam_size+ind_cross_size,fill::zeros);

//G_Sire的大矩阵	
	arma::uvec sire_pos1=regspace<uvec>(0,1,ind_sire_size-1);
	arma::uvec sire_pos2=regspace<uvec>(ind_sire_size,1,ind_sire_size+ind_cross_size-1);
	G_Sire.submat(sire_pos1,sire_pos1)=G_s;
	G_Sire.submat(sire_pos2,sire_pos2)=G_c_line_s;
	G_Sire.submat(sire_pos1,sire_pos2)=G_s_c_line_s;
	G_Sire.submat(sire_pos2,sire_pos1)=G_s_c_line_s.t();
	G_Sire.diag()=G_Sire.diag()+0.0001;
//G_Dam的大矩阵	
	arma::uvec dam_pos1=regspace<uvec>(0,1,ind_dam_size-1);
	arma::uvec dam_pos2=regspace<uvec>(ind_dam_size,1,ind_dam_size+ind_cross_size-1);
	G_Dam.submat(dam_pos1,dam_pos1)=G_d;
	G_Dam.submat(dam_pos2,dam_pos2)=G_c_line_d;
	G_Dam.submat(dam_pos1,dam_pos2)=G_d_c_line_d;
	G_Dam.submat(dam_pos2,dam_pos1)=G_d_c_line_d.t();
	G_Dam.diag()=G_Dam.diag()+0.0001;

    return List::create(Named("G_c_line_s") = G_c_line_s,_["G_c_line_d"] = G_c_line_d,
	                    Named("G_s") = G_s,Named("G_d") = G_d,
						Named("G_s_c_line_s") = G_s_c_line_s,Named("G_d_c_line_d") = G_d_c_line_d,
						Named("G_Sire") = G_Sire,Named("G_Dam") = G_Dam,
						Named("M_sire") = t_sire,Named("M_dam") = t_dam);
	
}




// [[Rcpp::export]]
List makeHA_partial_cpp(arma::Mat<double> & P_A, arma::Mat<double> & G_A, 
						int n_pure,           //基因型数据中，有多少个纯种父亲/母亲
						int n_cross,          //基因型数据中，有多少个杂种个体
						CharacterVector IND_geno,
						arma::uvec pos_A11,  //有系谱但是无基因型
						arma::uvec pos_A22,  //既有基因型又有系谱
						arma::uvec pos_geno, //既有基因型又有系谱的个体在 基因型数据中的位置
						arma::uvec pos_A,    //H矩阵的个体  在系谱中的位置 
						arma::uvec pos_H22,  //#既有基因型又有系谱的个体在  在H矩阵中的位置
						bool direct=false,
						bool inverse=true,
						double omega=0.05){	

				G_A=G_A.submat(pos_geno,pos_geno);
				arma::Mat<double> H_A,H_Ainv;
				arma::Mat<double> K(pos_A22.size(),pos_A22.size(),fill::ones);
				//构建K矩阵
				arma::uvec sire_pos1=regspace<uvec>(0,1,n_pure-1);
				arma::uvec sire_pos2=regspace<uvec>(n_pure,1,n_pure+n_cross-1);
				//K.submat(sire_pos1,sire_pos1)=K.submat(sire_pos1,sire_pos1)*1;
	            K.submat(sire_pos2,sire_pos2)=K.submat(sire_pos2,sire_pos2)*0.25;
	            K.submat(sire_pos1,sire_pos2)=K.submat(sire_pos1,sire_pos2)*0.5;
	            K.submat(sire_pos2,sire_pos1)=K.submat(sire_pos2,sire_pos1)*0.5;


				Rcout<<"Start constructing H_Additive relationship matrix......"<<endl; 						   
				
				arma::Mat<double> A22=P_A.submat(pos_A22,pos_A22);
				arma::Mat<double> A22_inv=arma::inv(A22);			

						
				arma::Mat<double> A11=P_A.submat(pos_A11,pos_A11);
				arma::Mat<double> A12=P_A.submat(pos_A11,pos_A22);				
                double diag_GA=arma::mean(G_A.diag());
				double diag_A22=arma::mean(A22.diag());
				double diag_K=arma::mean(K.diag());
				
                double ave_GA=mean(mean(G_A));
				double ave_A22=mean(mean(A22));
				double ave_K=mean(mean(K));
				
                Rcout<<"The mean of diagnal and all elements of A22 is : "<<diag_A22<<" and "<<ave_A22<<endl;
				Rcout<<"The mean of diagnal and all elements of G_A is : "<<diag_GA<<" and "<<ave_GA<<endl;
				Rcout<<"The mean of diagnal and all elements of K is : "<<diag_K<<" and "<<ave_K<<endl;


				double alpha=(ave_A22*diag_GA-ave_GA*diag_A22)/(ave_K*diag_GA-ave_GA*diag_K);	
				double beta=(ave_A22-ave_K*alpha)/ave_GA;			
				
                Rcout<<"Blending G_A to P_A , the adjusting parameter is: alpha= "<<alpha<<" , beta= "<<beta<<endl;
                G_A=alpha*K+beta*G_A;		
				Rcout<<"Adjusting G_A , the adjusting parameter omega is: omega= "<<omega<<endl;
				G_A=G_A*(1-omega)+A22*omega;
				int ind_n1=pos_A11.size(),ind_n2=pos_A22.size(),ind_n=pos_A11.size()+pos_A22.size();
				if(direct==true){
				arma::Mat<double> tmp=A12*A22_inv*G_A;	
				H_A=P_A.submat(pos_A,pos_A);
				H_A.submat(span(0,ind_n1-1),span(0,ind_n1-1))=A11+A12*A22_inv*(G_A-A22)*A22_inv*(A12.t());
				H_A.submat(span(0,ind_n1-1),span(ind_n1,ind_n-1))=tmp;
				H_A.submat(span(ind_n1,ind_n-1),span(0,ind_n1-1))=tmp.t();
				H_A.submat(span(ind_n1,ind_n-1),span(ind_n1,ind_n-1))=G_A;
				
				}
				if(inverse==true){
				
				arma::Mat<double> G_Ainv;

				G_Ainv=arma::inv(G_A);

				arma::Mat<double> P_Ainv=inv(P_A);

				H_Ainv=P_Ainv.submat(pos_A,pos_A);
				H_Ainv.submat(pos_H22,pos_H22)=H_Ainv.submat(pos_H22,pos_H22)+G_Ainv-A22_inv;				
				}			
				
				//cout.precision(5);
				//cout.setf(ios::fixed);

				//H_A=H_A.print(cout, "H_A:");
				//H_Ainv.print(cout, "H_Ainv:");
				
				H_Ainv=round(H_Ainv*100000)/100000;
				H_A=round(H_A*100000)/100000;
				
				return List::create(Named("H") =H_A,_["Hinv"] = H_Ainv);
}






//计算杂交的A矩阵
// [[Rcpp::export]]
List makeAinv_old_partial_cpp(arma::mat Pedigree, 
					   std::vector<std::string> IND_name,  //重命名前的个体号，与Pedigree一一对应
					   IntegerVector record_pos){   //前三列为个体号，父亲，母亲，第四列为品种(1,2,未知品种为0)


	int n_ind=Pedigree.n_rows;

	int Sire_id,Dam_id;
	double A1_svalue,A1_dvalue,A2_svalue,A2_dvalue,A12_svalue,A12_dvalue;
	arma::mat A1(n_ind,n_ind),A2(n_ind,n_ind); //品种的亲缘矩阵
	arma::mat D1_inv(n_ind,n_ind),D2_inv(n_ind,n_ind); 
	arma::mat T1_inv(n_ind,n_ind,fill::eye),T2_inv(n_ind,n_ind,fill::eye); 

    arma::mat Breed_matrix(n_ind,2);	
	A1.fill(0);
	A2.fill(0);
	D1_inv.fill(0);
	D2_inv.fill(0);
	Breed_matrix.fill(0);

	
	arma::vec Breed=Pedigree.col(3);
	arma::vec Animal=Pedigree.col(0);
	arma::vec Sire=Pedigree.col(1);
	arma::vec Dam=Pedigree.col(2);

	
	for(int i=0; i < n_ind; i++){
		
		Sire_id=Sire[i];
		Dam_id=Dam[i];
		
		if(Sire_id==0&&Dam_id==0){  //父母均未知,
		
			Breed_matrix(i,Breed[i]-1)=1; //父母未知的情况下,此时品种必须是已知的
			
			if(Breed[i]==1){			
				A1(i,i)=1;
				D1_inv(i,i)=1;
	
			}else if(Breed[i]==2){				
				A2(i,i)=1;
				D2_inv(i,i)=1;
			}	
		}else{                     //父母均已知		
			Breed_matrix(i,0)=(Breed_matrix(Sire_id-1,0)+Breed_matrix(Dam_id-1,0))/2; //父母均已知,计算品种1的比例
			Breed_matrix(i,1)=(Breed_matrix(Sire_id-1,1)+Breed_matrix(Dam_id-1,1))/2; //父母均已知,计算品种1的比例
			
			A1(i,i)=Breed_matrix(i,0)+(A1(Sire_id-1,Dam_id-1))/2;	
			A2(i,i)=Breed_matrix(i,1)+(A2(Sire_id-1,Dam_id-1))/2;	

			if(Breed_matrix(i,0)==1){ //品种1			
				D1_inv(i,i)=1-(A1(Sire_id-1,Sire_id-1)+A1(Dam_id-1,Dam_id-1))/4;
				T1_inv(i,Sire_id-1)=-0.5;
				T1_inv(i,Dam_id-1)=-0.5;	
		
				
			}else if(Breed_matrix(i,0)==0){	//品种2			
				D2_inv(i,i)=1-(A2(Sire_id-1,Sire_id-1)+A2(Dam_id-1,Dam_id-1))/4;
				T2_inv(i,Sire_id-1)=-0.5;
				T2_inv(i,Dam_id-1)=-0.5;
				
			}else {	 //杂种数据
				D1_inv(i,i)=Breed_matrix(i,0)-A1(Sire_id-1,Sire_id-1)/4-A1(Dam_id-1,Dam_id-1)/4;
				D2_inv(i,i)=Breed_matrix(i,1)-A2(Dam_id-1,Dam_id-1)/4-A2(Sire_id-1,Sire_id-1)/4;
				T1_inv(i,Sire_id-1)=-0.5;
				T2_inv(i,Dam_id-1)=-0.5;
				T1_inv(i,Dam_id-1)=-0.5;
				T2_inv(i,Sire_id-1)=-0.5;
			}


				for(int j=0;j<i;j++){
				
					A1_svalue=A1(j,Sire_id-1);	
					A2_svalue=A2(j,Sire_id-1);	

					A1_dvalue=A1(j,Dam_id-1);	
					A2_dvalue=A2(j,Dam_id-1);	
					
				A1(j,i)=A1(i,j)=(A1_svalue+A1_dvalue)/2;	
				A2(j,i)=A2(i,j)=(A2_svalue+A2_dvalue)/2;				
				}

		}
	}


	
	//删除 父母中有一个缺失的个体
	if(record_pos.size()>0){
		
		for(int t=0;t<record_pos.size();t++){
			Breed_matrix(record_pos[t],0)=0;
			Breed_matrix(record_pos[t],1)=0;
		}
	}
	
	 std::vector<std::string>  id1=IND_name;
	 std::vector<std::string>  id2=IND_name;

	arma::vec breed1=Breed_matrix.col(0);
	arma::vec breed2=Breed_matrix.col(1);

	arma::uvec row_condition1=arma::find(breed1!=0);  //A1矩阵中为零的行
	arma::uvec row_condition2=arma::find(breed2!=0);  //A2矩阵中为零的行
	
	D1_inv=D1_inv.submat(row_condition1,row_condition1);
	T1_inv=T1_inv.submat(row_condition1,row_condition1);
	A1=A1.submat(row_condition1,row_condition1);
	
	D2_inv=D2_inv.submat(row_condition2,row_condition2);
	T2_inv=T2_inv.submat(row_condition2,row_condition2);
	A2=A2.submat(row_condition2,row_condition2);
	
	for(int k=breed1.size()-1;k>-1;k--){
		if(breed1[k]==0){
		id1.erase(id1.begin()+k);		
		}
		
		if(breed2[k]==0){
		id2.erase(id2.begin()+k);		
		}		
	}
	D1_inv.diag()=1/D1_inv.diag();

	
	
	return List::create(Named("A1_inv") = T1_inv.t()*D1_inv*T1_inv,
						Named("A2_inv") = T2_inv.t()*D2_inv*T2_inv,
						Named("id_Breed1_A")=id1,
						Named("id_Breed2_A")=id2,
						Named("A1")=A1,
						Named("A2")=A2); 
}


//include segreation
// [[Rcpp::export]]
List makeA_partial_cpp1(arma::mat Pedigree){   //前三列为个体号，父亲，母亲，第四列为品种(1,2,未知品种为0)


	int n_ind=Pedigree.n_rows;

	int Sire_id,Dam_id;
	double A1_svalue,A1_dvalue,A2_svalue,A2_dvalue,A3_svalue,A3_dvalue;
	arma::mat A1(n_ind,n_ind),A2(n_ind,n_ind),A3(n_ind,n_ind); //品种的亲缘矩阵

    arma::mat Breed_matrix(n_ind,2);	
	A1.fill(0);
	A2.fill(0);
	A3.fill(0);
	Breed_matrix.fill(0);

	
	arma::vec Breed=Pedigree.col(3);
	arma::vec Animal=Pedigree.col(0);
	arma::vec Sire=Pedigree.col(1);
	arma::vec Dam=Pedigree.col(2);

	
	for(int i=0; i < n_ind; i++){
		Sire_id=Sire[i];
		Dam_id=Dam[i];
		
		if(Sire_id==0&&Dam_id==0){  //父母均未知,
		
			Breed_matrix(i,Breed[i]-1)=1; //父母未知的情况下,此时品种必须是已知的
			
			if(Breed[i]==1){			
				A1(i,i)=1;	
	
			}else if(Breed[i]==2){				
				A2(i,i)=1;					
			}	
		}else{                     //父母均已知
			
			Breed_matrix(i,0)=(Breed_matrix(Sire_id-1,0)+Breed_matrix(Dam_id-1,0))/2; //父母均已知,计算品种1的比例
			
			Breed_matrix(i,1)=(Breed_matrix(Sire_id-1,1)+Breed_matrix(Dam_id-1,1))/2; //父母均已知,计算品种1的比例
			
			A1(i,i)=Breed_matrix(i,0)+(A1(Sire_id-1,Dam_id-1))/2;	
			A2(i,i)=Breed_matrix(i,1)+(A2(Sire_id-1,Dam_id-1))/2;	
			A3(i,i)=2*(Breed_matrix(Sire_id-1,0)*Breed_matrix(Sire_id-1,1)+Breed_matrix(Dam_id-1,0)*Breed_matrix(Dam_id-1,1))+(A3(Sire_id-1,Dam_id-1))/2;
			
				for(int j=0;j<i;j++){
				
					A1_svalue=A1(j,Sire_id-1);	
					A2_svalue=A2(j,Sire_id-1);	
					A3_svalue=A3(j,Sire_id-1);
					
					A1_dvalue=A1(j,Dam_id-1);	
					A2_dvalue=A2(j,Dam_id-1);	
					A3_dvalue=A3(j,Dam_id-1);
				
				A1(j,i)=A1(i,j)=(A1_svalue+A1_dvalue)/2;	
				A2(j,i)=A2(i,j)=(A2_svalue+A2_dvalue)/2;
				A3(j,i)=A3(i,j)=(A3_svalue+A3_dvalue)/2;				
				}
			
		}
	}
	


	
	return List::create(Named("Breed1_A") = A1,
						Named("Breed2_A") = A2,
						Named("Breed3_A") = A3,						
						Named("Breed")=Breed_matrix); 
}
