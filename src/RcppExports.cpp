// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// allele_tracing_cpp
arma::mat allele_tracing_cpp(arma::mat& sire_hap, arma::mat& dam_hap, arma::mat& cross_hap, arma::mat& chr_pos, arma::mat& block_pos, CharacterMatrix Father_cross_ped, CharacterMatrix Mother_cross_ped, CharacterMatrix Unknow_cross_ped, std::vector<std::string> sire_hap_ind, std::vector<std::string> dam_hap_ind, std::vector<std::string> cross_hap_ind, int hap_win, int cpu_cores);
RcppExport SEXP _cBar2_allele_tracing_cpp(SEXP sire_hapSEXP, SEXP dam_hapSEXP, SEXP cross_hapSEXP, SEXP chr_posSEXP, SEXP block_posSEXP, SEXP Father_cross_pedSEXP, SEXP Mother_cross_pedSEXP, SEXP Unknow_cross_pedSEXP, SEXP sire_hap_indSEXP, SEXP dam_hap_indSEXP, SEXP cross_hap_indSEXP, SEXP hap_winSEXP, SEXP cpu_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type sire_hap(sire_hapSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type dam_hap(dam_hapSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type cross_hap(cross_hapSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type chr_pos(chr_posSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type block_pos(block_posSEXP);
    Rcpp::traits::input_parameter< CharacterMatrix >::type Father_cross_ped(Father_cross_pedSEXP);
    Rcpp::traits::input_parameter< CharacterMatrix >::type Mother_cross_ped(Mother_cross_pedSEXP);
    Rcpp::traits::input_parameter< CharacterMatrix >::type Unknow_cross_ped(Unknow_cross_pedSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type sire_hap_ind(sire_hap_indSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type dam_hap_ind(dam_hap_indSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type cross_hap_ind(cross_hap_indSEXP);
    Rcpp::traits::input_parameter< int >::type hap_win(hap_winSEXP);
    Rcpp::traits::input_parameter< int >::type cpu_cores(cpu_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(allele_tracing_cpp(sire_hap, dam_hap, cross_hap, chr_pos, block_pos, Father_cross_ped, Mother_cross_ped, Unknow_cross_ped, sire_hap_ind, dam_hap_ind, cross_hap_ind, hap_win, cpu_cores));
    return rcpp_result_gen;
END_RCPP
}
// makeA_partial_cpp
List makeA_partial_cpp(arma::mat Pedigree, std::vector<std::string> IND_name, IntegerVector record_pos, bool full_rank);
RcppExport SEXP _cBar2_makeA_partial_cpp(SEXP PedigreeSEXP, SEXP IND_nameSEXP, SEXP record_posSEXP, SEXP full_rankSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Pedigree(PedigreeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type IND_name(IND_nameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type record_pos(record_posSEXP);
    Rcpp::traits::input_parameter< bool >::type full_rank(full_rankSEXP);
    rcpp_result_gen = Rcpp::wrap(makeA_partial_cpp(Pedigree, IND_name, record_pos, full_rank));
    return rcpp_result_gen;
END_RCPP
}
// makeAinv_partial_cpp
List makeAinv_partial_cpp(arma::mat Pedigree, std::vector<std::string> IND_name, IntegerVector record_pos, bool full_rank);
RcppExport SEXP _cBar2_makeAinv_partial_cpp(SEXP PedigreeSEXP, SEXP IND_nameSEXP, SEXP record_posSEXP, SEXP full_rankSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Pedigree(PedigreeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type IND_name(IND_nameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type record_pos(record_posSEXP);
    Rcpp::traits::input_parameter< bool >::type full_rank(full_rankSEXP);
    rcpp_result_gen = Rcpp::wrap(makeAinv_partial_cpp(Pedigree, IND_name, record_pos, full_rank));
    return rcpp_result_gen;
END_RCPP
}
// G_matrix_partial_cpp
List G_matrix_partial_cpp(arma::Mat<int>& cross_hap, arma::Mat<int>& cross_boa, arma::Mat<int>& sire_hap, arma::Mat<int>& dam_hap);
RcppExport SEXP _cBar2_G_matrix_partial_cpp(SEXP cross_hapSEXP, SEXP cross_boaSEXP, SEXP sire_hapSEXP, SEXP dam_hapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int>& >::type cross_hap(cross_hapSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int>& >::type cross_boa(cross_boaSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int>& >::type sire_hap(sire_hapSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int>& >::type dam_hap(dam_hapSEXP);
    rcpp_result_gen = Rcpp::wrap(G_matrix_partial_cpp(cross_hap, cross_boa, sire_hap, dam_hap));
    return rcpp_result_gen;
END_RCPP
}
// makeHA_partial_cpp
List makeHA_partial_cpp(arma::Mat<double>& P_A, arma::Mat<double>& G_A, int n_pure, int n_cross, CharacterVector IND_geno, arma::uvec pos_A11, arma::uvec pos_A22, arma::uvec pos_geno, arma::uvec pos_A, arma::uvec pos_H22, bool direct, bool inverse, double omega);
RcppExport SEXP _cBar2_makeHA_partial_cpp(SEXP P_ASEXP, SEXP G_ASEXP, SEXP n_pureSEXP, SEXP n_crossSEXP, SEXP IND_genoSEXP, SEXP pos_A11SEXP, SEXP pos_A22SEXP, SEXP pos_genoSEXP, SEXP pos_ASEXP, SEXP pos_H22SEXP, SEXP directSEXP, SEXP inverseSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<double>& >::type P_A(P_ASEXP);
    Rcpp::traits::input_parameter< arma::Mat<double>& >::type G_A(G_ASEXP);
    Rcpp::traits::input_parameter< int >::type n_pure(n_pureSEXP);
    Rcpp::traits::input_parameter< int >::type n_cross(n_crossSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type IND_geno(IND_genoSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos_A11(pos_A11SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos_A22(pos_A22SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos_geno(pos_genoSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos_A(pos_ASEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos_H22(pos_H22SEXP);
    Rcpp::traits::input_parameter< bool >::type direct(directSEXP);
    Rcpp::traits::input_parameter< bool >::type inverse(inverseSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(makeHA_partial_cpp(P_A, G_A, n_pure, n_cross, IND_geno, pos_A11, pos_A22, pos_geno, pos_A, pos_H22, direct, inverse, omega));
    return rcpp_result_gen;
END_RCPP
}
// makeAinv_old_partial_cpp
List makeAinv_old_partial_cpp(arma::mat Pedigree, std::vector<std::string> IND_name, IntegerVector record_pos);
RcppExport SEXP _cBar2_makeAinv_old_partial_cpp(SEXP PedigreeSEXP, SEXP IND_nameSEXP, SEXP record_posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Pedigree(PedigreeSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type IND_name(IND_nameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type record_pos(record_posSEXP);
    rcpp_result_gen = Rcpp::wrap(makeAinv_old_partial_cpp(Pedigree, IND_name, record_pos));
    return rcpp_result_gen;
END_RCPP
}
// makeA_partial_cpp1
List makeA_partial_cpp1(arma::mat Pedigree);
RcppExport SEXP _cBar2_makeA_partial_cpp1(SEXP PedigreeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Pedigree(PedigreeSEXP);
    rcpp_result_gen = Rcpp::wrap(makeA_partial_cpp1(Pedigree));
    return rcpp_result_gen;
END_RCPP
}
// cumulativeSum
void cumulativeSum(const std::vector<int> input, std::vector<int>& result);
RcppExport SEXP _cBar2_cumulativeSum(SEXP inputSEXP, SEXP resultSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int> >::type input(inputSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type result(resultSEXP);
    cumulativeSum(input, result);
    return R_NilValue;
END_RCPP
}
// paste_vec_short
std::string paste_vec_short(arma::Col<int>& data);
RcppExport SEXP _cBar2_paste_vec_short(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Col<int>& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(paste_vec_short(data));
    return rcpp_result_gen;
END_RCPP
}
// get_haplotype_set_short
std::vector<int> get_haplotype_set_short(arma::Mat<int>& data_haplotype);
RcppExport SEXP _cBar2_get_haplotype_set_short(SEXP data_haplotypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int>& >::type data_haplotype(data_haplotypeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_haplotype_set_short(data_haplotype));
    return rcpp_result_gen;
END_RCPP
}
// allele_get_haplotype_set_short
void allele_get_haplotype_set_short(arma::Mat<int>& data_haplotype, Rcpp::List& haplotype_allele, int i_pos);
RcppExport SEXP _cBar2_allele_get_haplotype_set_short(SEXP data_haplotypeSEXP, SEXP haplotype_alleleSEXP, SEXP i_posSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int>& >::type data_haplotype(data_haplotypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type haplotype_allele(haplotype_alleleSEXP);
    Rcpp::traits::input_parameter< int >::type i_pos(i_posSEXP);
    allele_get_haplotype_set_short(data_haplotype, haplotype_allele, i_pos);
    return R_NilValue;
END_RCPP
}
// define_block_window_kb_cpp
List define_block_window_kb_cpp(std::vector<int> block1, std::vector<int> block2, std::vector<int> tmp_data_map_3);
RcppExport SEXP _cBar2_define_block_window_kb_cpp(SEXP block1SEXP, SEXP block2SEXP, SEXP tmp_data_map_3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type block1(block1SEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type block2(block2SEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type tmp_data_map_3(tmp_data_map_3SEXP);
    rcpp_result_gen = Rcpp::wrap(define_block_window_kb_cpp(block1, block2, tmp_data_map_3));
    return rcpp_result_gen;
END_RCPP
}
// single_base_factor_cpp
std::string single_base_factor_cpp(CharacterVector data, std::string miss_base);
RcppExport SEXP _cBar2_single_base_factor_cpp(SEXP dataSEXP, SEXP miss_baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::string >::type miss_base(miss_baseSEXP);
    rcpp_result_gen = Rcpp::wrap(single_base_factor_cpp(data, miss_base));
    return rcpp_result_gen;
END_RCPP
}
// pair_base_factor_cpp
std::string pair_base_factor_cpp(CharacterVector data, std::string miss_base);
RcppExport SEXP _cBar2_pair_base_factor_cpp(SEXP dataSEXP, SEXP miss_baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::string >::type miss_base(miss_baseSEXP);
    rcpp_result_gen = Rcpp::wrap(pair_base_factor_cpp(data, miss_base));
    return rcpp_result_gen;
END_RCPP
}
// get_blupf90_allele_string_phased_vcf
std::string get_blupf90_allele_string_phased_vcf(arma::Col<int> allele_string_col1, arma::Col<int> allele_string_col2, int max_length, std::vector<int> cumsum_haplo_type_num, int n_snp);
RcppExport SEXP _cBar2_get_blupf90_allele_string_phased_vcf(SEXP allele_string_col1SEXP, SEXP allele_string_col2SEXP, SEXP max_lengthSEXP, SEXP cumsum_haplo_type_numSEXP, SEXP n_snpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Col<int> >::type allele_string_col1(allele_string_col1SEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type allele_string_col2(allele_string_col2SEXP);
    Rcpp::traits::input_parameter< int >::type max_length(max_lengthSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type cumsum_haplo_type_num(cumsum_haplo_type_numSEXP);
    Rcpp::traits::input_parameter< int >::type n_snp(n_snpSEXP);
    rcpp_result_gen = Rcpp::wrap(get_blupf90_allele_string_phased_vcf(allele_string_col1, allele_string_col2, max_length, cumsum_haplo_type_num, n_snp));
    return rcpp_result_gen;
END_RCPP
}
// get_blupf90_allele_string_numeric
std::string get_blupf90_allele_string_numeric(arma::Row<int> allele_string_row, int max_length);
RcppExport SEXP _cBar2_get_blupf90_allele_string_numeric(SEXP allele_string_rowSEXP, SEXP max_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Row<int> >::type allele_string_row(allele_string_rowSEXP);
    Rcpp::traits::input_parameter< int >::type max_length(max_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(get_blupf90_allele_string_numeric(allele_string_row, max_length));
    return rcpp_result_gen;
END_RCPP
}
// make_bigmemory_object_cpp
SEXP make_bigmemory_object_cpp(int nrow, int ncol, std::string file_name, std::string file_path, std::string type);
RcppExport SEXP _cBar2_make_bigmemory_object_cpp(SEXP nrowSEXP, SEXP ncolSEXP, SEXP file_nameSEXP, SEXP file_pathSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< std::string >::type file_name(file_nameSEXP);
    Rcpp::traits::input_parameter< std::string >::type file_path(file_pathSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(make_bigmemory_object_cpp(nrow, ncol, file_name, file_path, type));
    return rcpp_result_gen;
END_RCPP
}
// make_bigmemory_object_address_cpp
SEXP make_bigmemory_object_address_cpp(int nrow, int ncol, std::string file_name, std::string file_path, std::string type);
RcppExport SEXP _cBar2_make_bigmemory_object_address_cpp(SEXP nrowSEXP, SEXP ncolSEXP, SEXP file_nameSEXP, SEXP file_pathSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< std::string >::type file_name(file_nameSEXP);
    Rcpp::traits::input_parameter< std::string >::type file_path(file_pathSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(make_bigmemory_object_address_cpp(nrow, ncol, file_name, file_path, type));
    return rcpp_result_gen;
END_RCPP
}
// NumericMatrix_to_arma
arma::Mat<int> NumericMatrix_to_arma(NumericMatrix& data_numeric);
RcppExport SEXP _cBar2_NumericMatrix_to_arma(SEXP data_numericSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type data_numeric(data_numericSEXP);
    rcpp_result_gen = Rcpp::wrap(NumericMatrix_to_arma(data_numeric));
    return rcpp_result_gen;
END_RCPP
}
// DataFrame_to_arma
arma::Mat<int> DataFrame_to_arma(DataFrame& data_numeric);
RcppExport SEXP _cBar2_DataFrame_to_arma(SEXP data_numericSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type data_numeric(data_numericSEXP);
    rcpp_result_gen = Rcpp::wrap(DataFrame_to_arma(data_numeric));
    return rcpp_result_gen;
END_RCPP
}
// get_offspring_generation_cpp
IntegerVector get_offspring_generation_cpp(DataFrame ped, CharacterVector IND_base);
RcppExport SEXP _cBar2_get_offspring_generation_cpp(SEXP pedSEXP, SEXP IND_baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type ped(pedSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type IND_base(IND_baseSEXP);
    rcpp_result_gen = Rcpp::wrap(get_offspring_generation_cpp(ped, IND_base));
    return rcpp_result_gen;
END_RCPP
}
// single_pedigree_cpp
List single_pedigree_cpp(CharacterMatrix ped);
RcppExport SEXP _cBar2_single_pedigree_cpp(SEXP pedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterMatrix >::type ped(pedSEXP);
    rcpp_result_gen = Rcpp::wrap(single_pedigree_cpp(ped));
    return rcpp_result_gen;
END_RCPP
}
// get_blupf90_allele_string_unphased_haplotype
std::string get_blupf90_allele_string_unphased_haplotype(arma::Col<int> allele_string_col1, arma::Col<int> allele_string_col2, int max_length);
RcppExport SEXP _cBar2_get_blupf90_allele_string_unphased_haplotype(SEXP allele_string_col1SEXP, SEXP allele_string_col2SEXP, SEXP max_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Col<int> >::type allele_string_col1(allele_string_col1SEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type allele_string_col2(allele_string_col2SEXP);
    Rcpp::traits::input_parameter< int >::type max_length(max_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(get_blupf90_allele_string_unphased_haplotype(allele_string_col1, allele_string_col2, max_length));
    return rcpp_result_gen;
END_RCPP
}
// union_cpp
CharacterVector union_cpp(CharacterVector X, CharacterVector Y);
RcppExport SEXP _cBar2_union_cpp(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(union_cpp(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// delete_bigmemory_file_cpp
void delete_bigmemory_file_cpp(std::string matrix_type, std::string bigmemory_data_name, std::string bigmemory_data_path, bool message);
RcppExport SEXP _cBar2_delete_bigmemory_file_cpp(SEXP matrix_typeSEXP, SEXP bigmemory_data_nameSEXP, SEXP bigmemory_data_pathSEXP, SEXP messageSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type matrix_type(matrix_typeSEXP);
    Rcpp::traits::input_parameter< std::string >::type bigmemory_data_name(bigmemory_data_nameSEXP);
    Rcpp::traits::input_parameter< std::string >::type bigmemory_data_path(bigmemory_data_pathSEXP);
    Rcpp::traits::input_parameter< bool >::type message(messageSEXP);
    delete_bigmemory_file_cpp(matrix_type, bigmemory_data_name, bigmemory_data_path, message);
    return R_NilValue;
END_RCPP
}
// get_allele
arma::Mat<int> get_allele(arma::Mat<int> Pedigree);
RcppExport SEXP _cBar2_get_allele(SEXP PedigreeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type Pedigree(PedigreeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_allele(Pedigree));
    return rcpp_result_gen;
END_RCPP
}
// matrix_col3_old
arma::Mat<double> matrix_col3_old(arma::Mat<double>& G, arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);
RcppExport SEXP _cBar2_matrix_col3_old(SEXP GSEXP, SEXP IND_genoSEXP, SEXP detSEXP, SEXP cpu_coresSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<double>& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type IND_geno(IND_genoSEXP);
    Rcpp::traits::input_parameter< bool >::type det(detSEXP);
    Rcpp::traits::input_parameter< int >::type cpu_cores(cpu_coresSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_col3_old(G, IND_geno, det, cpu_cores, threshold));
    return rcpp_result_gen;
END_RCPP
}
// matrix_col3_memory_old
SEXP matrix_col3_memory_old(SEXP pBigMat, std::string bigmemory_data_name, std::string bigmemory_data_path, arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);
RcppExport SEXP _cBar2_matrix_col3_memory_old(SEXP pBigMatSEXP, SEXP bigmemory_data_nameSEXP, SEXP bigmemory_data_pathSEXP, SEXP IND_genoSEXP, SEXP detSEXP, SEXP cpu_coresSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< std::string >::type bigmemory_data_name(bigmemory_data_nameSEXP);
    Rcpp::traits::input_parameter< std::string >::type bigmemory_data_path(bigmemory_data_pathSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type IND_geno(IND_genoSEXP);
    Rcpp::traits::input_parameter< bool >::type det(detSEXP);
    Rcpp::traits::input_parameter< int >::type cpu_cores(cpu_coresSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_col3_memory_old(pBigMat, bigmemory_data_name, bigmemory_data_path, IND_geno, det, cpu_cores, threshold));
    return rcpp_result_gen;
END_RCPP
}
// matrix_col3_memory_alt
SEXP matrix_col3_memory_alt(arma::Mat<double> G, std::string bigmemory_data_name, std::string bigmemory_data_path, arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);
RcppExport SEXP _cBar2_matrix_col3_memory_alt(SEXP GSEXP, SEXP bigmemory_data_nameSEXP, SEXP bigmemory_data_pathSEXP, SEXP IND_genoSEXP, SEXP detSEXP, SEXP cpu_coresSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<double> >::type G(GSEXP);
    Rcpp::traits::input_parameter< std::string >::type bigmemory_data_name(bigmemory_data_nameSEXP);
    Rcpp::traits::input_parameter< std::string >::type bigmemory_data_path(bigmemory_data_pathSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type IND_geno(IND_genoSEXP);
    Rcpp::traits::input_parameter< bool >::type det(detSEXP);
    Rcpp::traits::input_parameter< int >::type cpu_cores(cpu_coresSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_col3_memory_alt(G, bigmemory_data_name, bigmemory_data_path, IND_geno, det, cpu_cores, threshold));
    return rcpp_result_gen;
END_RCPP
}
// matrix_col3
arma::Mat<double> matrix_col3(arma::Mat<double>& G, arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);
RcppExport SEXP _cBar2_matrix_col3(SEXP GSEXP, SEXP IND_genoSEXP, SEXP detSEXP, SEXP cpu_coresSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<double>& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type IND_geno(IND_genoSEXP);
    Rcpp::traits::input_parameter< bool >::type det(detSEXP);
    Rcpp::traits::input_parameter< int >::type cpu_cores(cpu_coresSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_col3(G, IND_geno, det, cpu_cores, threshold));
    return rcpp_result_gen;
END_RCPP
}
// matrix_col3_memory
SEXP matrix_col3_memory(SEXP pBigMat, std::string bigmemory_data_name, std::string bigmemory_data_path, arma::Col<int> IND_geno, bool det, int cpu_cores, double threshold);
RcppExport SEXP _cBar2_matrix_col3_memory(SEXP pBigMatSEXP, SEXP bigmemory_data_nameSEXP, SEXP bigmemory_data_pathSEXP, SEXP IND_genoSEXP, SEXP detSEXP, SEXP cpu_coresSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< std::string >::type bigmemory_data_name(bigmemory_data_nameSEXP);
    Rcpp::traits::input_parameter< std::string >::type bigmemory_data_path(bigmemory_data_pathSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type IND_geno(IND_genoSEXP);
    Rcpp::traits::input_parameter< bool >::type det(detSEXP);
    Rcpp::traits::input_parameter< int >::type cpu_cores(cpu_coresSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_col3_memory(pBigMat, bigmemory_data_name, bigmemory_data_path, IND_geno, det, cpu_cores, threshold));
    return rcpp_result_gen;
END_RCPP
}
// cal_breed_pro
arma::mat cal_breed_pro(arma::mat Pedigree);
RcppExport SEXP _cBar2_cal_breed_pro(SEXP PedigreeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Pedigree(PedigreeSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_breed_pro(Pedigree));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cBar2_allele_tracing_cpp", (DL_FUNC) &_cBar2_allele_tracing_cpp, 13},
    {"_cBar2_makeA_partial_cpp", (DL_FUNC) &_cBar2_makeA_partial_cpp, 4},
    {"_cBar2_makeAinv_partial_cpp", (DL_FUNC) &_cBar2_makeAinv_partial_cpp, 4},
    {"_cBar2_G_matrix_partial_cpp", (DL_FUNC) &_cBar2_G_matrix_partial_cpp, 4},
    {"_cBar2_makeHA_partial_cpp", (DL_FUNC) &_cBar2_makeHA_partial_cpp, 13},
    {"_cBar2_makeAinv_old_partial_cpp", (DL_FUNC) &_cBar2_makeAinv_old_partial_cpp, 3},
    {"_cBar2_makeA_partial_cpp1", (DL_FUNC) &_cBar2_makeA_partial_cpp1, 1},
    {"_cBar2_cumulativeSum", (DL_FUNC) &_cBar2_cumulativeSum, 2},
    {"_cBar2_paste_vec_short", (DL_FUNC) &_cBar2_paste_vec_short, 1},
    {"_cBar2_get_haplotype_set_short", (DL_FUNC) &_cBar2_get_haplotype_set_short, 1},
    {"_cBar2_allele_get_haplotype_set_short", (DL_FUNC) &_cBar2_allele_get_haplotype_set_short, 3},
    {"_cBar2_define_block_window_kb_cpp", (DL_FUNC) &_cBar2_define_block_window_kb_cpp, 3},
    {"_cBar2_single_base_factor_cpp", (DL_FUNC) &_cBar2_single_base_factor_cpp, 2},
    {"_cBar2_pair_base_factor_cpp", (DL_FUNC) &_cBar2_pair_base_factor_cpp, 2},
    {"_cBar2_get_blupf90_allele_string_phased_vcf", (DL_FUNC) &_cBar2_get_blupf90_allele_string_phased_vcf, 5},
    {"_cBar2_get_blupf90_allele_string_numeric", (DL_FUNC) &_cBar2_get_blupf90_allele_string_numeric, 2},
    {"_cBar2_make_bigmemory_object_cpp", (DL_FUNC) &_cBar2_make_bigmemory_object_cpp, 5},
    {"_cBar2_make_bigmemory_object_address_cpp", (DL_FUNC) &_cBar2_make_bigmemory_object_address_cpp, 5},
    {"_cBar2_NumericMatrix_to_arma", (DL_FUNC) &_cBar2_NumericMatrix_to_arma, 1},
    {"_cBar2_DataFrame_to_arma", (DL_FUNC) &_cBar2_DataFrame_to_arma, 1},
    {"_cBar2_get_offspring_generation_cpp", (DL_FUNC) &_cBar2_get_offspring_generation_cpp, 2},
    {"_cBar2_single_pedigree_cpp", (DL_FUNC) &_cBar2_single_pedigree_cpp, 1},
    {"_cBar2_get_blupf90_allele_string_unphased_haplotype", (DL_FUNC) &_cBar2_get_blupf90_allele_string_unphased_haplotype, 3},
    {"_cBar2_union_cpp", (DL_FUNC) &_cBar2_union_cpp, 2},
    {"_cBar2_delete_bigmemory_file_cpp", (DL_FUNC) &_cBar2_delete_bigmemory_file_cpp, 4},
    {"_cBar2_get_allele", (DL_FUNC) &_cBar2_get_allele, 1},
    {"_cBar2_matrix_col3_old", (DL_FUNC) &_cBar2_matrix_col3_old, 5},
    {"_cBar2_matrix_col3_memory_old", (DL_FUNC) &_cBar2_matrix_col3_memory_old, 7},
    {"_cBar2_matrix_col3_memory_alt", (DL_FUNC) &_cBar2_matrix_col3_memory_alt, 7},
    {"_cBar2_matrix_col3", (DL_FUNC) &_cBar2_matrix_col3, 5},
    {"_cBar2_matrix_col3_memory", (DL_FUNC) &_cBar2_matrix_col3_memory, 7},
    {"_cBar2_cal_breed_pro", (DL_FUNC) &_cBar2_cal_breed_pro, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_cBar2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
