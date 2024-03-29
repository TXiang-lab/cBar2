# cBar2
Cross breeding analysis with R( two-way cross)
### OVERVIEW
<img src="https://user-images.githubusercontent.com/52848770/175520723-dcf2b3b1-29bd-4106-b08a-0f51a70d83c4.png"  align="right" height="250" style="zoom:50%  ;" />

`cBar2` is an useful and powerful tool for handling two-way cross breeding in animal and plant breeding.  In addition, it provides a  powerful function to do allele tracing for this case. 

By applying this package, user can complete  the genetic evaluation with  three different types of partial relationship matrix(pedigree, genomic, and single-step partial relationship matrix). 

😊 Good Luck Charlie !   If you have suggestion or question, please contact: [quanshun1994@gmail.com](mailto:quanshun1994@gmail.com) !

## GETTING STARTED

### 🙊Installation

`cBar2` links to R packages `Rcpp`, `RcppArmadillo` , `data.table` ,  `bigmemory`  .  These dependencies should be installed before installing `cBar2`.  

```R
install.packages(c("Rcpp", "RcppArmadillo","data.table","bigmemory","RcppProgress","BH"))
```


**👉 Note: In the analysis of DMU  and BLUPF90 , we recommend you to use package blupADC, it's super easy and effective !!!** 

#### Install cBar2

```R
devtools::install_github("TXiang-lab/cBar2")
```

After installed successfully, the `cBar2` package can be loaded by typing

``` {.r}
library(cBar2)
```

**Note**: In terms of the relationship matrix construction, we highly recommend Microsoft R Open(faster than traditional R many times)

### 🙊Features

-   Feature 1. Allele tracing 
-   Feature 2. Pedigree partial relationship matrix construction 
-   Feature 3. Genomic partial relationship matrix construction 
-   Feature 4. Single-step genomic partial relationship matrix construction 

## Usage

`cBar2` provides several datasets objects inside this package.

#### Feature 1. Allele tracing 

``` R
library(cBar2)
offspring_boa=allele_tracing(input_pedigree=ped,
                             hap_win=50,
                             cpu_cores=16,
			     haplotype_hap,
                             haplotype_map,
                             haplotype_sample)
```

#### Feature 2. Pedigree partial relationship matrix construction 

``` R
library(cBar2)
A_partial=makeA_partial(ped,
			output_matrix_type="col_three")
Ainv_partial=makeAinv_partial(ped,
			output_matrix_type="col_three")                      
```

#### Feature 3. Genomic partial relationship matrix construction 
Note: need to perform allele tracing in advance 

``` R
library(cBar2)
G_parital=makeGA_partial(ped,
                         offspring_boa,#obtained by allele tracing function 
                         haplotype_hap,
                         haplotype_sample,
			 output_matrix_type="col_three")
			 
Ginv_parital=makeGAinv_partial(ped,
                               offspring_boa, #obtained by allele tracing function 
                               haplotype_hap,
                               haplotype_sample,
			       output_matrix_type="col_three")
```

#### Feature 4. Single-step genomic partial relationship matrix construction 
Note: need to perform allele tracing in advance 
``` R
library(cBar2)
H_parital=makeHA_partial(ped,
                         offspring_boa,#obtained by allele tracing function 
                         haplotype_hap,
                         haplotype_sample,
			 output_matrix_type="col_three")
Hinv_parital=makeHAinv_partial(ped,
                               offspring_boa, #obtained by allele tracing function 
                               haplotype_hap,
                               haplotype_sample,
			       output_matrix_type="col_three")
```


