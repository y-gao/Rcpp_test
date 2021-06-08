//rcpp_get_cl_means.cpp

// [[Rcpp::depends(beachmat)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include "beachmat3/beachmat.h"
#include <RcppParallel.h>

#include <cmath>
#include <algorithm>
#include <unordered_map>
using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::plugins("cpp11")]]

//////////////////////////////////////////////////////////////////////////////////
// method1: iterate over only the non-zero elements in a column in a sparse matrix
//////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
Rcpp::NumericMatrix rcpp_sparse_get_cl_means(Rcpp::RObject mat, Rcpp::IntegerVector cl){	
	/* Set up hash map for column cluster lookup: mat_cluster_vector is the column number in the res matrix for whole dataset, 
	     cluster_col_map is for each cluster
	   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
	   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/

	std::vector<int> mat_cluster_vector(cl.length(), 0);// mat_cluster_vector has same length as cl
	std::unordered_map<int, int> cluster_col_map;//cluster_col_map is hash map used to find mat_cluster_vector
	int col_id = 0;//initialize 
	CharacterVector levs = cl.attr("levels");//initialize levs
	CharacterVector colnames_vec;//initialize colnames_vec
	for(int i = 0; i < cl.length(); i++){//loop through cl
	    // search cl[i] in hashmap cluster_col_map, cluster_col_map.find(cl[i]) returns a pointer to cl[i], if it does not exist, it will point to the end
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end()) {
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(levs[cl[i] - 1]);
		} else {
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}	
	
	
    /* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for(int i = 0; i < cl.length(); i++){
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}
	
	/* Loop though sparse matrix and access non-zero entry
	   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
	   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	// read mat as sparse Matrix block and build one pointer to the sparse matrix
	// (col, row, val): 00000000000000100000010 -> C++ sparse matrix (0000000000 0000001 00000010)
	auto ptr = beachmat::read_lin_sparse_block(mat);// ptr is pointer to mat
	NumericMatrix res(ptr->get_nrow(), unique_clusters);//initialize res matrix to store result
	
	std::vector<double> workspace_x(ptr->get_nrow());//initialize a vector workspace_x with length= number of rows:for value
	
    std::vector<int> workspace_i(ptr->get_nrow());// initialize a vector workspace_i with length= number of rows:for index
	for (int col_i = 0; col_i < ptr->get_ncol(); col_i++){
		// loop over col_i-th column of mat matrix
		auto indices = ptr->get_col(col_i, workspace_x.data(), workspace_i.data());//returns a pointer to a vector of data values for each requested column col_i,a pointer to a structure including non-zero row value and non-zero row index
		auto xptr = indices.x; // pointer to value of all non-zero entry in each column
        auto iptr = indices.i; // pointer to row number of all non-zero entry in each column
        auto nnzero = indices.n; // number of non-zero entry
		
		// loop through non zero in each column
		for (int row_i = 0; row_i < nnzero; row_i++){
			res(*(iptr+row_i), mat_cluster_vector[col_i]) += *(xptr+row_i);
		}
	}
	
	
	// loop over columns of res matrix, and divide by number of cells in each cluster
	 for (int i = 0; i < res.ncol(); i++){
	 	NumericMatrix::Column col = res( _ , i);
	 	col = col / cluster_id_number[i];
	    }
		
	colnames(res) = colnames_vec;
	return res;	
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// method2: iterate over only the non-zero elements in a column in a sparse matrix, using RcppParallel
//////////////////////////////////////////////////////////////////////////////////////////////////////


/*To use parallelFor, create a Worker object that defines an operator() which is called by the parallel scheduler.*/

struct ColumnMeanSparse : public Worker
{
  // source matrix and maps
  beachmat::lin_sparse_matrix* matrix_ptr;
    
  std::vector<int> mat_cluster_vector;
  
  // destination matrix
  RMatrix<double> output;
  
  // initialize with source, cl_id and destination
  ColumnMeanSparse(beachmat::lin_sparse_matrix* matrix_ptr, 
             std::vector<int> mat_cluster_vector,
             NumericMatrix output) 
    : matrix_ptr(matrix_ptr), mat_cluster_vector(mat_cluster_vector), output(output) {}
  
  // get the column sum of the given cols
  void operator()(std::size_t begin, std::size_t end) {
	std::vector<double> workspace_x(matrix_ptr->get_nrow());
    std::vector<int> workspace_i(matrix_ptr->get_nrow());
    for (int col_i = begin; col_i < end; col_i++){	  
		// loop over col_i-th column of mat matrix
		auto indices = matrix_ptr->get_col(col_i, workspace_x.data(), workspace_i.data());
		auto xptr = indices.x; // row values
        auto iptr = indices.i; // row indices
        auto nnzero = indices.n;
		
		for (int row_i = 0; row_i < nnzero; row_i++){
			output(*(iptr+row_i), mat_cluster_vector[col_i]) += *(xptr+row_i);
		}
    }
  }
};


/*The function calls the ColumnMeanSparse*/
//[[Rcpp::export]]
NumericMatrix rcpp_sparse_get_cl_means_parallel(Rcpp::RObject mat, IntegerVector cl){	
  /* Set up hash map for column number lookup: cluster_col_map is for each cluster
   mat_cluster_vector is a vector containing the column number for each column in mat corresponding to res
   levs is the level of the cluster id: "1" "2" "3" "4" "5" "6" "7" "8" "9" "23"
   colnames_vec is the ordered cluster id: "1" "23" "2" "3" "9" "5" "6" "4" "7" "8"*/
    std::vector<int> mat_cluster_vector(cl.length(), 0);
	std::unordered_map<int, int> cluster_col_map;
	int col_id = 0;
	CharacterVector levs = cl.attr("levels");
	CharacterVector colnames_vec;
	for(int i = 0; i < cl.length(); i++){
		if (cluster_col_map.find(cl[i]) == cluster_col_map.end()) {
			cluster_col_map[cl[i]] = col_id;
			mat_cluster_vector[i] = col_id;
			col_id++;
			colnames_vec.push_back(levs[cl[i] - 1]);
		} else {
			mat_cluster_vector[i] = cluster_col_map[cl[i]];
		}
	}	
	
    /* Count cells in each cluster 
	   unique_clusters = number of clusters = 10 
	   cluster_id_number = number of cells in each cluster correspond to colnames_vec */
	int unique_clusters = cluster_col_map.size();
	std::vector<int> cluster_id_number(unique_clusters, 0);
	for(int i = 0; i < cl.length(); i++){
		cluster_id_number[cluster_col_map[cl[i]]]++;
	}
  
  /* Loop though sparse matrix and access non-zero entry
   res is the matrix to store sum valuse in each cluster, dim=number of genes*number of clusters
   loop through each column of the sparse matrix mat. For each non-zero entry, find its cluster and column in res, add into the value*/
	auto ptr = beachmat::read_lin_sparse_block(mat);
	NumericMatrix res(ptr->get_nrow(), unique_clusters);
  
  // column mean sparse functor (pass input and output matrixes)
  auto matrix_ptr = ptr.get();
  ColumnMeanSparse column_mean_sparse(matrix_ptr, mat_cluster_vector, res);
  
  // call parallelFor to do the work
  parallelFor(0, ptr->get_ncol(), column_mean_sparse);
  
  // loop over columns of res matrix, and divide by number of cells in each cluster
  for (int i = 0; i < res.ncol(); i++){
    NumericMatrix::Column col = res( _ , i);
    col = col / cluster_id_number[i];
  }
  
  colnames(res) = colnames_vec;
  
  return res;	
}
