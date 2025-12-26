#include <Rcpp.h>
using namespace Rcpp;


//' @title Generate all non-empty subsets (Rcpp)
//' @name generatesubsetsC
//' @description
//' Generate all non-empty subsets of an input vector using C++.
//'
//' @param set An integer vector.
//' @return A list containing all non-empty subsets.
//'
//' @export
// [[Rcpp::export]]
List generate_subsetsC(IntegerVector set) {
 int n = set.size();
 int total = (1 << n);   // 2^n
 List subsets;
 
 for (int mask = 1; mask < total; ++mask) {
   IntegerVector subset;
   for (int i = 0; i < n; ++i) {
     if (mask & (1 << i)) {
       subset.push_back(set[i]);
     }
   }
   subsets.push_back(subset);
 }
   
 return subsets;
}

//' @title Generate nested subsets (Rcpp)
//' @name generatenestedsubsetsC
//' @description
//' Generate nested subsets \{1\}, \{1,2\}, ..., \{1,...,n\} using C++.
//'
//' @param set An integer vector.
//' @return A list of nested subsets.
//'
//' @export
// [[Rcpp::export]]
List generate_nested_subsetsC(IntegerVector set) {
 int n = set.size();
 List subsets(n);
 
 for (int i = 0; i < n; ++i) {
   subsets[i] = set[Range(0, i)];
 }
 
 return subsets;
}