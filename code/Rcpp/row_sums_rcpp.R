
test_mat <- matrix(1:4, nrow = 2)
(r_base <- rowSums(test_mat))

src <- '// Convert SEXP to NumericMatrix
Rcpp::NumericMatrix x(x_in); // x_in comes from R

// Proceed with algorithm like normal
Rcpp::NumericVector row_sums(x.nrow());

for (int i = 0; i < x.nrow(); ++i) {
for(int j = 0; j < x.ncol(); ++j) {
row_sums[i] += x(i,j);
}
}

return row_sums;     // returns the value to R
'




library("inline")

# Set the inputs of function
input_sig = signature(x_in = "numeric")


# Compile function and assign name
row_sum_inline <- cxxfunction(sig = input_sig,
                              body = src, plugin="Rcpp")

row_sum_inline(test_mat)














##############################################################################################################
























