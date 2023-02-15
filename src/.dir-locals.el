;((nil . ((compile-command . "R -e \"Rcpp::compileAttributes('..')\" && R CMD INSTALL .. && R -e 'devtools::test()'"))))
((nil . ((compile-command . "R -e \"Rcpp::compileAttributes('..')\" && R CMD INSTALL .. && R --vanilla < ../tests/testthat/test-CRAN-line-search.R"))))
