# singlecell
This is a workflow that I used it to analize singlecell data.

### problem 
 when you flow these codes, your may find this error: 
 Error in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  : 
  function 'as_cholmod_sparse' not provided by package 'Matrix'
when using `scDbiFinder()` to filter out doublets, if you prefer not to change the version of the `Matrix` package,first reinstall `irlba` package：
```
install.packages("irlba", type="source")
```
and then this problem would be solve.
