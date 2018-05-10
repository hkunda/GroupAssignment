
GroupAssignment
===============

Documentation
-------------

Installation
------------

``` r
devtools::install_github("hkunda/GroupAssignment")
```

### Inputs

@rankings \[NumericMatrix\]: a matrix of rankings for each student. The ith row corresponds to student i, and the jth column corresponds to the student in position j on any student's ranking list. If there are N students, then each student should have an ID in the range \[1,N\], which is used both as the row index (to refer to a student's rankings) and entries of the matrix. For example, if rankings\[i,j\] = k, then student i has placed student k in position j on i's ranking.

@students \[int\]: the number of students. This should be equal to the number of rows in @rankings.

@leaders \[int\]: the number of leaders to pick; equivalently, the number of groups to form.

@minGroupSize \[int\]: the minimum number of students in a group. NOTE: this value does include the group leader. So if minGroupSize = 4, then every group must consist of 1 leader and at least 3 other students.

@maxGroupSize \[int\]: the maximum number of students in a group. NOTE: this value does include the group leader. So if maxGroupSize = 6, then every group must consist of 1 leader and at most 5 other students.

### Outputs

Returns a NumericVector of length students+1. The ith element of the output is the id of the student leading the group that student i is in (for 0 &lt;= i &lt;= N). The N+1st element of the output is the total cost of the optimal matching.

### Examples

``` r
library(igraph)
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library(GroupAssignment)

set.seed(712)

x <- sample_smallworld(1, 20, 4, .1)
d <- distances(x) + 1
d <- max(d) - d + 1

ans <- optimalAssignment(d, leaders = 3, minGroupSize = 2, maxGroupSize = 10)
ans
```

    ##  [1]   1   4   4   4   4   4   4   4   4   4   1   4  13  13  13   1  13
    ## [18]  13  13   1 198
