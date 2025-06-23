### fast.svd.R  (2006-04-24)
###
###    Efficient Computation of the Singular Value Decomposition
###
### Copyright 2003-06 Korbinian Strimmer
###
###
### This file is part of the `corpcor' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

# private functions



# svd that retains only positive singular values
positive.svd = function(m,  trim = TRUE, tol = NULL)
{
  s = svd(m)
  if (!trim)
    return(s)
  
  if (is.null(tol))
    tol = max(dim(m)) * max(s$d) * .Machine$double.eps
  
  Positive = s$d > tol
  return(list(
    d = s$d[Positive],
    u = s$u[, Positive, drop = FALSE],
    v = s$v[, Positive, drop = FALSE]
  ))
}

# fast computation of svd(m) if n << p
# (n are the rows, p are columns)
utils$svd_small_nr <-
nsmall.svd <-
  function(m, trim = FALSE, tol = NULL, n = nrow(m)) 
{
  B = as.matrix(m %*% t(m), n, n)     # nxn matrix
  s = svd(B, nv = 0)    # of which svd is easy..
  if (!trim) {
    d = pmax(sqrt(s$d), .Machine$double.eps)
    u = as.matrix(s$u, n)
    return(list(
      d = d,
      u = u,
      v = crossprod(m, u) %*% diag(1 / d, length(d))
    ))
  }
  # determine rank of B  (= rank of m)
  if (is.null(tol))
    tol = dim(B)[1] * max(s$d) * .Machine$double.eps
  Positive = s$d > tol
  
  # positive singular values of m
  d = sqrt(s$d[Positive])
  
  # corresponding orthogonal basis vectors
  u = s$u[, Positive, drop = FALSE]
  v = crossprod(m, u) %*% diag(1 / d, nrow = length(d))
  
  return(list(d = d, u = u, v = v))
}

# fast computation of svd(m) if n >> p
# (n are the rows, p are columns)
utils$svd_small_nc <- 
  psmall.svd <-
  function(m, trim = FALSE, tol = NULL, p = ncol(m))
{
  B = as.matrix(crossprod(m), p, p)   # pxp matrix
  s = svd(B, nu = 0)    # of which svd is easy..
  if (!trim) {
    d = pmax(sqrt(s$d), .Machine$double.eps)
    v = as.matrix(s$v, nrow = p)
    return(list(
      d = d,
      v = v,
      u = m %*% t(t(v) / d) #m %*% v %*% diag(1 / d, length(d))
    ))
  }
  # determine rank of B  (= rank of m)
  if (is.null(tol))
    tol = dim(B)[1] * max(s$d) * .Machine$double.eps
  Positive = s$d > tol
  
  # positive singular values of m
  d = sqrt(s$d[Positive])
  
  # corresponding orthogonal basis vectors
  v = s$v[, Positive, drop = FALSE]
  u = m %*% v %*% diag(1 / d, nrow = length(d))
  
  return(list(d = d, u = u, v = v))
}


# public functions

# fast computation of svd(m)

# note that the signs of the columns vectors in u and v
# may be different from that given by svd()

# note that also only positive singular values are returned

fast.svd <-
  utils$fast.svd <-
  function(m, tol = NULL, trim = FALSE)
  {
    n = dim(m)[1]
    p = dim(m)[2]
    
    
    EDGE.RATIO = 2 # use standard SVD if matrix almost square
    if (n > EDGE.RATIO * p)
    {
      return(psmall.svd(m, trim, tol, p))
    }
    else if (EDGE.RATIO * n < p)
    {
      return(nsmall.svd(m, trim, tol, n))
    }
    else
      # if p and n are approximately the same
    {
      return(positive.svd(m, trim, tol))
    }
  }


utils$svdopt <-
  function(mat,
           k = NULL,
           nr = nrow(mat),
           nc = ncol(mat),
           rthin = nc > 2 * nr,
           cthin = nr > 2 * nc,
           trim = FALSE,
           tol = NULL)
  {
    if (is.null(k)) {
#      tryCatch({
        if (rthin)
          return(nsmall.svd(mat, trim, tol, nr))
        if (cthin)
          return(psmall.svd(mat, trim, tol, nc))
        return(positive.svd(mat, trim, tol))
#      }, error = function(e)
#        svd(mat))
    }
    if(k == min(nr, nc)) return(svd(mat))
#    tryCatch({
      if (rthin || cthin || k > 5)
        return(svds(mat, k)) # Rspectra
      return(irlba(mat, k)) # irlba
      
#    }, error = function(e)
#      utils$svd_simple(mat, k))
  }
