suvC <-
function(u,v,irow,pcol){
  dd=dim(u)
  nnrow=as.integer(dd[1])
  nncol=as.integer(nrow(v))
  nrank=dd[2]
  storage.mode(u)="double"
  storage.mode(v)="double"
  storage.mode(irow)="integer"
  storage.mode(pcol)="integer"
  nomega=as.integer(length(irow))
  .Fortran("suvC",
           nnrow,nncol,nrank,u,v,irow,pcol,nomega,
           r=double(nomega),
           PACKAGE="softImpute"
           )$r
}



sparse_prod <-
  function(n,m,r,H,sp,si,sx){
    storage.mode(H)="double"
    storage.mode(sx)="double"
    storage.mode(si)="integer"
    storage.mode(sp)="integer"
    storage.mode(n)="integer"
    storage.mode(m)="integer"
    storage.mode(r)="integer"
    result = double(r)
    
    .Fortran("sparse_prod",
             n,m,r,H,sp,si,sx
    )$result
  }


sparse_prod <- function(n, m, r, H, sp, si, sx) {
  # Ensure correct data types
  storage.mode(H) = "double"
  storage.mode(sx) = "double"
  storage.mode(si) = "integer"
  storage.mode(sp) = "integer"
  storage.mode(n) = "integer"
  storage.mode(m) = "integer"
  storage.mode(r) = "integer"
  
  # Preallocate the result vector
  result = double(r)
  
  # Call the Fortran subroutine
  .Fortran("sparse_prod",
           as.integer(n), as.integer(m), as.integer(r),
           as.double(H), as.integer(sp), as.integer(si),
           as.double(sx), result = as.double(result)
  )$result
}
