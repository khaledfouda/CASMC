Frob = function(Uold, Dsqold, Vold, U, Dsq, V) {
        denom = sum(Dsqold ^ 2)
        utu = Dsq * (t(U) %*% Uold)
        vtv = Dsqold * (t(Vold) %*% V)
        uvprod = sum(diag(utu %*% vtv))
        num = denom + sum(Dsq ^ 2) - 2 * uvprod
        num / max(denom, 1e-9)
}
