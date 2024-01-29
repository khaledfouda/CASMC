subroutine suvC_H (nrow, ncol, nrank, u, v, h, irow, pcol, nomega, r)
    integer nrow, ncol, nrank, nomega, irow(nomega), pcol(ncol+1)
    double precision u(nrow,nrank), v(ncol,nrank), h(nrow,nrow), r(nomega)
    double precision hu(nrow,nrank)
    integer ii, ni, jstart, jend, i, j, k

    ! Compute the intermediate matrix hu = h * u
    do i = 1, nrow
        do k = 1, nrank
            hu(i, k) = 0.0
            do j = 1, nrow
                hu(i, k) = hu(i, k) + h(i, j) * u(j, k)
            end do
        end do
    end do

    ! Compute the specified dot products
    ni = 0
    do 10 j = 1, ncol
        jstart = pcol(j) + 1
        jend = pcol(j + 1)
        if (jstart > jend) goto 10
        do 20 ii = jstart, jend
            ni = ni + 1
            i = irow(ii) + 1
            r(ni) = DOT_PRODUCT(hu(i,:), v(j,:))
    20   continue
    10   continue

    return
end
