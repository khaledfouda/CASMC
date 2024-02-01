subroutine sparse_prod(n, m, r, H, sp, si, sx, result)
    implicit none

    integer, intent(in) :: n, m, r
    integer, intent(in) :: sp(m+1), si(r)
    double precision, intent(in) :: H(n,n), sx(r)
    double precision, intent(out) :: result(r)

    integer :: j, jstart, jend, ii, i, counter, hrow
    double precision, dimension(n) :: h_sub_row

    counter = 1
    do j = 1, m
        jstart = sp(j) + 1
        jend = sp(j + 1)

        if (jstart <= jend) then
            do ii = jstart, jend
                hrow = si(ii) + 1

                ! Extract the sub-row from H corresponding to hrow
                h_sub_row = H(hrow, :)

                ! Compute the dot product for the specific elements
                result(counter) = 0.0
                do i = jstart, jend
                    result(counter) = result(counter) + h_sub_row(si(i) + 1) * sx(i)
                end do
                
                counter = counter + 1
            end do
        end if
    end do
end subroutine sparse_prod
