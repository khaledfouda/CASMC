subroutine sparse_prod(n, m, r, H, sp, si, sx, result)
    implicit none

    integer, intent(in) :: n, m, r
    integer, intent(in) :: sp(m+1), si(r)
    double precision, intent(in) :: H(n,n), sx(r)
    double precision, intent(out) :: result(r)

    integer :: j, jstart, jend, ii, counter, hrow
    double precision, dimension(n) :: h_sub_row, sx_sub

    counter = 1
    do j = 1, m
        jstart = sp(j) 
        jend = sp(j + 1) -1

        if (jstart <= jend) then
            sx_sub = sx(jstart:jend)
            do ii = jstart, jend
                hrow = si(ii) 
                call extract_row(H, hrow, n, h_sub_row)
                result(counter) = dot_product(h_sub_row, sx_sub)
                counter = counter + 1
            end do
        end if
    end do
end subroutine sparse_prod

subroutine extract_row(H, row, n, h_sub_row)
    implicit none

    integer, intent(in) :: row, n
    double precision, intent(in) :: H(n, n)
    double precision, intent(out), dimension(n) :: h_sub_row

    integer :: i

    do i = 1, n
        h_sub_row(i) = H(row, i)
    end do
end subroutine extract_row
