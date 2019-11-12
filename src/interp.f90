! Copyright or © or Copr. Actimar/IFREMER (contributor(s) : Stephane Raynaud) (2010-2019)
!
! raynaud@actimar.fr
!
!
! This software is a computer program whose purpose is to provide
! utilities for handling oceanographic and atmospheric data,
! with the ultimate goal of validating the MARS model from IFREMER.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and,  more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!
! =============================================================================
! ================================== 1D =======================================
! =============================================================================

subroutine interp1d(vari, yi, varo, yo, method, nx, nyi, nyo, extrap)
    ! Interpolation along the second axis (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yi: input y axis
    ! - yo: output y axis
    ! - method: 0 = nearest, 1 = linear, 2 = cubic, 3 = hermit
    ! - extrap: 0 = do not extrapolate, 1 = top, -1 = bottom, 2 = both
    !
    ! See: http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,method
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yi(nyi), yo(nyo)
    real(kind=8), intent(out) :: varo(nx,nyo)
    integer, intent(in), optional :: extrap

    ! Internal
    integer :: iyi,iyo
    real(kind=8) :: dy0,dy1,yi0,yi1,mu,m1
    real(kind=8) :: tension,bias,a0,a1,a2,a3 ! Hermit
    real(kind=8),allocatable :: vc0(:),vc1(:)
    logical :: bmask(nx,nyi)

    ! Initialisation
    m1 = -1d0
    varo = sqrt(m1)
    if(method>1) allocate(vc0(nx),vc1(nx))
    bias = 0.
    tension = 0.
    bmask = isnan(vari)

    ! Loop on input grid
    iyo = 1
    do iyi = 1, nyi-1

        yi0 = yi(iyi)
        yi1 = yi(iyi+1)

        if (yi1 < yo(1)) cycle
        if (yi0 > yo(nyo)) exit

        ! Loop on output grid
        do while (iyo <= nyo )

             ! Out of interval
            if (yo(iyo) < yi0) then
                iyo = iyo + 1
                cycle
            endif
            if (yo(iyo) > yi1) exit

            ! Distances to neighbours
            dy0 = yo(iyo) - yi0
            dy1 = yi1 - yo(iyo)

            ! Interpolation
            if (dy0==0.) then
                varo(:,iyo) = vari(:, iyi)
            else if (dy1==0.) then
                varo(:,iyo) = vari(:, iyi+1)

            elseif (method==0) then

                ! Nearest neighbour
                if (dy0 < dy1) then
                    varo(:,iyo) = vari(:, iyi)
                else
                    varo(:,iyo) = vari(:, iyi+1)
                endif

            elseif (method==1)then

                ! Linear
                varo(:,iyo) = &
                &    (vari(:, iyi)*dy1 + vari(:, iyi+1) * dy0) / &
                &    (dy0 + dy1)

            else

                ! Cubic and Hermit
                !
                if (iyi==1)then ! y0
                    vc0 = 2*vari(:, iyi) - vari(:, iyi+1)
                else
                    vc0 = vari(:, iyi-1)
                endif
                if (iyi==nyi-1)then ! y3
                    vc1 = 2*vari(:, iyi+1) - vari(:, iyi)
                else
                    vc1 = vari(:, iyi+2)
                endif
                mu = dy0 / (dy0+dy1)

                if (method==2)then

                    ! Cubic
                    !   mu2 = mu*mu;
                    !   a0 = y3 - y2 - y0 + y1;
                    !   a1 = y0 - y1 - a0;
                    !   a2 = y2 - y0;
                    !   a3 = y1;
                    !   return (a0*mu*mu2+a1*mu2+a2*mu+a3);

                    varo(:,iyo) = vc1 - vari(:, iyi+1) - vc0 + vari(:, iyi) !a0
                    varo(:,iyo) = mu**3*varo(:,iyo) + mu**2*(vc0-vari(:, iyi)-varo(:,iyo)) ! a0*mu^3 + a1*mu
                    varo(:,iyo) = varo(:,iyo) + mu*(vari(:, iyi+1)-vc0)
                    varo(:,iyo) = varo(:,iyo) + vari(:, iyi)

                else

                    ! Hermit
                    !   mu2 = mu * mu;
                    !   mu3 = mu2 * mu;
                    !   a0 =  2*mu3 - 3*mu2 + 1;
                    !   a1 =    mu3 - 2*mu2 + mu;
                    !   a2 =    mu3 -   mu2;
                    !   a3 = -2*mu3 + 3*mu2;
                    !   m0  = (y1-y0)*(1+bias)*(1-tension)/2;
                    !   m0 += (y2-y1)*(1-bias)*(1-tension)/2;
                    !   m1  = (y2-y1)*(1+bias)*(1-tension)/2;
                    !   m1 += (y3-y2)*(1-bias)*(1-tension)/2;
                    !   return(a0*y1 + a1*m0 + a2*m1 + a3*y2);
                    a0 = 2*mu**3 - 3*mu**2 + 1
                    a1 =    mu**3 - 2*mu**2 + mu
                    a2 =    mu**3 -   mu**2
                    a3 = -2*mu**3 + 3*mu**2
                    varo(:,iyo) = a0*vari(:, iyi)
                    varo(:,iyo) = varo(:,iyo) + a1*( &
                    &    (vari(:,iyi)-vc0)           *(1+bias)*(1-tension)/2 + &
                    &    (vari(:, iyi+1)-vari(:,iyi))*(1-bias)*(1-tension)/2)
                    varo(:,iyo) = varo(:,iyo) + a2*(&
                    &    (vari(:, iyi+1)-vari(:,iyi))*(1+bias)*(1-tension)/2 + &
                    &    (vc1-vari(:, iyi+1))        *(1-bias)*(1-tension)/2)
                    varo(:,iyo) = varo(:,iyo) + a3*vari(:, iyi+1)
                endif


            endif
            iyo = iyo + 1

        end do
    end do

    ! Extrapolation with nearest
    !if(method==0 .and.present(extrap).and.extrap/=0)then
    if(present(extrap).and.extrap/=0)then
        if((extrap==-1 .or. extrap==2) .and. (yo(1)<yi(1)))then
            do iyo=1,nyo
                if(yi(1)<yo(iyo))exit
                varo(:,iyo) = vari(:,1)
            enddo
        endif
        if((extrap==1 .or. extrap==2) .and. (yo(nyo)>yi(nyi)))then
            do iyo=nyo,1,-1
                if(yi(nyi)>yo(iyo))exit
                varo(:,iyo) = vari(:,nyi)
            enddo
        endif
    endif


end subroutine interp1d

subroutine interp1dx(vari, yi, varo, yo, method, nx, nxb, nyi, nyo, extrap)
    ! Interpolation along the second axis (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yi: input y axis
    ! - yo: output y axis
    ! - method: 0 = nearest, 1 = linear, 2 = cubic, 3 = hermit
    ! - extrap: 0 = do not extrapolate, 1 = top, -1 = bottom, 2 = both
    !
    ! See: http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,method,nxb
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nyo)
    real(kind=8), intent(out) :: varo(nx,nyo)
    integer, intent(in), optional :: extrap

    ! Internal
    integer :: iyi,iyo,ib
    real(kind=8) :: dy0(nxb),dy1(nxb)
    real(kind=8) :: tension,bias !
    real(kind=8) :: m1, mu(nxb)
    real(kind=8),allocatable :: vc0(:),vc1(:),a0(:),a1(:),a2(:),a3(:)
    integer :: ix0,ix1
    logical :: bitv(nxb), bmask(nx,nyi)


    ! Initialisation
    bmask = isnan(vari)
    m1 = -1d0
    varo = sqrt(m1)
    if(method>1) allocate(vc0(nxb),vc1(nxb))
    if(method==3)allocate(a0(nxb),a1(nxb),a2(nxb),a3(nxb))
    bias = 0.
    tension = 0.

    ! Loop on input grid
    do iyi = 1, nyi-1

        ! Loop on output grid
        do iyo = 1,nyo

            dy0 = yo(iyo)-yi(:,iyi)
            dy1 = yi(:,iyi+1)-yo(iyo)
            bitv = yo(iyo)>=yi(:,iyi).and.yo(iyo)<=yi(:,iyi+1)
            mu = dy0/(dy0+dy1)

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1

                if (method==0) then

                    ! Nearest
                    where(bitv)
                        where(dy0<dy1)
                            varo(ix0:ix1,iyo) = vari(ix0:ix1,iyi)
                        elsewhere
                            varo(ix0:ix1,iyo) = vari(ix0:ix1,iyi+1)
                        end where
                    end where

                elseif(method==1)then

                    ! Linear
                    where(bitv)
                        varo(ix0:ix1,iyo) = &
                        &    (vari(ix0:ix1, iyi)*dy1 + vari(ix0:ix1, iyi+1)*dy0) / &
                        &    (dy0+dy1)
                    end where

                else

                    ! Extrapolations
                    if (iyi==1)then ! y0
                        vc0 = 2*vari(ix0:ix1, iyi)-vari(ix0:ix1, iyi+1)
                    else
                        vc0 = vari(ix0:ix1, iyi-1)
                    endif
                    if (iyi==nyi-1)then ! y3
                        vc1 = 2*vari(ix0:ix1, iyi+1)-vari(ix0:ix1, iyi)
                    else
                        vc1 = vari(ix0:ix1, iyi+2)
                    endif

                    if (method==2)then

                        ! Cubic
                        where(bitv)
                            varo(ix0:ix1,iyo) = vc1 - vari(ix0:ix1, iyi+1) - vc0 + vari(ix0:ix1, iyi)
                            varo(ix0:ix1,iyo) = mu**3*varo(ix0:ix1,iyo) + mu**2*(vc0-vari(ix0:ix1, iyi)-varo(ix0:ix1,iyo))
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + mu*(vari(ix0:ix1, iyi+1)-vc0)
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1, iyi)
                        end where

                    else

                        ! Hermit
                        a0 = 2*mu**3 - 3*mu**2 + 1
                        a1 =    mu**3 - 2*mu**2 + mu
                        a2 =    mu**3 -   mu**2
                        a3 = -2*mu**3 + 3*mu**2
                        where(bitv)
                            varo(ix0:ix1,iyo) = a0*vari(ix0:ix1, iyi)
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + a1*( &
                            &    (vari(ix0:ix1,iyi)-vc0)           *(1+bias)*(1-tension)/2 + &
                            &    (vari(ix0:ix1, iyi+1)-vari(ix0:ix1,iyi))*(1-bias)*(1-tension)/2)
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + a2*(&
                            &    (vari(ix0:ix1, iyi+1)-vari(ix0:ix1,iyi))*(1+bias)*(1-tension)/2 + &
                            &    (vc1-vari(ix0:ix1, iyi+1))        *(1-bias)*(1-tension)/2)
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + a3*vari(ix0:ix1, iyi+1)
                        end where

                    endif
                endif
            end do
        end do
    end do

    ! Extrapolation with nearest
    !if(method==0 .and.present(extrap).and.extrap/=0)then
    if(present(extrap).and.extrap/=0)then
        do iyo = 1,nyo
            if(extrap==-1 .or. extrap==2)then
                bitv = yo(iyo)<yi(:,1) ! below
                do ib = 1, nx/nxb
                    ix0 = 1+(ib-1)*nxb
                    ix1 = ix0+nxb-1
                    where(bitv) varo(ix0:ix1,iyo) = vari(ix0:ix1, 1)
                enddo
            endif
            if(extrap==1 .or. extrap==2)then
                bitv = yo(iyo)>yi(:,nyi) ! above
                do ib = 1, nx/nxb
                    ix0 = 1+(ib-1)*nxb
                    ix1 = ix0+nxb-1
                    where(bitv) varo(ix0:ix1,iyo) = vari(ix0:ix1, nyi)
                enddo
            endif
        enddo
    endif

    ! Deallocations
    if(method>1) deallocate(vc0,vc1)
    if(method==3)deallocate(a0,a1,a2,a3)

end subroutine interp1dx

subroutine interp1dxx(vari, yi, varo, yo, method, nx, nxb, nyi, nyo, extrap)
    ! Interpolation along the second axis (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yi: input y axis
    ! - yo: output y axis
    ! - method: 0 = nearest, 1 = linear, 2 = cubic, 3 = hermit
    ! - extrap: 0 = do not extrapolate, 1 = top, -1 = bottom, 2 = both
    !
    ! See: http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,method,nxb
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nxb,nyo)
    real(kind=8), intent(out) :: varo(nx,nyo)
    integer, intent(in), optional :: extrap

    ! Internal
    integer :: iyi,iyo,ib
    real(kind=8) :: dy0(nxb),dy1(nxb)
    real(kind=8) :: tension, bias, m1 !
    real(kind=8) :: mu(nxb)
    real(kind=8),allocatable :: vc0(:),vc1(:),a0(:),a1(:),a2(:),a3(:)
    integer :: ix0,ix1
    logical :: bitv(nxb)


    ! Initialisation
    if(method>1) allocate(vc0(nxb),vc1(nxb))
    if(method==3)allocate(a0(nxb),a1(nxb),a2(nxb),a3(nxb))
    bias = 0d0
    tension = 0d0
    m1 = -1d0
    varo = sqrt(m1)

    ! Loop on input grid
    do iyi = 1, nyi-1

        ! Loop on output grid
        do iyo = 1,nyo

            dy0 = yo(:,iyo)-yi(:,iyi)
            dy1 = yi(:,iyi+1)-yo(:,iyo)
            bitv = yo(:,iyo)>=yi(:,iyi).and.yo(:,iyo)<=yi(:,iyi+1)
            mu = dy0/(dy0+dy1)

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1

                if (method==0) then

                    ! Nearest
                    where(bitv)
                        where(dy0<dy1)
                            varo(ix0:ix1,iyo) = vari(ix0:ix1,iyi)
                        elsewhere
                            varo(ix0:ix1,iyo) = vari(ix0:ix1,iyi+1)
                        end where
                    end where

                elseif(method==1)then

                    ! Linear
                    where(bitv)
                        varo(ix0:ix1,iyo) = &
                        &    (vari(ix0:ix1, iyi)*dy1 + vari(ix0:ix1, iyi+1)*dy0) / &
                        &    (dy0+dy1)
                    end where

                else

                    ! Extrapolations
                    if (iyi==1)then ! y0
                        vc0 = 2*vari(ix0:ix1, iyi)-vari(ix0:ix1, iyi+1)
                    else
                        vc0 = vari(ix0:ix1, iyi-1)
                    endif
                    if (iyi==nyi-1)then ! y3
                        vc1 = 2*vari(ix0:ix1, iyi+1)-vari(ix0:ix1, iyi)
                    else
                        vc1 = vari(ix0:ix1, iyi+2)
                    endif

                    if (method==2)then

                        ! Cubic
                        where(bitv)
                            varo(ix0:ix1,iyo) = vc1 - vari(ix0:ix1, iyi+1) - vc0 + vari(ix0:ix1, iyi)
                            varo(ix0:ix1,iyo) = mu**3*varo(ix0:ix1,iyo) + mu**2*(vc0-vari(ix0:ix1, iyi)-varo(ix0:ix1,iyo))
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + mu*(vari(ix0:ix1, iyi+1)-vc0)
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1, iyi)
                        end where

                    else

                        ! Hermit
                        a0 = 2*mu**3 - 3*mu**2 + 1
                        a1 =    mu**3 - 2*mu**2 + mu
                        a2 =    mu**3 -   mu**2
                        a3 = -2*mu**3 + 3*mu**2
                        where(bitv)
                            varo(ix0:ix1,iyo) = a0*vari(ix0:ix1, iyi)
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + a1*( &
                            &    (vari(ix0:ix1,iyi)-vc0)           *(1+bias)*(1-tension)/2 + &
                            &    (vari(ix0:ix1, iyi+1)-vari(ix0:ix1,iyi))*(1-bias)*(1-tension)/2)
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + a2*(&
                            &    (vari(ix0:ix1, iyi+1)-vari(ix0:ix1,iyi))*(1+bias)*(1-tension)/2 + &
                            &    (vc1-vari(ix0:ix1, iyi+1))        *(1-bias)*(1-tension)/2)
                            varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + a3*vari(ix0:ix1, iyi+1)
                        end where

                    endif
                endif
            end do
        end do
    end do

    ! Extrapolation with nearest
    !if(method==0 .and.present(extrap).and.extrap/=0)then
    if(present(extrap).and.extrap/=0)then
        do iyo = 1,nyo
            if(extrap==-1 .or. extrap==2)then
                bitv = yo(:,iyo)<yi(:,1) ! below
                do ib = 1, nx/nxb
                    ix0 = 1+(ib-1)*nxb
                    ix1 = ix0+nxb-1
                    where(bitv) varo(ix0:ix1,iyo) = vari(ix0:ix1, 1)
                enddo
            endif
            if(extrap==1 .or. extrap==2)then
                bitv = yo(:,iyo)>yi(:,nyi) ! above
                do ib = 1, nx/nxb
                    ix0 = 1+(ib-1)*nxb
                    ix1 = ix0+nxb-1
                    where(bitv) varo(ix0:ix1,iyo) = vari(ix0:ix1, nyi)
                enddo
            endif
        enddo
    endif

    ! Deallocations
    if(method>1) deallocate(vc0,vc1)
    if(method==3)deallocate(a0,a1,a2,a3)

end subroutine interp1dxx

subroutine extrap1d(vari, varo, extrap, nx, ny)
    ! Extrapolate valid data to the top and/or bottom
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - extrap: 0 = do not extrapolate, 1 = top, -1 = bottom, 2 = both
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: vari(nx,ny)
    real(kind=8), intent(out) :: varo(nx,ny)
    integer, intent(in), optional :: extrap

    ! Internal
    integer :: ix,iy, jj(ny), iymin, iymax
    logical :: valid(ny)

    ! Initialisation
    varo = vari
    if(extrap==0)return
    jj = (/(iy, iy=1,ny)/)

    ! Loop on extra dim
    !$OMP PARALLEL DO PRIVATE(ix,iymin,iymax,valid)
    !$& SHARED(vari,jj,varo,extrap)
    do ix = 1, nx

        valid = isnan(vari(ix,:))

        if(any(valid))then

            if(extrap==-1 .or. extrap==2)then
                iymin = minval(jj, mask=valid)
                if(iymin>1) varo(ix,1:iymin-1) = varo(ix,iymin)
            endif

            if(extrap==1 .or. extrap==2)then
                iymax = maxval(jj, mask=valid)
                if(iymax<ny) varo(ix,iymax+1:ny) = varo(ix,iymax)
            endif

        endif

    enddo
    !$OMP END PARALLEL DO

end subroutine extrap1d

! =============================================================================

subroutine cellave1d(vari, yib, varo, yob, conserv, nx, nyi, nyo, extrap)
    ! Remapping along the second axis (y)

    implicit none

    ! Extrernal
    integer, intent(in) :: nx, nyi, nyo
    integer, intent(in) :: conserv
    real(kind=8),intent(in) :: vari(nx,nyi)
    real(kind=8),intent(out) :: varo(nx,nyo)
    real(kind=8),intent(in) :: yib(nyi+1),yob(nyo+1)
    integer, intent(in), optional :: extrap

    ! Local
    integer :: iyi,iyo
    real(kind=8) :: yib0, yib1, wo(nx), dyio, m1


    ! Init
    varo = 0d0
    m1 = -1d0

    ! Loop on output cells
    iyi = 1
    do iyo = 1, nyo

        if(yob(iyo+1)==yob(iyo))cycle

        ! Loop on input cells
        wo = 0.
        do while (iyi<=nyi)

            ! Current input bounds
            yib0 = yib(iyi)
            yib1 = yib(iyi+1)
            if(iyi==1 .and. yib0>yob(iyo) .and. present(extrap) .and. &
                & (extrap==-1 .or. extrap==2)) yib0 = yob(iyo)
            if(iyi==nyi .and. yib1<yob(iyo+1) .and. present(extrap) .and. &
                & (extrap==1 .or. extrap==2)) yib1 = yob(iyo+1)

            ! No intersection
            if(yib0>yob(iyo+1)) exit
            if(yib1<yob(iyo))then
                iyi = iyi + 1
                cycle
            endif

            ! Contribution of intersection
            dyio = min(yib1,yob(iyo+1))-max(yib0,yob(iyo))
            if(conserv==1.and.yib0/=yib1) dyio = dyio / (yob(iyo+1)-yob(iyo))
            where(.not.isnan(vari(:,iyi)))
                wo = wo + dyio
                varo(:,iyo) = varo(:,iyo) + vari(:,iyi)*dyio
            endwhere

            ! Next input cell?
            if(yib1>=yob(iyo+1)) exit
            iyi = iyi + 1

        enddo

        ! Normalize
        if(conserv==1)where(wo/=0.)wo = 1.
        varo(:,iyo) = merge(varo(:,iyo)/wo, sqrt(m1), wo/=0.)

    enddo

end subroutine cellave1d

subroutine cellave1dx(vari, yib, varo, yob, conserv, nx, nxb, nyi, nyo, extrap)
    ! Remapping along the second axis (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yib: input y axis bounds
    ! - yob: output y axis bounds
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    integer, intent(in) ::  conserv
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yib(nxb,nyi+1), yob(nyo+1)
    real(kind=8), intent(out) :: varo(nx,nyo)
    integer, intent(in), optional :: extrap

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: wo(nx),dyo,yib0(nxb),yib1(nxb), m1
    integer :: ix0,ix1
    logical :: mapi(nxb,4)

    ! Initialisation
    m1 = -1d0
    varo = 0.

    ! Loop on output cells
    do iyo = 1,nyo

        dyo = yob(iyo+1) - yob(iyo)
        if(dyo==0)cycle
        if(conserv==0)dyo = 1.

        ! Loop on input cells
        wo = 0.
        do iyi = 1, nyi

            ! Current input bounds
            yib0 = yib(:, iyi)
            yib1 = yib(:, iyi+1)
            if(iyi==1 .and. present(extrap) .and. &
                    & (extrap==-1 .or. extrap==2)) &
                & where(yib0>yob(iyo)) yib0 = yob(iyo)
            if(iyi==nyi .and. present(extrap) .and. &
                    & (extrap==1 .or. extrap==2)) &
                & where(yib1<yob(iyo+1)) yib1 = yob(iyo+1)

            ! No intersection
            if(all(yib0>=yob(iyo+1)))exit
            if(all(yib1<yob(iyo)))cycle

            ! Conditional arrays
            mapi(:,1) = yib0>=yob(iyo).and.yib1<=yob(iyo+1)
            mapi(:,2) = yib0< yob(iyo).and.yib1> yob(iyo+1)
            mapi(:,3) = yib0< yob(iyo).and.yib1>yob(iyo)
            mapi(:,4) = yib0< yob(iyo+1).and.yib1> yob(iyo+1)

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(.not.isnan(vari(ix0:ix1,iyi)))
                    where(mapi(:,1))

                        ! Input inside
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(yib1-yib0) / dyo
                        wo(ix0:ix1) = wo(ix0:ix1) + yib1-yib0

                    elsewhere(mapi(:,2))

                        ! Output inside
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(yob(iyo+1)-yob(iyo)) / dyo
                        wo(ix0:ix1) = wo(ix0:ix1) + (yob(iyo+1)-yob(iyo))

                    elsewhere(mapi(:,3))

                        ! Input partly below
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(yib1-yob(iyo)) / dyo
                        wo(ix0:ix1) = wo(ix0:ix1) + yib1-yob(iyo)

                    elsewhere(mapi(:,4))

                        ! Input partly above
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(yob(iyo+1)-yib0) / dyo
                        wo(ix0:ix1) = wo(ix0:ix1) + (yob(iyo+1)-yib0)

                    endwhere
                endwhere
            enddo
        enddo

        ! Normalize
        if(conserv==1)where(wo/=0.)wo = 1.
        where(wo/=0.)
            varo(:,iyo) =  varo(:,iyo)/wo
        elsewhere
            varo(:,iyo) = sqrt(m1)
        endwhere
    enddo

end subroutine cellave1dx

subroutine cellave1dxx(vari, yib, varo, yob, conserv, nx, nxb, nyi, nyo, extrap)
    ! Remapping between two variable axes in space (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yib: input y axis bounds
    ! - yob: output y axis bounds
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    integer,intent(in) ::  conserv
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yib(nxb,nyi+1), yob(nxb,nyo+1)
    real(kind=8), intent(out) :: varo(nx,nyo)
    integer, intent(in), optional :: extrap

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: wo(nx),dyo(nxb),yib0(nxb),yib1(nxb),m1
    integer :: ix0,ix1
    logical :: mapi(nxb,4),dyook(nxb)



    ! Initialisation
    m1 = -1d0
    varo = 0.

    ! Loop on output grid
    do iyo = 1,nyo

        dyo = yob(:, iyo+1) - yob(:, iyo)
        dyook = dyo /= 0.
        if(conserv==0)dyo = 1.
        wo = 0.

        ! Loop on input grid
        do iyi = 1, nyi

            ! Current input bounds
            yib0 = yib(:, iyi)
            yib1 = yib(:, iyi+1)
            if(iyi==1 .and. present(extrap) .and. &
                    & (extrap==-1 .or. extrap==2)) &
                & where(yib0>yob(:, iyo)) yib0 = yob(:, iyo)
            if(iyi==nyi .and. present(extrap) .and. &
                    & (extrap==1 .or. extrap==2)) &
                & where(yib1<yob(:, iyo+1)) yib1 = yob(:, iyo+1)

            ! No intersection
            if(all(yib0>=yob(:, iyo+1)))exit
            if(all(yib1<yob(:, iyo)))cycle

            ! Conditional arrays
            mapi(:,1) = yib0>=yob(:, iyo).and.yib1<=yob(:, iyo+1) ! Inside
            mapi(:,2) = yib0< yob(:, iyo).and.yib1> yob(:, iyo+1) ! Embed
            mapi(:,3) = yib0< yob(:, iyo).and.yib1>yob(:, iyo) ! Below
            mapi(:,4) = yib0<yob(:, iyo+1).and.yib1> yob(:, iyo+1) ! Above

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(.not.isnan(vari(ix0:ix1,iyi)) .and. dyook)
                    where(mapi(:,1))

                        ! Input inside
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(yib1-yib0) / dyo
                        wo(ix0:ix1) = wo(ix0:ix1) + yib1-yib0

                    elsewhere(mapi(:,2))

                        ! Output inside
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(yob(:, iyo+1)-yob(:, iyo)) / dyo
                        wo(ix0:ix1) = wo(ix0:ix1) + yob(:, iyo+1)-yob(:, iyo)

                    elsewhere(mapi(:,3))

                        ! Input partly below
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(yib1-yob(:, iyo)) / dyo
                        wo(ix0:ix1) = wo(ix0:ix1) + yib1-yob(:, iyo)

                    elsewhere(mapi(:,4))

                        ! Input partly above
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(yob(:, iyo+1)-yib0) / dyo
                        wo(ix0:ix1) = wo(ix0:ix1) + yob(:, iyo+1)-yib0

                    endwhere
                endwhere
            enddo
        enddo

        ! Normalize
        if(conserv==1)where(wo>0d0)wo = 1.
        where(wo>0d0)
            varo(:,iyo) =  varo(:,iyo)/wo
        elsewhere
            varo(:,iyo) = sqrt(m1)
        endwhere
    enddo

end subroutine cellave1dxx


subroutine cellerr1d(vari, yi, varo, yo, yob, errm, errl, erro, nx, nxb, nyi, nyo)
    ! Error-based weighthed local mooving average interpolation between
    ! two variables axes in space (y) with interpolation error estimate
    !
    ! :Params:
    !
    !   - vari: input variable
    !   - varo: output variable
    !   - yi: input y axis
    !   - yob: output y axis bounds
    !   - errm: input measurement error
    !   - errl: input lag error per yi units
    !   - erro: output error
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    real(kind=8), intent(in) :: vari(nx,nyi), errm(nx,nyi)
    real(kind=8), intent(in) :: yi(nyi), yo(nyo), yob(nyo+1), errl(nxb)
    real(kind=8), intent(out) :: varo(nx,nyo), erro(nx,nyo)

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: zerrl(nxb),m1
    integer :: ix0,ix1
    logical :: goodi(nx,nyi), berrl(nxb)

    ! Masks
    goodi = .not. isnan(vari)
    berrl = isnan(errl)
    m1 = -1d0

    ! Initialisation
    varo = 0d0
    erro = 0d0

    ! Loop on output grid
    iyi = 1
    do iyo = 1,nyo

        ! Loop on input grid
        do while (iyi<=nyi)

            ! Not useful
            if(yi(iyi)>yob(iyo+1)) exit
            if(yi(iyi)<yob(iyo).or..not.any(goodi(:,iyi)))then
                iyi = iyi +1
                cycle
            endif

            ! Quadratic lag error
            where(berrl)
                zerrl = 0d0
            elsewhere(errl<0)
                zerrl = errl**2 * abs(yi(iyi)-yo(iyo))
            elsewhere
                zerrl = errl * abs(yi(iyi)-yo(iyo))
            endwhere

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(goodi(ix0:ix1,iyi))

                    varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + &
                        & vari(ix0:ix1,iyi) / (errm(ix0:ix1, iyi)**2 + zerrl)
                    erro(ix0:ix1,iyo) = erro(ix0:ix1,iyo) + &
                        & 1d0 / (errm(ix0:ix1, iyi)**2 + zerrl)

                endwhere
            enddo

            iyi = iyi + 1

        enddo

        ! Normalize
        where(erro(:,iyo)>0d0)
            varo(:,iyo) =  varo(:,iyo)/erro(:,iyo)
            erro(:,iyo) = 1d0/sqrt(erro(:,iyo))
        elsewhere
            varo(:,iyo) = sqrt(m1)
            erro(:,iyo) = sqrt(m1)
        endwhere

    enddo


end subroutine cellerr1d



subroutine cellerr1dx(vari, yi, varo, yo, yob, errm, errl, erro, nx, nxb, nyi, nyo)
    ! Error-based weighthed local mooving average interpolation between
    ! two variables axes in space (y) with interpolation error estimate
    !
    ! :Params:
    !
    !   - vari: input variable
    !   - varo: output variable
    !   - yi: input y axis
    !   - yob: output y axis bounds
    !   - errm: input measurement error
    !   - errl: input lag error per yi units
    !   - erro: output error
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    real(kind=8), intent(in) :: vari(nx,nyi), errm(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nyo), yob(nyo+1), errl(nxb)
    real(kind=8), intent(out) :: varo(nx,nyo), erro(nx,nyo)

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: zerrl(nxb), m1
    integer :: ix0,ix1
    logical :: goodi(nxb), berrl(nxb), bvari(nx,nyi)

    ! Masks
    bvari = isnan(vari)
    berrl = isnan(errl)
    m1 = -1d0

    ! Initialisation
    varo = 0d0
    erro = 0d0

    ! Loop on output grid
    do iyo = 1,nyo

        ! Loop on input grid
        do iyi = 1, nyi

            ! Valid?
            goodi(:) = .not. bvari(:,iyi) .and. &
                & yi(:, iyi)>=yob(iyo) .and. yi(:, iyi)<=yob(iyo+1)
            if(.not.any(goodi))cycle

            ! Quadratic lag error
            where(berrl)
                zerrl = 0d0
            elsewhere(errl<0)
                zerrl = errl**2 * abs(yi(:, iyi)-yo(iyo))
            elsewhere
                zerrl = errl * abs(yi(:, iyi)-yo(iyo))
            endwhere

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(goodi(ix0:ix1))

                    varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + &
                        & vari(ix0:ix1,iyi) / (errm(ix0:ix1, iyi)**2 + zerrl)
                    erro(ix0:ix1,iyo) = erro(ix0:ix1,iyo) + &
                        & 1d0 / (errm(:, iyi)**2 + zerrl)

                endwhere
            enddo
        enddo

        ! Normalize
        where(erro(:,iyo)>0d0)
            varo(:,iyo) =  varo(:,iyo)/erro(:,iyo)
            erro(:,iyo) = 1d0/sqrt(erro(:,iyo))
        elsewhere
            varo(:,iyo) = sqrt(m1)
            erro(:,iyo) = sqrt(m1)
        endwhere
    enddo


end subroutine cellerr1dx



subroutine cellerr1dxx(vari, yi, varo, yo, yob, errm, errl, erro, nx, nxb, nyi, nyo)
    ! Error-based weighthed local mooving average interpolation between
    ! two variables axes in space (y) with interpolation error estimate
    !
    ! :Params:
    !
    !   - vari: input variable
    !   - varo: output variable
    !   - yi: input y axis
    !   - yo: output y axis bounds
    !   - errm: input measurement error
    !   - errl: input lag error per yi units
    !   - erro: output error
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    real(kind=8), intent(in) :: vari(nx,nyi), errm(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nxb,nyo), yob(nxb,nyo+1), errl(nxb)
    real(kind=8), intent(out) :: varo(nx,nyo), erro(nx,nyo)

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: zerrl(nxb), m1
    integer :: ix0,ix1
    logical :: goodi(nxb), berrl(nxb), bvari(nx,nyi)

    ! Masks
    bvari = isnan(vari)
    berrl = isnan(errl)
    m1 = -1d0

    ! Initialisation
    varo = 0d0
    erro = 0d0

    ! Loop on output grid
    do iyo = 1,nyo

        ! Loop on input grid
        do iyi = 1, nyi

            ! Valid?
            goodi(:) = .not. bvari(:,iyi) .and. &
                & yi(:, iyi)>=yob(:, iyo) .and. yi(:, iyi)<=yob(:, iyo+1)
            if(.not.any(goodi))cycle

            ! Quadratic lag error
            where(berrl)
                zerrl = 0d0
            elsewhere(errl<0)
                zerrl = errl**2 * abs(yi(:, iyi)-yo(:, iyo))
            elsewhere
                zerrl = errl * abs(yi(:, iyi)-yo(:, iyo))
            endwhere

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(goodi(ix0:ix1))

                    varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + &
                        & vari(ix0:ix1,iyi) / (errm(ix0:ix1, iyi)**2 + zerrl)
                    erro(ix0:ix1,iyo) = erro(ix0:ix1,iyo) + &
                        & 1d0 / (errm(:, iyi)**2 + zerrl)

                endwhere
            enddo
        enddo

        ! Normalize
        where(erro(:,iyo)>0d0)
            varo(:,iyo) =  varo(:,iyo)/erro(:,iyo)
            erro(:,iyo) = 1d0/sqrt(erro(:,iyo))
        elsewhere
            varo(:,iyo) = sqrt(m1)
            erro(:,iyo) = sqrt(m1)
        endwhere
    enddo


end subroutine cellerr1dxx


! =============================================================================
! ======================================= 2D ==================================
! =============================================================================

subroutine nearest2d(vari, xxi, yyi, varo, xxo, yyo, nb, nogeo, nxi, nyi, nxo, nyo, nz)
    ! Nearest neighbour interpolation of two curvilinear grids

    implicit none

    ! Parameters
    integer,intent(in) :: nxi, nyi, nxo, nyo, nz,nb
    real(kind=8), intent(in) :: vari(nz,nyi,nxi)
    real(kind=8), intent(in) :: xxi(nyi,nxi),yyi(nyi,nxi),xxo(nyo,nxo),yyo(nyo,nxo)
    real(kind=8),intent(out) :: varo(nz,nyo,nxo)
    integer, intent(in), optional :: nogeo

    ! Local variables
    integer :: ixo,iyo, znb, ixlastline,iylastline,znb2,ixlast,iylast
    integer :: ixmin,ixmax,iymin,iymax, imin,jmin
    logical :: geo

    geo = .not. present(nogeo) .or. nogeo==0

    ! Loop on output points
    if(nb==0)then

        ! Scan all input points everytime
        do ixo = 1, nxo
            do iyo = 1, nyo
                call closest2d(xxi, yyi, xxo(iyo,ixo), yyo(iyo,ixo), &
                    & nxi, nyi, imin, jmin, nogeo)
                varo(:,iyo,ixo) = vari(:, jmin, imin)
            enddo
!             exit
        enddo

     else
        if(nb<0)then
            znb = 10
        else
            znb = max(2,nb)
        endif
        znb2 = znb/2
        ixlast = 1
        iylast = 1
        ixlastline = 1
        iylastline = 1
        do ixo = 1, nxo
            do iyo = 1, nyo

                ! Try a small block
                ixmin = max(1, ixlast-znb2)
                ixmax = min(nxi, ixlast+znb2)
                iymin = max(1, iylast-znb2)
                iymax = min(nyi, iylast+znb2)
                call closest2d(xxi(iymin:iymax, ixmin:ixmax), &
                    & yyi(iymin:iymax, ixmin:ixmax), xxo(iyo, ixo), yyo(iyo, ixo), &
                    & ixmax-ixmin+1, iymax-iymin+1, imin, jmin, nogeo)
                imin = imin+ixmin-1
                jmin = jmin+iymin-1

                ! Fall on bounds so use full block
                if((imin==ixmin.and.ixmin/=1).or.(imin==ixmax.and.ixmax/=nxi).or.&
                    & (jmin==iymin.and.iymin/=1).or.(jmin==iymax.and.iymax/=nyi))&
                    & call closest2d(xxi, yyi, xxo(iyo, ixo), yyo(iyo, ixo), &
                    &   nxi, nyi, imin, jmin, nogeo)

                ! Store value
                varo(:,iyo,ixo) = vari(:, jmin, imin)

                ! Update min/max positions
                if(ixo==nxo)then
                    ixlast = ixlastline
                    iylast = iylastline
                else
                    ixlast = imin
                    iylast = jmin
                    if(ixo==1)then
                        ixlastline = ixlast
                        iylastline = iylast
                    endif
                endif
            enddo

         enddo
    endif
end subroutine nearest2d

subroutine closest2d(xxi, yyi, xo, yo, nxi, nyi, i, j, nogeo)
    ! Find indices of closest point on 2D axes
    implicit none
    real(kind=8),intent(in) :: xxi(nyi,nxi),yyi(nyi,nxi),xo,yo
    integer, intent(in):: nyi, nxi
    integer,intent(out) :: i,j
    integer,intent(in) :: nogeo
    real(kind=8) :: dx(nyi,nxi), mindist, dist
    integer:: ij(2), it, jt

    if(nogeo==1)then ! like meters
        dx = xxi-xo
        ij = minloc(dx**2+(yyi-yo)**2)
        i = ij(2)
        j = ij(1)
    else ! degrees
        mindist = 3.15d0
        do it=1,nxi
            do jt=1,nyi
                call haversine(xo, yo, xxi(jt, it), yyi(jt, it), dist)
                if(dist<mindist)then
                    i = it
                    j = jt
                    mindist = dist
                endif
            end do
        end do
    endif

end subroutine closest2d

subroutine haversine(lon0, lat0, lon1, lat1, dist)
    ! Haversine distance on a unit sphere
    real(kind=8), intent(in) :: lon0, lat0, lon1, lat1
    real(kind=8), intent(out) :: dist
    real(kind=8) :: pi
    deg2rad = acos(-1d0) / 180d0
    dist = 2d0 * asin(sqrt(sin(deg2rad*(lat0-lat1)*0.5d0)**2 + &
        & cos(deg2rad*lat0) * cos(deg2rad*lat1) * sin(deg2rad*(lon0-lon1)*0.5d0)**2))
end subroutine haversine

! =============================================================================

subroutine linear2d(vari, xi, yi, varo, xo,  yo, nogeo, nxi, nyi, nxo, nyo, nz)
    ! Simple bilinear interpolation between two regular grids
    !

    implicit none

    ! Parameters
    integer,intent(in) :: nxi,nyi,nxo,nyo,nz
    real(kind=8),intent(in) :: vari(nz,nyi,nxi)
    real(kind=8),intent(in) :: xi(nxi),yi(nyi),xo(nxo),yo(nyo)
    integer,intent(in), optional :: nogeo
    real(kind=8),intent(out) :: varo(nz,nyo,nxo)

    ! Local variables
    real(kind=8) :: dxi,dyi,fx,fy,fdxi,m1
    integer :: ixi,iyi,ixo,iyo
    logical :: geo


    ! Geographic
    geo = .not. present(nogeo) .or. nogeo==0

    ! Missing
    m1 = -1.
    varo = sqrt(m1)

    ! Loop on input x
    iyo = 1
    do iyi = 1, nyi-1

        ! Must overlap
        if (yi(iyi+1)<yo(1))cycle
        if (yi(iyi)>yo(nyo))exit

        ! Cell height
        dyi = yi(iyi+1)-yi(iyi)

        ! Loop on output y
        do while (iyo<=nyo.and.yo(iyo)<=yi(iyi+1))

            ! Still not inside
            if(yo(iyo)<yi(iyi))then
                iyo = iyo+1
                cycle
            endif

            ! Y weight
            fy = (yo(iyo)-yi(iyi))/dyi

            ! Loop on input x
            ixo = 1
            do ixi = 1, nxi-1

                ! Must overlap
                if (xi(ixi+1)<xo(1))cycle
                if (xi(ixi)>xo(nxo))exit

                ! Cell width
                dxi = xi(ixi+1)-xi(ixi)
                if(geo.and.dxi>180.)dxi=360.-dxi

                ! Loop on output x
                do while (ixo<=nxo.and.xo(ixo)<=xi(ixi+1))

                    ! Still not inside
                    if(xo(ixo)<xi(ixi))then
                        ixo = ixo+1
                        cycle
                    endif

                    ! X weight
                    fdxi = xo(ixo)-xi(ixi)
                    if(geo.and.fdxi>180.)fdxi=360.-fdxi
                    fx = fdxi/dxi

                    ! Interpolation
                    varo(:,iyo,ixo) = &
                        & vari(:,  iyi  ,ixi)*(1-fx)*(1-fy) + &
                        & vari(:,iyi  ,ixi+1)*fx    *(1-fy) + &
                        & vari(:,iyi+1,ixi)*(1-fx)*fy     + &
                        & vari(:,iyi+1,ixi+1)*fx    *fy

                    ! Next output x
                    ixo = ixo+1
                enddo
            enddo

            ! Next output y
            iyo = iyo+1
        enddo
    enddo
end subroutine linear2d

subroutine dstwgt2d(vari, xi, yi, varo, xo,  yo, nogeo, nxi, nyi, nxo, nyo, nz)
    ! Simple distance weight interpolation between two regular grids
    ! It does not interpolate missing values

    implicit none

    ! Parameters
    integer,intent(in) :: nxi,nyi,nxo,nyo,nz
    real(kind=8),intent(in) :: vari(nz,nyi,nxi)
    real(kind=8),intent(in) :: xi(nxi),yi(nyi),xo(nxo),yo(nyo)
    integer,intent(in), optional :: nogeo
    real(kind=8),intent(out) :: varo(nz,nyo,nxo)

    ! Local variables
    real(kind=8) :: small,m1
    real(kind=8) :: dxi,dyi,dd(4),ww(nz,4),dx0,dx1,dy0,dy1,wsum(nz),vv(nz,4),mvs(nz)
    integer :: ixi,iyi,ixo,iyo,i4
    logical :: geo

    geo  = .not.present(nogeo) .or. nogeo==0


    ! Missing
    m1 = -1.
    varo = sqrt(m1)
    vv = sqrt(m1)
    small = epsilon(1d0)*1.1

    ! Loop on input x
    iyo = 1
    do iyi = 1, nyi-1

        ! Must overlap
        if (yi(iyi+1)<yo(1))cycle
        if (yi(iyi)>yo(nyo))exit

        ! Cell height
        dyi = yi(iyi+1)-yi(iyi)

        ! Loop on output y
         do while (iyo<=nyo.and.yo(iyo)<yi(iyi+1))

            ! Still not inside
            if(yo(iyo)<yi(iyi))then
                iyo = iyo+1
                cycle
            endif

            ! Y axis
            dy0 = yo(iyo)-yi(iyi)
            dy1 = yi(iyi+1)-yo(iyo)

            ! Loop on input x
            ixo = 1
            do ixi = 1, nxi-1

                ! Must overlap
                if (xi(ixi+1)<xo(1))cycle
                if (xi(ixi)>xo(nxo))exit

                ! Cell width
                dxi = xi(ixi+1)-xi(ixi)
                if(geo.and.dxi>180.)dxi=360.-dxi

                ! Loop on output x
                do while (ixo<=nxo.and.xo(ixo)<xi(ixi+1))

                    ! Still not inside
                    if(xo(ixo)<xi(ixi))then
                        ixo = ixo+1
                        cycle
                    endif

                    ! X axis
                    dx0 = xo(ixo)-xi(ixi)
                    dx1 = xi(ixi+1)-xo(ixo)
                    if(geo)then
                        if(dx0>180.)dx0=360.-dx0
                        if(dx1>180.)dx1=360.-dx1
                        dx0 = dx0*cos(yo(iyo)*3.14159d0/180.)
                        dx1 = dx1*cos(yo(iyo)*3.14159d0/180.)
                    endif

                    ! Distances and weights
                    dd = (/sqrt(dx0**2+dy0**2),sqrt(dx1**2+dy0**2),&
                    &      sqrt(dx0**2+dy1**2),sqrt(dx1**2+dy1**2)/)
                    ww = 0.
                    vv(:,1:2) = vari(:,iyi,ixi:ixi+1)
                    vv(:,3:4) = vari(:,iyi+1,ixi:ixi+1)
                    where(dd==0.)dd=tiny(1.)
                    do i4=1,4
                        ww(:,i4) = 1./dd(i4)
                    enddo
                    where(isnan(vv))
                        ww = 0.
                        vv = 0.
                    endwhere
                    wsum = sum(ww,dim=2)!,mask=vv/=mv)
                    mvs = 0.
                    where(wsum==0.)
                        wsum = 1.
                        mvs = sqrt(m1)
                    endwhere

                    ! Interpolation
                     varo(:,iyo,ixo) = sum(ww*vv, dim=2)/wsum + mvs

                    ! Next output x
                    ixo = ixo+1
                enddo
            enddo

            ! Next output y
            iyo = iyo+1
        enddo
    enddo
end subroutine dstwgt2d



subroutine mlinear2d (vari, xi,  yi,  varo, xo,  yo, ext,  nxi, nyi, no,  nogeo)
    ! Binilear interpolation of random point to a regular grid
    ! with minimal missing data handling

    implicit none

    ! Parameters
    integer, intent(in) :: nxi, nyi, no
    real(kind=8),    intent(in) :: xi(nxi), yi(nyi), vari(nyi, nxi), xo(no), yo(no)
    logical, intent(in) :: ext
    integer, intent(in), optional :: nogeo
    real(kind=8),    intent(out) :: varo(no)

    ! Local parameters
    integer :: ix, iy, io, ix1, iy1
    real(kind=8) :: ww(2,2), zz(2,2), dx, dy, fx, fy, pi, zxi(nxi),zxo(no), m1
    logical ::  ms(2,2), ma(2,2),geo


    ! Inits
    pi = 3.14159d0
    m1 = -1.
    varo = sqrt(m1)

    ! Geo
    geo = .not. present(nogeo) .or. nogeo==0
    if(geo)then
        zxi = modulo(zxi,360.)
        zxo = modulo(zxo,360.)
    else
        zxi = xi
        zxo = xo
    endif

    ! Loop on output points
    do io = 1, no

        ! Check bounds
        if (zxo(io)<zxi(1).or.zxo(io)>zxi(nxi)) cycle
        if (yo(io)<yi(1).or.yo(io)>yi(nyi)) cycle

        ! Find indices for current position
        ix = count(zxi<=zxo(io))
        iy = count(yi<=yo(io))

        ! Get bloc of data
        ix1 = min(ix+1,nxi)
        iy1 = min(iy+1,nyi)
        zz(1,1) = vari(iy,ix)
        zz(1,2) = vari(iy,ix1)
        zz(2,1) = vari(iy1,ix)
        zz(2,2) = vari(iy1,ix1)

        ! Check that we have all needed points
        ms = isnan(zz)
        if(all(ms).or.(.not.ext.and.any(ms)))cycle

        ! Cell size
        dx = abs(zxi(ix1)-zxi(ix))
        if(ix==ix1)then
            dx = 1.
        elseif(geo)then
            if(dx>180.)dx = 360.-dx
            dx = dx*cos((yi(iy)+yi(iy1))*.5)
        endif
        dy = yi(iy1)-yi(iy)
        if(iy==iy1)dy = 1.

        ! Guess weights
        fx = zxo(io)-zxi(ix)
        if(geo) fx = fx*cos((yi(iy)+yi(iy1))*.5)
        fy = yo(io)-yi(iy)
        ww(1,1) = (dy-fy)*(dx-fx)
        ww(1,2) = (dy-fy)*fx
        ww(2,1) = fy*(dx-fx)
        ww(2,2) = fy*fx

        ! Average
        ma = ww/=0..and..not.ms
        if(all(.not.ma))cycle
        varo(io) = sum(zz*ww,mask=ma)/sum(ww,mask=ma)

    enddo

end subroutine mlinear2d

! =============================================================================

subroutine nearest2dto1d(xi,yi,zi,xo,yo,zo,nxi,nyi,no,nz)
    ! nearest neighbour interpolation of gridded data to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xi(nxi), yi(nyi), xo(no), yo(no)
    real(kind=8),intent(in) :: zi(nz,nyi,nxi)
    real(kind=8),intent(out) :: zo(nz,no)
    real(kind=8) :: dx(nxi)

    integer :: io,i,j
    real(kind=8) :: m1

    m1 = -1.
    zo = sqrt(m1)
    !$OMP PARALLEL DO PRIVATE(io,i,j,dx)
    !$& SHARED(xi,yi,zi,xo,yo,zo,nxi,nyi,no)
    do io = 1, no
!        if(xo(io)<xi(1).or.xo(io)>xi(nxi))cycle
!        if(yo(io)<yi(1).or.yo(io)>yi(nyi))cycle
        dx = abs(xi-xo(io))
        where(dx>180.)dx=360.-dx
        i = minloc(dx,dim=1)
        j = minloc(abs(yi-yo(io)),dim=1)

        zo(:,io) = zi(:,j,i)
    enddo
    !$OMP END PARALLEL DO
end subroutine nearest2dto1d

subroutine linear2dto1d(xi,yi,zi,xo,yo,zo,nxi,nyi,no,nz)
    ! bilinear interpolation of gridded data to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xi(nxi), yi(nyi), xo(no), yo(no)
    real(kind=8),intent(in) :: zi(nz,nyi,nxi)
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j
    real(kind=8) :: a,b,m1

    m1 = -1.
    zo = sqrt(m1)
    !$OMP PARALLEL DO PRIVATE(io,i,j,a,b)
    !$& SHARED(xi,yi,zi,xo,yo,zo,nxi,nyi,no)
    do io = 1, no
        if(xo(io)>=xi(1).and.xo(io)<=xi(nxi).and.&
            & yo(io)>=yi(1).and.yo(io)<=yi(nyi))then

            if(xi(nxi)==xo(io))then
                i = nxi-1
            else
                i = minloc(xi,dim=1,mask=xi>xo(io))-1
            endif
            if(yi(nyi)==yo(io))then
                j = nyi-1
            else
                j = minloc(yi,dim=1,mask=yi>yo(io))-1
            endif

            a = abs(xo(io)-xi(i))
            if(a>180.)a=360.-180.

            a = a/(xi(i+1)-xi(i))
            b = (yo(io)-yi(j))/(yi(j+1)-yi(j))

            zo(:, io) = (1-b)*(1-a)*zi(:,j,  i) + &
                &       (1-b)*a*zi(:,j,  i+1) + &
                &        b*(1-a)*zi(:,j+1,i) + &
                &        b*a*zi(:,j+1,i+1)

        endif

    enddo
    !$OMP END PARALLEL DO
end subroutine linear2dto1d



subroutine dstwgt2dto1d(xi,yi,zi,xo,yo,zo,nxi,nyi,no,nz)
    ! Distance weight interpolation of gridded data to random positions
    !
    ! Distances are computed with the four corners of a cell and are relative
    ! to the cell sizes.

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xi(nxi), yi(nyi), xo(no), yo(no)
    real(kind=8),intent(in) :: zi(nz,nyi,nxi)
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j,i4
    real(kind=8) :: dx0,dx1,dy0,dy1,dd(4),vv(nz,4),wsum(nz),ww(nz,4),m1
    logical :: bmask(nz,nyi,nxi),bb(nz,4)

    m1 = -1
    zo = sqrt(m1)
    bmask = isnan(zi)
    !$OMP PARALLEL DO PRIVATE(io,i,j,dx0,dx1,dy0,dy1,vv,bb,dd,ww,wsum,i4)
    !$& SHARED(xi,yi,zi,xo,yo,zo,nxi,nyi,no,bmask)
    do io = 1, no

        if(xo(io)>=xi(1).and.xo(io)<=xi(nxi).and.&
            & yo(io)>=yi(1).and.yo(io)<=yi(nyi))then

            ! Indices
            if(xi(nxi)==xo(io))then
                i = nxi-1
            else
                i = minloc(xi,dim=1,mask=xi>xo(io))-1
            endif
            if(yi(nyi)==yo(io))then
                j = nyi-1
            else
                j = minloc(yi,dim=1,mask=yi>yo(io))-1
            endif

            ! Distances
            dx0 = (xo(io)-xi(i))/(xi(i+1)-xi(i))
            dx1 = 1d0-dx0
            dy0 = (yo(io)-yi(i))/(yi(i+1)-yi(i))
            dy1 = 1d0-dy0
            dd = (/sqrt(dx0**2+dy0**2),sqrt(dx0**2+dy1**2),&
            &      sqrt(dx1**2+dy0**2),sqrt(dx1**2+dy1**2)/)
            where(dd==0d0)dd=tiny(0d0)

            ! Weights
            vv = reshape(zi(:,j:j+1,i:i+1),(/nz,4/))
            bb = reshape(bmask(:,j:j+1,i:i+1),(/nz,4/))
            ww = 0d0
            do i4=1,4
                ww(:,i4) = 1d0/dd(i4)
            enddo
            where(bb)
                ww=0d0
                vv=0d0
            endwhere
            wsum = merge(1d0, sum(ww,dim=2), all(bb,dim=2))

            ! Interpolation
            zo(:,io) = merge(sqrt(m1), sum(ww*vv, dim=2)/wsum, all(bb,dim=2))

        endif

    enddo
    !$OMP END PARALLEL DO

end subroutine dstwgt2dto1d

subroutine linear4dto1d(xi,yi,zi,ti,vi,xo,yo,zo,to,vo,nxi,nyi,nzi,nti,no)
    ! nearest neighbour interpolation of gridded data to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,nzi,nti,no
    real(kind=8),intent(in) :: xi(nxi), yi(nyi), zi(nzi), ti(nti)
    real(kind=8),intent(in) :: xo(no), yo(no), zo(no), to(no)
    real(kind=8),intent(in) :: vi(nti,nzi,nyi,nxi)
    real(kind=8),intent(out) :: vo(no)

    real(kind=8) :: a , b, c, d
    real(kind=8) :: ximin, yimin, zimin, timin, ximax, yimax, zimax, timax, m1
    logical :: bmask(nti,nzi,nyi,nxi)
    integer :: io, i, j, k, l, ii, jj, kk, ll, npi, npj, npk, npl

    m1 = -1.
    vo = sqrt(m1)
    ximin = minval(xi)
    ximax = maxval(xi)
    yimin = minval(yi)
    yimax = maxval(yi)
    zimin = minval(zi)
    zimax = maxval(zi)
    timin = minval(ti)
    timax = maxval(ti)

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(io,i,j,k,l,ii,jj,kk,ll,a,b,c,d)&
    !$OMP PRIVATE(npi,npj,npl,npk)
    do io = 1, no
        if(       (nxi==1 .or. (xo(io)>=ximin.and.xo(io)<=ximax)) .and.&
                & (nyi==1 .or. (yo(io)>=yimin.and.yo(io)<=yimax)) .and.&
                & (nzi==1 .or. (zo(io)>=zimin.and.zo(io)<=zimax)) .and.&
                & (nti==1 .or. (to(io)>=timin.and.to(io)<=timax)))then

            ! Weights
            ! - X
            if(nxi==1)then
                i = 1
                a = 0d0
                npi = 1
            else if(xi(nxi)==xo(io))then
                i = nxi
                npi = 1
                a = 0d0
            else
                i = minloc(xi,dim=1,mask=xi>xo(io))-1
                npi = 2
                a = abs(xo(io)-xi(i))
                if(a>180d0)a=a-180d0
                a = a/(xi(i+1)-xi(i))
            endif
            ! - Y
            if(nyi==1)then
                j = 1
                b = 0d0
                npj = 1
            else if(yi(nyi)==yo(io))then
                j = nyi
                b  = 0d0
                npj = 1
            else
                j = minloc(yi, dim=1, mask=yi>yo(io)) - 1
                b = (yo(io)-yi(j))/(yi(j+1)-yi(j))
                npj = 2
            endif
            ! - Z
            if(nzi==1)then
                k = 1
                c = 0d0
                npk = 1
            else if(zi(nzi)==zo(io))then
                k = nzi
                c = 0d0
                npk = 1
            else
                k = minloc(zi,dim=1,mask=zi>zo(io))-1
                if(zi(k+1)==zi(k))then
                    c = 0d0
                    npk = 1
                else
                    c = (zo(io)-zi(k))/(zi(k+1)-zi(k))
                    npk = 2
                endif
            endif
            ! - T
            if(nti==1)then
                l = 1
                d = 0d0
                npl = 1
            else if(ti(nti)==to(io))then
                l = nti
                d = 0d0
                npl = 1
            else
                l = minloc(ti, dim=1, mask=ti>to(io))-1
                if(ti(l+1)==ti(l))then
                    d = 0d0
                    npl = 1
                else
                    d = (to(io)-ti(l))/(ti(l+1)-ti(l))
                    npl = 2
                endif
            endif

            ! Interpolate
            vo(io) = 0d0
            do ll=0,npl-1
                do kk=0,npk-1
                    do jj=0,npj-1
                        do ii=0,npi-1
                            vo(io) = vo(io) +vi(l+ll, k+kk, j+jj, i+ii) * &
                                & ((1-a) * (1-ii) + a * ii)* &
                                & ((1-b) * (1-jj) + b * jj)* &
                                & ((1-c) * (1-kk) + c * kk)* &
                                & ((1-d) * (1-ll) + d * ll)
                        enddo
                    enddo
                enddo
            enddo
        endif
    enddo
    !$OMP END PARALLEL DO

end subroutine linear4dto1d

subroutine linear4dto1dx(xi,yi,zi,ti,vi,xo,yo,zo,to,vo,nex,nxi,nyi,nzi,nti,no)
    ! nearest neighbour interpolation of gridded data to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,nzi,nti,no,nex
    real(kind=8),intent(in) :: xi(nxi), yi(nyi), zi(nzi), ti(nti)
    real(kind=8),intent(in) :: xo(no), yo(no), zo(no), to(no)
    real(kind=8),intent(in) :: vi(nex,nti,nzi,nyi,nxi)
    real(kind=8),intent(out) :: vo(nex,no)

    real(kind=8) :: a , b, c, d, m1
    logical :: bmask(nex,nti,nzi,nyi,nxi)
    integer :: io, i, j, k, l, ii, jj, kk, ll, npi, npj, npk, npl

    m1 = -1
    vo = sqrt(m1)

    !$OMP PARALLEL DO PRIVATE(io,i,j,k,l,ii,jj,kk,ll)
    !$& SHARED(xi,yi,zi,ti,vi,xo,yo,zo,vo,nei,nxi,nyi,nzi,nti,no)
    do io = 1, no
        if((nxi==1 .or. (xo(io)>=xi(1).and.xo(io)<=xi(nxi))) .and.&
                & (nyi==1 .or. (yo(io)>=yi(1).and.yo(io)<=yi(nyi))) .and.&
                & (nzi==1 .or. (zo(io)>=zi(1).and.zo(io)<=zi(nzi))) .and.&
                & (nti==1 .or. (to(io)>=ti(1).and.to(io)<=ti(nti))))then

            ! Weights
            ! - X
            if(nxi==1)then
                i = 1
                a = 0d0
                npi = 1
            else if(xi(nxi)==xo(io))then
                i = nxi
                npi = 1
                a = 0d0
            else
                i = minloc(xi,dim=1,mask=xi>xo(io))-1
                npi = 2
                a = abs(xo(io)-xi(i))
                if(a>180d0)a=360d0-180d0
                a = a/(xi(i+1)-xi(i))
            endif
            ! - Y
            if(nyi==1)then
                j = 1
                b = 0d0
                npj = 1
            else if(yi(nyi)==yo(io))then
                j = nyi
                b  = 0d0
                npj = 1
            else
                j = minloc(yi, dim=1, mask=yi>yo(io)) - 1
                b = (yo(io)-yi(j))/(yi(j+1)-yi(j))
                npj = 2
            endif
            ! - Z
            if(nzi==1)then
                k = 1
                c = 0d0
                npk = 1
            else if(zi(nzi)==zo(io))then
                k = nzi
                c = 0d0
                npk = 1
            else
                k = minloc(zi,dim=1,mask=zi>zo(io))-1
                if(zi(k+1)==zi(k))then
                    c = 0d0
                    npk = 1
                else
                    c = (zo(io)-zi(k))/(zi(k+1)-zi(k))
                    npk = 2
                endif
            endif
            ! - T
            if(nti==1)then
                l = 1
                d = 0d0
                npl = 1
            else if(ti(nti)==to(io))then
                l = nti
                d = 0d0
                npl = 1
            else
                l = minloc(ti, dim=1, mask=ti>to(io))-1
                if(ti(l+1)==ti(l))then
                    d = 0d0
                    npl = 1
                else
                    d = (to(io)-ti(l))/(ti(l+1)-ti(l))
                    npl = 2
                endif
            endif

            ! Interpolate
            vo(:,io) = 0d0
            do ii=0,npi-1
                do jj=0,npj-1
                    do kk=0,npk-1
                        do ll=0,npl-1
                            vo(:,io) = vo(:,io) +vi(:,l+ll, k+kk, j+jj, i+ii) * &
                                & ((1-a) * (1-ii) + a * ii)* &
                                & ((1-b) * (1-jj) + b * jj)* &
                                & ((1-c) * (1-kk) + c * kk)* &
                                & ((1-d) * (1-ll) + d * ll)
                        enddo
                    enddo
                enddo
            enddo

        endif
    enddo
    !$OMP END PARALLEL DO

end subroutine linear4dto1dx

subroutine linear4dto1dxx(xxi,yyi,zzi,ti,vi,xo,yo,zo,to,vo,&
        nxi,nyi,nyix,nxiy, nyiz,nxiz,nzi, nti,ntiz, no,nex,nexz)
    ! linear interpolation of gridded data to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,nyix,nxiy, nyiz,nxiz,nzi, nti,ntiz, no,nex,nexz
    real(kind=8),intent(in) :: xxi(nyix,nxi), yyi(nyi,nxiy), zzi(nexz,ntiz,nzi,nyiz,nxiz), ti(nti)
    real(kind=8),intent(in) :: xo(no), yo(no), zo(no), to(no)
    real(kind=8),intent(in) :: vi(nex,nti,nzi,nyi,nxi)
    real(kind=8),intent(out) :: vo(nex,no)

    real(kind=8) :: a , b, c(nexz), d, p, q, zi(nexz,nzi), az, bz, dz, m1
    real(kind=8) :: ximin, yimin, zimin, timin, ximax, yimax, zimax, timax
    logical :: bmask(nex,nti,nzi,nyi,nxi), curved
    integer :: io, i, j, k(nexz), l, ii, jj, kk, ll, npi, npj, npk(nexz), npl, &
        & npiz, npjz, nplz, iz, jz, lz, iez, ieb, ie0, ie1

    m1 = -1d0
    vo = sqrt(m1)
    bmask = isnan(vi)
    ximin = minval(xxi)
    ximax = maxval(xxi)
    yimin = minval(yyi)
    yimax = maxval(yyi)
    zimin = minval(zzi)
    zimax = maxval(zzi)
    timin = minval(ti)
    timax = maxval(ti)
    curved = nyix/=1
    if(curved .and. nxi/=nxiy .and. nyi/=nyix)stop "linear4dto1: Invalid curved dimensions"
    if(nxiz/=1 .and. nxiz/=nxi)stop "linear4dto1: Invalid nxiz dimension"
    if(nyiz/=1 .and. nyiz/=nyi)stop "linear4dto1: Invalid nyiz dimension"
    if(ntiz/=1 .and. ntiz/=nti)stop "linear4dto1: Invalid ntiz dimension"


    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(io,i,j,k,l,a,b,c,d,ii,jj,kk,ll,p,q,npi,npj,npk,npl) &
    !$OMP PRIVATE(npiz,npjz,nplz,iz,jz,lz,az,bz,dz,ieb,ie0,iez,zi)
    do io = 1, no
        if(       (nxi==1 .or. (xo(io)>=ximin.and.xo(io)<=ximax)) .and.&
                & (nyi==1 .or. (yo(io)>=yimin.and.yo(io)<=yimax)) .and.&
                & (nzi==1 .or. (zo(io)>=zimin.and.zo(io)<=zimax)) .and.&
                & (nti==1 .or. (to(io)>=timin.and.to(io)<=timax)))then

            ! Weights
            if(curved)then
                call curv2rel_single(xxi, yyi, xo(io), yo(io), p, q, nxi, nyi)
                if(p<1 .or. p>nxi .or. q<1 .or. q>nyi) continue
                i = int(p)
                j = int(q)
                a = p - i
                b = q - j
                npi = 2
                npj = 2

            else
                ! - X
                if(nxi==1)then
                    i = 1
                    a = 0d0
                    npi = 1
                else if(xxi(nxi,1)==xo(io))then
                    i = nxi
                    npi = 1
                    a = 0d0
                else
                    i = minloc(xxi(1,:), dim=1, mask=xxi(1,:)>xo(io))-1
                    npi = 2
                    a = xo(io)-xxi(1,i)
                    if(abs(a)>180d0)a=a-180d0 ! FIXME: linear4dto1dxx: abs(a)>180d0
                    a = a/(xxi(1,i+1)-xxi(1,i))
                endif

                ! - Y
                if(nyi==1)then
                    j = 1
                    b = 0d0
                    npj = 1
                else if(yyi(1,nyi)==yo(io))then
                    j = nyi
                    b  = 0d0
                    npj = 1
                else
                    j = minloc(yyi(:,1), dim=1, mask=yyi(:,1)>yo(io)) - 1
                    b = (yo(io)-yyi(j,1))/(yyi(j+1,1)-yyi(j,1))
                    npj = 2
                endif

            endif

            ! - T
            if(nti==1)then
                l = 1
                d = 0d0
                npl = 1
            else if(ti(nti)==to(io))then
                l = nti
                d = 0d0
                npl = 1
            else
                l = minloc(ti, dim=1, mask=ti>to(io))-1
                if(ti(l+1)==ti(l))then
                    d = 0d0
                    npl = 1
                else
                    d = (to(io)-ti(l))/(ti(l+1)-ti(l))
                    npl = 2
                endif
            endif

            ! - Z
            if(nzi==1)then
                k = 1
                c = 0d0
                npk = 1
            else

                ! Local zi

                if(nxiz==1)then
                    npiz = 1
                    az = 0
                    iz = 1
                else
                    npiz = npi
                    az = a
                    iz = i
                endif

                if(nyiz==1)then
                    npjz = 1
                    bz = 0
                    jz = 1
                else
                    npjz = npj
                    bz = b
                    jz = j
                endif

                if(ntiz==1)then
                    nplz = 1
                    dz = 0
                    lz = 1
                else
                    nplz = npl
                    dz = d
                    lz = l
                endif

                zi = 0d0
                do ll=0,nplz-1
                    do jj=0,npjz-1
                        do ii=0,npiz-1
                            zi = zi + zzi(:, lz+ll, :, jz+jj, iz+ii) * &
                                & ((1-az) * (1-ii) + az * ii)* &
                                & ((1-bz) * (1-jj) + bz * jj)* &
                                & ((1-dz) * (1-ll) + dz * ll)
                        enddo
                    enddo
                enddo

                ! Normal stuff (c(nexz),zi(nexz,nzi),k(nexz)
                do iez = 1, nexz ! extra dim
                    if(zi(iez,nzi)==zo(io))then
                        k(iez) = nzi
                        c(iez) = 0d0
                        npk(iez) = 1
                    else
                        k(iez) = minloc(zi(iez,:), dim=1, mask=zi(iez,:)>zo(io))-1
                        if(zi(iez,k(iez)+1)==zi(iez,k(iez)))then
                            c(iez) = 0d0
                            npk(iez) = 1
                        else
                            c(iez) = (zo(io)-zi(iez,k(iez))) / (zi(iez,k(iez)+1)-zi(iez,k(iez)))
                            npk(iez) = 2
                        endif
                    endif
                enddo


            endif

            ! Interpolate
            vo(:,io) = 0d0
            do ieb = 1, nex/nexz

                ie0 = 1+(ieb-1)*nexz

                do iez=0,nexz-1
                    if(.not. any(bmask(ie0:ie0+iez,l:l+npl-1, k(iez+1):k(iez+1)+npk(iez+1)-1, j:j+npj-1, i:i+npi-1)))then
                        do ll=0,npl-1
                            do kk=0,npk(iez+1)-1
                                do jj=0,npj-1
                                    do ii=0,npi-1
                                        vo(ie0+iez,io) = vo(ie0+iez,io) +vi(ie0+iez,l+ll, k(iez+1)+kk, j+jj, i+ii) * &
                                            & ((1-a) * (1-ii) + a * ii)* &
                                            & ((1-b) * (1-jj) + b * jj)* &
                                            & ((1-c(iez+1)) * (1-kk) + c(iez+1) * kk)* &
                                            & ((1-d) * (1-ll) + d * ll)
                                    enddo
                                enddo
                            enddo
                        enddo

                    endif
                enddo

            end do

        endif
    enddo
    !$OMP END PARALLEL DO

end subroutine linear4dto1dxx


subroutine nearest4dto1dxx(xxi,yyi,zzi,ti,vi,xo,yo,zo,to,vo,&
        nxi,nyi,nyix,nxiy, nyiz,nxiz,nzi, nti,ntiz, no,nex,nexz)
    ! linear interpolation of gridded data to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,nyix,nxiy, nyiz,nxiz,nzi, nti,ntiz, no,nex,nexz
    real(kind=8),intent(in) :: xxi(nyix,nxi), yyi(nyi,nxiy), zzi(nexz,ntiz,nzi,nyiz,nxiz), ti(nti)
    real(kind=8),intent(in) :: xo(no), yo(no), zo(no), to(no)
    real(kind=8),intent(in) :: vi(nex,nti,nzi,nyi,nxi)
    real(kind=8),intent(out) :: vo(nex,no)

    real(kind=8) :: p, q, zi(nexz,nzi), dx(nxi),m1
    logical ::  curved
    integer :: io, i, j, k(nexz), l, iz, jz, lz, iez, ieb, ie0, ie1

    m1 = -1d0
    vo = sqrt(m1)
    curved = nyix/=1
    if(curved .and. nxi/=nxiy .and. nyi/=nyix)stop "linear4dto1: Invalid curved dimensions"
    if(nxiz/=1 .and. nxiz/=nxi)stop "linear4dto1: Invalid nxiz dimension"
    if(nyiz/=1 .and. nyiz/=nyi)stop "linear4dto1: Invalid nyiz dimension"
    if(ntiz/=1 .and. ntiz/=nti)stop "linear4dto1: Invalid ntiz dimension"

!    print*,'xxi',xxi
!    print*,'yyi',yyi
!    print*,'ti',ti
!    print*,'zzi',zzi


    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(io,i,j,k,l,p,q) &
    !$OMP PRIVATE(iz,jz,lz,ieb,ie0,iez,zi)
    do io = 1, no

            ! Weights
            if(curved)then
!print*,'curved'
                call curv2rel_single(xxi, yyi, xo(io), yo(io), p, q, nxi, nyi)
                p = max(p, 1d0)
                p = min(p, dble(nxi))
                q = max(q, 1d0)
                q = min(q, dble(nyi))
                i = int(p)
                j = int(q)
                if(p-i>.5)i = i + .5
                if(q-j>.5)j = j + .5

            else
                ! - X
                if(nxi==1)then
                    i = 1
                else
                    dx = abs(xo(io)-xxi(1,:))
                    where(dx>180.)dx=360.-dx
                    i = minloc(dx, dim=1)
                endif

                ! - Y
                if(nyi==1)then
                    j = 1
                else
                    j = minloc(abs(yo(io)-yyi(:,1)), dim=1)
                endif

            endif

            ! - T
!            print*,'nti',nti
            if(nti==1)then
                l = 1
            else
                l = minloc(abs(to(io)-ti), dim=1)
            endif

            ! - Z
            if(nzi==1)then
                k = 1
            else

                ! Local zi

                if(nxiz==1)then
                    iz = 1
                else
                    iz = i
                endif

                if(nyiz==1)then
                    jz = 1
                else
                    jz = j
                endif

                if(ntiz==1)then
                    lz = 1
                else
                    lz = l
                endif

                zi = zzi(:, lz, :, jz, iz)
!                print*,'zi',zi

                ! Normal stuff zi(nexz,nzi),k(nexz)
                do iez = 1, nexz ! extra dim
                    k(iez) = minloc(abs(zo(io)-zi(iez,:)), dim=1)
                enddo


            endif

!            print*,'X i',i
!            print*,'Y j',j
!            print*,'Z k',k
!            print*,'T l',l

            ! Interpolate
            vo(:,io) = 0d0
!            print*,'nex/nexz',nex/nexz
            do ieb = 1, nex/nexz
!                print*,'ieb',ieb

                ie0 = 1+(ieb-1)*nexz
!                ie1 = ie0+nexz-1

                do iez=0,nexz-1
                    vo(ie0+iez,io) = vi(ie0+iez,l, k(iez+1), j, i)
                enddo

            end do


    enddo
    !$OMP END PARALLEL DO

end subroutine nearest4dto1dxx


subroutine nearest4dto1d(xi,yi,zi,ti,vi,xo,yo,zo,to,vo,nxi,nyi,nzi,nti,no)
    ! nearest neighbour interpolation of gridded data to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,nzi,nti,no
    real(kind=8),intent(in) :: xi(nxi), yi(nyi), zi(nzi), ti(nti)
    real(kind=8),intent(in) :: xo(no), yo(no), zo(no), to(no)
    real(kind=8),intent(in) :: vi(nti,nyi,nyi,nxi)
    real(kind=8),intent(out) :: vo(no)

    real(kind=8) :: dx(nxi),m1
    integer :: io,i,j,k,l

    m1 = -1d0
    vo = sqrt(m1)
    k = 1
    l = 1
    !$OMP PARALLEL DO PRIVATE(io,i,j,k,l,dx)
    !$& SHARED(xi,yi,zi,ti,vi,xo,yo,zo,nxi,nyi,nzi,nti,no)
    do io = 1, no
!        if(xo(io)<xi(1).or.xo(io)>xi(nxi))cycle
!        if(yo(io)<yi(1).or.yo(io)>yi(nyi))cycle
        dx = abs(xi-xo(io))
        where(dx>180.)dx=360.-dx
        i = minloc(dx,dim=1)
        j = minloc(abs(yi-yo(io)),dim=1)
        if(nzi/=1)k = minloc(abs(zi-zo(io)),dim=1)
        if(nti/=1)l = minloc(abs(ti-to(io)),dim=1)

        vo(io) = vi(l,k,j,i)
    enddo
    !$OMP END PARALLEL DO

end subroutine nearest4dto1d



subroutine linept(x,y,x1,x2,y1,y2,xc,yc)
    ! Coordinates of line point closest to target point
    implicit none
    real(kind=8),intent(in) :: x,y,x1,x2,y1,y2
    real(kind=8), intent(out) :: xc, yc
    real(kind=8) :: dx, dy
    dy = y2-y1
    dx = x2-x1
    xc = x1 + (dx**2*(x-x1) + dx*dy*(y-y1))/(dx**2+dy**2)
    yc = y + (x-xc)*dy/dx
end subroutine linept

subroutine linepts(xx,yy,x1,x2,y1,y2,xxc,yyc,np)
    ! Coordinates of line points closest to target points
    implicit none
    integer, intent(in) :: np
    real(kind=8),intent(in) :: xx(np),yy(np),x1,x2,y1,y2
    real(kind=8), intent(out) :: xxc(np), yyc(np)
    integer :: i
    do i=1,np
        call linept(xx(i),yy(i),x1,x2,y1,y2,xxc(i),yyc(i))
    enddo
end subroutine linepts

subroutine lineptss(xx,yy,xx1,xx2,yy1,yy2,xxc,yyc,np)
    ! Coordinates of lines points closest to target points
    implicit none
    integer, intent(in) :: np
    real(kind=8),intent(in) :: xx(np),yy(np),xx1(np),xx2(np),yy1(np),yy2(np)
    real(kind=8), intent(out) :: xxc(np), yyc(np)
    integer :: i
    do i=1,np
        call linept(xx(i),yy(i),xx1(np),xx2(np),yy1(np),yy2(np),xxc(i),yyc(i))
    enddo
end subroutine lineptss

subroutine dstpt2line(x,y,x1,x2,y1,y2,d)
    ! Distance from a point to a line
    real(kind=8),intent(in) :: x,y,x1,x2,y1,y2
    real(kind=8), intent(out) :: d
    real(kind=8) :: xc,yc
    call linept(x,y,x1,x2,y1,y2,xc,yc)
    d = sqrt((xc-x)**2+(yc-y)**2)
end subroutine dstpt2line

subroutine dstpts2line(x,y,x1,x2,y1,y2,d,np)
    ! Distance from points to a line
    real(kind=8),intent(in) :: x(np),y(np),x1,x2,y1,y2
    real(kind=8), intent(out) :: d(np)
    real(kind=8) :: xc(np),yc(np)
    call linepts(x,y,x1,x2,y1,y2,xc,yc,np)
    d = sqrt((xc-x)**2+(yc-y)**2)
end subroutine dstpts2line

subroutine dstpts2lines(x,y,x1,x2,y1,y2,d,np)
    ! Distance from points to lines
    real(kind=8),intent(in),dimension(np) :: x,y,x1,x2,y1,y2
    real(kind=8), intent(out) :: d(np)
    real(kind=8) :: xc(np),yc(np)
    call lineptss(x,y,x1,x2,y1,y2,xc,yc,np)
    d = sqrt((xc-x)**2+(yc-y)**2)
end subroutine dstpts2lines

subroutine curv2rect(x1,x2,x3,x4,y1,y2,y3,y4,x,y,p,q)
    ! Coordinate transform from curvilinear to rectangular
    !
    ! Cell shape:
    !
    !   2 - 3
    !   |   |
    !   1 - 4
    !
    !:Source: http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19890018062_1989018062.pdf

    implicit none

    real(kind=8), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4,x,y
    real(kind=8), intent(out) :: p,q
    real(kind=8) :: p1,p2 ,q1,q2, AA, BB, CC, DD, a,b,c,d,e,f,sDD,xx,yy,small

    small = epsilon(1d0)*2

    ! Coefs
    a = x4 -x1
    b = x2 -x1
    c = x3-x4-x2 +x1
    d = y4 -y1
    e = y2 -y1
    f = y3-y4-y2 +y1

    ! Solve A*p**2 + B*p + C = 0
    yy = y-y1
    xx = x-x1
    AA = c*d - a*f
    BB = -c*yy + b*d + xx*f - a*e
    CC = -yy*b + e*xx
    if(abs(AA)<small)then
        p1 = -CC/BB
        p2 = p1
    else
        DD = BB**2 - 4*AA*CC
        sDD = sqrt(DD)
        p1 = (-BB-sDD)/(2*AA)
        p2 = (-BB+sDD)/(2*AA)
    endif

    ! Get q from p
    if(abs(b+c*p1)>small)then
        q1 = (xx-a*p1)/(b+c*p1)
    else
        q1 = (yy-d*p1)/(e+f*p1)
    endif

    ! Select point closest to center
    if(p1<0d0 .or. p1>1d0 .or. q1<0d0 .or. q1>1d0)then
        if(abs(b+c*p2)>small)then
            q2 = (xx-a*p2)/(b+c*p2)
        else
            q2 = (yy-d*p2)/(e+f*p2)
        endif
        p = p2
        q = q2
    else
        p = p1
        q = q1
    endif

end subroutine curv2rect

subroutine curv2rectss(xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,xx,yy,pp,qq,np)
    ! Same as curv2rel but for a definite and same number of quadrangles and points
    implicit none

    integer, intent(in) :: np
    real(kind=8), intent(in), dimension(np) :: xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,xx,yy
    real(kind=8), intent(out), dimension(np) :: pp,qq
    integer :: i
    !$OMP PARALLEL DO PRIVATE(i)
    !$& SHARED(xx1,xx2,xx3,xx4,yy1,yy2,yy3,yy4,xx,yy,pp,qq,np)
    do i=1,np
        call curv2rect(xx1(i),xx2(i),xx3(i),xx4(i), &
            & yy1(i),yy2(i),yy3(i),yy4(i),xx(i),yy(i),pp(i),qq(i))
    enddo
    !$OMP END PARALLEL DO
end subroutine curv2rectss

subroutine curv2rel(xxi, yyi, xo, yo, p, q, nxi, nyi, no)
    ! Convert a series of absolute coordinates to coordinates relative to
    ! a curved grid

    implicit none

    integer,intent(in) :: nxi,nyi,no
    real(kind=8),intent(in) :: xxi(nyi,nxi), yyi(nyi,nxi), xo(no), yo(no)
    real(kind=8),intent(out) :: p(no), q(no)

    integer :: io!,i,j,ic,jc
!    real(kind=8) :: a,b
!    logical :: binside

    p = -1d0
    q = -1d0

    !$OMP PARALLEL DO PRIVATE(io)
    !$& SHARED(xxi,yyi,xo,yo,nxi,nyi,no,p,q)
    do io = 1, no

        call curv2rel_single(xxi, yyi, xo(io), yo(io), p(io), q(io), nxi, nyi)

    enddo
    !$OMP END PARALLEL DO

end subroutine curv2rel

subroutine curv2rel_single(xxi, yyi, xo, yo, p, q, nxi, nyi)
    ! curv2rel for a single output point

    implicit none

    integer,intent(in) :: nxi,nyi
    real(kind=8),intent(in) :: xxi(nyi,nxi), yyi(nyi,nxi), xo, yo
    real(kind=8),intent(out) :: p, q

    integer :: i,j,ic,jc
    real(kind=8) :: a,b

    p = -1d0
    q = -1d0

    ! Find the closest corner
    call closest2d(xxi,yyi,xo,yo,nxi,nyi,ic,jc,1)

    ! Curvilinear to rectangular
    main: do i=max(ic-1,1), min(ic,nxi-1)
        do j = max(jc-1,1), min(jc,nyi-1)

            ! Get relative position
            call curv2rect(xxi(j,i),xxi(j+1,i),xxi(j+1,i+1),xxi(j,i+1), &
                         & yyi(j,i),yyi(j+1,i),yyi(j+1,i+1),yyi(j,i+1), &
                         & xo, yo, a, b)

            ! Store absolute indices
            if(a>=0d0-tiny(0d0) .and. a<=1d0+tiny(0d0) &
                & .and. b>=0d0-tiny(0d0) .and. b<=1d0+tiny(0d0))then
                p = dble(i) + a
                q = dble(j) + b
                return
            endif

        enddo
    enddo main

end subroutine curv2rel_single


subroutine nearest2dto1dc_reduc(p,q,zzi,zo,nxi,nyi,no,nz)
    ! Nearest interpolation of gridded data with 2D AXES to random positions
    ! This version takes relative positions with respect to output grid

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: p(no),q(no),zzi(nz,nyi,nxi)
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j
    real(kind=8) :: a,b,m1

    m1 = -1d0
    zo = sqrt(m1)

    !$OMP PARALLEL DO PRIVATE(io,i,j,a,b)
    !$& SHARED(p,q,zzi,zo,nxi,nyi,no)
    do io = 1, no

        if(p(io)>0d0 .and. q(io)>0d0)then

            ! Cell
            a = mod(p(io),1d0)
            b = mod(q(io),1d0)
            i = int(p(io)-a)
            j = int(q(io)-b)

            ! Interpolation
            if(a>0.5)i = i+1
            if(b>0.5)j = j+1
            zo(:, io) = zzi(:,j, i)
        endif

    enddo
    !$OMP END PARALLEL DO

end subroutine nearest2dto1dc_reduc

subroutine nearest2dto1dc(xxi,yyi,zzi,xo,yo,zo,nxi,nyi,no,nz)
    ! nearest interpolation of gridded data with 2D AXES to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xxi(nyi,nxi), yyi(nyi,nxi), xo(no), yo(no)
    real(kind=8),intent(in) :: zzi(nz,nyi,nxi)
    real(kind=8),intent(out) :: zo(nz,no)

    real(kind=8) :: p(no), q(no)

    ! Relative positions
    call curv2rel(xxi, yyi, xo, yo, p, q, nxi, nyi, no)

    ! Interpolation
    call nearest2dto1dc_reduc(p, q, zzi, zo, nxi, nyi, no, nz)

end subroutine nearest2dto1dc

subroutine linear2dto1dc_reduc(p,q,zzi,zo,nxi,nyi,no,nz)
    ! Bilinear interpolation of gridded data with 2D AXES to random positions
    ! This version takes relative positions with respect to output grid

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: p(no),q(no),zzi(nz,nyi,nxi)
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j
    real(kind=8) :: a,b,m1

    m1 = -1d0
    zo = sqrt(m1)

    !$OMP PARALLEL DO PRIVATE(io,i,j,a,b)
    !$& SHARED(p,q,zzi,zo,nxi,nyi,no)
    do io = 1, no

        if(p(io)>0d0 .and. q(io)>0d0)then

            ! Cell
            a = mod(p(io),1d0)
            b = mod(q(io),1d0)
            i = int(p(io)-a)
            j = int(q(io)-b)


            ! Interpolation
            zo(:, io) = (1-b)*(1-a)*zzi(:,j,  i) + &
            &       (1-b)*a*zzi(:,j,  i+1) + &
            &        b*(1-a)*zzi(:,j+1,i) + &
            &        b*a*zzi(:,j+1,i+1)

        endif

    enddo
    !$OMP END PARALLEL DO

end subroutine linear2dto1dc_reduc


subroutine linear2dto1dc(xxi,yyi,zzi,xo,yo,zo,nxi,nyi,no,nz)
    ! Bilinear interpolation of gridded data with 2D AXES to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xxi(nyi,nxi), yyi(nyi,nxi), xo(no), yo(no)
    real(kind=8),intent(in) :: zzi(nz,nyi,nxi)
    real(kind=8),intent(out) :: zo(nz,no)

    real(kind=8) :: p(no), q(no)

    ! Relative positions
    call curv2rel(xxi, yyi, xo, yo, p, q, nxi, nyi, no)

    ! Interpolation
    call linear2dto1dc_reduc(p, q, zzi, zo, nxi, nyi, no, nz)

end subroutine linear2dto1dc


subroutine dstwgt2dto1dc_reduc(p,q,zzi,zo,nxi,nyi,no,nz)
    ! Bilinear interpolation of gridded data with 2D AXES to random positions
    ! This version takes relative positions with respect to output grid

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: p(no),q(no),zzi(nz,nyi,nxi)
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j,i4
    real(kind=8) :: a,b, vv(nz,4),dx0,dx1,dy0,dy1,ww(nz,4),wsum(nz),dd(4),m1
    logical :: bmask(nz,nyi,nxi),bb(nz,4)

    m1 = -1d0
    zo = sqrt(m1)
    bmask = isnan(zzi)

    !$OMP PARALLEL DO PRIVATE(io,i,j,a,b,vv,bb,dd,ww,wsum,i4)
    !$& SHARED(p,q,zzi,zo,bmask,nxi,nyi,no)
    do io = 1, no

        if(p(io)>0d0 .and. q(io)>0d0)then

            ! Cell
            a = mod(p(io),1d0)
            b = mod(q(io),1d0)
            i = int(p(io)-a)
            j = int(q(io)-b)

            ! Distances
            dx0 = a
            dx1 = 1d0-a
            dy0 = b
            dy1 = 1d0-b
            dd = (/sqrt(dx0**2+dy0**2),sqrt(dx0**2+dy1**2),&
            &      sqrt(dx1**2+dy0**2),sqrt(dx1**2+dy1**2)/)

            ! Data
            vv = reshape(zzi(:,j:j+1,i:i+1),(/nz,4/))

            ! On a point?
            if(any(dd==0d0))then

                zo(:,io) = vv(:, minloc(dd,dim=1))

            else

                ! Weights
                bb = reshape(bmask(:,j:j+1,i:i+1),(/nz,4/))
                ww = 0d0
                do i4=1,4
                    ww(:,i4) = 1d0/dd(i4)
                enddo
                where(bb)
                    ww=0d0
                    vv=0d0
                endwhere
                wsum = merge(1d0, sum(ww,dim=2), all(bb,dim=2))

                ! Interpolation
                zo(:,io) = merge(sqrt(m1), sum(ww*vv, dim=2)/wsum, all(bb,dim=2))

            endif

        endif

    enddo
    !$OMP END PARALLEL DO

end subroutine dstwgt2dto1dc_reduc


subroutine dstwgt2dto1dc(xxi,yyi,zzi,xo,yo,zo,nxi,nyi,no,nz)
    ! Distance weight interpolation of gridded data with 2D AXES to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xxi(nyi,nxi), yyi(nyi,nxi), xo(no), yo(no)
    real(kind=8),intent(in) :: zzi(nz,nyi,nxi)
    real(kind=8),intent(out) :: zo(nz,no)

    real(kind=8) :: p(no), q(no)

    ! Relative positions
    call curv2rel(xxi, yyi, xo, yo, p, q, nxi, nyi, no)

    ! Interpolation
    call dstwgt2dto1dc_reduc(p, q, zzi, zo, nxi, nyi, no, nz)

end subroutine dstwgt2dto1dc


