! Copyright or © or Copr. Actimar (contributor(s) : Stephane Raynaud) (2010)
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

subroutine interp1d(vari, yi, varo, yo, mv, method, nx, nyi, nyo, extrap)
    ! Interpolation along the second axis (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yi: input y axis
    ! - yo: output y axis
    ! - mv: missing value (used only for initialisation)
    ! - method: 0 = nearest, 1 = linear, 2 = cubic, 3 = hermit
    ! - extrap: 0 = do not extrapolate, 1 = top, -1 = bottom, 2 = both
    !
    ! See: http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,method
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yi(nyi), yo(nyo)
    real(kind=8), intent(in) :: mv
    real(kind=8), intent(out) :: varo(nx,nyo)
    integer, intent(in), optional :: extrap

    ! Internal
    integer :: iyi,iyo
    real(kind=8) :: dy0,dy1,yi0,yi1,mu
    real(kind=8) :: tension,bias,a0,a1,a2,a3 ! Hermit
    real(kind=8) :: zyi(nyi),zyo(nyo)
    real(kind=8),allocatable :: vc0(:),vc1(:)
    logical :: bmask(nx,nyi)

    ! Treat only monotonically increasing grids
    zyi = yi
    zyo = yo
    if (yi(nyi)<yi(1))zyi = -yi
    if (nyo>1.and.yo(nyo)<yo(1))zyo = -yo

    ! Initialisation
    varo = mv
    if(method>1) allocate(vc0(nx),vc1(nx))
    bias = 0.
    tension = 0.
    bmask = abs(vari-mv)<=abs(epsilon(0d0)*1.1*mv)

    ! Loop on input grid
    iyo = 1
    do iyi = 1, nyi-1

        yi0 = zyi(iyi)
        yi1 = zyi(iyi+1)

        if (yi1<zyo(1)) cycle
        if (yi0>zyo(nyo)) exit

        ! Loop on output grid
        do while (iyo <= nyo )

             ! Out of interval
            if (zyo(iyo) < yi0) then
                iyo = iyo + 1
                cycle
            endif
            if (zyo(iyo) > yi1) exit

            ! Distances to neighbours
            dy0 = zyo(iyo)-yi0
            dy1 = yi1-zyo(iyo)

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
                &    (vari(:, iyi)*dy1 + vari(:, iyi+1)*dy0) / &
                &    (dy0+dy1)
                varo(:,iyo) = merge(mv, varo(:,iyo), any(bmask(:, iyi:iyi+1),dim=2) )

            else

                ! Cubic and Hermit
                !
                if (iyi==1)then ! y0
                    vc0 = 2*vari(:, iyi)-vari(:, iyi+1)
                else
                    vc0 = vari(:, iyi-1)
                endif
                if (iyi==nyi-1)then ! y3
                    vc1 = 2*vari(:, iyi+1)-vari(:, iyi)
                else
                    vc1 = vari(:, iyi+2)
                endif
                mu = dy0/(dy0+dy1)

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

                ! Mask
                varo(:,iyo) = merge(mv, varo(:,iyo), &
                    & any(bmask(:, max(iyi-1,1):min(iyi+2,nyi)), dim=2) )

            endif
            iyo = iyo + 1

        end do
    end do

    ! Extrapolation with nearest
    !if(method==0 .and.present(extrap).and.extrap/=0)then
    if(present(extrap).and.extrap/=0)then
        if((extrap==-1 .or. extrap==2) .and. (zyo(1)<zyi(1)))then
            do iyo=1,nyo
                if(zyi(1)<zyo(iyo))exit
                varo(:,iyo) = vari(:,1)
            enddo
        endif
        if((extrap==1 .or. extrap==2) .and. (zyo(nyo)>zyi(nyi)))then
            do iyo=nyo,1,-1
                if(zyi(nyi)>zyo(iyo))exit
                varo(:,iyo) = vari(:,nyi)
            enddo
        endif
    endif


end subroutine interp1d

subroutine interp1dx(vari, yi, varo, yo, mv, method, nx, nxb, nyi, nyo, extrap)
    ! Interpolation along the second axis (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yi: input y axis
    ! - yo: output y axis
    ! - mv: missing value (used only for initialisation)
    ! - method: 0 = nearest, 1 = linear, 2 = cubic, 3 = hermit
    ! - extrap: 0 = do not extrapolate, 1 = top, -1 = bottom, 2 = both
    !
    ! See: http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,method,nxb
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nyo)
    real(kind=8), intent(in) :: mv
    real(kind=8), intent(out) :: varo(nx,nyo)
    integer, intent(in), optional :: extrap

    ! Internal
    integer :: iyi,iyo,ib
    real(kind=8) :: dy0(nxb),dy1(nxb)
    real(kind=8) :: tension,bias !
    real(kind=8) :: zyi(nxb,nyi),zyo(nyo),mu(nxb)
    real(kind=8),allocatable :: vc0(:),vc1(:),a0(:),a1(:),a2(:),a3(:)
    integer :: ix0,ix1
    logical :: bitv(nxb), bmask(nx,nyi)

    ! Treat only monotonically increasing grids
    bmask = abs(vari-mv)<=(epsilon(0d0)*1.1*mv)
    do ib = 1,nxb
        if(all(bmask(ib,:)))cycle
        if (yi(ib,nyi)<yi(ib,1)) then
            zyi = -yi
        else
            zyi = yi
        endif
        exit
    enddo
    if (yo(nyo)<yo(1)) then
        zyo = -yo
    else
        zyo = yo
    endif

    ! Initialisation
    varo = mv
    if(method>1) allocate(vc0(nxb),vc1(nxb))
    if(method==3)allocate(a0(nxb),a1(nxb),a2(nxb),a3(nxb))
    bias = 0.
    tension = 0.

    ! Loop on input grid
    do iyi = 1, nyi-1

        ! Loop on output grid
        do iyo = 1,nyo

            dy0 = zyo(iyo)-zyi(:,iyi)
            dy1 = zyi(:,iyi+1)-zyo(iyo)
            bitv = zyo(iyo)>=zyi(:,iyi).and.zyo(iyo)<=zyi(:,iyi+1)
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
                        varo(ix0:ix1,iyo) = merge(mv, varo(ix0:ix1,iyo), &
                            & any(bmask(ix0:ix1, iyi:iyi+1),dim=2) )
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
                            varo(ix0:ix1,iyo) = merge(mv, varo(ix0:ix1,iyo), &
                                & any(bmask(ix0:ix1, max(iyi-1,1):min(iyi+2,nyi)), dim=2) )
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
                            varo(ix0:ix1,iyo) = merge(mv, varo(ix0:ix1,iyo), &
                                & any(bmask(ix0:ix1, max(iyi-1,1):min(iyi+2,nyi)), dim=2) )
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
                bitv = zyo(iyo)<zyi(:,1) ! below
                do ib = 1, nx/nxb
                    ix0 = 1+(ib-1)*nxb
                    ix1 = ix0+nxb-1
                    where(bitv) varo(ix0:ix1,iyo) = vari(ix0:ix1, 1)
                enddo
            endif
            if(extrap==1 .or. extrap==2)then
                bitv = zyo(iyo)>zyi(:,nyi) ! above
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

subroutine interp1dxx(vari, yi, varo, yo, mv, method, nx, nxb, nyi, nyo, extrap)
    ! Interpolation along the second axis (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yi: input y axis
    ! - yo: output y axis
    ! - mv: missing value (used only for initialisation)
    ! - method: 0 = nearest, 1 = linear, 2 = cubic, 3 = hermit
    ! - extrap: 0 = do not extrapolate, 1 = top, -1 = bottom, 2 = both
    !
    ! See: http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,method,nxb
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nxb,nyo)
    real(kind=8), intent(in) :: mv
    real(kind=8), intent(out) :: varo(nx,nyo)
    integer, intent(in), optional :: extrap

    ! Internal
    integer :: iyi,iyo,ib
    real(kind=8) :: dy0(nxb),dy1(nxb)
    real(kind=8) :: tension,bias !
    real(kind=8) :: zyi(nxb,nyi),zyo(nxb,nyo),mu(nxb)
    real(kind=8),allocatable :: vc0(:),vc1(:),a0(:),a1(:),a2(:),a3(:)
    integer :: ix0,ix1
    logical :: bitv(nxb), bmask(nx,nyi)

    ! Treat only monotonically increasing grids
    bmask = abs(vari-mv)<=abs(epsilon(0d0)*1.1*mv)
    do ib = 1,nxb
        if(all(bmask(ib,:)))cycle
        if (yi(ib,nyi)<yi(ib,1)) then
            zyi = -yi
        else
            zyi = yi
        endif
        if (yo(ib,nyo)<yo(ib,1)) then
            zyo = -yo
        else
            zyo = yo
        endif
        exit
    enddo
!    if (maxvalyo(nyo)<yo(1)) then
!        zyo = -yo
!    else
!        zyo = yo
!    endif

    ! Initialisation
    varo = mv
    if(method>1) allocate(vc0(nxb),vc1(nxb))
    if(method==3)allocate(a0(nxb),a1(nxb),a2(nxb),a3(nxb))
    bias = 0.
    tension = 0.

    ! Loop on input grid
    do iyi = 1, nyi-1

        ! Loop on output grid
        do iyo = 1,nyo

            dy0 = zyo(:,iyo)-zyi(:,iyi)
            dy1 = zyi(:,iyi+1)-zyo(:,iyo)
            bitv = zyo(:,iyo)>=zyi(:,iyi).and.zyo(:,iyo)<=zyi(:,iyi+1)
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
                        varo(ix0:ix1,iyo) = merge(mv, varo(ix0:ix1,iyo), &
                            & any(bmask(ix0:ix1, iyi:iyi+1),dim=2) )
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
                            varo(ix0:ix1,iyo) = merge(mv, varo(ix0:ix1,iyo), &
                                & any(bmask(ix0:ix1, max(iyi-1,1):min(iyi+2,nyi)), dim=2) )
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
                            varo(ix0:ix1,iyo) = merge(mv, varo(ix0:ix1,iyo), &
                                & any(bmask(ix0:ix1, max(iyi-1,1):min(iyi+2,nyi)), dim=2) )
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
                bitv = zyo(:,iyo)<zyi(:,1) ! below
                do ib = 1, nx/nxb
                    ix0 = 1+(ib-1)*nxb
                    ix1 = ix0+nxb-1
                    where(bitv) varo(ix0:ix1,iyo) = vari(ix0:ix1, 1)
                enddo
            endif
            if(extrap==1 .or. extrap==2)then
                bitv = zyo(:,iyo)>zyi(:,nyi) ! above
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


subroutine extrap1d(vari, varo, mv, extrap, nx, ny)
    ! Extrapolate valid data to the top and/or bottom
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - mv: missing value (used only for initialisation)
    ! - extrap: 0 = do not extrapolate, 1 = top, -1 = bottom, 2 = both
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: vari(nx,ny)
    real(kind=8), intent(in) :: mv
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

        valid = abs(vari(ix,:)-mv)>tiny(1d0)

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

subroutine remap1d(vari, yi, varo, yo, mv, conserv, nx, nyi, nyo, yib, yob,extrap)
    ! Remapping along the second axis (y)

    implicit none

    ! Extrernal
    integer, intent(in) :: nx, nyi, nyo
    integer, intent(in) :: conserv
    real(kind=8),intent(in) :: vari(nx,nyi),mv
    real(kind=8),intent(in) :: yi(nyi),yo(nyo)
    real(kind=8),intent(out) :: varo(nx,nyo)
    real(kind=8),intent(in), optional :: yib(nyi+1),yob(nyo+1)
    integer, intent(in), optional :: extrap

    ! Local
    integer :: iyi,iyo
    real(kind=8) :: zyib(nyi+1),zyob(nyo+1),wo(nx),dyio

    ! Bounds
    if(present(yib).and. .not.all(yib==0d0))then
        zyib = yib
    else
        zyib(2:nyi-1) = (yi(2:nyi)+yi(1:nyi-1))*.5
        zyib(1)     = yi(1)  +(yi(1)  -yi(2))    *.5
        zyib(nyi+1) = yi(nyi)+(yi(nyi)-yi(nyi-1))*.5
    endif
    if(present(yob).and. .not.all(yob==0d0))then
        zyob = yob
    else
        zyob(2:nyo-1) = (yo(2:nyo)+yo(1:nyo-1))*.5
        zyob(1)     = yo(1)  +(yo(1)  -yo(2))    *.5
        zyob(nyo+1) = yo(nyo)+(yo(nyo)-yo(nyo-1))*.5
    endif
    if (zyib(2)<zyib(1)) zyib = -zyib
    if (zyob(2)<zyob(1)) zyob = -zyob

    ! Extrapolation with nearest bound
    if(present(extrap).and.extrap/=0)then
       if(extrap==-1 .or. extrap==2)then
           if (zyib(1) > zyob(1)) zyib(1)=zyob(1) !below
       endif
       if(extrap==1 .or. extrap==2)then
           if (zyib(nyi+1) < zyob(nyo+1)) zyib(nyi+1)=zyob(nyo+1) !above
       endif
    endif
    ! Init
    varo = 0d0

    iyi = 1
    do iyo = 1, nyo

        wo = 0.
        do while (iyi<=nyi)

            ! No intersection
            if(zyib(iyi)>zyob(iyo+1)) exit
            if(zyib(iyi+1)<zyob(iyo))then
                iyi = iyi +1
                cycle
            endif

            ! Contribution of intersection
            dyio = min(zyib(iyi+1),zyob(iyo+1))-max(zyib(iyi),zyob(iyo))
            if(conserv==1.and.zyib(iyi)/=zyib(iyi+1)) dyio = dyio / (zyib(iyi+1)-zyib(iyi))
            where(vari(:,iyi)/=mv)
                wo = wo + dyio
                varo(:,iyo) = varo(:,iyo) + vari(:,iyi)*dyio
            endwhere

            ! Next input cell?
            if(zyib(iyi+1)>zyob(iyo+1)) exit
            iyi = iyi +1

        enddo

        ! Normalize
        if(conserv==1)where(wo/=0.)wo = 1.
        varo(:,iyo) = merge(varo(:,iyo)/wo, mv, wo/=0.)
!        where(wo/=0.)
!            varo(:,iyo) = varo(:,iyo)/wo
!        elsewhere
!            varo(:,iyo) = mv
!        endwhere

    enddo

end subroutine remap1d

subroutine remap1dx(vari, yi, varo, yo, mv, conserv, nx, nxb, nyi, nyo, yib, yob,extrap)
    ! Remapping along the second axis (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yi: input y axis
    ! - yo: output y axis
    ! - mv: missing value (used only for initialisation)
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    integer, intent(in) ::  conserv
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nyo)
    real(kind=8), intent(in) :: mv
    real(kind=8), intent(out) :: varo(nx,nyo)
    real(kind=8), intent(in), optional :: yib(nxb,nyi+1), yob(nyo+1)
    integer, intent(in), optional :: extrap

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: zyib(nxb,nyi+1),zyob(nyo+1),wo(nx),dyi(nxb)
    integer :: ix0,ix1
    logical :: mapi(nxb,4)

    ! Bounds
    if(present(yib).and. .not.all(yib==0d0))then
        zyib = yib
    else
        zyib(:,2:nyi-1) = (yi(:,2:nyi)+yi(:,1:nyi-1))*.5
        zyib(:,1)     = yi(:,1)  +(yi(:,1)  -yi(:,2))    *.5
        zyib(:,nyi+1) = yi(:,nyi)+(yi(:,nyi)-yi(:,nyi-1))*.5
    endif
    if(present(yob).and. .not.all(yob==0d0))then
        zyob = yob
    else
        zyob(2:nyo-1) = (yo(2:nyo)+yo(1:nyo-1))*.5
        zyob(1)     = yo(1)  +(yo(1)  -yo(2))    *.5
        zyob(nyo+1) = yo(nyo)+(yo(nyo)-yo(nyo-1))*.5
    endif
    if (zyib(1,nyi)<zyib(1,1)) zyib = -zyib
    if (zyob(nyo)<zyob(1)) zyob = -zyob

    ! Extrapolation with nearest bound
    if(present(extrap).and.extrap/=0)then
       if(extrap==-1 .or. extrap==2)then
           where (zyib(:,1) > zyob(1)) zyib(:,1)=zyob(1)
       endif
       if(extrap==1 .or. extrap==2)then
           where (zyib(:,nyi+1) < zyob(nyo+1)) zyib(:,nyi+1)=zyob(nyo+1)
       endif
    endif

    ! Initialisation
    varo = 0.
    if(conserv==0)dyi = 1.

    ! Loop on output grid
    do iyo = 1,nyo

        wo = 0.

        ! Loop on input grid
        do iyi = 1, nyi

            ! No intersection
            if(all(zyib(:, iyi)>=zyob(iyo+1)))exit
            if(all(zyib(:, iyi+1)<zyob(iyo)))cycle

            ! Conditional arrays
            mapi(:,1) = zyib(:, iyi)>=zyob(iyo).and.zyib(:, iyi+1)<=zyob(iyo+1)
            mapi(:,2) = zyib(:, iyi)< zyob(iyo).and.zyib(:, iyi+1)> zyob(iyo+1)
            mapi(:,3) = zyib(:, iyi)< zyob(iyo).and.zyib(:, iyi+1)>zyob(iyo)
            mapi(:,4) = zyib(:, iyi)< zyob(iyo+1).and.zyib(:, iyi+1)> zyob(iyo+1)

            ! Conservative
            if(conserv==1) dyi = zyib(:, iyi+1)-zyib(:, iyi)

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(vari(ix0:ix1,iyi)/=mv)
                    where(mapi(:,1))

                        ! Input inside
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(zyib(:, iyi+1)-zyib(:, iyi))/dyi
                        wo(ix0:ix1) = wo(ix0:ix1) + zyib(:, iyi+1)-zyib(:, iyi)

                    elsewhere(mapi(:,2))

                        ! Output inside
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(zyob(iyo+1)-zyob(iyo))/dyi
                        wo(ix0:ix1) = wo(ix0:ix1) + (zyob(iyo+1)-zyob(iyo))

                    elsewhere(mapi(:,3))

                        ! Input partly below
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(zyib(:, iyi+1)-zyob(iyo))/dyi
                        wo(ix0:ix1) = wo(ix0:ix1) + zyib(:, iyi+1)-zyob(iyo)

                    elsewhere(mapi(:,4))

                        ! Input partly above
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(zyob(iyo+1)-zyib(:, iyi))/dyi
                        wo(ix0:ix1) = wo(ix0:ix1) + (zyob(iyo+1)-zyib(:, iyi))

                    endwhere
                endwhere
            enddo
        enddo

        ! Normalize
        if(conserv==1)where(wo/=0.)wo = 1.
        where(wo/=0.)
            varo(:,iyo) =  varo(:,iyo)/wo
        elsewhere
            varo(:,iyo) = mv
        endwhere
    enddo

end subroutine remap1dx

subroutine remap1dxx(vari, yi, varo, yo, mv, conserv, nx, nxb, nyi, nyo, yib, yob,extrap)
    ! Remapping between two variable axes in space (y)
    !
    ! - vari: input variable
    ! - varo: output variable
    ! - yi: input y axis
    ! - yo: output y axis
    ! - mv: missing value (used only for initialisation)
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    integer,intent(in) ::  conserv
    real(kind=8), intent(in) :: vari(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nxb,nyo)
    real(kind=8), intent(in) :: mv
    real(kind=8), intent(out) :: varo(nx,nyo)
    real(kind=8), intent(in), optional :: yib(nxb,nyi+1), yob(nxb,nyo+1)
    integer, intent(in), optional :: extrap

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: zyib(nxb,nyi+1),zyob(nxb,nyo+1),wo(nx),dyi(nxb)
    integer :: ix0,ix1
    logical :: mapi(nxb,4)

    ! Bounds
    if(present(yib).and. .not.all(yib==0d0))then
        zyib = yib
    else
        zyib(:,2:nyi-1) = (yi(:,2:nyi)+yi(:,1:nyi-1))*.5
        zyib(:,1)     = yi(:,1)  +(yi(:,1)  -yi(:,2))    *.5
        zyib(:,nyi+1) = yi(:,nyi)+(yi(:,nyi)-yi(:,nyi-1))*.5
    endif
    if(present(yob).and. .not.all(yob==0d0))then
        zyob = yob
    else
        zyob(:,2:nyo-1) = (yo(:,2:nyo)+yo(:,1:nyo-1))*.5
        zyob(:,1)     = yo(:,1)  +(yo(:,1)  -yo(:,2))    *.5
        zyob(:,nyo+1) = yo(:,nyo)+(yo(:,nyo)-yo(:,nyo-1))*.5
    endif
    if (zyib(1,nyi)<zyib(1,1)) zyib = -zyib
    if (zyob(1,nyo)<zyob(1,1)) zyob = -zyob

    ! Extrapolation with nearest bound
    if(present(extrap).and.extrap/=0)then
       if(extrap==-1 .or. extrap==2)then
           where (zyib(:,1) > zyob(:,1)) zyib(:,1)=zyob(:,1)
       endif
       if(extrap==1 .or. extrap==2)then
           where (zyib(:,nyi+1) < zyob(:,nyo+1)) zyib(:,nyi+1)=zyob(:,nyo+1)
       endif
    endif

    ! Initialisation
    varo = 0.
    if(conserv==0)dyi = 1.

    ! Loop on output grid
    do iyo = 1,nyo

        wo = 0.

        ! Loop on input grid
        do iyi = 1, nyi

            ! No intersection
            if(all(zyib(:, iyi)>=zyob(:, iyo+1)))exit
            if(all(zyib(:, iyi+1)<zyob(:, iyo)))cycle

            ! Conditional arrays
            mapi(:,1) = zyib(:, iyi)>=zyob(:, iyo).and.zyib(:, iyi+1)<=zyob(:, iyo+1) ! Inside
            mapi(:,2) = zyib(:, iyi)< zyob(:, iyo).and.zyib(:, iyi+1)> zyob(:, iyo+1) ! Embed
            mapi(:,3) = zyib(:, iyi)< zyob(:, iyo).and.zyib(:, iyi+1)>zyob(:, iyo) ! Below
            mapi(:,4) = zyib(:, iyi)<zyob(:, iyo+1).and.zyib(:, iyi+1)> zyob(:, iyo+1) ! Above

            ! Conservative
            if(conserv==1) dyi = zyib(:, iyi+1)-zyib(:, iyi)

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(vari(ix0:ix1,iyi)/=mv)
                    where(mapi(:,1))

                        ! Input inside
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(zyib(:, iyi+1)-zyib(:, iyi))/dyi
                        wo(ix0:ix1) = wo(ix0:ix1) + zyib(:, iyi+1)-zyib(:, iyi)

                    elsewhere(mapi(:,2))

                        ! Output inside
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(zyob(:, iyo+1)-zyob(:, iyo))/dyi
                        wo(ix0:ix1) = wo(ix0:ix1) + zyob(:, iyo+1)-zyob(:, iyo)

                    elsewhere(mapi(:,3))

                        ! Input partly below
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(zyib(:, iyi+1)-zyob(:, iyo))/dyi
                        wo(ix0:ix1) = wo(ix0:ix1) + zyib(:, iyi+1)-zyob(:, iyo)

                    elsewhere(mapi(:,4))

                        ! Input partly above
                        varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + vari(ix0:ix1,iyi)&
                        &    *(zyob(:, iyo+1)-zyib(:, iyi))/dyi
                        wo(ix0:ix1) = wo(ix0:ix1) + zyob(:, iyo+1)-zyib(:, iyi)

                    endwhere
                endwhere
            enddo
        enddo

        ! Normalize
        if(conserv==1)where(wo>0d0)wo = 1.
        where(wo>0d0)
            varo(:,iyo) =  varo(:,iyo)/wo
        elsewhere
            varo(:,iyo) = mv
        endwhere
    enddo

end subroutine remap1dxx


subroutine cellerr1d(vari, yi, varo, yo, mv, errm, errl, erro, yob, nx, nxb, nyi, nyo)
    ! Error-based weighthed local mooving average interpolation between
    ! two variables axes in space (y) with interpolation error estimate
    !
    ! :Params:
    !
    !   - vari: input variable
    !   - varo: output variable
    !   - yi: input y axis
    !   - yo: output y axis
    !   - mv: missing value (used only for initialisation)
    !   - errm: input measurement error
    !   - errl: input lag error per yi units
    !   - erro: output error
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    real(kind=8), intent(in) :: vari(nx,nyi), errm(nx,nyi)
    real(kind=8), intent(in) :: yi(nyi), yo(nyo), errl(nxb)
    real(kind=8), intent(in) :: mv
    real(kind=8), intent(out) :: varo(nx,nyo), erro(nx,nyo)
    real(kind=8), intent(in), optional :: yob(nyo+1)

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: zyob(nyo+1),zyi(nyi), zyo(nyo), zerrl(nxb)
    integer :: ix0,ix1
    logical :: goodi(nx,nyi), berrl(nxb)

    ! Masks
    goodi = abs(vari-mv)>abs(epsilon(1d0)*mv)
    berrl = abs(errl-mv)<=abs(epsilon(1d0)*mv)

   ! Coordinates
    if (yi(nyi)<yi(1)) then
        zyi = -yi
    else
        zyi = yi
    endif
     if (yo(nyo)<yo(1)) then
        zyo = -yo
    else
        zyo = yo
    endif
    if(present(yob).and. .not.all(yob==0d0))then
        zyob = yob
    else
        if(nyo==1)then
            zyob(1) = min(yo(1), minval(yi))
            zyob(2) = max(yo(2), maxval(yi))
        else
            zyob(2:nyo-1) = (yo(2:nyo)+yo(1:nyo-1))*.5
            zyob(1)     = yo(1)  +(yo(1)  -yo(2))    *.5
            zyob(nyo+1) = yo(nyo)+(yo(nyo)-yo(nyo-1))*.5
        endif
    endif
    if (zyob(nyo)<zyob(1)) zyob = -zyob


    ! Initialisation
    varo = 0d0
    erro = 0d0

    ! Loop on output grid
    iyi = 1
    do iyo = 1,nyo

        ! Loop on input grid
        do while (iyi<=nyi)

            ! Not useful
            if(zyi(iyi)>zyob(iyo+1)) exit
            if(zyi(iyi)<zyob(iyo).or..not.any(goodi(:,iyi)))then
                iyi = iyi +1
                cycle
            endif

            ! Quadratic lag error
            where(berrl)
                zerrl = 0d0
            elsewhere(errl<0)
                zerrl = errl**2 * abs(zyi(iyi)-zyo(iyo))
            elsewhere
                zerrl = errl * abs(zyi(iyi)-zyo(iyo))
            endwhere

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(goodi(ix0:ix1,iyi))

                    varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + &
                        & vari(ix0:ix1,iyi) / (errm(:, iyi)**2 + zerrl)
                    erro(ix0:ix1,iyo) = erro(ix0:ix1,iyo) + &
                        & 1d0 / (errm(:, iyi)**2 + zerrl)

                endwhere
            enddo

            iyi = iyi + 1

        enddo

        ! Normalize
        where(erro(:,iyo)>0d0)
            varo(:,iyo) =  varo(:,iyo)/erro(:,iyo)
            erro(:,iyo) = 1d0/sqrt(erro(:,iyo))
        elsewhere
            varo(:,iyo) = mv
            erro(:,iyo) = mv
        endwhere

    enddo


end subroutine cellerr1d



subroutine cellerr1dx(vari, yi, varo, yo, mv, errm, errl, erro, yob, nx, nxb, nyi, nyo)
    ! Error-based weighthed local mooving average interpolation between
    ! two variables axes in space (y) with interpolation error estimate
    !
    ! :Params:
    !
    !   - vari: input variable
    !   - varo: output variable
    !   - yi: input y axis
    !   - yo: output y axis
    !   - mv: missing value (used only for initialisation)
    !   - errm: input measurement error
    !   - errl: input lag error per yi units
    !   - erro: output error
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    real(kind=8), intent(in) :: vari(nx,nyi), errm(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nyo), errl(nxb)
    real(kind=8), intent(in) :: mv
    real(kind=8), intent(out) :: varo(nx,nyo), erro(nx,nyo)
    real(kind=8), intent(in), optional :: yob(nyo+1)

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: zyob(nyo+1),zyi(nxb,nyi), zyo(nyo), zerrl(nxb)
    integer :: ix0,ix1
    logical :: goodi(nxb), berrl(nxb), bvari(nx,nyi)

    ! Masks
    bvari = abs(vari-mv)<=abs(epsilon(1d0)*mv)
    berrl = abs(errl-mv)<=abs(epsilon(1d0)*mv)

   ! Coordinates
    do ib = 1,nxb
        if(all(bvari(ib,:)))cycle
        if (yi(ib,nyi)<yi(ib,1)) then
            zyi = -yi
        else
            zyi = yi
        endif
       exit
    enddo
     if (yo(nyo)<yo(1)) then
        zyo = -yo
    else
        zyo = yo
    endif
    if(present(yob).and. .not.all(yob==0d0))then
        zyob = yob
    else
        if(nyo==1)then
            zyob(1) = min(yo(1), minval(yi))
            zyob(2) = max(yo(2), maxval(yi))
        else
            zyob(2:nyo-1) = (yo(2:nyo)+yo(1:nyo-1))*.5
            zyob(1)     = yo(1)  +(yo(1)  -yo(2))    *.5
            zyob(nyo+1) = yo(nyo)+(yo(nyo)-yo(nyo-1))*.5
        endif
    endif
    if (zyob(nyo)<zyob(1)) zyob = -zyob


    ! Initialisation
    varo = 0d0
    erro = 0d0

    ! Loop on output grid
    do iyo = 1,nyo

        ! Loop on input grid
        do iyi = 1, nyi

            ! Valid?
            goodi(:) = .not. bvari(:,iyi) .and. &
                & zyi(:, iyi)>=zyob(iyo) .and. zyi(:, iyi)<=zyob(iyo+1)
            if(.not.any(goodi))cycle

            ! Quadratic lag error
            where(berrl)
                zerrl = 0d0
            elsewhere(errl<0)
                zerrl = errl**2 * abs(zyi(:, iyi)-zyo(iyo))
            elsewhere
                zerrl = errl * abs(zyi(:, iyi)-zyo(iyo))
            endwhere

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(goodi(ix0:ix1))

                    varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + &
                        & vari(ix0:ix1,iyi) / (errm(:, iyi)**2 + zerrl)
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
            varo(:,iyo) = mv
            erro(:,iyo) = mv
        endwhere
    enddo


end subroutine cellerr1dx



subroutine cellerr1dxx(vari, yi, varo, yo, mv, errm, errl, erro, yob, nx, nxb, nyi, nyo)
    ! Error-based weighthed local mooving average interpolation between
    ! two variables axes in space (y) with interpolation error estimate
    !
    ! :Params:
    !
    !   - vari: input variable
    !   - varo: output variable
    !   - yi: input y axis
    !   - yo: output y axis
    !   - mv: missing value (used only for initialisation)
    !   - errm: input measurement error
    !   - errl: input lag error per yi units
    !   - erro: output error
    !

    implicit none

    ! Extrernal
    integer, intent(in) :: nx,nyi,nyo,nxb
    real(kind=8), intent(in) :: vari(nx,nyi), errm(nx,nyi)
    real(kind=8), intent(in) :: yi(nxb,nyi), yo(nxb,nyo), errl(nxb)
    real(kind=8), intent(in) :: mv
    real(kind=8), intent(out) :: varo(nx,nyo), erro(nx,nyo)
    real(kind=8), intent(in), optional :: yob(nxb,nyo+1)

    ! Local
    integer :: iyi,iyo,ib
    real(kind=8) :: zyob(nxb,nyo+1),zyi(nxb,nyi), zyo(nxb,nyo), zerrl(nxb)
    integer :: ix0,ix1
    logical :: goodi(nxb), berrl(nxb), bvari(nx,nyi)

    ! Masks
    bvari = abs(vari-mv)<=abs(epsilon(1d0)*mv)
    berrl = abs(errl-mv)<=abs(epsilon(1d0)*mv)

   ! Coordinates
    do ib = 1,nxb
        if(all(bvari(ib,:)))cycle ! FIXME: bad bmask
        if (yi(ib,nyi)<yi(ib,1)) then
            zyi = -yi
        else
            zyi = yi
        endif
         if (yo(ib,nyo)<yo(ib,1)) then
            zyo = -yo
        else
            zyo = yo
        endif
       exit
    enddo
    if(present(yob).and. .not.all(yob==0d0))then
        zyob = yob
    else
        if(nyo==1)then
            zyob(:,1) = merge(yo(:,1), yi(:,1), yo(:,1)<yi(:,1))
            zyob(:,2) = merge(yo(:,1), yi(:,nyi), yo(:,1)>yi(:,nyi))
        else
            zyob(:,2:nyo-1) = (yo(:,2:nyo)+yo(:,1:nyo-1))*.5
            zyob(:,1)     = yo(:,1)  +(yo(:,1)  -yo(:,2))    *.5
            zyob(:,nyo+1) = yo(:,nyo)+(yo(:,nyo)-yo(:,nyo-1))*.5
        endif
    endif
    if (zyob(1,nyo)<zyob(1,1)) zyob = -zyob


    ! Initialisation
    varo = 0d0
    erro = 0d0

    ! Loop on output grid
    do iyo = 1,nyo

        ! Loop on input grid
        do iyi = 1, nyi

            ! Valid?
            goodi(:) = .not. bvari(:,iyi) .and. &
                & zyi(:, iyi)>=zyob(:, iyo) .and. zyi(:, iyi)<=zyob(:, iyo+1)
            if(.not.any(goodi))cycle

            ! Quadratic lag error
            where(berrl)
                zerrl = 0d0
            elsewhere(errl<0)
                zerrl = errl**2 * abs(zyi(:, iyi)-zyo(:, iyo))
            elsewhere
                zerrl = errl * abs(zyi(:, iyi)-zyo(:, iyo))
            endwhere

            ! Loop on blocks
            do ib = 1, nx/nxb

                ! Block
                ix0 = 1+(ib-1)*nxb
                ix1 = ix0+nxb-1
                where(goodi(ix0:ix1))

                    varo(ix0:ix1,iyo) = varo(ix0:ix1,iyo) + &
                        & vari(ix0:ix1,iyi) / (errm(:, iyi)**2 + zerrl)
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
            varo(:,iyo) = mv
            erro(:,iyo) = mv
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
!    real(kind=8) :: georectif(nyi,nxi)
!    real(kind=8) :: zxxi(nyi,nxi), zxxo(nyo,nxo)
!     real(kind=8),allocatable :: zxib(:,:) ,zyib(:,:)

!    interface
!        function closest2d(xxi,yyi,xo,yo,georectif,geo)
!        real(kind=8),intent(in) :: xxi(:,:),yyi(:,:),xo,yo,georectif(:,:)
!        real(kind=8) :: closest2d(2)
!        logical,intent(in) :: geo
!        end function closest2d
!    end interface

!    ! Latitude rectification
!    if(geo)then
!        georectif = cos(yyi*3.14159d0/180.)
!        zxxi = xxi
!        zxxo = xxo
!!         zxxi = modulo(xxi,360.)
!!         zxxo = modulo(xxo,360.)
!!         where(zxi<zxi(1))zxi=zxi+360.
!!         where(zxo<zxo(1))zxo=zxo+360.
!    else
!        georectif = 1.
!        zxxi = xxi
!        zxxo = xxo
!    endif

    geo = .not. present(nogeo) .or. nogeo==0

    ! Loop on output points
    if(nb==0)then

        ! Scan all input points everytime
        do iyo = 1, nyo
            do ixo = 1, nxo
                call closest2d(xxi,yyi,xxo(iyo,ixo),yyo(iyo,ixo),nxi,nyi,imin,jmin,.not. geo)
                varo(:,iyo,ixo) = vari(:,jmin,imin)
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
        do ixo = 1, nxo
            do iyo = 1, nyo

                ! Try a small block
                ixmin = max(1,ixlast-znb2)
                ixmax = min(nxi,ixlast+znb2)
                iymin = max(1,iylast-znb2)
                iymax = min(nyi,iylast+znb2)
                call closest2d(xxi(iymin:iymax,ixmin:ixmax), &
                    & yyi(iymin:iymax,ixmin:ixmax),xxo(iyo,ixo),yyo(iyo,ixo),nxi,nyi,imin,jmin,.not. geo)
                imin = imin+ixmin-1
                jmin = jmin+iymin-1

                ! Fall on bounds so use full block
                if((imin==ixmin.and.ixmin/=1).or.(imin==ixmax.and.ixmax/=nxi).or.&
                    & (jmin==iymin.and.iymin/=1).or.(jmin==iymax.and.iymax/=nyi))&
                    & call closest2d(xxi,yyi,xxo(iyo,ixo),yyo(iyo,ixo),nxi,nyi,imin,jmin,.not. geo)

                ! Store value
                varo(:,iyo,ixo) = vari(:,jmin,imin)

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

subroutine closest2d(xxi,yyi,xo,yo,nxi,nyi,i,j,nogeo)
    ! Find indices of closest point on 2D axes
    implicit none
    real(kind=8),intent(in) :: xxi(nyi,nxi),yyi(nyi,nxi),xo,yo
    integer, intent(in):: nyi,nxi
    integer,intent(out) :: i,j
    logical,intent(in) :: nogeo
    real(kind=8) :: dx(nyi,nxi)
    integer:: ij(2)

    if(nogeo)then
        dx = xxi-xo
    else
        dx = abs(xxi-xo)
        where(dx>180d0)dx = 360d0-dx
        dx = dx * cos(yyi*3.14159d0/180d0)
    endif
    ij = minloc(dx**2+(yyi-yo)**2)
    i = ij(2)
    j = ij(1)

end subroutine closest2d

! =============================================================================

subroutine bilin    (vari, xi, yi, varo, xo,  yo, mv, nogeo, nxi, nyi, nxo, nyo, nz)
    ! Simple bilinear interpolation between two regular grids
    !
    ! See also :f:subr:`mixt2d` as an alternative method

    implicit none

    ! Parameters
    integer,intent(in) :: nxi,nyi,nxo,nyo,nz
    real(kind=8),intent(in) :: vari(nz,nyi,nxi),mv
    real(kind=8),intent(in) :: xi(nxi),yi(nyi),xo(nxo),yo(nyo)
    integer,intent(in), optional :: nogeo
    real(kind=8),intent(out) :: varo(nz,nyo,nxo)

    ! Local variables
    real(kind=8) :: zxi(nxi),zyi(nyi),zxo(nxo),zyo(nyo)
    real(kind=8) :: dxi,dyi,fx,fy,fdxi
    integer :: ixi,iyi,ixo,iyo
    logical :: geo, bmaski(nz,nyi,nxi)

    ! Monotonically increasing
    zxi = xi
    zyi = yi
    zxo = xo
    zyo = yo
    if (xi(nxi)<xi(1))zxi = -xi
    if (yi(nyi)<yi(1))zyi = -yi
    if (xo(nxo)<xo(1))zxo = -xo
    if (yo(nyo)<yo(1))zyo = -yo

    ! Geographic
!     if(geo)then
!         zxi = modulo(zxi,360.)
!         zxo = modulo(zxo,360.)
!         where(zxi<zxi(1))zxi=zxi+360.
!         where(zxo<zxo(1))zxo=zxo+360.
!     endif
    geo = .not. present(nogeo) .or. nogeo==0

    ! Missing
    varo = mv
    bmaski = abs(vari-mv)<abs(epsilon(1d0)*mv)

    ! Loop on input x
    iyo = 1
    do iyi = 1, nyi-1

        ! Must overlap
        if (zyi(iyi+1)<zyo(1))cycle
        if (zyi(iyi)>zyo(nyo))exit

        ! Cell height
        dyi = zyi(iyi+1)-zyi(iyi)

        ! Loop on output y
        do while (iyo<=nyo.and.zyo(iyo)<=zyi(iyi+1))

            ! Still not inside
            if(zyo(iyo)<zyi(iyi))then
                iyo = iyo+1
                cycle
            endif

            ! Y weight
            fy = (zyo(iyo)-zyi(iyi))/dyi

            ! Loop on input x
            ixo = 1
            do ixi = 1, nxi-1

                ! Must overlap
                if (zxi(ixi+1)<zxo(1))cycle
                if (zxi(ixi)>zxo(nxo))exit

                ! Cell width
                dxi = zxi(ixi+1)-zxi(ixi)
                if(geo.and.dxi>180.)dxi=360.-dxi

                ! Loop on output x
                do while (ixo<=nxo.and.zxo(ixo)<=zxi(ixi+1))

                    ! Still not inside
                    if(zxo(ixo)<zxi(ixi))then
                        ixo = ixo+1
                        cycle
                    endif

                    ! X weight
                    fdxi = zxo(ixo)-zxi(ixi)
                    if(geo.and.fdxi>180.)fdxi=360.-fdxi
                    fx = fdxi/dxi

                    ! Interpolation
                    varo(:,iyo,ixo) = &
                        & vari(:,  iyi  ,ixi)*(1-fx)*(1-fy) + &
                        & vari(:,iyi  ,ixi+1)*fx    *(1-fy) + &
                        & vari(:,iyi+1,ixi)*(1-fx)*fy     + &
                        & vari(:,iyi+1,ixi+1)*fx    *fy

                    ! Masking
                    varo(:,iyo,ixo) = merge(mv, varo(:,iyo,ixo), &
                        & bmaski(:,iyi,ixi).or.bmaski(:,iyi,ixi+1).or. &
                        & bmaski(:,iyi+1,ixi).or.bmaski(:,iyi+1,ixi+1))

                    ! Next output x
                    ixo = ixo+1
                enddo
            enddo

            ! Next output y
            iyo = iyo+1
        enddo
    enddo
end subroutine bilin

subroutine dstwgt   (vari, xi, yi, varo, xo,  yo, mv, nogeo, nxi, nyi, nxo, nyo, nz)
    ! Simple distance weight interpolation between two regular grids
    ! It does not interpolate missing values

    implicit none

    ! Parameters
    integer,intent(in) :: nxi,nyi,nxo,nyo,nz
    real(kind=8),intent(in) :: vari(nz,nyi,nxi),mv
    real(kind=8),intent(in) :: xi(nxi),yi(nyi),xo(nxo),yo(nyo)
    integer,intent(in), optional :: nogeo
    real(kind=8),intent(out) :: varo(nz,nyo,nxo)

    ! Local variables
    real(kind=8) :: zxi(nxi),zyi(nyi),zxo(nxo),zyo(nyo),small
    real(kind=8) :: dxi,dyi,dd(4),ww(nz,4),dx0,dx1,dy0,dy1,wsum(nz),vv(nz,4),mvs(nz)
    integer :: ixi,iyi,ixo,iyo,i4
    logical :: geo

    geo  = .not.present(nogeo) .or. nogeo==0

    ! Monotonically increasing
    zxi = xi
    zyi = yi
    zxo = xo
    zyo = yo
    if (xi(nxi)<xi(1))zxi = -xi
    if (yi(nyi)<yi(1))zyi = -yi
    if (xo(nxo)<xo(1))zxo = -xo
    if (yo(nyo)<yo(1))zyo = -yo

    ! Geographic
!     if(geo)then
!         zxi = modulo(zxi,360.)
!         zxo = modulo(zxo,360.)
!         where(zxi<zxi(1))zxi=zxi+360.
!         where(zxo<zxo(1))zxo=zxo+360.
!     endif

    varo = mv
    vv = mv
    small = epsilon(1d0)*1.1

    ! Loop on input x
    iyo = 1
    do iyi = 1, nyi-1

        ! Must overlap
        if (zyi(iyi+1)<zyo(1))cycle
        if (zyi(iyi)>zyo(nyo))exit

        ! Cell height
        dyi = zyi(iyi+1)-zyi(iyi)

        ! Loop on output y
         do while (iyo<=nyo.and.zyo(iyo)<zyi(iyi+1))

            ! Still not inside
            if(zyo(iyo)<zyi(iyi))then
                iyo = iyo+1
                cycle
            endif

            ! Y axis
            dy0 = zyo(iyo)-zyi(iyi)
            dy1 = zyi(iyi+1)-zyo(iyo)

            ! Loop on input x
            ixo = 1
            do ixi = 1, nxi-1

                ! Must overlap
                if (zxi(ixi+1)<zxo(1))cycle
                if (zxi(ixi)>zxo(nxo))exit

                ! Cell width
                dxi = zxi(ixi+1)-zxi(ixi)
                if(geo.and.dxi>180.)dxi=360.-dxi

                ! Loop on output x
                do while (ixo<=nxo.and.zxo(ixo)<zxi(ixi+1))

                    ! Still not inside
                    if(zxo(ixo)<zxi(ixi))then
                        ixo = ixo+1
                        cycle
                    endif

                    ! X axis
                    dx0 = zxo(ixo)-zxi(ixi)
                    dx1 = zxi(ixi+1)-zxo(ixo)
                    if(geo)then
                        if(dx0>180.)dx0=360.-dx0
                        if(dx1>180.)dx1=360.-dx1
                        dx0 = dx0*cos(zyo(iyo)*3.14159d0/180.)
                        dx1 = dx1*cos(zyo(iyo)*3.14159d0/180.)
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
                    where(abs((vv-mv))<mv*small)ww=0.
                    wsum = sum(ww,dim=2)!,mask=vv/=mv)
                    mvs = 0.
                    where(wsum==0.)
                        wsum = 1.
                        mvs = mv
                    endwhere

                    ! Interpolation
                     varo(:,iyo,ixo) = sum(ww*vv, dim=2)/wsum+mvs

                    ! Next output x
                    ixo = ixo+1
                enddo
            enddo

            ! Next output y
            iyo = iyo+1
        enddo
    enddo
end



subroutine mbilin2d (vari, xi,  yi,  varo, xo,  yo, mv, ext,  nxi, nyi, no,  nogeo)
    ! Binilear interpolation of random point to a regular grid
    ! with minimal missing data handling

    implicit none

    ! Parameters
    integer, intent(in) :: nxi, nyi, no
    real(kind=8),    intent(in) :: xi(nxi), yi(nyi), vari(nyi, nxi), xo(no), yo(no)
    real(kind=8),    intent(in) :: mv
    logical, intent(in) :: ext
    integer, intent(in), optional :: nogeo
    real(kind=8),    intent(out) :: varo(no)

    ! Local parameters
    integer :: ix, iy, io, ix1, iy1
    real(kind=8) :: ww(2,2), zz(2,2), dx, dy, fx, fy, pi, zxi(nxi),zxo(no)
    logical ::  ms(2,2), ma(2,2),geo


    ! Inits
    pi = 3.14159d0
    varo = mv

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
        ms = zz==mv
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

end subroutine mbilin2d

! =============================================================================

subroutine nearest2dto1d(xi,yi,zi,xo,yo,zo,mv,nxi,nyi,no,nz)
    ! nearest neighbour interpolation of gridded data to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xi(nxi), yi(nyi), xo(no), yo(no)
    real(kind=8),intent(in) :: zi(nz,nyi,nxi),mv
    real(kind=8),intent(out) :: zo(nz,no)
    real(kind=8) :: dx(nxi)

    integer :: io,i,j

    zo = mv
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

subroutine bilin2dto1d(xi,yi,zi,xo,yo,zo,mv,nxi,nyi,no,nz)
    ! bilinear interpolation of gridded data to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xi(nxi), yi(nyi), xo(no), yo(no)
    real(kind=8),intent(in) :: zi(nz,nyi,nxi),mv
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j
    real(kind=8) :: a,b
    logical :: bmask(nz,nyi,nxi)

    zo = mv
    bmask = abs(zi-mv)<abs(mv*epsilon(0d0)*1.1)
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

            where(bmask(:,j,i).or.bmask(:,j,i+1).or.bmask(:,j+1,i).or.bmask(:,j+1,i+1))&
                & zo(:, io)=mv
        endif

    enddo
    !$OMP END PARALLEL DO
end subroutine bilin2dto1d

subroutine dstwgt2dto1d(xi,yi,zi,xo,yo,zo,mv,nxi,nyi,no,nz)
    ! Distance weight interpolation of gridded data to random positions
    !
    ! Distances are computed with the four corners of a cell and are relative
    ! to the cell sizes.

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xi(nxi), yi(nyi), xo(no), yo(no)
    real(kind=8),intent(in) :: zi(nz,nyi,nxi),mv
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j,i4
    real(kind=8) :: dx0,dx1,dy0,dy1,dd(4),vv(nz,4),wsum(nz),ww(nz,4)
    logical :: bmask(nz,nyi,nxi),bb(nz,4)

    zo = mv
    bmask = abs(zi-mv)<abs(mv*epsilon(1d0)*1.1)
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
            where(bb)ww=0d0
            wsum = merge(1d0, sum(ww,dim=2), all(bb,dim=2))

            ! Interpolation
            zo(:,io) = merge(mv, sum(ww*vv, dim=2)/wsum, all(bb,dim=2))

        endif

    enddo
    !$OMP END PARALLEL DO

end subroutine dstwgt2dto1d



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

!function dstpt2line(x,y,x1,x2,y1,y2)
!    ! Distance from a point to a line
!    real(kind=8),intent(in) :: x,y,x1,x2,y1,y2
!    real(kind=8) :: xy(2), dstpt2dt
!    xy = linept(x,y,x1,x2,y1,y2)
!    dstpt2line = sqrt((xy(1)-x)**2+(xy(2)-y)**2)
!end function dstpt2line

subroutine curv2rect(x1,x2,x3,x4,y1,y2,y3,y4,x,y,p,q)
    ! Coordinate transform from curvilinear to rectangular
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

    integer :: io,i,j,ic,jc
    real(kind=8) :: a,b
    logical :: binside

    p = -1d0
    q = -1d0

    !$OMP PARALLEL DO PRIVATE(io,i,j,ic,jc,a,b,binside)
    !$& SHARED(xxi,yyi,xo,yo,nxi,nyi,no,p,q)
    do io = 1, no

        ! Find the closest corner
        call closest2d(xxi,yyi,xo(io),yo(io),nxi,nyi,ic,jc,.true.)

        ! Curvilinear to rectangular
        binside = .false.
        main: do i=max(ic-1,1), min(ic,nxi-1)
            do j = max(jc-1,1), min(jc,nyi-1)

                ! Get relative position
                call curv2rect(xxi(j,i),xxi(j+1,i),xxi(j+1,i+1),xxi(j,i+1), &
                             & yyi(j,i),yyi(j+1,i),yyi(j+1,i+1),yyi(j,i+1), &
                             & xo(io), yo(io), a, b)

                ! Store absolute indices
                binside = a>=0d0-tiny(0d0) .and. a<=1d0+tiny(0d0) &
                    & .and. b>=0d0-tiny(0d0) .and. b<=1d0+tiny(0d0)
                if(binside)then
                    p(io) = dble(i) + a
                    q(io) = dble(j) + b
                    exit main
                endif

            enddo
!            if(binside)exit
        enddo main

    enddo
    !$OMP END PARALLEL DO

end subroutine curv2rel

subroutine nearest2dto1dc_reduc(p,q,zzi,zo,mv,nxi,nyi,no,nz)
    ! Nearest interpolation of gridded data with 2D AXES to random positions
    ! This version takes relative positions with respect to output grid

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: p(no),q(no),zzi(nz,nyi,nxi),mv
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j
    real(kind=8) :: a,b

    zo = mv

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

subroutine nearest2dto1dc(xxi,yyi,zzi,xo,yo,zo,mv,nxi,nyi,no,nz)
    ! nearest interpolation of gridded data with 2D AXES to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xxi(nyi,nxi), yyi(nyi,nxi), xo(no), yo(no)
    real(kind=8),intent(in) :: zzi(nz,nyi,nxi),mv
    real(kind=8),intent(out) :: zo(nz,no)

    real(kind=8) :: p(no), q(no)

    ! Relative positions
    call curv2rel(xxi, yyi, xo, yo, p, q, nxi, nyi, no)

    ! Interpolation
    call nearest2dto1dc_reduc(p, q, zzi, zo, mv, nxi, nyi, no, nz)

end subroutine nearest2dto1dc

subroutine bilin2dto1dc_reduc(p,q,zzi,zo,mv,nxi,nyi,no,nz)
    ! Bilinear interpolation of gridded data with 2D AXES to random positions
    ! This version takes relative positions with respect to output grid

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: p(no),q(no),zzi(nz,nyi,nxi),mv
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j
    real(kind=8) :: a,b
    logical :: bmask(nz,nyi,nxi)

    zo = mv
    bmask = abs(zzi-mv)<abs(epsilon(0d0)*1.1*mv)

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

            ! Mask
            zo(:, io) = merge(mv, zo(:, io), &
                & any(reshape(bmask(:, j:j+1, i:i+1),(/nz,4/)), dim=2))


        endif

    enddo
    !$OMP END PARALLEL DO

end subroutine bilin2dto1dc_reduc


subroutine bilin2dto1dc(xxi,yyi,zzi,xo,yo,zo,mv,nxi,nyi,no,nz)
    ! Bilinear interpolation of gridded data with 2D AXES to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xxi(nyi,nxi), yyi(nyi,nxi), xo(no), yo(no)
    real(kind=8),intent(in) :: zzi(nz,nyi,nxi),mv
    real(kind=8),intent(out) :: zo(nz,no)

    real(kind=8) :: p(no), q(no)

    ! Relative positions
    call curv2rel(xxi, yyi, xo, yo, p, q, nxi, nyi, no)

    ! Interpolation
    call bilin2dto1dc_reduc(p, q, zzi, zo, mv, nxi, nyi, no, nz)

end subroutine bilin2dto1dc


subroutine dstwgt2dto1dc_reduc(p,q,zzi,zo,mv,nxi,nyi,no,nz)
    ! Bilinear interpolation of gridded data with 2D AXES to random positions
    ! This version takes relative positions with respect to output grid

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: p(no),q(no),zzi(nz,nyi,nxi),mv
    real(kind=8),intent(out) :: zo(nz,no)

    integer :: io,i,j,i4
    real(kind=8) :: a,b, vv(nz,4),dx0,dx1,dy0,dy1,ww(nz,4),wsum(nz),dd(4)
    logical :: bmask(nz,nyi,nxi),bb(nz,4)

    zo = mv
    bmask = abs(zzi-mv)<abs(epsilon(1d0)*1.1*mv)

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
                where(bb)ww=0d0
                wsum = merge(1d0, sum(ww,dim=2), all(bb,dim=2))

                ! Interpolation
                zo(:,io) = merge(mv, sum(ww*vv, dim=2)/wsum, all(bb,dim=2))

            endif

        endif

    enddo
    !$OMP END PARALLEL DO

end subroutine dstwgt2dto1dc_reduc


subroutine dstwgt2dto1dc(xxi,yyi,zzi,xo,yo,zo,mv,nxi,nyi,no,nz)
    ! Distance weight interpolation of gridded data with 2D AXES to random positions

    implicit none

    integer,intent(in) :: nxi,nyi,no,nz
    real(kind=8),intent(in) :: xxi(nyi,nxi), yyi(nyi,nxi), xo(no), yo(no)
    real(kind=8),intent(in) :: zzi(nz,nyi,nxi),mv
    real(kind=8),intent(out) :: zo(nz,no)

    real(kind=8) :: p(no), q(no)

    ! Relative positions
    call curv2rel(xxi, yyi, xo, yo, p, q, nxi, nyi, no)

    ! Interpolation
    call dstwgt2dto1dc_reduc(p, q, zzi, zo, mv, nxi, nyi, no, nz)

end subroutine dstwgt2dto1dc


! =============================================================================

subroutine mixt2dx   (vari, xi,  yi,  varo, xo,  yo,  mv, ext,       nxi,nyi,nxo,nyo,nz)
    ! Extension of avgext2d to third dimension (simple loop)

    implicit none

    ! Parameters
    integer, intent(in) :: nxi,nyi,nxo,nyo,nz
    real(kind=8),intent(in) :: vari(nxi,nyi,nz),mv
    real(kind=8),intent(in) :: xi(nxi),yi(nyi),xo(nxo),yo(nyo)
    logical,intent(in) :: ext
    real(kind=8),intent(out) :: varo(nxo,nyo,nz)

    ! Local variables
    integer :: iz,ixy
    real(kind=8) :: xymin,xymax

    do iz = 1, nz
        ! Inter- and extra-polation
        call mixt2d(nxi,nyi,xi,yi,vari(:,:,iz),nxo,nyo,xo,yo,varo(:,:,iz))
        ! No extrapolation
        if(.not.ext.and.iz==nz)then
            xymin = minval(xi)
            xymax = maxval(xi)
            do ixy = 1, nxo
                if(xo(ixy)<xymin.or.xo(ixy)>xymax)varo(ixy,:,:) = mv
            enddo
            xymin = minval(yi)
            xymax = maxval(yi)
            do ixy = 1, nyo
                if(yo(ixy)<xymin.or.yo(ixy)>xymax)varo(:,ixy,:) = mv
            enddo
        endif
    enddo

end subroutine mixt2dx

SUBROUTINE mixt2d(IDIM,JDIM,XX,YY,FF,MDIM,NDIM,RR,SS,GG)

    implicit none
    integer :: idim, jdim, mdim, ndim
    real(kind=8), intent(in) :: XX(IDIM),YY(JDIM),RR(MDIM),SS(NDIM),FF(IDIM,JDIM)
    real(kind=8), intent(out) :: GG(MDIM,NDIM)
!  ROUTINE TO INTERPOLATE FROM THE 2-DIMENSIONAL ARRAY FF DEFINED ON THE
!  GRID XX X YY ONTO THE NEW GRID RR X SS FOR ARRAY GG.
!
!  :ASSUMPTIONS:
!               1) LINEAR INTERPOLATION AND EXTRAPOLATION ARE USED TO
!                  EXPAND AND EXTENDTFF OVER THE INPUT GRID.
!               2) GG IS OBTAINED BY AREA AVERAGING FF OVER THE CELL CE
!                  ON (M,N).
!               3) NO ACCOUNT IS TAKEN OF SPHERICITY - I.E. THE INPUT A
!                  OUTPUT GRIDS ARE ASSUMED TO BE RECTILINEAR.
!               4) THE OPERATION IS NOT INVERTABLE, NOR WILL IT YIELD
!                  THE IDENTITY OPERATION IF THE INPUT AND OUTPUT GRIDS
!                  ARE IDENTICAL.  THIS SUITS FINE TO COARSE INTERPOLAT
!                  BUT MAY NOT ALWAYS SUIT COARSE TO FINE EXTRAPOLATION
!
!
!  :Method: SCAN THE OBJECT GRID FOR N=1,N
!          LOOP THROUGH INPUT GRID FOR ALL STRIPS ENCOMPASSING N
!          SCAN THE OBJECT GRID FOR M=1,M FOR EACH
!          LOOP THROUGH INPUT GRID FOR ALL X-CELLS AROUND M
!          DO AREA AVERAGE.
!
!
    real(kind=8) :: s1, s2, dy, s3, s4, dsy, ds, x1, x2, r1, r2, r3, r4, &
        & dr, drx, FACR, facs, y1, y2, dx
    integer :: n, m, j, jj, i, ii
      IF(XX(IDIM).LT.RR(1).OR.XX(1).GT.RR(MDIM).OR. &
     &     YY(JDIM).LT.SS(1).OR.YY(1).GT.SS(NDIM)) THEN
         PRINT *,'GRIDS GIVEN TO INTERP_2 DO NOT OVER-LAP  '
           STOP 'FORCE EXIT'
      ENDIF

      DO 1 N=1,NDIM
      DO 1 M=1,MDIM
    1 GG(M,N)=0.0

!  INCREMENT YY AXIS TILL Y2 OVERLAPS S1 AT LEAST

      S1=1.5*SS(1)-0.5*SS(2)
      DO 88 JJ=1,JDIM
          J=JJ
         IF(YY(J+1).GT.S1) GOTO 89
   88 CONTINUE
   89 DY=YY(J+1)-YY(J)
      Y1=YY(J)
      IF(J.EQ.1) Y1=Y1-1.E10*DY
      Y2=YY(J+1)

      DO 8 N=1,NDIM
      IF(N.EQ.1) THEN
          S1=1.5*SS(N)-0.5*SS(N+1)
          S2=0.5*(SS(N)+SS(N+1))
      GOTO 2
      ENDIF
      IF(N.EQ.NDIM) THEN
          S1=0.5*(SS(N)+SS(N-1))
          S2=1.5*SS(N)-0.5*SS(N-1)
          GOTO 2
      ENDIF
      S1=0.5*(SS(N)+SS(N-1))
      S2=0.5*(SS(N)+SS(N+1))
    2 CONTINUE

! J LOOP: THE N CELL SHOULD BE INTERSECTING THE J CELL. AT EDGES OR
!         WHEN THE J CELL STRADDLES N, DO INTEGRATION TO S2.

      S3=MAX(Y1,S1)
      S4=MIN(Y2,S2)
      DS=S4-S3
      DSY=0.5*(S4+S3)-YY(J)

!  MOVE ALONG XX AXIS UNTIL X2 AT LEAST OVERLAPS R1

      R1=1.5*RR(1)-0.5*RR(2)
      DO 66 II=1,IDIM
      I=II
         IF(XX(I+1).GT.R1) GOTO 67
   66 CONTINUE
   67 DX=XX(I+1)-XX(I)
      X1=XX(I)
      IF(I.EQ.1) X1=X1-1.E10*DX
      X2=XX(I+1)
!
      DO 6 M=1,MDIM
      IF(M.EQ.1) THEN
      R1=1.5*RR(M)-0.5*RR(M+1)
      R2=0.5*(RR(M)+RR(M+1))
      GOTO 4
      ENDIF
      IF(M.EQ.MDIM) THEN
      R1=0.5*(RR(M)+RR(M-1))
      R2=1.5*RR(M)-0.5*RR(M-1)
      GOTO 4
       ENDIF
      R1=0.5*(RR(M)+RR(M-1))
      R2=0.5*(RR(M)+RR(M+1))
    4 CONTINUE

! I LOOP


      R3=MAX(X1,R1)
      R4=MIN(X2,R2)
      DR=R4-R3
      DRX=0.5*(R4+R3)-XX(I)

! DO AVERAGE OVER THIS (PART) CELL

      FACR=DRX/DX
      FACS=DSY/DY
      GG(M,N)=GG(M,N) + DR*DS*(FF(I,J)*(1.-FACR)*(1.-FACS) + &
     &                           FF(I+1,J)*FACR*(1.-FACS) + &
     &                           FF(I,J+1)*FACS*(1.-FACR) + &
     &                           FF(I+1,J+1)*FACR*FACS )

!  REPEAT FOR FURTHER CELL SECTIONS NEEDED FOR THIS CELL, M,N

      IF(X2.GE.R2) GOTO 6
         I=I+1
         DX=XX(I+1)-XX(I)
         X1=XX(I)
         X2=XX(I+1)
         IF(I+1.EQ.IDIM) X2=X2+1.E10*DX
         GOTO 4
    6 CONTINUE
!  REPEAT FOR FUTHER STRIPS NEEDED FOR THIS STRIP N

      IF(Y2.GE.S2) GOTO 8
         J=J+1
         DY=YY(J+1)-YY(J)
         Y1=YY(J)
         Y2=YY(J+1)
         IF(J+1.EQ.JDIM) Y2=Y2+1.E10*DY
         GOTO 2
    8 CONTINUE

!  DIVIDE EACH SUM BY THE AREA OF THE CELL

      DO 14 N=1,NDIM
       IF(N.EQ.1)THEN
         DS=SS(N+1)-SS(N)
        ELSE
         IF(N.EQ.NDIM) THEN
        DS=SS(N)-SS(N-1)
        ELSE
         DS=0.5*(SS(N+1)-SS(N-1))
       ENDIF
       ENDIF
         DO 12 M=1,MDIM
       IF(M.EQ.1)THEN
         DR=RR(M+1)-RR(M)
        ELSE
         IF(M.EQ.MDIM) THEN
        DR=RR(M)-RR(M-1)
        ELSE
         DR=0.5*(RR(M+1)-RR(M-1))
       ENDIF
       ENDIF
            GG(M,N)=GG(M,N)/DR/DS
   12   CONTINUE
   14 CONTINUE


      RETURN
END subroutine mixt2d





! =============================================================================

subroutine cargen(xi,yi,zi,xo,yo,zo,mv,npt,nx,ny)
    ! The famous cargen interpolator


  implicit none
  ! declaration lecture hxhy
  integer,intent(in):: nx,ny ,npt
  REAL*8,intent(out):: zo(nx,ny)
  REAL*8,intent(in):: xi(npt),yi(npt),zi(npt)
  REAL*8,intent(in):: xo(nx),yo(ny),mv
  integer imin,imax,jmin,jmax
  integer npint
  real*8 xk(npt),yk(npt),zk(npt)
  real*8 x0,y0,z0,pourcent
  real*4 ytemp,taillytemp,fitemp

  !    pour fichier head
  REAL*8 finac,fisac,geac,gwac,dfiac
  REAL*8 finac2,fisac2,geac2,gwac2
  REAL*8 dgac,creuxdon
  real*8 del2,creux,deltafi,deltag,epsfi,epsg
  integer nbpts,nbx,nby,ix,jy
  real*8 fibas,fihaut,fibasp,fihautp
  real*8 gouest,gest,gouestp,gestp
  real*8 fial,gal,fi,g,taillx,epsig,epsifi,tailly
  real*8 minig,maxig,minifi,maxifi

  integer i,j,l
  integer k,nk
  real*8 rad
  data rad/57.29578/

  !    autre
!   integer modeinterp
  real*8 lag,cova,factdiv
!   character*3 testnbcarreaux

  !==========================================================================
  !==========================================================================
  !                    A REMPLIR EN DUR PAR L'UTILISATEUR CONFIRME
  !==========================================================================
  !==========================================================================
  !===== Parametre d'interpolation
  !==========================================================================
  !parametres pour le nombre de gros carreaux (par d�faut mettre creuxdon=0.1
  !et factdiv=2.5), modifi� principalement factdiv.
  !si testnbcarreaux='oui' le programme s'arrete apr�s avoir donn� le quadrillage
  creuxdon=0.1  ! 0<creuxdon<0.9
  factdiv=2.5
!   testnbcarreaux='non'
  !Methode d'interpolation
!   modeinterp=1
  !1. Mod�le lin�aire ifremer (� prendre par defaut)
  !2. Mod�le exponentiel (� develloper pour parametrer lag et cova)
  ! tester avec cova=1 et trouver la gamme de valeur pour la quelle lag
  ! n'a plus d'influence sur hmax et hymax
  ! Puis caler cova par rapport � la bathy max de la zone
  lag=10.0   !distance entre 2 points bathy
  cova=0.25  !covariance

  !=========================================================================
  !=========================================================================
  !================                           ==============================
  !================         PROGRAMME         ==============================
  !================                           ==============================
  !=========================================================================
  !=========================================================================
  !----------------------------------------------------------------------
  ! Lecture fichier head
  !----------------------------------------------------------------------
  imin=1
  imax = nx
  jmin=1
  jmax = ny
  dgac = sum(xo(2:nx)-xo(1:nx-1))/float(nx-1)
  dfiac = sum(yo(2:ny)-yo(1:ny-1))/float(ny-1)
  gwac = xo(2)
  geac = xo(nx)
  fisac = yo(2)
  finac = yo(ny)


  !----------------------------------------------------------------------
  ! Lecture bathy
  !----------------------------------------------------------------------
  nbpts = npt

  zo = mv

  !-------------------------------------------------------
  !         INTERPOLATION
  !-------------------------------------------------------

  del2=(dgac/rad)**2
  creux=creuxdon
  nbx=sqrt(nbpts*0.001)/(1.0-min(creux,0.9))/factdiv
  nbx=max(1,nbx)
  nby=nbx
  deltafi=(finac-fisac)/nby
  deltag=(gwac-geac)/nbx
  !-----modif ben-----
  epsfi=0.5*deltafi
  epsg=0.5*deltag
  fisac2=fisac-abs(epsfi)
  finac2=finac+abs(epsfi)
  gwac2=gwac-abs(epsg)
  geac2=geac+abs(epsg)
  nbx=nbx+1
  nby=nby+1
  !-------------------
  deltafi=(finac2-fisac2)/nby
  deltag=(gwac2-geac2)/nbx
  epsfi=0.2*deltafi
  epsg=epsfi/cos(0.5*(finac2+fisac2)/rad)

  !
  !----------------------------------------------------------
  !  on decoupe en gros carreaux
  !----------------------------------------------------------
  !

  minig=99.0
  maxig=-99.0
  minifi=99.0
  maxifi=-99.0


  do jy=1,nby

     pourcent=jy/nby*100.0
     !    write(*,*)pourcent
     fibas =fisac2+(jy-1)*deltafi
     fihaut=fibas+deltafi
     fibasp=fibas-epsfi
     fihautp=fihaut+epsfi
     do ix=1,nbx
        gouest=gwac2-(ix-1)*deltag
        gest=gouest-deltag
        !-----modif ben-----
        !       gouestp=gouest+epsg
        !       gestp=gest-epsg
        gouestp=gouest-epsg
        gestp=gest+epsg
        !-------------------
        npint=0
!         write(*,*) '-----------------------------'
!         write(*,*) 'carreau :',jy,ix
        !    write(*,*) fibas,fihaut,fibasp,fihautp
        !    write(*,*) gouest,gest,gouestp,gestp
        !
        !----------------------------------------------------------
        ! on recherche dans le gros carreau les points
        ! parmi les sondes qui sont dans le gros carreaux.
        !----------------------------------------------------------
        !

        do l=1,nbpts
           fial=yi(l)
           gal=xi(l)
           if(gal.ge.gouestp.and.gal.le.gestp&
                &.and.fial.ge.fibasp.and.fial.le.fihautp)then
              npint=npint+1
              xk(npint)=-gal/rad
              fitemp=dble(fial)
              ytemp=alog(tan(fitemp*0.0087266463+0.7853981))
              yk(npint)=dble(ytemp)
              zk(npint)=zi(l)
           endif
        end do
        !    if (npint.eq.0) stop
        !
        taillx=dgac/rad
        epsig=dgac/100.0
        epsifi=dfiac/100.0
        tailly=taillx
        !
        !----------------------------------------------------------
        !  on balaie la matrice de calcul et on regarde les mailles
        !  qui sont dans ces gros carreaux
        !----------------------------------------------------------
        !

        do j=jmin,jmax
          fi=yo(j)
          do i=imin,imax

                 g=xo(i)

                 fitemp=dble(fi)
                 ytemp=alog(tan(fitemp*0.0087266463+0.7853981))
                 y0=dble(ytemp)
                 fitemp=dble(fi+dfiac)
                 taillytemp=alog(tan(fitemp*0.0087266463+0.7853981))-ytemp
                 tailly=dble(taillytemp)

                 if (minig.ge.(gouest-epsig)) then
                    minig=(gouest-epsig)
                 endif
                 if (maxig.le.(gest+epsig)) then
                    maxig=(gest+epsig)
                 endif
                 if (minifi.ge.(fibas-epsifi)) then
                    minifi=(fibas-epsifi)
                 endif
                 if (maxifi.le.(fihaut+epsifi)) then
                    maxifi=(fihaut+epsifi)
                 endif


                 if(g.ge.(gouest-epsig).and.&
                      &       g.le.(gest+epsig).and.&
                      &      fi.ge.(fibas-epsifi).and.&
                      &      fi.le.(fihaut+epsifi)) then
                    x0=-g/rad
                    !
                    !--------------------------------------------------------------------
                    ! s il y a au moins trois sondes dans ce gros carreau on va pouvoir
                    ! faire un krigeage
                    !--------------------------------------------------------------------
                    !
                    if(npint.gt.3)then
!                        if(modeinterp.eq.1) then
                          call zkrig(x0,y0,z0,xk,yk,zk,npint,npt)
!                        elseif(modeinterp.eq.2) then
!                           call zkriben(x0,y0,z0,xk,yk,zk,npint,lag,cova,npt)
!                        endif

                       !            if(z0.lt.hmin) z0=-77.7
                       !
                       !--------------------------------------------------------------------
                       ! ce -99.9 est-il la terre ?
                       !--------------------------------------------------------------------
                       !
                    elseif(npint.eq.0)then
                       z0=mv
                    elseif(npint.eq.1)then
                       z0=mv
                    else
                       z0=0.5*(zk(1)+zk(2))
                    endif


                   zo(i,j)=z0

                 endif


           enddo !i
        enddo !j

     enddo !jy
  enddo !jx

end subroutine cargen

subroutine zkrig(x0,y0,z0,x,y,z,n,npt)

  implicit none

  integer, intent(in) :: npt,n
  integer :: ordre(20),m!,nbsondmx
  real*8 :: yvois(20),zvois(20),dist(20),xvois(20)
  real*8,intent(in) :: x(npt),y(npt),z(npt)
  real*8,intent(in) :: x0,y0
  real*8,intent(out) :: z0
  real*8 :: a2,b2,a1,b1,denom,d0i,bid
  real*8 :: xp(3),yp(3),zp(3),xv,yv
  integer :: i,ifois,ia,ic,ib,it,id,j,iv,l
  interface
    function plan(xpl,ypl,zpl,x1,y1)
      real*8,intent(in):: xpl(3),ypl(3),zpl(3),x1,y1
      real*8 :: plan
    end function plan
  end interface

   m=min0(12,n)
  do j=1,m
     ordre(j)=j
     dist(j)=1.e+20
  end do
  do i=1,n
     d0i=(x0-x(i))*(x0-x(i))+(y0-y(i))*(y0-y(i))
     do j=1,m
        if(d0i-dist(j).lt.0.0) then
           do l=m,j+1,-1
              ordre(l)=ordre(l-1)
              dist(l)=dist(l-1)
           end do
           ordre(j)=i
           dist(j)=d0i
           go to 10
        endif
     end do
10   continue
  end do
  !
  do j=1,m
     i=ordre(j)
     zvois(j)=z(i)
     xvois(j)=x(i)
     yvois(j)=y(i)
  enddo

  ifois=0
500 ifois=ifois+1
  ia=0
  ib=0
  ic=0
  id=0
  denom=xvois(1)-x0
  if(denom.eq.0.0)  denom=denom+0.000000001
  a1=(yvois(1)-y0)/denom
  b1=y0-a1*x0
  if(a1.eq.0.0) a1=a1+1.0e-20
  a2=-1.0/a1
  b2=y0-a2*x0
  xp(1)=xvois(1)
  yp(1)=yvois(1)
  zp(1)=zvois(1)
  it=1
  !
  do iv=2,m
     xv=xvois(iv)
     yv=yvois(iv)
     if(yvois(1)-a2*xvois(1)-b2.gt.0.0) then
        if(yv-a2*xv-b2.le.0.0) then
           !
           if(yv-a1*xv-b1.gt.0.0) then
              ib=ib+1
              if(ib.lt.2) then
                 it=it+1
                 xp(it)=xv
                 yp(it)=yv
                 zp(it)=zvois(iv)
                 if(it.ge.3) then
                    z0=plan(xp,yp,zp,x0,y0)
                    go to 800
                 end if
              end if
           else
              ia=ia+1
              if(ia.lt.2) then
                 it=it+1
                 xp(it)=xv
                 yp(it)=yv
                 zp(it)=zvois(iv)
                 if(it.ge.3) then
                    z0=plan(xp,yp,zp,x0,y0)
                    go to 800
                 end if
              end if
           end if
           !
        end if
     end if
  end do
  !
  if(ifois.lt.2) then
     bid=xvois(1)
     xvois(1)=xvois(2)
     xvois(2)=bid
     bid=yvois(1)
     yvois(1)=yvois(2)
     yvois(2)=bid
     bid=zvois(1)
     zvois(1)=zvois(2)
     zvois(2)=bid
     go to 500
  elseif(ifois.eq.2) then
     bid=xvois(1)
     xvois(1)=xvois(3)
     xvois(3)=bid
     bid=yvois(1)
     yvois(1)=yvois(3)
     yvois(3)=bid
     bid=zvois(1)
     zvois(1)=zvois(3)
     zvois(3)=bid
     go to 500
  endif
  !
  z0=z(ordre(1))



  !
800 continue


  !
end subroutine zkrig

function plan(xpl,ypl,zpl,x1,y1)
  ! Planar interpolation with 3 points
  implicit none

  real*8,intent(in):: xpl(3),ypl(3),zpl(3),x1,y1
  real*8 :: plan
  real*8 p,q,r,s,t,u,a,d,denom,b

  !
  !    write(*,*) (xpl(j),j=1,3)
  !    write(*,*) (ypl(j),j=1,3)
  !    write(*,*) (zpl(j),j=1,3)
  p=zpl(1)-zpl(2)
  q=xpl(1)-xpl(2)
  r=ypl(1)-ypl(2)
  s=zpl(2)-zpl(3)
  t=xpl(2)-xpl(3)
  u=ypl(2)-ypl(3)
  denom=r*t-u*q
  !
  if(denom.ne.0.0.and.(q.ne.0.0.or.t.ne.0.0)) then
     b=(s*q-p*t)/denom
     if(q.ne.0.0) then
        a=-(p+b*r)/q
     elseif(t.ne.0.0) then
        a=-(s+b*u)/t
     endif
     d=-zpl(1)-a*xpl(1)-b*ypl(1)
     plan=-a*x1-b*y1-d
  else
     write(*,*)'interpolation impossible en x0,y0, avec (x,y)...'
     write(*,7)x1,y1,xpl(1),ypl(1),xpl(2),ypl(2),xpl(3),ypl(3)
     plan=zpl(1)
  endif
  !
7 format(2f8.3,2x,6f10.4)
end function plan


