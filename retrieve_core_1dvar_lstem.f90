

!!   The subroutine returns the retrieved emissivity  
!! Jiheng Hu, 2024/1/21


subroutine retrieve_core_1dvar(nlevel,nchannel,imonth,incident,longitude,latitude,&
		PLEVEL,PHALF,TBobs,Tatm,QWatm,LST,T2m,Snowc,Smc,SfcPress,&
		EmissAnalyt,Emiss1st,Em1DVar,TbTune,LstTune,nIters)


integer,					intent(IN)	:: nlevel
integer,					intent(IN)	:: nchannel
integer,					intent(IN)	:: imonth
real,						intent(IN)	:: incident
real,  						intent(IN)	:: longitude
real,  						intent(IN)	:: latitude

real, 	dimension(nlevel),  intent(IN)	:: PLEVEL
real, 	dimension(nlevel),  intent(IN)	:: PHALF
real, 	dimension(nchannel),intent(IN) 	:: TBobs
real, 	dimension(nlevel), 	intent(IN)	:: Tatm,QWatm
real, 						intent(IN)	:: LST,T2m,Snowc,Smc,SfcPress

real, 	dimension(nchannel),intent(OUT)	:: EmissAnalyt,Em1DVar,Emiss1st,TbTune
real,						intent(OUT)	:: LstTune
INTEGER,					intent(OUT) :: nIters


INTEGER, parameter :: nG=10, nR=10
real, 	dimension(nchannel)				:: Emss_K,LST_K
real, 	dimension(nchannel,nlevel)		:: Ta_K,Qw_K
real, 	dimension(nchannel,nG)			:: Ja

!!! covariance
integer,	parameter					:: nEOF=8
real, 		dimension(nchannel)         :: EY

!! TB and EMissivity
real, dimension(nchannel,nchannel)  	:: LambdaTB  !! NEDT*NEDT

real, dimension(nG)      				:: Xb, Xg, DX
REAL, dimension(nR,nR)					:: Lambda

REAL*8, dimension(nchannel,nchannel) :: Q,invM
REAL, dimension(nR,nchannel) :: G
REAL,    DIMENSION(nchannel)      :: B

integer :: file_unit

INTEGER ::  nLayEff
INTEGER ::  iter,ich,i,j
LOGICAL ::  CvgceReached
REAL, parameter ::  ChiSqThresh=1.0


!! departures 
REAL, DIMENSION(nchannel) :: dY2

EY=(/0.77,0.78,0.63,0.60,0.51,0.41,0.42,0.32,0.31/)  

! 0.25*0.25
Lambda(1,:) =(/12.773 ,0.,0.,0.,0.,0.,0.,0.,0.,0./)
Lambda(2,:) =(/0., 0.0625,0.,0.,0.,0.,0.,0.,0.,0./)
Lambda(3,:) =(/0., 0.,0.0625,0.,0.,0.,0.,0.,0.,0./)
Lambda(4,:) =(/0., 0.,0.,0.0625,0.,0.,0.,0.,0.,0./)
Lambda(5,:) =(/0., 0.,0.,0.,0.0625,0.,0.,0.,0.,0./)
Lambda(6,:) =(/0., 0.,0.,0.,0.,0.0625,0.,0.,0.,0./)
Lambda(7,:) =(/0., 0.,0.,0.,0.,0.,0.0625,0.,0.,0./)
Lambda(8,:) =(/0., 0.,0.,0.,0.,0.,0.,0.0625,0.,0./)
Lambda(9,:) =(/0., 0.,0.,0.,0.,0.,0.,0.,0.0625,0./)
Lambda(10,:)=(/0., 0.,0.,0.,0.,0.,0.,0.,0.,0.0625/)

! NEDT*NEDT
LambdaTB(1,:)=(/0.77,0.,0.,0.,0.,0.,0.,0.,0./)
LambdaTB(2,:)=(/0.,0.78,0.,0.,0.,0.,0.,0.,0./)
LambdaTB(3,:)=(/0.,0.,0.63,0.,0.,0.,0.,0.,0./)
LambdaTB(4,:)=(/0.,0.,0.,0.60,0.,0.,0.,0.,0./)
LambdaTB(5,:)=(/0.,0.,0.,0.,0.51,0.,0.,0.,0./)
LambdaTB(6,:)=(/0.,0.,0.,0.,0.,0.41,0.,0.,0./)
LambdaTB(7,:)=(/0.,0.,0.,0.,0.,0.,0.42,0.,0./)
LambdaTB(8,:)=(/0.,0.,0.,0.,0.,0.,0.,0.32,0./)
LambdaTB(9,:)=(/0.,0.,0.,0.,0.,0.,0.,0.,0.31/)


!! Determine the layer number:  
CALL DeterminNlayEff(nlevel,PLEVEL,SfcPress,nLayEff)
  
!!! perform a inital clear sky forward modeling to get the emissivity 1st guess, and analytical emissivity solution.
CALL rttov_fwd_clearsky(imonth,nLayEff,nchannel,incident,longitude,latitude,plevel(1:nLayEff),phalf(1:nLayEff),&
  				QWatm(1:nLayEff),Tatm(1:nLayEff),LST,SfcPress,T2m,Snowc,Smc,TBobs,EmissAnalyt,Emiss1st)
 
LstTune	=	LST
TbTune	=	TBobs
nIters	=	0				

Xb=0.0
Xb(1)=LstTune
Xb(2:10)=EmissAnalyt

Ja=0.0

DO WHILE (.True.)
	nIters=nIters+1

	! input emissin, Tatm, QWatm,
	! output jacobian matrix	
	CALL rttov_fwd_jacobian(nLayEff,nchannel,incident,plevel(1:nLayEff),&
				  Xg(nlevel+1:nlevel+nLayEff),Xg(1:nLayEff),Xg(nlevel+nlevel+1),T2m,Xg(nlevel+nlevel+2:nlevel+nlevel+10),&
				  TbTune,Emss_K, Ta_K(:,1:nLayEff),Qw_K(:,1:nLayEff),Ja(:,,1))
				  
	do ich=1,nchannel
		Ja(1+ich, 1+ich)=Emss_K(ich)
	end do
		
	CALL Convgce(nchannel,TBobs,TbTune,EY,ChiSq,CvgceReached,ChiSqThresh,dY2)
	! print*, nIters, dY2

	IF (CvgceReached) GOTO 221  

    !---Compute the Levenberg-Marquardt optimal solution 
    Q = MATMUL(dble(Ja),MATMUL(dble(Lambda),TRANSPOSE(dble(Ja))))
    Q = dble(LambdaTB) + Q
    ! Q = matinv_dbl(Q)
	CALL matinv_dbl(Q,invM,nchannel)
    G = REAL(MATMUL(dble(Lambda),MATMUL(TRANSPOSE(dble(Ja)),invM)))
    B=TBobs-TbTune+matmul(Ja,DX)
	
    DX = matmul(G,B)  !! update dx	

    !---Compute geophysical vector 
	Xg = Xg + DX

	print*, Xg(nlevel+nlevel+2:nlevel+nlevel+10)
	print*, DX(nlevel+nlevel+2:nlevel+nlevel+10)
	
	IF (nIters.ge.20) GOTO 221 
END DO

stop

221 continue		
! return EmissAnalyt,Emiss1st,Em1DVar,LstTune,nIters
stop 






end subroutine retrieve_core_1dvar

SUBROUTINE DeterminNlayEff(nLay,pres_lay,SfcPress,nLayEff)
  INTEGER               :: nLay,nLayEff,i
  REAL                  :: SfcPress
  REAL,    DIMENSION(nLay) :: pres_lay
  nLayEff=0
  Do i=1,nLay
     IF (pres_lay(i).le.SfcPress .and. SfcPress/pres_lay(i).gt.1.01) THEN
        nLayEff = nLayEff + 1
     ENDIF
  ENDDO
  RETURN
END SUBROUTINE DeterminNlayEff

SUBROUTINE DisabLayBelowSfc(nLay,nLayEff,SaAtm)
  REAL, DIMENSION(nLay,nLay) :: SaAtm
  REAl                 :: SfcPress
  INTEGER              :: nLayEff,i,nLay,iattempt
  !---Disable retrieval of below-sfc layers by setting covar to high value and decorrel.
  !! TOA(1) -> srf(nLayEff) -> bottom (nLay)
  !! decorrel.
  SaAtm(nLayEff+1:nLay,:)  =0.
  SaAtm(:,nLayEff+1:nLay)  =0.
  !! covar to high value
  Do i=nLayEff+1,nLay
     SaAtm(i,i) =0.   
	 !! The source code of MIRS set it to 0., shouldnt we set it to a large number for large uncertainty?
  ENDDO
	
  RETURN
END SUBROUTINE DisabLayBelowSfc


SUBROUTINE ProjCov(nR,nG,Ustar,Sa,Lambda)
  !---Input/Output variables
  INTEGER              :: nR,nG
  REAL, DIMENSION(nG,nR) :: Ustar
  REAl, DIMENSION(nG,nG) :: Sa
  REAL, DIMENSION(nR,nR) :: Lambda
  !---Local variables
  REAl, DIMENSION(nR,nG) :: UstarT
  REAL, DIMENSION(nG,nR) :: A
  UstarT = TRANSPOSE(Ustar(1:nG,1:nR))
  A      = MATMUL(Sa(1:nG,1:nG),Ustar(1:nG,1:nR))
  Lambda(1:nR,1:nR)=MATMUL(UstarT,A)
  RETURN
END SUBROUTINE ProjCov

SUBROUTINE Convgce(nchan,Ym,Y,EY,ChiSq,CvgceReached,ChiSqThresh,dY2)
  INTEGER               :: nchan,nch,ichan
  REAL,  DIMENSION(nchan) :: Ym,Y,EY,dY2
  ! INTEGER, DIMENSION(:) :: ChanSel
  REAL                  :: ChiSq,ChiSqThresh
  LOGICAL               :: CvgceReached

  ChiSq = 0.
  DO ichan=1,nchan
	 dY2(ichan) = 0.
     dY2(ichan)= (Ym(ichan)-Y(ichan))**2.
     ChiSq=ChiSq+dY2(ichan)/EY(ichan)
  ENDDO
  ChiSq=(ChiSq/nchan)
  CvgceReached=.FALSE.
  IF (ChiSq.le.ChiSqThresh) CvgceReached=.TRUE.
  RETURN
END SUBROUTINE Convgce

subroutine matinv_dbl(A,retM,nch)
    ! Invert matrix by Gauss method
    ! --------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n,nch
    REAL*8, intent(in),dimension(nch,nch) :: A
    REAL*8, dimension(size(A,1),size(A,2)) :: b
    REAL*8, intent(out), dimension(size(A,1),size(A,2)) :: retM
    REAL*8, dimension(size(A,1)) :: temp
    ! - - - Local Variables - - -
    REAL( SELECTED_REAL_KIND(15) ) :: c, d
    INTEGER :: i, j, k, m, imax(1), ipvt(size(A,1))
    ! - - - - - - - - - - - - - -
    b = A
    n=size(A,1)
    retM=A
    ipvt = (/ (i, i = 1, n) /)
    ! Go into loop- b, the inverse, is initialized equal to a above
    DO k = 1,n
       ! Find the largest value and position of that value
       imax = MAXLOC(ABS(b(k:n,k)))
       m = k-1+imax(1)
       !   sigular matrix check
       if(ABS(b(m,k)).LE.(1.D-40)) then
          !CALL ErrHandl(ErrorType,Err_SingularMatrx,'')
          retM(1,1) = -99999999.0
          return 
       ENDIF
       ! get the row beneath the current row if the current row will
       ! not compute
       IF (m .ne. k) THEN
          ipvt( (/m,k/) ) = ipvt( (/k,m/) )
          b((/m,k/),:) = b((/k,m/),:)
       END IF
       ! d is a coefficient - brings the pivot value to one and then is applied
       ! to the rest of the row
       d = 1/b(k,k)
       temp = b(:,k)
       DO j = 1, n
          c = b(k,j)*d
          b(:,j) = b(:,j)-temp*c
          b(k,j) = c
       END DO
       b(:,k) = temp*(-d)
       b(k,k) = d
    END DO
    retM(:,ipvt) = b 
end subroutine matinv_dbl