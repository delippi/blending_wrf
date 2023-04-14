program da_diagnose

!----------------------------------------------------------------------
! Purpose:  diagnose for WRFDA backgound/analysis blending scheme   
!
! Input     : fg               -- WRF model forecast (eg. wrfout_d01) 
!             fg_blend         -- Blended background for WRFDA 
!
! Output    : fg_blend_diag    -- Diagnosed for WRFDA Blended background
! 
! Compile:
!
!   pgf90 -o da_diagnose.exe -I$NETCDF/include -L$NETCDF/lib -lnetcdf \
!         da_diagnose.f90
!
! Usage:
!
!   da_diagnose.exe  [-fg filename] [-fgb filename] [-o outputfile] [-h]
!
!     -fg    Optional, WRFDA background from WRF forecast                 default - fg
!     -fgb   Optional, Blended background for WRFDA                       default - fg_blend
!     -o     Optional, Diagnosed for WRFDA Blended background             default - fg_blend_diag
!     -h     Show this help
!
! LF     2013-02-16
!----------------------------------------------------------------------

  use netcdf

  implicit none

  integer, parameter :: total_var = 13
  character (len=6), dimension(1:total_var) :: vNam

  integer :: i, j, n, t, ierr
  integer :: nLat_fg,  nLon_fg,  nLev_fg
  integer :: nLat_fgb, nLon_fgb, nLev_fgb
  integer :: nLon, nLat, nLev, nTim
  integer :: ii, jj, kk

  integer :: ncidfg, ncidfgb, ncidout
  integer :: nDims, dLen, Dimid
  integer :: varid, varid_fg, varid_fgb

  real, dimension(:,:,:,:), allocatable :: fg, fgb, var_out
  real, dimension(:,:,:), allocatable :: fg_12, fgb_12, diff_12
  real, dimension(:,:,:), allocatable :: fg_13, fgb_13, diff_13
  real, dimension(:,:), allocatable :: fg_07, fgb_07, diff_07
  real, dimension(:,:,:), allocatable :: fg_11, fg_10
  real, dimension(:,:), allocatable :: fg_09, fg_08, fg_06
  real, dimension(:), allocatable :: fg_05, fg_04, fg_03, fg_02, fg_01
  real, dimension(:,:,:,:), allocatable :: varr_fg,varr_fgb,varr_diff

  real, dimension(:,:,:),   allocatable :: qw,qvapor,PH1,PHB,t1,tw,PH
  real, dimension(:,:),     allocatable :: MU1,psw,MUB_2d,HGT_2d,MU_2d
  real, dimension(:),       allocatable :: DNW_1d,ZNU_1d,RDNW_1d,RDN_1d
  real                                  :: Dx_fg, Dx_fgb, Dx, ptop 

  integer                               :: dbg_lev

  integer, dimension(nf90_max_var_dims) :: vDimIDs
  integer, dimension(4)                 :: vdimsizes

  character (len = 19)   :: times_fg, times_fgb
  character (len = 255)  :: appname         = ""
  character (len = 255)  :: arg             = ""
  character (len = 255)  :: f_fg            = "fg"
  character (len = 255)  :: f_fgb           = "fg_blend"
  character (len = 255)  :: f_out           = "fg_blend_diag"
  character (len = 255)  :: errmsg          = ""

  LOGICAL :: file_exists

  integer iargc

  real, parameter :: R_ZERO = 1.0E-9


  vNam(1)="P_TOP"
  vNam(2)="RDN"
  vNam(3)="RDNW"
  vNam(4)="ZNU"
  vNam(5)="DNW"
  vNam(6)="HGT"
  vNam(7)="PSFC"
  vNam(8)="MUB"
  vNam(9)="MU"
  vNam(10)="PHB"
  vNam(11)="PH"
  vNam(12)="QVAPOR"
  vNam(13)="T"

  call getarg(0, appname)
  n=index(appname, '/', BACK=.true.)
  appname = trim(appname(n+1:))

  DO i = 1, iargc(), 2
    arg=""
    call getarg(i, arg)
    call lcase(arg)
    select case ( trim(arg) )
      case ("-fg")
        call getarg(i+1, arg)
        f_fg=trim(arg)
      case ("-fgb")
        call getarg(i+1, arg)
        f_fgb=trim(arg)
      case ("-o")
        call getarg(i+1, arg)
        f_out=trim(arg)
      case default
        call show_usage()
        call exit(0)
    end select
  END DO

  inquire(FILE=trim(f_fg), EXIST=file_exists)

  if ( .not. file_exists ) then
    Write(*,*) "\nError: "//trim(f_fg)//" not exists\n"
    call show_usage()
    call exit(-1)
  endif

  inquire(FILE=trim(f_fgb), EXIST=file_exists)

  if ( .not. file_exists ) then
    Write(*,*) "\nError: "//trim(f_fgb)//" not exists\n"
    call show_usage()
    call exit(-1)
  endif

  if ( trim(f_fg) .eq. trim(f_out) .or. trim(f_fgb) .eq. trim(f_out) ) then
    Write(*,*) "\nError: Output file same as input file\n"
    call exit(-1)
  endif

  ierr = nf90_open(trim(f_fg), NF90_NOWRITE, ncidfg)
  errmsg = trim(f_fg)
  if ( ierr /= nf90_noerr ) call nf90_handle_err(ierr, errmsg)

  ierr = nf90_open(trim(f_fgb), NF90_NOWRITE, ncidfgb)
  errmsg = trim(f_fgb)
  if ( ierr /= nf90_noerr ) call nf90_handle_err(ierr, errmsg)

  ierr = nf90_get_att(ncidfg, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION",  nLon_fg) 
  ierr = nf90_get_att(ncidfg, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION",nLat_fg)
  ierr = nf90_get_att(ncidfg, NF90_GLOBAL, "BOTTOM-TOP_GRID_DIMENSION", nLev_fg) 
  ierr = nf90_get_att(ncidfg, NF90_GLOBAL, "DX",                          Dx_fg) 

  ierr = nf90_get_att(ncidfgb, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION",  nLon_fgb)
  ierr = nf90_get_att(ncidfgb, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION",nLat_fgb)
  ierr = nf90_get_att(ncidfgb, NF90_GLOBAL, "BOTTOM-TOP_GRID_DIMENSION", nLev_fgb)
  ierr = nf90_get_att(ncidfgb, NF90_GLOBAL, "DX",                          Dx_fgb)

  if ( nLon_fg /= nLon_fgb .or. &
       nLat_fg /= nLat_fgb .or. &
       nLev_fg /= nLev_fgb ) then
    Write(*,fmt='(a)')      "Error : Domain size NOT match"
    Write(*,fmt='(a,3(I5.1,1X))') "  WRF forecast = ", nLon_fg,  nLat_fg,  nLev_fg
    Write(*,fmt='(a,3(I5.1,1X))') "  Host model   = ", nLon_fgb, nLat_fgb, nLev_fgb
    call exit(-1) 
  end if

  if ( (Dx_fg - Dx_fgb) > R_ZERO ) then
    Write(*,fmt='(a)')      "Error : Grid size NOT match"
    Write(*,fmt='(a,I5.1)') "  WRF forecast = ", Dx_fg
    Write(*,fmt='(a,I5.1)') "  Host model   = ", Dx_fgb
    call exit(-1)
  end if

  if ( ( Dx_fg - 1000.) > R_ZERO ) Dx_fg = Dx_fg / 1000. 

  Write(*,*) " Input "
  Write(*,*) "   WRF forecast                  : "//trim(f_fg)
  Write(*,*) "   Blended background for WRFDA  : "//trim(f_fgb)
  Write(*,*) "   Dx                            : ",Dx_fg
  Write(*,*) "Output "
  Write(*,*) "   Diagnosed background file     : "//trim(f_out)
  !STOP
  errmsg = ""

  call system("cp "//f_fgb//" "//f_out)
  ierr = nf90_open(f_out, NF90_WRITE, ncidout)
  errmsg= trim(f_out)
  if ( ierr /= nf90_noerr ) call nf90_handle_err(ierr, errmsg)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(12)), varid_fg)
    ierr = nf90_inq_varid(ncidfgb, trim(vNam(12)), varid_fgb)
    ierr = nf90_inquire_variable(ncidfgb, varid_fgb, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfgb, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(12), nLon, nLat, nLev, nTim
  if ( nTim /= 1 ) then
    Write(*,*) "\nError: Input file request Time=1\n"
    call exit(-1)
  endif
    allocate(  fg_12(nLon, nLat, nLev), stat=ierr)
    allocate( fgb_12(nLon, nLat, nLev), stat=ierr)
    allocate(diff_12(nLon, nLat, nLev))
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_12)
    ierr = nf90_get_var(ncidfgb, varid_fgb,fgb_12)

    diff_12(:,:,:)=fgb_12(:,:,:)-fg_12(:,:,:)
    !print*,diff_12(100,100,3)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(13)), varid_fg)
    ierr = nf90_inq_varid(ncidfgb, trim(vNam(13)), varid_fgb)
    ierr = nf90_inquire_variable(ncidfgb, varid_fgb, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfgb, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(13), nLon, nLat, nLev, nTim
    allocate(  fg_13(nLon, nLat, nLev), stat=ierr)
    allocate( fgb_13(nLon, nLat, nLev), stat=ierr)
    allocate(diff_13(nLon, nLat, nLev))
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_13)
    ierr = nf90_get_var(ncidfgb, varid_fgb,fgb_13)

    diff_13(:,:,:)=fgb_13(:,:,:)-fg_13(:,:,:)
    !print*,diff_13(100,100,3)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(7)), varid_fg)
    ierr = nf90_inq_varid(ncidfgb, trim(vNam(7)), varid_fgb)
    ierr = nf90_inquire_variable(ncidfgb, varid_fgb, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfgb, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(7), nLon, nLat, nLev, nTim
    allocate(  fg_07(nLon, nLat), stat=ierr)
    allocate( fgb_07(nLon, nLat), stat=ierr)
    allocate(diff_07(nLon, nLat))
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_07)
    ierr = nf90_get_var(ncidfgb, varid_fgb,fgb_07)

    diff_07(:,:)=fgb_07(:,:)-fg_07(:,:)
    !print*,diff_07(100,100)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(11)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
      ii = nLon
      jj = nLat
      kk = nLev
    !print*, vNam(11), nLon, nLat, nLev, nTim
    allocate(  fg_11(nLon, nLat, nLev), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_11)
    !print*,fg_11(100,100,3)
    !print*,ii,jj,kk

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(10)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(10), nLon, nLat, nLev, nTim
    allocate(  fg_10(nLon, nLat, nLev), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_10)
    !print*,fg_10(100,100,3)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(9)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(9), nLon, nLat, nLev, nTim
    allocate(  fg_09(nLon, nLat), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_09)
    !print*,fg_09(100,100)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(8)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(8), nLon, nLat, nLev, nTim
    allocate(  fg_08(nLon, nLat), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_08)
    !print*,fg_08(100,100)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(6)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(6), nLon, nLat, nLev, nTim
    allocate(  fg_06(nLon, nLat), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_06)
    !print*,fg_06(100,100)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(5)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(5), nLon, nLat, nLev, nTim
    allocate(  fg_05(nLon), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_05)
    !print*,fg_05(30)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(4)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(4), nLon, nLat, nLev, nTim
    allocate(  fg_04(nLon), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_04)
    !print*,fg_04(30)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(3)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(3), nLon, nLat, nLev, nTim
    allocate(  fg_03(nLon), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_03)
    !print*,fg_03(30)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(2)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(2), nLon, nLat, nLev, nTim
    allocate(  fg_02(nLon), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_02)
    !print*,fg_02(30)

    ierr = nf90_inq_varid(ncidfg,  trim(vNam(1)), varid_fg)
    ierr = nf90_inquire_variable(ncidfg, varid_fg, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidfg, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    !print*, vNam(1), nLon, nLat, nLev, nTim
    allocate(  fg_01(nLon), stat=ierr)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg_01)
    !print*,fg_01(1)

    allocate ( PH(ii, jj, kk) )
    allocate ( MU_2d(ii,  jj) )

    call balance(ii, jj, kk, &
                 diff_12,fg_12,fg_11,fg_10,fg_13,diff_13, &
                 fg_09,diff_07,fg_08,fg_06, &
                 fg_05,fg_04,fg_03,fg_02,fg_01, &
                 PH,MU_2d)

    ierr = nf90_inq_varid(ncidout, "PH", varid)
    ierr = nf90_inquire_variable(ncidout, varid, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidout, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)

    allocate(var_out(nLon, nLat, nLev, nTim), stat=ierr)
    ierr = nf90_get_var(ncidout, varid,    var_out)
      var_out(:,:,:,1) = PH(:,:,:)
    ierr = nf90_put_var(ncidout, varid, var_out)
    deallocate(var_out, stat=ierr)

    ierr = nf90_inq_varid(ncidout, "MU", varid)
    ierr = nf90_inquire_variable(ncidout, varid, ndims=nDims,dimids=vDimIDs)
    vdimsizes = 1
    do j=1, nDims
      ierr = nf90_inquire_dimension(ncidout, vDimIDs(j), len = dLen )
      vdimsizes(j) = dLen
    end do
    nLon = vdimsizes(1)
    nLat = vdimsizes(2)
    nLev = vdimsizes(3)
    nTim = vdimsizes(4)
    allocate(var_out(nLon, nLat, nLev, nTim), stat=ierr)
    ierr = nf90_get_var(ncidout, varid,    var_out)
      var_out(:,:,1,1) = MU_2d(:,:)
    ierr = nf90_put_var(ncidout, varid, var_out)
    deallocate(var_out, stat=ierr)

  ierr = nf90_close(ncidfg)
  ierr = nf90_close(ncidfgb)
  ierr = nf90_close(ncidout)

    deallocate(fg_12,  stat=ierr)
    deallocate(fgb_12, stat=ierr)
    deallocate(diff_12 )
    deallocate(fg_13,  stat=ierr)
    deallocate(fgb_13, stat=ierr)
    deallocate(diff_13 )
    deallocate(fg_07,  stat=ierr)
    deallocate(fgb_07, stat=ierr)
    deallocate(diff_07 )
    deallocate(fg_11,  stat=ierr)
    deallocate(fg_10,  stat=ierr)
    deallocate(fg_09,  stat=ierr)
    deallocate(fg_08,  stat=ierr)
    deallocate(fg_06,  stat=ierr)
    deallocate(fg_05,  stat=ierr)
    deallocate(fg_04,  stat=ierr)
    deallocate(fg_03,  stat=ierr)
    deallocate(fg_02,  stat=ierr)
    deallocate(fg_01,  stat=ierr)
    deallocate( PH )
    deallocate( MU_2d )

  Write(*,*) "Diagnose completed successfully"

contains
  subroutine show_usage()
     Write(*,*) 'Usage : '//trim(appname)// &
     '[-h] [-fg filename] [-fgb filename] [-o outputfile]'
     Write(*,*) "  -fg    Optional, first guest from last forecast, default - fg"
     Write(*,*) "  -fgb   Optional, Blended background for WRFDA,   default - fg_blend"
     Write(*,*) "  -o     Optional, Diagnosed background,           default - fg_blend_diag"
     Write(*,*) "  -h     Show this help"
  end subroutine show_usage

  subroutine nf90_handle_err(ierr, errmsg)
    integer,          intent(in) :: ierr
    character(len=*), intent(in) :: errmsg

    if(ierr /= nf90_noerr) then
      Write(*,*) trim(nf90_strerror(ierr))//" : "//trim(errmsg)
      call system("rm "//f_out)
      STOP
    end if
  end subroutine nf90_handle_err

  subroutine lcase(string)
  ! convert a word to lower case
    character (len=*) , intent(in out) :: string
    integer :: i,ic,nlen
    nlen = len(trim(string))
    do i=1,nlen
      ic = ichar(string(i:i))
      if (ic >= 65 .and. ic < 90) string(i:i) = char(ic+32)
    end do
  end subroutine lcase 

end program da_diagnose
