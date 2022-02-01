program da_blending

!----------------------------------------------------------------------
! Purpose:  WRFDA backgound/analysis blending scheme   
!
! Input     : fg          -- WRF model forecast (eg. wrfout_d01) 
!             bg          -- Host model forecast (eg. wrfinput_d01) 
!
! Output    : fg_blend    -- Blended background for WRFDA 
!           : filtered_fg -- Filtered WRF model forecast  (optional -debug 1)
!           : filtered_bg -- Filtered Host model forecast (optional -debug 1)
!           : filtered_diff -- Difference between filtered_bg and filtered_fg (optional -debug 1)
! Compile:
!
!   pgf90 -o da_blending.exe -I$NETCDF/include -L$NETCDF/lib -lnetcdf \
!         da_blending.f90 raymond.f
!
! Usage:
!
!   da_blending.exe  [-fg filename] [-bg filename] [-o outputfile] 
!                    [-Lx cut-off length-scale ] [-debug level] [-h]
!
!     -fg    Optional, WRFDA background from WRF forecast                 default - fg
!     -bg    Optional, Host model forecast interpolated into WRF grid     default - bg
!     -o     Optional, Blended bacground for WRFDA                        default - fg_blend
!     -Lx    Optional, Cut-off length-scale (km) for Raymont filter       default - 1200 km 
!     -debug Optional, debug level                      default - 0 ; > 0: output filtered files.
!     -h     Show this help
!
! jliu@ucar.edu     2012-09-15
! hlwang@ucar.edu   2013-01-31
!----------------------------------------------------------------------

  use netcdf

  implicit none

  character (len=6), dimension(1:13) :: vNam 

  integer :: i, j, n, t, ierr
  integer :: nLat_fg, nLon_fg, nLev_fg
  integer :: nLat_bg, nLon_bg, nLev_bg
  integer :: nLon, nLat, nLev, nTim,nbdy

  integer :: ncidfg, ncidbg, ncidout, ncidfg_filtered, ncidbg_filtered, nciddiff_filtered
  integer :: nDims, dLen, Dimid
  integer :: varid, varid_fg, varid_bg

  real, dimension(:,:,:,:), allocatable :: fg, bg, var_out, fg_filtered, bg_filtered,var_work
  real, dimension(:,:),     allocatable :: field, field_work
  real                                  :: eps, Lx 
  real                                  :: Dx_fg, Dx_bg, Dx 

  integer                               :: dbg_lev

  integer, dimension(nf90_max_var_dims) :: vDimIDs
  integer, dimension(4)                 :: vdimsizes

  character (len = 19)   :: times_fg, times_bg
  character (len = 255)  :: out_path        = ""
  character (len = 255)  :: appname         = ""
  character (len = 255)  :: arg             = ""
  character (len = 255)  :: f_fg            = "fg"
  character (len = 255)  :: f_bg            = "bg"
  character (len = 255)  :: f_out           = "fg_blend"
  character (len = 255)  :: f_fg_filtered   = "fg_filtered"
  character (len = 255)  :: f_bg_filtered   = "bg_filtered"
  character (len = 255)  :: f_diff_filtered = "diff_filtered"
  character (len = 255)  :: errmsg          = ""

  LOGICAL :: file_exists

  integer iargc

  real, parameter :: R_ZERO = 1.0E-9
  real            :: PI

  nbdy = 40

  PI = 4.D0*DATAN(1.D0)
  Lx = 0.D0
  dbg_lev = 0

  !Fields will be blended
  vNam(1)="U"
  vNam(2)="V"
  vNam(3)="T"
  vNam(4)="QVAPOR"
  vNam(5)="PH"
  vNam(6)="P"
  vNam(7)="MU"
  vNam(8)="U10"
  vNam(9)="V10"
  vNam(10)="T2"
  vNam(11)="Q2"
  vNam(12)="PSFC"
  vNam(13)="TH2"

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
      case ("-bg")
        call getarg(i+1, arg)
        f_bg=trim(arg)
      case ("-o")
        call getarg(i+1, arg)
        f_out=trim(arg)
      case ("-lx")
        call getarg(i+1, arg)
        if ( index(arg, '.') .eq. 0 )  arg=trim(arg)//".0"
        read(arg, '(F)') Lx
      case ("-debug")
        call getarg(i+1, arg)
        read(arg, '(I)') dbg_lev
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

  inquire(FILE=trim(f_bg), EXIST=file_exists)

  if ( .not. file_exists ) then
    Write(*,*) "\nError: "//trim(f_bg)//" not exists\n"
    call show_usage()
    call exit(-1)
  endif

  if ( trim(f_fg) .eq. trim(f_out) .or. trim(f_bg) .eq. trim(f_out) ) then
    Write(*,*) "\nError: Output file same as input file\n"
    call exit(-1)
  endif

  ierr = nf90_open(trim(f_fg), NF90_NOWRITE, ncidfg)
  errmsg = trim(f_fg)
  if ( ierr /= nf90_noerr ) call nf90_handle_err(ierr, errmsg)

  ierr = nf90_open(trim(f_bg), NF90_NOWRITE, ncidbg)
  errmsg = trim(f_bg)
  if ( ierr /= nf90_noerr ) call nf90_handle_err(ierr, errmsg)

  ierr = nf90_get_att(ncidfg, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION",  nLon_fg) 
  ierr = nf90_get_att(ncidfg, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION",nLat_fg)
  ierr = nf90_get_att(ncidfg, NF90_GLOBAL, "BOTTOM-TOP_GRID_DIMENSION", nLev_fg) 
  ierr = nf90_get_att(ncidfg, NF90_GLOBAL, "DX",                          Dx_fg) 

  ierr = nf90_get_att(ncidbg, NF90_GLOBAL, "WEST-EAST_GRID_DIMENSION",  nLon_bg)
  ierr = nf90_get_att(ncidbg, NF90_GLOBAL, "SOUTH-NORTH_GRID_DIMENSION",nLat_bg)
  ierr = nf90_get_att(ncidbg, NF90_GLOBAL, "BOTTOM-TOP_GRID_DIMENSION", nLev_bg)
  ierr = nf90_get_att(ncidbg, NF90_GLOBAL, "DX",                          Dx_bg)

  if ( nLon_fg /= nLon_bg .or. &
       nLat_fg /= nLat_bg .or. &
       nLev_fg /= nLev_bg ) then
    Write(*,fmt='(a)')      "Error : Domain size NOT match"
    Write(*,fmt='(a,3(I5.1,1X))') "  WRF forecast = ", nLon_fg, nLat_fg, nLev_fg
    Write(*,fmt='(a,3(I5.1,1X))') "  Host model   = ", nLon_bg, nLat_bg, nLev_bg
    call exit(-1) 
  end if

  if ( (Dx_fg - Dx_bg) > R_ZERO ) then
    Write(*,fmt='(a)')      "Error : Grid size NOT match"
    Write(*,fmt='(a,I5.1)') "  WRF forecast = ", Dx_fg
    Write(*,fmt='(a,I5.1)') "  Host model   = ", Dx_bg
    call exit(-1)
  end if

  if (  Lx * ( -1. ) >  R_ZERO ) then
    Write(*,fmt='(a)')      "Error : Negative Lx"
    call exit(-1)
  end if

  if ( ( Dx_fg - 1000.) > R_ZERO ) Dx_fg = Dx_fg / 1000. 

!  if ( ( Lx - 0. ) < R_ZERO ) Lx = 4 * Dx_fg
  if ( ( Lx - 0. ) < R_ZERO ) Lx = 1200 

  eps = (tan(PI*Dx_fg/Lx))**(-6)

  Write(*,*) " Input "
  Write(*,*) "   WRF forecast                  : "//trim(f_fg)
  Write(*,*) "   Host model analysis/forecast  : "//trim(f_bg)
  Write(*,*) "   Lx                            : ",Lx
  Write(*,*) "   Dx                            : ",Dx_fg
  Write(*,*) "   eps                           : ",eps
  Write(*,*) "Output "
  Write(*,*) "   Blended background file       : "//trim(f_out)
  !STOP
  errmsg = ""

  call system("cp "//f_fg//" "//f_out)
  ierr = nf90_open(f_out, NF90_WRITE, ncidout)
  errmsg= trim(f_out)
  if ( ierr /= nf90_noerr ) call nf90_handle_err(ierr, errmsg)

  if ( dbg_lev > 0 ) then
    n=index(f_out, '/', BACK=.true.)
    out_path = trim(f_out(1:n))
    !f_fg_filtered=trim(out_path)//"filtered_fg"
    !f_bg_filtered=trim(out_path)//"filtered_bg"
    f_diff_filtered=trim(out_path)//"filtered_diff"
    !call system("cp "//trim(f_fg)//" "//trim(f_fg_filtered))
    !call system("cp "//trim(f_bg)//" "//trim(f_bg_filtered))
    call system("cp "//trim(f_fg)//" "//trim(f_diff_filtered))
    !ierr = nf90_open(f_fg_filtered, NF90_WRITE, ncidfg_filtered)
    !ierr = nf90_open(f_bg_filtered, NF90_WRITE, ncidbg_filtered)
    ierr = nf90_open(f_diff_filtered, NF90_WRITE, nciddiff_filtered)
  endif

  n = ubound(vNam,1)
  do i=1,n 

    Write (*,*) "Blending backgrounds for "//trim(vNam(i))

    ierr = nf90_inq_varid(ncidout, trim(vNam(i)), varid)
    if ( ierr /= 0 ) cycle ! CSS added

    ierr = nf90_inq_varid(ncidfg, trim(vNam(i)), varid_fg)
    if ( ierr /= 0 ) cycle ! CSS added

    ierr = nf90_inq_varid(ncidbg, trim(vNam(i)), varid_bg)
    if ( ierr /= 0 ) cycle ! CSS added

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
    allocate(     fg(nLon, nLat, nLev, nTim), stat=ierr)
    allocate(     bg(nLon, nLat, nLev, nTim), stat=ierr)
    allocate(  field(nLon*nLat,nLev),         stat=ierr)

    allocate(var_work(nLon+nbdy, nLat+nbdy, nLev, nTim), stat=ierr)
    allocate(field_work((nLon+nbdy)*(nLat+nbdy), nLev), stat=ierr)

    var_work = 0.0
    field_work = 0.0

    if ( dbg_lev > 0 ) then
      allocate(fg_filtered(nLon, nLat, nLev, nTim), stat=ierr)
      allocate(bg_filtered(nLon, nLat, nLev, nTim), stat=ierr)
    endif

    ierr = nf90_get_var(ncidout, varid,    var_out)
    ierr = nf90_get_var(ncidfg,  varid_fg, fg)
    ierr = nf90_get_var(ncidbg,  varid_bg, bg)

    do t=1, nTim

      var_work(nbdy/2:nLon+nbdy/2,nbdy/2:nLat+nbdy/2,:,t) = bg(:,:,:,t)-fg(:,:,:,t)

      !field = reshape(fg(:,:,:,t), (/nLon*nLat,nLev/))
      field_work = reshape(var_work(:,:,:,t),(/(nLon+nbdy)*(nLat+nbdy),nLev/))
      call RAYMOND (field_work ,nLon+nbdy,nLat+nbdy,nLev,EPS)
      var_work(:,:,:,t) = reshape(field_work, (/nLon+nbdy, nLat+nbdy, nLev/))      
      var_out(:,:,:,t)  = var_work(nbdy/2:nLon+nbdy/2,nbdy/2:nLat+nbdy/2,:,t) 
      !field = reshape(var_out(:,:,:,t), (/nLon*nLat,nLev/))
      var_out(:,:,:,t) = var_out(:,:,:,t) + fg(:,:,:,t)

      !field = reshape(bg(:,:,:,t), (/nLon*nLat, nLev/))
      !call RAYMOND (field ,nLon,nLat,nLev,EPS)

      if ( dbg_lev > 0 ) then 
        !bg_filtered(:,:,:,t) = reshape(field, (/nLon, nLat, nLev/))
        !fg_filtered(:,:,:,t) = fg(:,:,:,t) - (var_out(:,:,:,t)-fg_filtered(:,:,:,t))
      end if
    end do

    ierr = nf90_put_var(ncidout, varid, var_out)

    deallocate(var_out, stat=ierr)
    deallocate(field,   stat=ierr)
    deallocate(fg,      stat=ierr)
    deallocate(bg,      stat=ierr)
    deallocate(var_work, stat=ierr)
    deallocate(field_work,   stat=ierr)

    if ( dbg_lev > 0 ) then

      !ierr = nf90_put_var(ncidfg_filtered,   varid_fg, fg_filtered)
      !ierr = nf90_put_var(ncidbg_filtered,   varid_bg, bg_filtered)
      ierr = nf90_put_var(nciddiff_filtered, varid_fg,var_out-fg )

      deallocate(fg_filtered, stat=ierr)
      deallocate(bg_filtered, stat=ierr)

    endif

  end do

  ierr = nf90_close(ncidfg)
  ierr = nf90_close(ncidbg)
  ierr = nf90_close(ncidout)

  if ( dbg_lev > 0 ) then
    !ierr = nf90_close(ncidfg_filtered)
    !ierr = nf90_close(ncidbg_filtered)
    ierr = nf90_close(nciddiff_filtered)
  endif

  Write(*,*) "Blending completed successfully"

contains
  subroutine show_usage()
     Write(*,*) 'Usage : '//trim(appname)// &
     '[-h] [-fg filename] [-bg filename] [-o outputfile] [-Lx cut-off length-scale]'
     Write(*,*) "  -fg    Optional, first guest from last forecast, default - fg"
     Write(*,*) "  -bg    Optional, host model analysis/forecast,   default - bg"
     Write(*,*) "  -o     Optional, blended first guest,            default - fg_blend"
     Write(*,*) "  -Lx    Optional, cut-off length-scale in km      default - 4 * Dx"
     Write(*,*) "  -debug Optional, debug level                     default - 0"
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

end program da_blending
