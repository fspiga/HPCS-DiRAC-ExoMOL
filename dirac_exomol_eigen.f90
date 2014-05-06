module d_module
 
!dec$ define blacs_  = 1
!dec$ define mpi_    = 0


  use accuracy
  use timer
  implicit none

  public FLReadInput

  !
!  integer(ik),parameter:: verbose = 4
  integer(ik),parameter:: verbose = 2
  integer, parameter :: trk        = selected_real_kind(12)
  character(len=cl)    :: matrix_file='matrix'
  !
  integer(ik),parameter ::  maxnprocs=1024
  integer(ik)           ::  iterout = 100
  !
  logical :: print_vectors = .false.
  !
  real(rk) :: walltime=1e7,time_factor = 0.95
  !
  real(trk):: rtime
  !
  logical :: gen_mat = .false.
  integer(ik) :: mat_len
#if 1
  integer(hik)	::	seed(1000000)
#endif
  !
  contains
  !
  subroutine FLReadInput(Jrot,gamma,factor,nroots,tol,sparse,eigensolver,chkpoint,zpe,memory,energy_thresh,coef_thresh)
    !
    use input
    !
    integer(ik),intent(out)            :: Jrot,gamma,nroots,chkpoint,eigensolver
    logical,intent(out)                :: sparse
    !
    integer                            :: sparse_
    !
    real(rk),intent(out)               :: factor,tol,memory,zpe,energy_thresh,coef_thresh
    character(len=cl)                  :: diagonalizer
    !
    character(len=cl) :: w
    !
    logical :: eof
    !
    !write(out,"('Read the input')")
    !
    sparse = 0
    !
    jrot = -1
    gamma = -1
    nroots = 1e6
    energy_thresh = 1.0d+05
    coef_thresh = 1.0d-16
    diagonalizer = "PDSYEVD" ! default
    factor = 1.0_rk
    tol = small_
    sparse = .false.
    walltime = safe_max
    zpe = -small_
    memory = 256.0_rk
    !
    call input_options(echo_lines=.false.,error_flag=1)
    !
    ! read the general input 
    !
    do
        call read_line(eof) ; if (eof) exit
        call readu(w)
        select case(w)
        case("STOP","FINISH","END")
          exit
        case("")
          print "(1x)"    !  Echo blank lines
          !
        case ("J","JROT")
          !
          call readi(jrot)
          !
        case ("GAMMA","SYM","SYMMETRY")
          !
          call readi(gamma)
          !
        case ("ZPE")
          !
          call readf(zpe)
          !
        case ("ENERGY_THRESH")
          !
          call readf(energy_thresh)
          !
        case ("COEF_THRESH")
          !
          call readf(coef_thresh)
          !
        case ("FACTOR")
          !
          call readf(factor)
          !
        case ("CHECKPOINT")
          !
          ! chkpoint :  0 - none, 1 - save, 2  - read
          !
          call readi(chkpoint)
          !
          if (chkpoint<0.or.chkpoint>2) then 
           !
           write (out,"(' Illegal checkpoint value = ',i)") chkpoint
           stop 'Illegal checkpoint value '
           !
          endif
          !
        case ("WALLCLOCK","MAXTIME","WALLTIME")
          !
          call readf(walltime)
          !
          ! convert to seconds
          !
          walltime = walltime*3600.0_rk
          !
        case ("MEMORY","MEM")
          !
          call readf(memory)
          !
        case ("NROOTS")
          !
          call readi(nroots)
          !
        case ("ITEROUT")
          !
          call readi(iterout)
          !
        case ("TOL","TOLARENCE")
          !
          call readf(tol)
          !
        case ("GEN-MAT")
        	call readi(mat_len)
        	gen_mat = .true.
        case ("SPARSE")
          !
          !call readi(sparse_)
          !
          sparse = .true.
          !
        case ("DIAGONALIZER")
          !
          call readu(diagonalizer)
          !
          select case ( trim(diagonalizer) )
            !
          case ('PDSYEVD')
            !
            eigensolver = 1
            !
          case ('PDSYEVX')
            !
            eigensolver = 2
            !
          case ('ELPA-1STAGE')
            !
            eigensolver = 3
            !
          case ('ELPA-2STAGE')
            !
            eigensolver = 4
            !
          case default
            !
            write(out,'("Uknown eigensolver = ",a)'), diagonalizer
            stop "Uknown eigensolver"
            !
          end select
          !
        case default
          call report ("Principal keyword "//trim(w)//" not recognized",.true.)
        end select
        !
    end do
    !
    ! write(out,"('...done!')")
    !
  end subroutine FLReadInput

#if 1
  integer(hik) function  rrr(i)
  	integer(ik)	::	i
  	seed(i) = mod((1103515245 * seed(i) + 12345),1099511627776);
  	rrr=seed(i)
  	return
  end function rrr
#endif
  

end module d_module


program dirac_exomol_eigen
  
  use accuracy
  use d_module
  use timer
#if defined(__ELPA)
  use elpa1
  use elpa2
#endif
  use mpi

  implicit none

  !integer(ik) :: verbose=6

    ! Related to matrix reading
    integer(ik) ::              jrot, gamma, chkptIO, info, i, j, i_loc, j_loc, iroot, idimen, jdimen
    character(cl) ::            jchar, symchar, filename, buf
    integer(hik)                matsize
    integer(ik)                 dimen_s, nelem, max_nelem
    integer(ik), allocatable :: bterm(:,:)
    real(rk), allocatable ::    a_temp(:), a_loc(:,:)   
    real(rk) ::                 thresh=30000.0_rk, zpe, energy_thresh, coef_thresh,rrr_value
    ! Related to diagonalization
    integer(ik)  ::     nroots, nsolv
    ! Parameters
    integer(ik)         lwork, liwork
    double precision    zero, t1, t2, work_(10)
    parameter           ( zero = 0.0d0 )
    integer(ik)         loc_r, loc_c, lda
    double precision    mone
    integer(ik)         maxprocs,global_i,global_j
    parameter           ( mone = -1.0d0, maxprocs = 512 )
    ! Local Scalars
    integer(ik)         context, iam, m, mycol, myrow, nb, npcol, nprocs, nprow, nz, trilwmin, proc_row, proc_col,iwork_(10),l_nrows, l_ncols, l_eigvec, iprow, ipcol
#if defined(__ELPA)
    integer(ik)         mpi_comm_rows, mpi_comm_cols
#endif
    ! Local Arrays
    integer(ik)                      desca(50), descz(50), irecl
    integer(ik), allocatable ::      iwork(:)
    double precision, allocatable :: work(:), w(:), z_loc(:,:), eigvec(:),local_vecs(:)
    !
    double precision :: vl,vu
    integer(ik)      :: il,iu,nvals,nvects,nvalsmax,nvals_,nvects_,chkpoint,eigensolver
    !
    real(rk)         :: gfactor,tol,memory_now,memory_max,memory
    !
    logical   :: sparse = .false.
    !
    integer(ik), external :: numroc, iceil, indxl2g
    !real(rk), external ::    MPI_Wtime, pdlamch
    real(rk), external ::    pdlamch
    !
    character(len=1)   :: range
    !
    real(rk)           :: abstol,orfac
    integer(ik),allocatable     ::  iclustr(:),ifail(:)
    real(rk),allocatable        ::  gap(:)
    !
    character(len=cl)    :: unitfname
    character(len=4)    :: str_row,str_col
    !
    !call setrteopts("ufmt_littleendian=10") ! on-the-fly Intel -> IBM binaries conversion
    !
    ! -------------------- !
    ! setup processor grid !
    ! -------------------- !
    !
    call blacs_pinfo(iam, nprocs)
    !
    if ( verbose>=4 ) write(out,"('iam = ',i4)") iam
    !
    do i=1,int( sqrt( dble(nprocs) ) + 1 )
      if(mod(nprocs,i) .eq. 0) nprow = i
    end do
    npcol = nprocs/nprow

    if (iam == 0) then
        write(out, "('BLACS topology: ',i8,' x ',i8,' (',i8,' processors over ',i8,' involved)')") nprow, npcol,nprow*npcol,nprocs
    endif
    !
    call blacs_get( -1, 0, context )
    call blacs_gridinit( context, 'r', nprow, npcol )
    call blacs_gridinfo( context, nprow, npcol, myrow, mycol )
    !
#if defined(__ELPA)
    call get_elpa_row_col_comms(MPI_COMM_WORLD, myrow, mycol, mpi_comm_rows, mpi_comm_cols)
#endif

    if (verbose>=4) then
      !
      do i=0,nprocs-1
        call blacs_barrier(context, 'a')
        if (iam .eq. i) then
          write(out,"(/'PE = ',i4,':',i4,' Grid-coord (',i4,',',i4,') PROW = ',i4,':',i4,' PCOL = ',i4,':',i4)") iam,nprocs,myrow,mycol,myrow,nprow,mycol,npcol
        endif
      enddo
    endif
    !
    call blacs_barrier(context, 'a')
    !
    ! --------------- !
    ! read input file !
    ! --------------- !
    !
    if (iam==0) then 
      !
      call FLReadInput(Jrot,gamma,gfactor,nroots,tol,sparse,eigensolver,chkpoint,zpe,memory,energy_thresh,coef_thresh)
      !
      if (verbose>=4) write (out,"(' iam = ',i8,' gfactor  = ',f20.8)") iam,gfactor
      !
      if (verbose>=4) write(out,"(' Only for iam = ',i,' jrot,gamma,nroots,chkpoint ',5i8)"), iam,jrot,gamma,nroots,chkpoint
      !
    endif
    !
    !dec$ if (blacs_ > 0)
      !
      call blacs_barrier(context,'A')
      !
      if (iam==0.and.verbose>=5) write(out,"(' destributing parameters...')") 
      !
      if (verbose>=6) write(out,"(' iam = ',i4)") iam
      !
      if ( (myrow.eq.0) .and. (mycol.eq.0) ) then
         call igebs2d(context, 'all', 'i-ring', 1, 1, eigensolver, 1 )
         call igebs2d(context, 'all', 'i-ring', 1, 1, gen_mat, 1 )
         call igebs2d(context, 'all', 'i-ring', 1, 1, mat_len, 1 )
         call igebs2d(context, 'all', 'i-ring', 1, 1, nroots, 1 )
         call igebs2d(context, 'all', 'i-ring', 1, 1, sparse, 1 )
         call igebs2d(context, 'all', 'i-ring', 1, 1, jrot , 1 )
         call igebs2d(context, 'all', 'i-ring', 1, 1, gamma, 1 )
         call dgebs2d(context, 'all', 'i-ring', 1, 1, walltime, 1 )
         call dgebs2d(context, 'all', 'i-ring', 1, 1, zpe, 1 )
         call dgebs2d(context, 'all', 'i-ring', 1, 1, energy_thresh, 1 )
         call dgebs2d(context, 'all', 'i-ring', 1, 1, coef_thresh, 1 )
         call dgebs2d(context, 'all', 'i-ring', 1, 1, tol, 1 )
         call dgebs2d(context, 'all', 'i-ring', 1, 1, gfactor, 1 )
         call dgebs2d(context, 'all', 'i-ring', 1, 1, memory, 1 )
         
    

      else
         call igebr2d(context, 'all', 'i-ring', 1, 1, eigensolver, 1, 0, 0 )
         call igebr2d(context, 'all', 'i-ring', 1, 1, gen_mat, 1 ,0 ,0)
         call igebr2d(context, 'all', 'i-ring', 1, 1, mat_len, 1, 0, 0 )
         call igebr2d(context, 'all', 'i-ring', 1, 1, nroots, 1, 0, 0 )
         call igebr2d(context, 'all', 'i-ring', 1, 1, sparse , 1, 0, 0 )
         call igebr2d(context, 'all', 'i-ring', 1, 1, jrot , 1, 0, 0 )
         call igebr2d(context, 'all', 'i-ring', 1, 1, gamma, 1, 0, 0 )
         call dgebr2d(context, 'all', 'i-ring', 1, 1, walltime, 1, 0, 0 )
         call dgebr2d(context, 'all', 'i-ring', 1, 1, zpe, 1, 0, 0 )
         call dgebr2d(context, 'all', 'i-ring', 1, 1, energy_thresh, 1, 0, 0 )
         call dgebr2d(context, 'all', 'i-ring', 1, 1, coef_thresh, 1, 0, 0 )
         call dgebr2d(context, 'all', 'i-ring', 1, 1, tol, 1, 0, 0 )
         call dgebr2d(context, 'all', 'i-ring', 1, 1, gfactor, 1, 0, 0 )
         call dgebr2d(context, 'all', 'i-ring', 1, 1, memory, 1, 0, 0 )
         

      endif
      !
      if (iam==0.and. verbose>=5) write(out,"(' ...done!')")
      !
      !if (sparse==1) blacs_init = .true.
      !
    !dec$ end if

    if(gen_mat == .true.) then
	    
	if (iam == 0) then
            write(out, "('Generating matrix of size        : ',i8)") mat_len
            write(out, "('Number of eigenstates to compute : ',i8)") nroots
        endif

        t1 = MPI_Wtime()
        zpe = 0.0
	    dimen_s = mat_len
	    chkptIO = 10+iam
	    
	    !
	    ! --------------------- !
	    ! define the block size !
	    ! --------------------- !
	    !
	    nb = min ( dimen_s/nprow, dimen_s/npcol )
	    nb = min ( nb, 64 )
	    nb = max(nb,1)    	

#if defined(__DEBUG)
            if (iam == 0) then
               write(out, "('NB : ',i3)") nb
            endif
#endif
	    
            loc_r = numroc(dimen_s,nb,myrow,0,nprow)
	    loc_c = numroc(dimen_s,nb,mycol,0,npcol)
	    lda = max (1,loc_r)
	    
	    call descinit( desca, dimen_s, dimen_s, nb, nb, 0, 0, context, lda, info)
   	    call descinit( descz, dimen_s, dimen_s, nb, nb, 0, 0, context, lda, info)
	    
#if defined(__DEBUG)
            write(out, "('I am',I4,', Local problem size: ',i6,' x ',i6)") iam, loc_r, loc_c
#endif

	    allocate(a_loc(loc_r,loc_c),z_loc(loc_r,loc_c),w(dimen_s),stat=info)
	    matsize = int(loc_r,kind=hik)*int(loc_c,kind=hik)

	    call ArrayStart(context,iam,'diag_scalapack:a_loc',info,size(a_loc),kind(a_loc),matsize)
	    call ArrayStart(context,iam,'diag_scalapack:z_loc',info,size(z_loc),kind(z_loc),matsize)   
	    call ArrayStart(context,iam,'diag_scalapack:w',info,size(w),kind(w))

            do i = 1, dimen_s
                seed(i) = 123456789-i-1;
            enddo
	   
            do j = 1,loc_c
                do i=1,loc_r
                    !
                    global_i = indxl2g( i, nb, myrow, 0, nprow )
                    global_j = indxl2g( j, nb, mycol, 0, npcol )
                    !

#if defined(__DEBUG)
                    call infog2l(global_i, global_j, desca, nprow, npcol, myrow, mycol, i_loc, j_loc, proc_row, proc_col)
                    if (( j_loc .ne. j ) .OR. (i_loc .ne. i )) then
                        write(out, "('(myrow: ',i4,', mycol: ',i4,') First set:',i4,' x ',i4,', Second set:',i4,' x ',i4)") myrow, mycol, global_i, global_j, i, j
                    endif
#endif

                    if(global_i >= global_j) a_loc(i,j) = real(rrr(global_i),rk)/real(1099511627776.0,rk)
                    if(global_i <  global_j) a_loc(i,j) = real(rrr(global_j),rk)/real(1099511627776.0,rk)
                    if(global_i == global_j) a_loc(i,j) = a_loc(i,j) + real(10.0,rk) + real(global_j,rk)

                enddo
            enddo
        !
      	call blacs_barrier(context, 'a')
        !
        t2 = MPI_Wtime()
        if (iam == 0) then
            write(out,'(/a,f12.6,a)') 'Time to Initialize the input matrix is',t2-t1,' sec'
        endif
        
     else    
	    	  
	    ! There is no reading from file...
        call blacs_abort( context, -1 )

    endif
    !
    ! -------------------------- !
    ! define remaining workspace !
    ! -------------------------- !
    !
    ! I will suggest to call the solvers and probe the sugegsted (hopefully minimum) space required
    !   by auxiliary *WORK data structures (NdFilipo)
    !
    select case (eigensolver)
      !
    case (1)
      !
      trilwmin = 3*dimen_s + max( nb*( loc_r+1 ), 3*nb)
      !
      lwork = max( 1 + 6*dimen_s + 2*loc_r*loc_c, trilwmin) + 2*dimen_s
      liwork = 7*dimen_s + 8*npcol + 2
      !
      allocate(work(lwork), iwork(liwork), stat=info)
      call ArrayStart(context,iam,'diag_scalapack:work',info,size(work),kind(work))
      call ArrayStart(context,iam,'diag_scalapack:iwork',info,size(iwork),kind(iwork))
      !
      if (iam == 0.and.verbose>=4) then
        write(out,"(/'lwork =  ', i16, ' liwork = ', i16, ' lda = ', i16, ' loc_r = ', i16, ' loc_c = ', i16)") lwork, liwork, lda, loc_r, loc_c
      endif
      !
      if (verbose>=4) call MemoryReport(context,iam,memory_now,memory_max)
      !
      ! -------------------------------------------------------------------------------------- !
      ! do any tests or evaluations here, before calling diagonalization subroutine, if needed !
      ! -------------------------------------------------------------------------------------- !
      !
      !lwork = -1
      !liwork = -1
      !call pdsyevd('V', 'U', dimen_s, a_loc, 1, 1, desca, w, z_loc, 1, 1, descz, work, lwork, iwork, liwork, info)
      !write(out,*) 'computed optimal lwork =  ', work(1), 'liwork = ', iwork(1)
      !call blacs_gridexit(context)
      !call blacs_exit(0)
      !stop
      !
      ! ---------------------------------------------- !
      ! call pdsyevd to compute the eigendecomposition !
      ! ---------------------------------------------- !
      !
      call blacs_barrier(context, 'a')
      t1 = MPI_Wtime()
      !
      !lwork = 10 ; liwork = 10
      !
#if defined(__ELPA)
      !
      if (iam == 0) then
#if defined(__2STAGE)
        write(out,"(/'Starting ELPA solve_evp_real_2stage...')")
#else
        write(out,"(/'Starting ELPA solve_evp_real...')")
#endif
      endif

      ! How many eigenvalues/eigenvectors ?
      ! nsolv=dimen_s == requests all eigenstates
      nsolv=nroots      
      !
#if defined(__2STAGE)
      call solve_evp_real_2stage(dimen_s, nsolv, a_loc, lda, w, z_loc, lda, nb, mpi_comm_rows, mpi_comm_cols, mpi_comm_world)
#else
      call solve_evp_real(dimen_s, nsolv, a_loc, lda, w, z_loc, lda, nb,mpi_comm_rows, mpi_comm_cols)
#endif
      info = 0
      !
#else
      !
      if (iam == 0) then
        write(out,"(/'Starting pdsyevd...')")
      endif

      call pdsyevd('V', 'L', dimen_s, a_loc, 1, 1, desca, w, z_loc, 1, 1, descz, work, lwork, iwork, liwork, info)
      !
#endif
      !
      nvals = dimen_s
      !
    case (2)
      !
      range = 'I'
      !
      VL = 0
      VU = 1e6
      !
      IL = 1
      IU = nroots
      !
      nvals  = dimen_s
      nvects = dimen_s
      lwork = -1 ; liwork = -1
      !
      abstol = 2.0*PDLAMCH(context,'U')
      orfac = -1e-4
      !
      allocate(iclustr(2*nprow*npcol),gap(nprow*npcol),ifail(dimen_s), stat=info)
      call ArrayStart(context,iam,'diag_scalapack:iclustr',info,size(iclustr),kind(iclustr))
      call ArrayStart(context,iam,'diag_scalapack:gap',info,size(gap),kind(gap))
      call ArrayStart(context,iam,'diag_scalapack:ifail',info,size(ifail),kind(ifail))
      !
      call pdsyevx('V', range, 'L', dimen_s, a_loc, 1, 1, desca, vl, vu, il,iu, abstol, nvals, nvects, w, orfac, z_loc, 1, 1, descz, work_, & 
                   lwork, iwork_, liwork, ifail, iclustr, gap,  info)
      !
      lwork = work_(1)*1.3
      liwork = iwork_(1)*1.3
      !
      if (verbose>=4 ) write(out,"(/'iam = ',i4,' lwork =  ', i16, ' liwork = ', i16, ' lda = ', i16)") iam,lwork, liwork, lda
      !
      allocate(work(lwork), iwork(liwork), stat=info)
      call ArrayStart(context,iam,'diag_scalapack:work',info,size(work),kind(work))
      call ArrayStart(context,iam,'diag_scalapack:iwork',info,size(iwork),kind(iwork))
      !
#if defined(__DEBUG)
      if (iam == 0) then
        write(out,"(/'lwork =  ', i16, 'liwork = ', i16, 'lda = ', i16, 'loc_r = ', i16, 'loc_c = ', i16)") lwork, liwork, lda, loc_r, loc_c
      endif
#endif
      !
      
      if (iam == 0) then
        write(out,"(/'Starting pdsyevx...')")
      endif

      call pdsyevx('V', range, 'L', dimen_s, a_loc, 1, 1, desca, vl, vu, il,iu, abstol, nvals, nvects, w, orfac, z_loc, 1, 1, descz, & 
                   work, lwork, iwork, liwork, ifail, iclustr, gap,  info)
      !
      if (iam == 0) then
        write(out,"(/'pdsyevx: nvals =  ', i16, 'nvetcs = ', i16)") nvals,nvects
      endif
      !
    case default
      !
      write(out,'("Uknown eigensolver = ",i)'), eigensolver
      stop "Uknown eigensolver = "
      !
    end select
    !
    call blacs_barrier(context, 'a')
    t2 = MPI_Wtime()
    !
    ! ------------------------------ !
    ! do tests on convergence if any !
    ! ------------------------------ !
    !
!    if (iam == 0 .and.verbose>=3) then
    if (iam == 0) then
!      write(out,"(/'info = ', i8)") info
      if (info .eq. 0) then
        write(out,"(/'Diagonalization finished successfully!')")
        write(out,'(/a,f12.6,a)') 'Time to diagonalize matrix is ',t2-t1,' sec'

        ! We can mimic the same timing granularity in ScaLAPACK too by patching the netlib code [NdFilippo]

#if defined(__ELPA)
        print *,'Time tridiag_real     :',time_evp_fwd
        print *,'Time solve_tridi      :',time_evp_solve
        print *,'Time trans_ev_real    :',time_evp_back
        print *,'Total time (sum above):',time_evp_back+time_evp_solve+time_evp_fwd
#endif

      else if (info .lt. 0) then
        write(out,"(/'Info is less than zero. Info is equal to ', i8)") info
      else
        write(out,"(/'Info is larger than zero. Info is equal to ', i8)") info
      endif
    endif
    !
    call blacs_barrier(context, 'a')
    t2 = MPI_Wtime()
    !   
    !if (iam == 0.and.verbose>=4) then
    !  write(out,'(/a,f12.6,a)') 'Time to print vectors with PDLAPRNT is ',t2-t1,' sec'
    !endif
    !
    deallocate(a_loc, z_loc, work, iwork, w)
    call ArrayStop(context,'diag_scalapack:work')
    call ArrayStop(context,'diag_scalapack:iwork')
    call ArrayStop(context,'diag_scalapack:a_loc')
    call ArrayStop(context,'diag_scalapack:z_loc')
    call ArrayStop(context,'diag_scalapack:w')
    !
    ! Exit BLACS
    !
    !if (verbose>=3) call MemoryReport(context,iam,memory_now,memory_max)
    if (iam == 0) then
       call MemoryReport(context,iam,memory_now,memory_max)
    endif
    !

    call blacs_gridexit(context)
    call blacs_exit(0)


    !
end program dirac_exomol_eigen

