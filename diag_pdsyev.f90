module d_module
 
!dec$ define blacs_  = 1
!dec$ define mpi_    = 0


  use accuracy
  use timer
  implicit none

  public FLReadInput

  !
  integer(ik),parameter:: verbose = 4
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
  integer(hik)	::	seed(1000000)
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
    write(out,"('Read the input')")
    !
    sparse = 0
    !
    jrot = -1
    gamma = -1
    nroots = 1e6
    energy_thresh = 1.0d+05
    coef_thresh = 1.0d-16
    diagonalizer = "PDSYEVD"
    factor = 1.0_rk
    tol = small_
    sparse = .false.
    walltime = safe_max
    zpe = -small_
    memory = 256.0_rk
    !
    call input_options(echo_lines=.true.,error_flag=1)
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
          case default
            !
            write(out,'("Uknown eigensolver = ",a)'), diagonalizer
            stop "Uknown eigensolver = "
            !
          end select
          !
        case default
          call report ("Principal keyword "//trim(w)//" not recognized",.true.)
        end select
        !
    end do
    !
    write(out,"('...done!')")
    !
  end subroutine FLReadInput
  
  integer(hik) function  rrr(i)
  	integer(ik)	::	i
  	seed(i) = mod((1103515245 * seed(i) + 12345),1099511627776);
  	rrr=seed(i)
  	return
  end function rrr
  
  

end module d_module


program diag_pdsyev
  
  use accuracy
  use d_module
  use timer
  implicit none

  !integer(ik) :: verbose=6

  !include 'mpif.h'

    ! Related to matrix reading
    integer(ik) ::              jrot, gamma, chkptIO, info, i, j, i_loc, j_loc, iroot, idimen, jdimen
    character(cl) ::            jchar, symchar, filename, buf
    integer(hik)                matsize
    integer(ik)                 dimen_s, nelem, max_nelem
    integer(ik), allocatable :: bterm(:,:)
    real(rk), allocatable ::    a_temp(:), a_loc(:,:)   
    real(rk) ::                 thresh=30000.0_rk, zpe, energy_thresh, coef_thresh
    ! Related to diagonalization
    integer(ik)  ::     nroots
    ! Parameters
    integer(ik)         lwork, liwork
    double precision    zero, t1, t2, work_(10)
    parameter           ( zero = 0.0d0 )
    integer(ik)         loc_r, loc_c, lda
    double precision    mone
    integer(ik)         maxprocs
    parameter           ( mone = -1.0d0, maxprocs = 512 )
    ! Local Scalars
    integer(ik)         context, iam, m, mycol, myrow, nb, npcol, nprocs, nprow, nz, trilwmin, proc_row, proc_col,iwork_(10),l_nrows, l_ncols, l_eigvec, iprow, ipcol
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
    integer(ik), external :: numroc, iceil
    real(rk), external ::    MPI_Wtime, pdlamch
    !real(rk), external ::    pdlamch
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
    !
    call blacs_get( -1, 0, context )
    call blacs_gridinit( context, 'r', nprow, npcol )
    call blacs_gridinfo( context, nprow, npcol, myrow, mycol )
    !
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
      if (iam==0.and.verbose>=5) write(out,"(' ...done!')")
      !
      !if (sparse==1) blacs_init = .true.
      !
    !dec$ end if

    if(gen_mat == .true.) then
	    if ( verbose>=4 ) write(out, "('Generating matrix of size ',i4)") mat_len
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
	    
			!
	    loc_r = numroc(dimen_s,nb,myrow,0,nprow)
	    loc_c = numroc(dimen_s,nb,mycol,0,npcol)
	    lda = max (1,loc_r)
	    
	    call descinit( desca, dimen_s, dimen_s, nb, nb, 0, 0, context, lda, info)
   	    call descinit( descz, dimen_s, dimen_s, nb, nb, 0, 0, context, lda, info)
	    
	    !
	    allocate(a_loc(loc_r,loc_c),z_loc(loc_r,loc_c),w(dimen_s),stat=info)
	    matsize = int(loc_r,kind=hik)*int(loc_c,kind=hik)
	    call ArrayStart(context,iam,'diag_scalapack:a_loc',info,size(a_loc),kind(a_loc),matsize)
	    call ArrayStart(context,iam,'diag_scalapack:z_loc',info,size(z_loc),kind(z_loc),matsize)   
	    call ArrayStart(context,iam,'diag_scalapack:w',info,size(w),kind(w))
	    
	    
	    do i = 1, dimen_s
	    	seed(i) = 123456789-i-1;
	    enddo
	    
	    
    	    do j = 1,dimen_s
    	    	do i=1,dimen_s
          	!
          	call infog2l(i, j, desca, nprow, npcol, myrow, mycol, i_loc, j_loc, proc_row, proc_col)
          	!
           	if (myrow == proc_row .and. mycol == proc_col) then
           		if(i>=j) a_loc(i_loc,j_loc) = real(rrr(i),rk)/real(1099511627776.0,rk)
           		if(i < j) a_loc(i_loc,j_loc) = real(rrr(j),rk)/real(1099511627776.0,rk)
           		if(i==j) a_loc(i_loc,j_loc) = a_loc(i_loc,j_loc) + real(10.0,rk) + real(j,rk)
           	endif
          !
          !if (j == i .and. a_temp(j) <= thresh) nroots = nroots + 1
          !
      	  enddo
      	  enddo
      	  
      	  call blacs_barrier(context, 'a')
        !
        
        
     else    
	    	  

	    !
	    ! ---------------------------------------------------- !
	    ! read headings and matrix dimensions from matrix file !
	    ! ---------------------------------------------------- !

	    if ( verbose>=4 ) write(out, '("j,gamma = ",2i4)') jrot,gamma

	    !
	    write(jchar, '(i4)') jrot
	    write(symchar, '(i4)') gamma
	    !filename = 'bterm'//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))
	    !


	    filename = 'matrix'//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'.chk'

	    if ( verbose>=4 ) write(out, '(a)') filename

	    !
	    chkptIO = 10+iam
	    !
	    open(chkptIO,form='unformatted',action='read',position='rewind',file=filename)
	    !
	    !read(chkptIO) buf(1:11)
	    !if (buf(1:11)/='Start bterm') then
	    !  write(out, '(/a)') 'error: buf(1:11)/=Start bterm'
	    !  stop
	    !endif
	    !
	    rewind(chkptIO)
	    !
	    if (sparse) then 
	      !
	      read(chkptIO) buf(1:11)
	      if (buf(1:11)/='Start bterm') then
		write(out, '(/a)') 'error: buf(1:11)/= Start bterm'
		stop
	      endif
	      !
	    endif
	    !
	    read(chkptIO) dimen_s
	    !
	    allocate(bterm(dimen_s,2),w(dimen_s),a_temp(dimen_s),stat=info)
	    call ArrayStart(context,iam,'diag_scalapack:bterm',info,size(bterm),kind(bterm))
	    call ArrayStart(context,iam,'diag_scalapack:w',info,size(w),kind(w))
	    call ArrayStart(context,iam,'diag_scalapack:a_temp',info,size(a_temp),kind(a_temp))
	    !
	    if (verbose>=3) call MemoryReport(context,iam,memory_now,memory_max)
	    !
	    !dec$ if (blacs_ > 0)
	      call dgsum2d( context,'A', 'i-ring',1,1,memory_now,1,0,0)
	      call dgsum2d( context,'A', 'i-ring',1,1,memory_max,1,0,0)
	      !
	      if (verbose>=3.and.iam==0) write (out,"(t2,'Total memory   = ',t47,f18.8,' Gb')") memory_now
	      if (verbose>=3.and.iam==0) write (out,"(t2,'Maximal memory = ',t47,f18.8,' Gb (',f16.1,')')") memory_max,memory_limit
	      !
	    !dec$ end if
	    !
	    if (sparse) then 
	      !
	      read(chkptIO) bterm(1:dimen_s,1:2)
	      !
	      !close(chkptIO)
	      !
	      !call blacs_barrier(context, 'a')
	      !
	      max_nelem = 0
	      !
	      do i=1,dimen_s
		nelem = bterm(i,2) - bterm(i,1) + 1
		!if (iam == 0) then
		!  write(out,'(/a,i4,a,i4,a,i4)') 'row: ', i, ' bterm(1): ', bterm(i,1), ' bterm(2): ', bterm(i,2)
		!endif
		if (nelem .gt. max_nelem) then
		  max_nelem = nelem
		endif
	      enddo
	      !
	    else
	      !
	      max_nelem = dimen_s
	      !
	    endif
	    !
	    if (iam == 0.and.verbose>=4) then
	      write(out,"(/'The matrix of the size ', i16, ' will be read. The maximum number of non-zero elements in a row is ', i16)") dimen_s, max_nelem 
	    endif
	    !
	    ! --------------------- !
	    ! define the block size !
	    ! --------------------- !
	    !
	    nb = min ( dimen_s/nprow, dimen_s/npcol )
	    nb = min ( nb, 64 )
	    nb = max(nb,1)
	    !
	    ! -------------------------------------------------- !
	    ! define local array size and initialize descriptors !
	    ! -------------------------------------------------- !
	    !
	    loc_r = numroc(dimen_s,nb,myrow,0,nprow)
	    loc_c = numroc(dimen_s,nb,mycol,0,npcol)
	    lda = max (1,loc_r)
	    !
	    allocate(a_loc(loc_r,loc_c),z_loc(loc_r,loc_c),stat=info)
	    matsize = int(loc_r,kind=hik)*int(loc_c,kind=hik)
	    call ArrayStart(context,iam,'diag_scalapack:a_loc',info,size(a_loc),kind(a_loc),matsize)
	    call ArrayStart(context,iam,'diag_scalapack:z_loc',info,size(z_loc),kind(z_loc),matsize) 
	    !
	    !
	    call descinit( desca, dimen_s, dimen_s, nb, nb, 0, 0, context, lda, info)
	    call descinit( descz, dimen_s, dimen_s, nb, nb, 0, 0, context, lda, info)
	    !
	    ! ---------------------------------------------- !
	    ! read the input matrix in a distributed fassion !
	    ! ---------------------------------------------- !
	    !
	    !filename = 'hmat'//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))
	    !
	    !open(chkptIO,form='unformatted',action='read',position='rewind',file=filename)
	    !
	    !call blacs_barrier(context, 'a')
	    t1 = MPI_Wtime()
	    !
	    !nroots = 0
	    !
	    read(chkptIO) buf(1:12)
	    if (buf(1:12)/='Start matrix') then
	      write(out, '(/a)') 'error: buf(1:12)/= Start matrix'
	      stop
	    endif
	    !
	    do i=1,dimen_s
	      !
	      if (sparse) then 
		!
		read(chkptIO) a_temp(bterm(i,1):bterm(i,2))
		!
		if (bterm(i,1) > 1) then
		  !
		  do j=1,bterm(i,1)-1
		    !
		    call infog2l(i, j, desca, nprow, npcol, myrow, mycol, i_loc, j_loc, proc_row, proc_col)
		    !
		    if (myrow == proc_row .and. mycol == proc_col) a_loc(i_loc,j_loc) = 0.0d0
		    !
		  enddo
		  !
		endif
		!
		do j=bterm(i,1),bterm(i,2)
		  !
		  call infog2l(i, j, desca, nprow, npcol, myrow, mycol, i_loc, j_loc, proc_row, proc_col)
		  !
		  if (j <= i) then
		    !
		    if (myrow == proc_row .and. mycol == proc_col) a_loc(i_loc,j_loc) = a_temp(j)
		    !
		  endif
		  !
		  !if (j == i .and. a_temp(j) <= thresh) nroots = nroots + 1
		  !
		enddo
		!
	      else
		!
		read(chkptIO) a_temp
		!
		do j=1,dimen_s
		  !
		  call infog2l(i, j, desca, nprow, npcol, myrow, mycol, i_loc, j_loc, proc_row, proc_col)
		  !
		  if (j <= i) then
		    !
		    if (myrow == proc_row .and. mycol == proc_col) a_loc(i_loc,j_loc) = a_temp(j)
		    !
		  endif
		  !
		  !if (j == i .and. a_temp(j) <= thresh) nroots = nroots + 1
		  !
		enddo
		!
	      endif
	      !
	    enddo
	    !
	    read(chkptIO) buf(1:10)
	    if (buf(1:10)/='End matrix') then
	      write(out, '(/a)') 'error: buf(1:10)/= End matrix'
	      stop
	    endif
	    !
	    call blacs_barrier(context, 'a')
	    t2 = MPI_Wtime()
	    !
	    if (iam == 0.and.verbose>=4) then
	      write(out,'(/a,f12.6,a)') 'Time to read and distribute matrix is ',t2-t1,' sec'
	    endif
	    !
	    close(chkptIO)
	    !
	    deallocate(a_temp, bterm)
	    call ArrayStop(context,'diag_scalapack:a_temp')
	    call ArrayStop(context,'diag_scalapack:bterm')
    endif
    !
    ! -------------------------- !
    ! define remaining workspace !
    ! -------------------------- !
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
      call pdsyevd('V', 'L', dimen_s, a_loc, 1, 1, desca, w, z_loc, 1, 1, descz, work, lwork, iwork, liwork, info)
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
      !if (iam == 0) then
      !  write(out,"(/'lwork =  ', i16, 'liwork = ', i16, 'lda = ', i16, 'loc_r = ', i16, 'loc_c = ', i16)") lwork, liwork, lda, loc_r, loc_c
      !endif
      !
      call pdsyevx('V', range, 'L', dimen_s, a_loc, 1, 1, desca, vl, vu, il,iu, abstol, nvals, nvects, w, orfac, z_loc, 1, 1, descz, & 
                   work, lwork, iwork, liwork, ifail, iclustr, gap,  info)
      !
      if (iam == 0.and.verbose>=4) then
        write(out,"(/'nvals =  ', i16, 'nvetcs = ', i16)") nvals,nvects
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
    if (iam == 0.and.verbose>=3) then
      write(out,"(/'info = ', i8)") info
      if (info .eq. 0) then
        write(out,"(/'Diagonalization finished successfully!')")
        write(out,'(/a,f12.6,a)') 'Time to diagonalize matrix is ',t2-t1,' sec'
      else if (info .lt. 0) then
        write(out,"(/'Info is less than zero. Info is equal to ', i8)") info
      else
        write(out,"(/'Info is larger than zero. Info is equal to ', i8)") info
      endif
    endif
    !
    ! ------------------------------------------- !
    ! Print or store eigenvalues and eigenvectors !
    ! ------------------------------------------- !
    !
    if (iam == 0 .and. info == 0) then
      !
      write(out,"(/'Computed eigenvalues:')")
      !
      do i=1,dimen_s
          !
          !nroots = nroots + 1
          !
          write(out,'(f18.8)') w(i)-zpe
          !
      enddo
      !
    endif
    !
    ! -------------------------------- !
    ! Print eigenvectors with PDLAPRNT !
    ! -------------------------------- !
    !
    
    if(gen_mat==.false.) then
	    
	    call blacs_barrier(context, 'a')
	    t1 = MPI_Wtime()
	    !
	    !
	    !open(12,form='unformatted',action='write',file='eigenvectors')
	    !
	    write(unitfname,"('eigensolution for j = ',i6,' sym = ',i4)") jrot,gamma
	    !
	    call IOStart(trim(unitfname),chkptIO)
	    !
	    !
	    write(str_row,'(i4)') myrow
	    write(str_col,'(i4)') mycol
	    !
	    !
	    write(jchar, '(i4)') jrot
	    write(symchar, '(i4)') gamma
	    filename = 'solution_'//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'_'//trim(adjustl(str_row))//'_'//trim(adjustl(str_col))//'.chk'


	    !--------------------------------------------------------------------------
	    !              PRINT SOLUTION
	    !
	    call blacs_barrier(context, 'a')
	    !
	    filename = 'solution_'//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'_'//trim(adjustl(str_row))//'_'//trim(adjustl(str_col))//'.chk'
	    !
	    write(out, '(/1x,a,1x,a,1x,a,1x,i3)') 'write eigenvectors into file=', trim(filename), 'iounit=', chkptIO
	    !
	    open(chkptIO, form='unformatted', action='write', status='unknown', position='rewind', file=filename, buffered='yes')
	    !
	    if (info/=0) then
	      write(out, '(/a,1x,a)') 'error while opening file=', trim(filename)
	      stop
	    endif
	    !
	    nvals = 0
	    do idimen=1, dimen_s
	      if (abs(w(idimen)-zpe)<=energy_thresh) then
		nvals = nvals + 1
	      else
		exit
	      endif
	    enddo
	    !
	    write(chkptIO) nprow, npcol
	    write(chkptIO) 'Start energies'
	    write(chkptIO) dimen_s, nvals
	    write(chkptIO) w(1:nvals)
	    write(chkptIO) 'End energies'
	    !
	    write(chkptIO) loc_r,loc_c
	    write(chkptIO) 'Start vectors'
	    !
	    do idimen=1, dimen_s
	      !
	      do jdimen=1, nvals
		!
		call infog2l(idimen, jdimen, descz, nprow, npcol, myrow, mycol, i, j, iprow, ipcol)
		!
		if (myrow==iprow .and. mycol==ipcol) then
		  !
		  if (abs(z_loc(i,j))>=coef_thresh) write(chkptIO) idimen, jdimen, z_loc(i,j)
		  !
		endif
		!
	      enddo
	      !
	    enddo
	    !
	    write(chkptIO) 'End vectors'
	    !
	    call blacs_barrier(context, 'a')
	    !
	    t2 = MPI_Wtime()
	    !
	    if (iam == 0.and.verbose>=4) then
	      write(out,'(/a,f12.6,a)') 'Time to write vectors: ',t2-t1,' sec'
	    endif
	    !
	    close(chkptIO)
	    !
	    write(out, '(1x,a,1x,i3)') 'done for iounit=', chkptIO
	    !
	    call blacs_barrier(context, 'a')
	    !
	    !--------------------------------------------------------------------------



	    !
	    !write(chkptIO,*) dimen_s, nprow, npcol, myrow, mycol, loc_r, loc_c
	    !
	    !!write(chkptIO,*) myrow, mycol, loc_r, loc_r,nroots
	    !
	    !write(chkptIO,*) 'Start vectors'
	    !write(chkptIO,*)  z_loc
	    !write(chkptIO,*) 'End vectors'

	    !
	    if (myrow==0.and.mycol==0) then 
	      !
	      filename = 'solution_'//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'.tmp'
	      open(chkptIO,  action='write', status='unknown', position='rewind', file=filename, iostat=info)
	      ! 
	      write(chkptIO,*) 'Start energies'
	      write(chkptIO,*) dimen_s,nvals
	      write(chkptIO,*) w(1:nvals)
	      write(chkptIO,*) 'End energies'
	      !
	      close(chkptIO)
	      !
	    endif

	    !
	    call blacs_barrier(context, 'a')
	    !
	    !filename = 'solution_'//trim(adjustl(jchar))//'_'//trim(adjustl(symchar))//'_'//trim(adjustl(str_row))//'_'//trim(adjustl(str_col))//'.tmp'

	    !write(out, '(/1x,a,1x,a,1x,a,1x,i3)') 'write eigenvectors into file=', trim(filename), 'iounit=', chkptIO
	    !!
	    !open(chkptIO,  action='write', status='unknown', position='rewind', file=filename, iostat=info)
	    !write(chkptIO,*) myrow, mycol, loc_r, loc_r,nroots
	    !write(chkptIO,*) 'Start vectors'
	    !write(chkptIO,*)  z_loc(:,1:nroots)
	    !write(chkptIO,*) 'End vectors'

	    !
	    !if (myrow==0.and.mycol==0) then 
	    !  ! 
	    !  !open(chkptIO,form='unformatted',action='write',position='rewind',status='replace',file=filename)
	    !  !
	    !  write(chkptIO,*) 'Start energies'
	    !  write(chkptIO,*) dimen_s,nvals
	    !  write(chkptIO,*) w(1:nvals)
	    !  write(chkptIO,*) 'End energies'
	    !  !
	    !  !write(chkptIO) 'Start vectors'
	    !  !
	    !endif
	    !
	    !close(chkptIO)


	    !call blacs_barrier(context, 'a')
	    !
	    !do j=1,nvals
	    !   !
	    !   do i=1,dimen_s
	    !     !
	    !     call infog2l(i, j, descz, nprow, npcol, myrow, mycol, i_loc, j_loc, proc_row, proc_col)
	    !       !
	    !       if (myrow == proc_row .and. mycol == proc_col) write(2000,"(2i8,2x,g26.18)") i,j,z_loc(i_loc,j_loc)
	    !       !
	    !   enddo
	    !   !
	    !   !
	    !enddo



	    !
	    !call pdlaprnt(dimen_s, dimen_s, z_loc, 1, 1, descz, 0, 0, 'A', chkptIO, work)
	    !
	    !call pdlaprnt_local(dimen_s, dimen_s, z_loc, 1, 1, descz, 0, 0, 'A', 12, work)
	    !
	    !close(12)
	    !


	    !nvals,nvects
	    !dec$ if (blacs_ > 0)
	    if ( iam==0 ) then
	       call igebs2d(context, 'all', 'i-ring', 1, 1,nvals, 1 )
	       call igebs2d(context, 'all', 'i-ring', 1, 1,nvects, 1 )
	    else
	       call igebr2d(context, 'all', 'i-ring', 1, 1,nvals_, 1, 0, 0 )
	       call igebr2d(context, 'all', 'i-ring', 1, 1,nvects_, 1, 0, 0 )
	       !
	       if (nvals/=nvals_.or.nvects/=nvects_) then 
		 !
		 write(out,"('nvals = ',i8,' or nvects from ',i4,' do not agree with values from iam = 0: ',i8)") nvals_,nvects_,iam,nvals,nvects
		 !
		 stop 'nvals do not agree'
		 !
	       endif
	       !
	    endif
	    !dec$ endif
	    !
    
    else
    	    t1 = MPI_Wtime()
    	    write(jchar, '(i4)') mat_len
    	    call blacs_barrier(context, 'a')
    	    
    	    allocate(local_vecs(1:dimen_s))
    	    call ArrayStart(context,iam,'diag_scalapack:local_vecs',info,size(local_vecs),kind(local_vecs))
    	    
	    !
	    filename = 'solution_gen-mat_N-'//trim(adjustl(jchar))//'_'//'.chk'
	    !
	    if(iam==0) write(out, '(/1x,a,1x,a,1x,a,1x,i3)') 'writing first 10 eigenvectors into file=', trim(filename), 'iounit=', chkptIO
	    !
	     if(iam==0) open(chkptIO, form='unformatted', action='write', status='unknown', position='rewind', file=filename, buffered='yes')
	    !
	    do idimen=1, dimen_s
	    	local_vecs = 0
	      !
	      do jdimen=1, min(dimen_s,10)
		!
		call infog2l(idimen, jdimen, descz, nprow, npcol, myrow, mycol, i, j, iprow, ipcol)
		!
		if (myrow==iprow .and. mycol==ipcol) then
		  !
		   local_vecs(idimen) = z_loc(i,j)
		  !
		endif
		!
	      enddo
	      !
	      
	      call dgsum2d(context, 'all', ' ', dimen_s, 1, local_vecs, -1, -1, 0 )
	      if(iam==0) write(chkptIO) local_vecs
	      
	      
	    enddo
	    
	    deallocate(local_vecs)
	    call ArrayStop(context,'diag_scalapack:local_vecs')
     
	    
	    
	    !
	    !
	    call blacs_barrier(context, 'a')
	    !
	    t2 = MPI_Wtime()
	    !
	    if (iam == 0.and.verbose>=4) then
	      write(out,'(/a,f12.6,a)') 'Time to write vectors: ',t2-t1,' sec'
	    endif
	    !
	     if(iam==0) close(chkptIO)
	    !
	    !write(out, '(1x,a,1x,i3)') 'done for iounit=', chkptIO
	    !
	    call blacs_barrier(context, 'a')
	    !
    	
     endif  
    
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
    if (verbose>=3) call MemoryReport(context,iam,memory_now,memory_max)
    !

    call blacs_gridexit(context)
    call blacs_exit(0)


    !
end program diag_pdsyev

