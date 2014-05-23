module d_module
    !
    use timer
    use accuracy
    !
    implicit none
    !
    integer(ik),parameter:: verbose = 2
    integer(ik) :: mat_len
    integer(hik), allocatable	::	seeded_array(:)
    character(len=cl)           :: diagonalizer, generator_str

contains
    ! Shifting pseudo-random function based on MOD
    integer(hik) function  rrr(i)
        integer(ik)	::	i
        seeded_array(i) = mod((1103515245 * seeded_array(i) + 12345),1099511627776);
        rrr=seeded_array(i)
        return
    end function rrr

end module d_module


program dirac_exomol_eigen
  
    use d_module
    use accuracy
    use timer

#if defined(__ELPA)
    use elpa1
    use elpa2
#endif

#if !defined(__IPM)
    use mpi
#endif

    implicit none

#if defined(__IPM)
    include 'mpif.h'
#endif

    !integer(ik) :: verbose=6

    ! Related to matrix reading
    integer(ik) ::              jrot, gamma, info, i, j, i_loc, j_loc, iroot, idimen, jdimen
    character(cl) ::            jchar, symchar, filename, buf
    integer(hik)                matsize
    integer(ik)                 dimen_s, nelem, max_nelem
    integer(ik), allocatable :: bterm(:,:)
    real(rk), allocatable ::    a_temp(:), a_loc(:,:), c_loc(:,:)
    real(rk) ::                 thresh=30000.0_rk, zpe, energy_thresh, coef_thresh,rrr_value
    ! Related to diagonalization
    integer(ik)  ::     nroots, nsolv
    ! Parameters
    integer(ik)         lwork, liwork
    double precision    zero, t1, t2, t3, t4, work_(10)
    parameter           ( zero = 0.0d0 )
    integer(ik)         loc_r, loc_c, lda
    double precision    mone
    integer(ik)         maxprocs,global_i,global_j
    parameter           ( mone = -1.0d0, maxprocs = 512 )
    ! Local Scalars
    integer(ik)         context, iam, m, mycol, myrow, nb, npcol, nprocs, nprow, nz, trilwmin, proc_row, proc_col,iwork_(10),l_nrows, l_ncols, l_eigvec, iprow, ipcol
    integer(ik)         mpi_comm_rows, mpi_comm_cols
    ! Local Arrays
    integer(ik)                      desca(9), descz(9), descc(9), my_seed(1), irecl
    integer(ik), allocatable ::      iwork(:)
    double precision, allocatable :: work(:), w(:), z_loc(:,:), eigvec(:),local_vecs(:)
    !
    double precision :: vl,vu
    integer(ik)      :: il,iu,nvals,nvects,nvalsmax,nvals_,nvects_,chkpoint,eigensolver,matrix_generator,dims(2),ndims
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
    INTEGER*4 :: iargc
    character*16 arg1, arg2, arg3, arg4
    !
    ! --------------- !
    ! read input file !
    ! --------------- !
    !
    mat_len = 10000                ! default
    nroots = 10000                 ! default
    diagonalizer  = "PDSYEVD"      ! default - name
    generator_str = "RANDOM-LOCAL" ! default - name
    !
    if (iargc() == 4) then
       call getarg(1, arg1)
       call getarg(2, arg2)
       call getarg(3, arg3)
       call getarg(4, arg4)
       read(arg1, *) mat_len
       read(arg2, *) nroots
       read(arg3, *) diagonalizer
       read(arg4, *) generator_str
    endif
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
            write(out,'("Unknown eigensolver = ",a)'), diagonalizer
            stop "Unknown eigensolver"
            !
    end select
    !
    select case ( trim(generator_str) )
        !
        case ('SYM-POSITIVE-O3')
            !
            matrix_generator = 1
            !
        case ('SYM-POSITIVE-O2')
            !
            matrix_generator = 2
            !
        case ('RANDOM-LOCAL')
            !
            matrix_generator = 3
            !
        case default
            !
            write(out,'("Unknown matrix generator = ",a)'), generator_str
            stop "Unknown generator"
            !
    end select
    !
    ! ------------------ !
    ! MPI initialization !
    ! ------------------ !
    !
    call MPI_INIT(info)
    !
    ! -------------------- !
    ! setup processor grid !
    ! -------------------- !
    !
    call blacs_pinfo(iam, nprocs)
    !
    ndims = 2
    dims  = 0
    CALL MPI_DIMS_CREATE( nprocs, ndims, dims, info)
#if 1
    nprow = dims(1)  !  cartesian direction 0
    npcol = dims(2)  !  cartesian direction 1
#else
    ! Reversed...
    nprow = dims(2)
    npcol = dims(1)
#endif
    !
    if (iam == 0) then
        write(out, "('BLACS topology: ',i8,' x ',i8,' (',i8,' processors over ',i8,' involved)')") nprow, npcol,nprow*npcol,nprocs
    endif
    !
    ! ---------------- !
    ! initialize BLACS !
    ! ---------------- !
    !
    call blacs_get( -1, 0, context )
    call blacs_gridinit( context, 'r', nprow, npcol )
    call blacs_gridinfo( context, nprow, npcol, myrow, mycol )
    !
#if defined(__ELPA)
    call get_elpa_row_col_comms(MPI_COMM_WORLD, myrow, mycol, mpi_comm_rows, mpi_comm_cols)
#endif

#if 0
    do i=0,nprocs-1
        call blacs_barrier(context, 'a')
        if (iam .eq. i) then
            write(out,"(/'PE = ',i4,':',i4,' Grid-coord (',i4,',',i4,') PROW = ',i4,':',i4,' PCOL = ',i4,':',i4)") iam,nprocs,myrow,mycol,myrow,nprow,mycol,npcol
        endif
    enddo
#endif
    !
    if (iam == 0) then
        write(out, "('Generating matrix of size        : ',i8)") mat_len
        write(out, "('Number of eigenstates to compute : ',i8)") nroots
    endif
    !
    t1 = MPI_Wtime()
    zpe = 0.0
    dimen_s = mat_len
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

    allocate(a_loc(loc_r,loc_c),stat=info)
    matsize = int(loc_r,kind=hik)*int(loc_c,kind=hik)

    call ArrayStart(context,iam,'diag_scalapack:a_loc',info,size(a_loc),kind(a_loc),matsize)

    call descinit( descc, dimen_s, dimen_s, nb, nb, 0, 0, context, lda, info)

    ! THEORY:
    ! Generate a dense n x n symmetric, positive definite matrix
    ! 1) A = rand(n,n); % generate a random n x n matrix using [0,1) values
    ! 2) A = A+A' ( O(n^2) complexity )
    !     _ or _
    !    A = A*A' ( O(n^3) complexity )
    ! 3) A = A + n*I ( A symmetric diagonally dominant matrix and symmetric positive definite)

    ! PRACTICE:

    ! Init random number generator... in a non entirely random way!
    i = 1
    call RANDOM_SEED(size = i)
    my_seed(1)=19830607
    call RANDOM_SEED(put=my_seed)
    !
    if (iam == 0) then
        write(out,"(/'Generating matrix using ',a,' ...')") trim(generator_str)
    endif
    !
    ! Select matrix generator
    select case (matrix_generator)
        !
        case (1)
            ! "SYM-POSITIVE-O3" - O(n^3) complexity
            !
            allocate(c_loc(loc_r,loc_c),stat=info)
            matsize = int(loc_r,kind=hik)*int(loc_c,kind=hik)
            call ArrayStart(context,iam,'aux:c_loc',info,size(c_loc),kind(c_loc),matsize)
            !
            if (iam == 0) then
                write(out,"(/'Fill randomly C...')")
            endif
            !
            do j = 1,loc_c
                do i=1,loc_r
                    !
                    call RANDOM_NUMBER (HARVEST=c_loc(i,j))
                    !
                    ! c = n*I
#if defined(__PDGEMM_C_TERM)
                    global_j = indxl2g( j, nb, mycol, 0, npcol )
                    global_i = indxl2g( i, nb, myrow, 0, nprow )
                    if(global_i == global_j) then
                        a_loc(i,j) = dimen_s*1.0d0
                    else
                        a_loc(i,j) = 0.0d0
                    endif
#endif
                enddo
            enddo
            !
            if (iam == 0) then
                write(out,"(/'Compute symmetric positive A (PDGEMM)... ')", advance='no')
            endif
            !
            call blacs_barrier(context, 'a')
            t3 = MPI_Wtime()
            !
#if defined(__PDGEMM_C_TERM)
            call PDGEMM('N', 'T', dimen_s, dimen_s, dimen_s, 1.0d0, c_loc, 1, 1, descc, c_loc, 1, 1, descc, 1.0d0, a_loc, 1, 1, desca)
#else
            call PDGEMM('N', 'T', dimen_s, dimen_s, dimen_s, 1.0d0, c_loc, 1, 1, descc, c_loc, 1, 1, descc, 0.0d0, a_loc, 1, 1, desca)
#endif
            !
            call blacs_barrier(context, 'a')
            t4 = MPI_Wtime()
            if (iam == 0) then
                write(out,'(f7.2,a)') t4-t3,' sec'
            endif
            !
            call ArrayStop(context,'aux:c_loc')
            deallocate (c_loc)
            !
        case (2)
            ! "SYM-POSITIVE-O2" - O(n^2) complexity
            !
            allocate(c_loc(loc_r,loc_c),stat=info)
            matsize = int(loc_r,kind=hik)*int(loc_c,kind=hik)
            call ArrayStart(context,iam,'aux:c_loc',info,size(c_loc),kind(c_loc),matsize)
            !
            if (iam == 0) then
                write(out,"(/'Fill randomly C...')")
            endif
            !
            do j = 1,loc_c
                do i=1,loc_r
                    !
                    call RANDOM_NUMBER (HARVEST=c_loc(i,j))
                enddo
            enddo
            !
            if (iam == 0) then
                write(out,"(/'Transpose C (PDTRAN)... ')", advance='no')
            endif
            !
            call blacs_barrier(context, 'a')
            t3 = MPI_Wtime()
            !
            call PDTRAN( dimen_s, dimen_s, 1.0d0, c_loc, 1, 1, descc, 0.0d0, a_loc, 1, 1, desca )
            !
            call blacs_barrier(context, 'a')
            t4 = MPI_Wtime()
            if (iam == 0) then
                write(out,'(f7.2,a)') t4-t3,' sec'
            endif
            !
            if (iam == 0) then
                write(out,"(/'Making A diagonal dominant...')")
            endif
            !
            do j = 1,loc_c
                global_j = indxl2g( j, nb, mycol, 0, npcol )
                do i=1,loc_r
                    global_i = indxl2g( i, nb, myrow, 0, nprow )
                    if(global_i == global_j) then
                        a_loc(i,j) = dimen_s*1.0d0
                    else
                        a_loc(i,j) = 0.0d0
                    endif
                enddo
            enddo
            !
            call ArrayStop(context,'aux:c_loc')
            deallocate (c_loc)
            !
        case (3)
            ! "RANDOM-LOCAL" - fast
            !
            allocate (seeded_array(dimen_s), stat=info)
            !
            do i = 1, dimen_s
                seeded_array(i) = 123456789-i-1;
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
                    !
                    if(global_i >= global_j) a_loc(i,j) = real(rrr(global_i),rk)/real(1099511627776.0,rk)
                    if(global_i <  global_j) a_loc(i,j) = real(rrr(global_j),rk)/real(1099511627776.0,rk)
                    if(global_i == global_j) a_loc(i,j) = a_loc(i,j) + real(10.0,rk) + real(global_j,rk)
                    !
                enddo
            enddo
            !
            deallocate (seeded_array)
            !
        case default
            !
            write(out,'("Unknown matrix generator = ",i2)'), matrix_generator
            stop "Unknown eigensolver = "
            !
    end select
    !
    call blacs_barrier(context, 'a')
    !
    t2 = MPI_Wtime()
    if (iam == 0) then
        write(out,'(/a,f12.6,a)') 'Time to Initialize the input matrix is',t2-t1,' sec'
    endif
    !
    ! --------------------------------- !
    ! allocate solution data structures !
    ! --------------------------------- !
    !
    allocate(z_loc(loc_r,loc_c),w(dimen_s),stat=info)
    matsize = int(loc_r,kind=hik)*int(loc_c,kind=hik)
    !
    call ArrayStart(context,iam,'diag_scalapack:z_loc',info,size(z_loc),kind(z_loc),matsize)
    call ArrayStart(context,iam,'diag_scalapack:w',info,size(w),kind(w))
    !
    !
    if (iam == 0) then
        write(out,"(/'Starting ',a,' ...')") trim(diagonalizer)
    endif
    !
    ! --------------------------- !
    ! initialize the eigensolver  !
    ! --------------------------- !
    !
#if defined(__ELPA) && defined(__ELPA_TIMING)
    elpa_print_times = .true.
#endif
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
#if defined(__DEBUG)
            if (iam == 0) then
                write(out,"(/'lwork =  ', i16, ' liwork = ', i16, ' lda = ', i16, ' loc_r = ', i16, ' loc_c = ', i16)") lwork, liwork, lda, loc_r, loc_c
            endif
#endif
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
#if defined(__ELPA)
        case (3)
            !
            ! How many eigenvalues/eigenvectors ?
            ! nsolv=dimen_s == requests all eigenstates
            nsolv=nroots
            !
        case (4)
            !
            ! How many eigenvalues/eigenvectors ?
            ! nsolv=dimen_s == requests all eigenstates
            nsolv=nroots
            !
#endif
        case default
            !
            write(out,'("Uknown eigensolver = ",i2)'), eigensolver
            stop "Uknown eigensolver = "
            !
    end select
    !
    ! ------------------------ !
    ! print memory footprint   !
    ! ------------------------ !
    !
    if (iam == 0) then
        call MemoryReport(context,iam,memory_now,memory_max)
    endif
    !
    ! --------------------- !
    ! call the eigensolver  !
    ! --------------------- !
    !
    call blacs_barrier(context, 'a')
    t1 = MPI_Wtime()
    !
#if defined(__IPM)
    call mpi_pcontrol( 1,"solver"//char(0))
#endif
    !
    select case (eigensolver)
          !
        case (1)
            !
            call pdsyevd('V', 'L', dimen_s, a_loc, 1, 1, desca, w, z_loc, 1, 1, descz, work, lwork, iwork, liwork, info)
            !
            ! ... what is this?
            nvals = dimen_s
            !
        case (2)
            !
            !
            call pdsyevx('V', range, 'L', dimen_s, a_loc, 1, 1, desca, vl, vu, il,iu, abstol, nvals, nvects, w, orfac, z_loc, 1, 1, descz, &
                work, lwork, iwork, liwork, ifail, iclustr, gap,  info)
            !
#if defined(__ELPA)
        case (3)
            !
            call solve_evp_real(dimen_s, nsolv, a_loc, lda, w, z_loc, lda, nb,mpi_comm_rows, mpi_comm_cols)
            !
            ! How check if ELPA ended successfully? ... boh!
            info = 0
            !
        case (4)
            !
            call solve_evp_real_2stage(dimen_s, nsolv, a_loc, lda, w, z_loc, lda, nb, mpi_comm_rows, mpi_comm_cols, mpi_comm_world)
            !
            ! How check if ELPA ended successfully? ... boh!
            info = 0
            !
#endif
        case default
            !
            write(out,'("Uknown eigensolver = ",i2)'), eigensolver
            stop "Uknown eigensolver = "
            !
    end select
    !
#if defined(__IPM)
    call mpi_pcontrol( -1,"solver"//char(0))
#endif
   ! 
   call blacs_barrier(context, 'a')
    t2 = MPI_Wtime()
    !
    ! ------------------------------ !
    ! do tests on convergence if any !
    ! ------------------------------ !
    !
    if (iam == 0) then
        write(out,'(/a,f12.6,a)') 'Time to diagonalize matrix is ',t2-t1,'sec'
        if (info .eq. 0) then
            write(out,"(/'Diagonalization finished successfully!')")
        else
            write(out,"(/'Something went bananas... info = ', i6)") info
        endif
    endif
    !
    call blacs_barrier(context, 'a')
    !   
    if (eigensolver .le. 2 ) then
        deallocate(work, iwork)
        call ArrayStop(context,'diag_scalapack:work')
        call ArrayStop(context,'diag_scalapack:iwork')
    endif

    deallocate(a_loc, z_loc, w)
    call ArrayStop(context,'diag_scalapack:a_loc')
    call ArrayStop(context,'diag_scalapack:z_loc')
    call ArrayStop(context,'diag_scalapack:w')
    !
    ! Exit BLACS
    !
    call blacs_gridexit(context)
    call blacs_exit(0)
    !
end program dirac_exomol_eigen

