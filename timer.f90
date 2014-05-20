module timer
    !
    use accuracy
    implicit none
    private
    public ArrayStart, ArrayStop, MemoryReport, memory_limit
    !
    integer, parameter :: trk        = selected_real_kind(12)
    integer, parameter :: table_size = 1000 ! Max number of entries to track
    integer, parameter :: tarray_size = 10000 ! Max number of entries to track
    integer, parameter :: name_len   =   40 ! Max length of timer name
    !
    real(rk)           :: maxmemory   =   0 ! Maximal memory allocated
    real(rk)           :: memory_limit =   3.8 ! Gb
    !
    type tim
        logical                 :: used       ! Slot used?
        logical                 :: active     ! Currently active?
        character(len=name_len) :: name       ! Timer name
        real(trk)               :: calls      ! Number of times the timer was invoked
                                              !
                                              ! All times below are in seconds
                                              !
        real(trk)               :: real_time  ! Total real time on this timer
        real(trk)               :: cpu_time   ! ditto for CPU time
        real(trk)               :: real_kids  ! Real time spent in nested timers
        real(trk)               :: cpu_kids   ! ditto for CPU time
        real(trk)               :: real_start ! For active timers, time of activation
        real(trk)               :: cpu_start  ! ditto for CPI time
        integer(ik)             :: stack_p    ! For active timers, position in the stack
    end type tim
    !
    !  Timers and sundry data
    !
    type(tim), target :: t_table (table_size) ! All of our timers
    integer(ik)       :: t_nested(table_size) ! Stack of currently active timers
    integer(ik)       :: t_appear(table_size) ! Appearance order for the timers
    integer(ik)       :: t_count  = 0         ! Number of defined timers
    integer(ik)       :: t_active = 0         ! Number of currently active timers
    real(trk)         :: prog_start           ! Timebase for the real time


    ! -----  io-units
    type tio_unit
        logical                 :: used       ! Slot used?
        logical                 :: active     ! Currently active?
        character(len=name_len) :: name       ! Timer name
        integer(ik)             :: stack_p    ! For active units, position in the stack
        integer(ik)             :: slot     ! Number of times the unit was invoked
    end type tio_unit


    ! -----  io-units
    type tarray_unit
        logical                 :: used       ! Slot used?
        logical                 :: active     ! Currently active?
        character(len=cl)       :: name       ! Timer name
        integer(ik)             :: stack_p    ! For active arrays, position in the stack
        integer(ik)             :: slot       ! Number of times the unit was invoked
        real(trk)             :: size       ! total size
    end type tarray_unit


    type(tarray_unit), target :: array_table (tarray_size) ! All of our array uints
    integer(ik)       :: array_appear(tarray_size) ! Appearance order for the i/o-unit
    integer(ik)       :: array_count =  0         ! Number of defined io units (we will start counting from unit=10)
    integer(ik)       :: array_active = 0         ! Number of currently active units
    real(trk)         :: array_size   = 0         ! Number of currently active units


    type(tio_unit), target :: io_table (table_size) ! All of our i/o uints
    integer(ik)       :: io_appear(table_size) ! Appearance order for the i/o-unit
    integer(ik)       :: io_count =  0         ! Number of defined io units (we will start counting from unit=10)
    integer(ik)       :: io_active = 0         ! Number of currently active units


!
contains
    !
    !  Linear insertion with linear search (Algorithm L)
    !
    function insert_item(name) result(pos)
        character(len=*), intent(in) :: name
        integer(ik)                  :: pos
        !
        pos = string_hash(name)
        search: do
            if (.not.t_table(pos)%used) then
                !
                ! This is a new key, insert it
                !
                t_count = t_count + 1
                if (t_count>=table_size/5) then
                    write (out,"('Too many timers. Increase table_size in "// &
                        "timer.f90 to at least ',i5)") t_count*5
                    stop 'timer%insert_item'
                end if
                t_appear(t_count)      = pos
                t_table(pos)%used      = .true.
                t_table(pos)%active    = .false.
                t_table(pos)%name      = name
                t_table(pos)%calls     = 0
                t_table(pos)%real_time = 0
                t_table(pos)%cpu_time  = 0
                t_table(pos)%real_kids = 0
                t_table(pos)%cpu_kids  = 0
                exit search
            end if
            if (t_table(pos)%name==name) then
                !
                ! This is an existing key, simply return the location
                !
                exit search
            end if
            pos = 1 + modulo(pos-2,table_size)
        end do search
      !
    end function insert_item
    !
    integer function string_hash(str) result(h)
        character(len=*), intent(in) :: str
        !
        integer :: i, chr, g
        integer :: mask
        data mask/Z"1FFFFFF"/
        !
        !    This hash assumes at least 29-bit integers. It is supposedly
        !    documented in Aho, Sethi, and Ullman, pp. 434-438
        !
        h = 0
        do i=1,len_trim(str)
            chr = ichar(str(i:i))
            h   = ishft(h,  4) + chr
            g   = ishft(h,-24)
            h   = iand(ieor(h,g),mask)
        end do
        h = 1 + modulo(h,table_size)
    end function string_hash
    !

    !
    !  Start new array.
    !
    subroutine ArrayStart(comm,iam,name,alloc,isize,ikind,hsize)
        integer(ik), intent(in)      :: comm,iam   ! communicator
        character(len=*), intent(in) ::  name  ! Unit name
        integer(ik), intent(in)      :: alloc  ! allocation error
        integer(ik), intent(in)      :: isize  ! Unit size
        integer(ik), intent(in)      :: ikind  ! Unit kind
        integer(hik), intent(in),optional   :: hsize  ! Unit size
        !
        integer(ik)          :: pos  ! Unit position
        type(tarray_unit), pointer   :: t    ! Current io_unit (for convenience)
        logical                   :: ifopen
        real(rk)                  :: size_,memory_now,memory_
        !
        !  One-time initialization
        !

        if (alloc/=0) then
            write(out,"(/' Error ',i8,' trying to allocate array ',a)") alloc,name
            !write(out,"( ' Total memory ',f14.2,' array size =  ',f18.2,' Gb')") array_size,(real(ikind)*real(isize))/real(1024**3)
            !call MemoryReport(comm,iam,memory_now,memory_)
            stop 'ArrayStart - allocation error'
        end if
        !
        if (array_count==0) then
            array_table(:)%used = .false.
        end if
        !
        pos =  insert_arrayunit(comm,name)
        t   => array_table(pos)
        !
        size_ = (ikind*real(isize))/real(1024**3) ! size in GByte
        if (present(hsize)) size_ = (ikind*real(hsize))/real(1024**3)
        !
#if defined(__MEM_DEBUG)
        write (out,"(/' Size of array ',a,' contributed of ',f20.8,' GByte')") name, size_
#endif
        !
        if (size_<0.0_trk) then
            write (out,"(' Size of array ',a,' is negative ',f20.8)") name,size_
            size_ = 1
            call MemoryReport(comm,iam,memory_now,memory_)
            !
            stop 'ArrayStart - negative size'
        end if
        !
        if (.not.t%active) then
            !write (out,"('ArrayStart: ArrayStart ',a,' is already active')") trim(name)
            !stop 'ArrayStart - nested ArrayStart'
            !
            !  Push the new array to the array stack
            !
            array_active = array_active + 1
            t%active     = .true.
          !
        end if
        !
        array_size   = array_size + size_
#if defined(__MEM_DEBUG)
        write (out,"(/' Cumulated memory : ',f20.8,' GByte')") array_size
#endif
        !
        !maxmemory = max(array_size,maxmemory)
        !
#if 0
        if (array_size>memory_limit) then
            !
            write(out,"('Run out of memory')")
            call MemoryReport(comm,iam,memory_now,memory_)
            stop 'Run out of memory'
          !
        endif
#endif
        !
        if (trim(array_table(pos)%name)==name) then
            !
            t%size = t%size + size_
            !
        end if
      !
      !slot  = t%slot
      !
    end subroutine ArrayStart

    !
    !  End an array-unit
    !
    subroutine ArrayStop(comm,name)
        integer(ik), intent(in)      :: comm   ! communicator
        character(len=*), intent(in) :: name  ! Unit name
        !
        integer(ik)        :: pos        ! unit position
        type(tarray_unit), pointer :: t          ! Current unit (for convenience)
        real(rk)           :: mem
        !
        pos =  insert_arrayunit(comm,name)
        t   => array_table(pos)
        !
        if (.not.t%active) then
            write (out,"('ArrayStop: array ',a,' is not running')") trim(name)
            stop 'ArrayStop - inactive array counter'
        end if
        !
        array_size   = array_size - t%size
        t%size = 0
        mem = array_size
        !
        !if (array_size<=0) then
        !
        !  Mark as inactive
        !
        t%active    = .false.
        ! 
        !  Pop the timer from stack, and update counts for the parent
        !
        array_active = array_active - 1
        !
      !endif 
      !
    end subroutine ArrayStop

    !
    !  Linear insertion with linear search (Algorithm L)
    !
    function insert_arrayunit(comm,name) result(pos)
        integer(ik), intent(in)      :: comm   ! communicator
        character(len=*), intent(in)  :: name
        integer(ik)                  :: pos
        !
        pos = string_hash(name)
        search: do
            if (.not.array_table(pos)%used) then
                !
                ! This is a new key, insert it
                !
                array_count = array_count + 1
                if (array_count>=tarray_size/4) then
                    write (out,"('Too many array_units. Increase tarray_size in "// &
                        "timer.f90 to at least ',i5)") array_count*4+1
                    stop 'timer%insert_item'
                end if
                array_appear(array_count)      = pos
                array_table(pos)%used      = .true.
                array_table(pos)%active    = .false.
                array_table(pos)%name      = name
                array_table(pos)%slot      = array_count
                array_table(pos)%size      = 0
                exit search
            endif
            if (trim(array_table(pos)%name)==name) then
                !
                ! This is an existing key, simply return the location
                !
                !array_table(pos)%size = array_table(pos)%size + size_
                !
                exit search
              !
            end if
            pos = 1 + modulo(pos-2,tarray_size)
        end do search
      !
    end function insert_arrayunit

    !
    !  Produce timing report
    !
    subroutine MemoryReport(comm,iam,memory_now,memory_max)
        integer(ik), intent(in)      :: comm,iam   ! communicator
        real(trk)             :: mem
        real(trk),intent(out) :: memory_now,memory_max
        real(trk)          :: mem_threshold
        integer(ik)        :: ord
        integer(ik)        :: pos, kid_pos
        type(tarray_unit), pointer :: t, k
        character(len=1)   :: active
        integer(ik)        :: omitted
        !
        memory_now  = array_size
        !
        maxmemory = max(memory_now,maxmemory)
        !
        memory_max = maxmemory
        !
#if defined(__MEM_DETAILED)
        write (out,"(/t2,'Detailed memory report:')")
#endif

        omitted = 0
        scan: do ord=1,array_count
            pos = array_appear(ord)
            t => array_table(pos)
            if (.not.t%used) then
                write (out,"('Array ',i4,' in slot ',i5,' is defined but unused?!')") ord, pos
                stop 'MemoryReport - logic error'
            end if
            !
#if defined(__MEM_DETAILED)
            write (out,"(t5,'(Processor: ',i6,') Array = ',a25,', Total memory = ',f7.4,' GByte')") iam, TRIM(t%name), t%size
#endif
            !
            ! Calculate active-array corrections
            !
            mem= 0
            mem = mem   + t%size
            !
            !  Output needed?
            !
            if (mem<mem_threshold.or.mem<1e-7.or..not.t%active) then
                omitted = omitted + 1
                cycle scan
            end if
        end do scan
 
        write (out,"(/t2,'(Processor: ',I6,') Total memory = ',f7.4,' GByte ( ',f5.2,' GByte )')") iam, memory_now, memory_limit

#if 0
        if (omitted>0) then
            write (out,"(/' (',i9,' arrays contributing less than 1% are not shown)')") &
                omitted
        end if
        write (out,"()")
#endif
    end subroutine MemoryReport

end module timer
