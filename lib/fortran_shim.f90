module shims
  use, intrinsic :: iso_c_binding
  ! ---- Replace these with your actual modules ----
  use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
  use type_crystalstructure,          only: lo_crystalstructure
  use lo_memtracker, only: lo_mem_helper
  ! ------------------------------------------------
  implicit none
contains

  ! --- tiny helper: convert C string -> allocatable Fortran string ---
  pure integer function cstr_len(c) result(n)
    character(kind=c_char), intent(in) :: c(*)
    integer :: i
    n = 0
    do i = 1, huge(1)
      if (c(i) == c_null_char) exit
      n = n + 1
    end do
  end function cstr_len

  subroutine cstr_to_fstr(c, f)
    character(kind=c_char), intent(in)  :: c(*)
    character(len=:),       allocatable, intent(out) :: f
    integer :: n, i
    n = cstr_len(c)
    allocate(character(len=n) :: f)
    do i = 1, n
      f(i:i) = achar(iachar(c(i)))
    end do
  end subroutine cstr_to_fstr

  ! =========================
  ! MEM helper: new / destroy
  ! =========================
  subroutine lo_mem_new_ptr(out_ptr) bind(C, name="lo_mem_new_ptr")
    type(c_ptr),  intent(out) :: out_ptr
    type(lo_mem_helper), pointer :: mem

    out_ptr = c_null_ptr
    allocate(mem, stat=status); if (status /= 0) return
    call mem%init()
    out_ptr = c_loc(mem)
  end subroutine lo_mem_new_ptr

  subroutine lo_mem_destroy_ptr(ptr) bind(C, name="lo_mem_destroy_ptr")
    type(c_ptr), value :: ptr
    type(lo_mem_helper), pointer :: mem
    if (c_associated(ptr)) then
      call c_f_pointer(ptr, mem)
      if (associated(mem)) deallocate(mem)
    end if
  end subroutine lo_mem_destroy_ptr

  ! =======================================
  ! Crystal: readfromfile -> pointer / dtor
  ! =======================================
  subroutine lo_crystal_from_file_ptr(c_fname, mem_ptr, verbosity, out_ptr, status) &
      bind(C, name="lo_crystal_from_file_ptr")
    character(kind=c_char), intent(in) :: c_fname(*)
    type(c_ptr),            value      :: mem_ptr
    integer(c_int),         value      :: verbosity
    type(c_ptr),            intent(out):: out_ptr
    integer(c_int),         intent(out):: status

    type(lo_crystalstructure), pointer :: p
    type(lo_mem_helper),      pointer :: mem
    character(len=:), allocatable :: fname

    out_ptr = c_null_ptr; status = 0

    if (.not. c_associated(mem_ptr)) then
      status = -1; return
    end if
    call c_f_pointer(mem_ptr, mem)

    call cstr_to_fstr(c_fname, fname)

    allocate(p, stat=status); if (status /= 0) return
    call p%readfromfile(p, fname, mem, verbosity)

    out_ptr = c_loc(p)
  end subroutine lo_crystal_from_file_ptr

  subroutine lo_crystal_destroy_ptr(ptr) bind(C, name="lo_crystal_destroy_ptr")
    type(c_ptr), value :: ptr
    type(lo_crystalstructure), pointer :: p
    if (c_associated(ptr)) then
      call c_f_pointer(ptr, p)
      if (associated(p)) deallocate(p)
    end if
  end subroutine lo_crystal_destroy_ptr

  ! ====================================================
  ! FC2: readfromfile (needs crystal + mem) / destructor
  ! ====================================================
  subroutine lo_fc2_from_file_ptr(c_fname, crys_ptr, mem_ptr, verbosity, out_ptr, status) &
      bind(C, name="lo_fc2_from_file_ptr")
    character(kind=c_char), intent(in) :: c_fname(*)
    type(c_ptr),            value      :: crys_ptr, mem_ptr
    integer(c_int),         value      :: verbosity
    type(c_ptr),            intent(out):: out_ptr
    integer(c_int),         intent(out):: status

    type(lo_forceconstant_secondorder), pointer :: fc
    type(lo_crystalstructure),          pointer :: p
    type(lo_mem_helper),                pointer :: mem
    character(len=:), allocatable :: fname

    out_ptr = c_null_ptr; status = 0

    if (.not. c_associated(crys_ptr)) then; status = -1; return; end if
    if (.not. c_associated(mem_ptr))  then; status = -1; return; end if
    call c_f_pointer(crys_ptr, p)
    call c_f_pointer(mem_ptr,  mem)

    call cstr_to_fstr(c_fname, fname)

    allocate(fc, stat=status); if (status /= 0) return
    call fc%readfromfile(fc, p, fname, mem, verbosity)

    out_ptr = c_loc(fc)
  end subroutine lo_fc2_from_file_ptr

  subroutine lo_fc2_destroy_ptr(ptr) bind(C, name="lo_fc2_destroy_ptr")
    type(c_ptr), value :: ptr
    type(lo_forceconstant_secondorder), pointer :: fc
    if (c_associated(ptr)) then
      call c_f_pointer(ptr, fc)
      if (associated(fc)) deallocate(fc)
    end if
  end subroutine lo_fc2_destroy_ptr


end module shims
