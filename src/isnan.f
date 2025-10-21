      module my_is_nan_mod
        !poor man's replacement for the intrinsic module
        !do not use if the compiler supports the F2003 version 
        !make sure not to use fast math optimizations when compiling
        use iso_fortran_env
        !
        !the common extension isnan() may actually fail with fast math optimizations
        !
        interface my_is_nan
          module procedure my_is_nan_real32
          module procedure my_is_nan_real64
        end interface
        !
        contains
          logical function my_is_nan_real32(x) result(res)
            real(real32), intent(in) :: x
            res = x /= x
          end function
          !
          logical elemental function my_is_nan_real64(x) result(res)
            real(real64), intent(in) :: x
            res = x /= x
          end function
      end module