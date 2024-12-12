!    Copyright (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

module diagnostic_output
contains
subroutine set_debug_level(level)
  ! Temporarily set the verbosity of the program.
  use fldebug_parameters
  implicit none
  integer, intent(in) :: level

  current_debug_level=level

end subroutine set_debug_level

subroutine set_global_debug_level(level) bind(c)
  ! Set the global verbosity of the program.
  use, intrinsic :: iso_c_binding, only : c_int
  use fldebug_parameters
  implicit none
  integer(c_int), intent(in), value :: level

  global_debug_level = level
  current_debug_level = global_debug_level
end subroutine set_global_debug_level

subroutine reset_debug_level
  ! Temporarily set the verbosity of the program.
  use fldebug_parameters
  implicit none

  current_debug_level=global_debug_level

end subroutine reset_debug_level

function debug_level()
  ! Simply return the current debug level. This makes the debug level
  ! effectively global.
  use fldebug_parameters
  implicit none
  integer :: debug_level

  debug_level=current_debug_level

end function debug_level

function get_global_debug_level() bind(c)
  ! Simply return the global debug level. This makes the debug level
  ! effectively global.
  use, intrinsic :: iso_c_binding, only : c_int
  use fldebug_parameters
  implicit none
  integer(c_int) :: get_global_debug_level

  get_global_debug_level = global_debug_level
end function get_global_debug_level
end module diagnostic_output
