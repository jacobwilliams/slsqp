!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/22/2015
!  license: BSD
!
!  Numeric kind definitions.

    module slsqp_kinds

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    private

    integer,parameter,public :: wp = real64  !! Using "double precision" real kinds

    end module slsqp_kinds
!*****************************************************************************************
