module test_simple
  use funit
  implicit none

contains
  !@test
  subroutine test_assert_true_and_false()
#line 8 "/home/ccaveayl/projects/u-dales/u-dales/unit_tests/test_simple.pf"
  call assertTrue(1 == 1, &
 & location=SourceLocation( &
 & 'test_simple.pf', &
 & 8) )
  if (anyExceptions()) return
#line 9 "/home/ccaveayl/projects/u-dales/u-dales/unit_tests/test_simple.pf"
#line 9 "/home/ccaveayl/projects/u-dales/u-dales/unit_tests/test_simple.pf"
  call assertFalse(1 == 2, &
 & location=SourceLocation( &
 & 'test_simple.pf', &
 & 9) )
  if (anyExceptions()) return
#line 10 "/home/ccaveayl/projects/u-dales/u-dales/unit_tests/test_simple.pf"
  end subroutine test_assert_true_and_false

  !@test
  subroutine test_dqsatdT()
    use initfac, only: dqsatdT
#line 15 "/home/ccaveayl/projects/u-dales/u-dales/unit_tests/test_simple.pf"
  call assertEqual(0.1384832710e-2, dqsatdT(300.0), &
 & location=SourceLocation( &
 & 'test_simple.pf', &
 & 15) )
  if (anyExceptions()) return
#line 16 "/home/ccaveayl/projects/u-dales/u-dales/unit_tests/test_simple.pf"
  end subroutine
end module test_simple

module Wraptest_simple
   use FUnit
   use test_simple
   implicit none
   private

contains


end module Wraptest_simple

function test_simple_suite() result(suite)
   use FUnit
   use test_simple
   use Wraptest_simple
   implicit none
   type (TestSuite) :: suite

   class (Test), allocatable :: t

   suite = TestSuite('test_simple_suite')

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_assert_true_and_false', &
      test_assert_true_and_false))
   call suite%addTest(t)

   if(allocated(t)) deallocate(t)
   allocate(t, source=TestMethod('test_dqsatdT', &
      test_dqsatdT))
   call suite%addTest(t)


end function test_simple_suite

