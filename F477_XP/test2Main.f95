
!! program test2Main
subroutine test2Main
    use mGeoidUn

    integer*8, parameter :: CLen = 3
    integer(kind = 4) :: Norder = 30, vecLen = CLen
    real*8 :: VecFLAT(CLen), VecFLON(CLen), VecU(CLen), dElapseT, &
    outputM(CLen,CLen)

    data  VecFLAT/ 0.0, 90.0, -90.0 /, VecFLON/0.0, 180.0, -180.0/

    call init_const1(Norder)
    call f477(dElapseT, VecFLAT, VecFLON , VecU, vecLen, Norder)

    outputM(:,1) = VecFLAT
    outputM(:,2) = VecFLON
    outputM(:,3) = VecU

    open(20, file='OUTF477.DAT', status='unknown')

    WRITE(20, 101) transpose(outputM)

101  FORMAT(2F14.7,F10.3)

    close(20)

    stop

end
