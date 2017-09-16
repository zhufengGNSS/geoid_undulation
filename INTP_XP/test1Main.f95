
program test1Main
    use mGeoidUn
    implicit none
    real*8 tElapsedT, VecFlat(3), VecFlon(3), VecU(3)
    integer*8 :: order = 30

    data VecFlat /0.0, 90.0, -90.0/
    data VecFlon / 0.0, 180.0, -180.0/

    call init_const1(order)
        write(*,*) 'MAXN：', MAXN
    call f477(tElapsedT, VecFlat, VecFlon, VecU)
        write(*,*) 'MAXN：', MAXN
    write(*,*) '时间消耗统计：', tElapsedT, ' 秒'

end program test1Main
