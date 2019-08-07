

  subroutine kram3(aa,b,x)
  real aa(3,3), a(3,3),b(3),x(3),d,det
  call detr3(aa,DET)
  do i=1,3
     a=aa
     do J=1,3
        a(j,i)=b(j)
     end do
     call detr3(a,d)
     x(i)=d/det
  end do
  return
  end



     SUBROUTINE DETR3(A,D)
     REAL A(3,3),d
     INTEGER num(3,3,3), i1,i2,i3
     num=0
     num(1,2,3)=1
     num(1,3,2)=-1
     num(2,1,3)=-1
     num(2,3,1)=1
     num(3,1,2)=1
     num(3,2,1)=-1
     D=0.
     DO I1=1,3
        DO I2=1,3
           DO I3=1,3
              D=D+A(1,I1)*A(2,I2)*A(3,I3)*num(i1,i2,i3)
           end do
        end do
     end do
     return
     end

