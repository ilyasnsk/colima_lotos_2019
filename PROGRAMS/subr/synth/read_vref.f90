subroutine read_vref(ar,md)
character*8 ar,md
common/refmod/nref,zref(600),vref(600,2)

open(1,file='../../../DATA/'//ar//'/'//md//'/ref_syn.dat')
read(1,*,end=81)vpvs
iref=0
82 continue
    read(1,*,end=81)z,vp,vs
    iref=iref+1
    zref(iref)=z
    vref(iref,1)=vp
    if(vpvs.lt.0.000001) then
        vref(iref,2)=vs
    else
        vref(iref,2)=vref(iref,1)/vpvs
    end if
    goto 82
81 close(1)
nref=iref
write(*,*)' nref=',nref

return
end