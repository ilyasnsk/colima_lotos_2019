USE DFPORT

call prepare_noise()


open(1,file='data_in/rays.dat')
open(11,file='data_out/rays.dat')

nray=0
nsrc=0
np=0
ns=0
! Read the sources:
992	continue
    read(1,*,end=991)fini,tini,zold,nkrat
    write(11,*)fini,tini,zold,nkrat

    nsrc=nsrc+1

    do i=1,nkrat
        nline=nline+1
        read(1,*)ips,ist,tobs
        err=our_noise(ips)
        write(11,*)ips,ist,tobs+err
        !write(*,*)ips,ist,tobs,err,tobs+err
        if(ips.eq.1) np=np+1
        if(ips.eq.2) ns=ns+1
        nray=nray+1
    end do
    goto 992

991 close(1)
close(11)

write(*,*)' nsrc=',nsrc,' nray=',nray
write(*,*)' np=',np,' ns=',ns
stop
end