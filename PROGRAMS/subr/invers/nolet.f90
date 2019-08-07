! ----------------------------------------------------------
	subroutine pstomo(m,n,x,u,v,w,itmax,nz,aaa,ncolrow,ncol)
!
!   subroutine to solve the linear tomographic problem Ax=u using the 
!   lsqr algorithm.
!
!   Input: m is the number of data (rows), n the number of unknowns (columns),
!      u contains the data (is overwritten), itmax is the number of iterations.
!   Output: x is the solution; intermediate results are written to a diskfile
!      named <tomo.int> which must not already exist before the call.
!   Scratch: arrays v(n) and w(n)
!   Subroutines: routines avpu and atupv must be supplied by the user.
!      avpu(m,n,u,v) computes u=u+A*v for given input u,v (overwrites u)
!      atupv(m,n,u,v) computes v=v+A(transpose)*u for given u,v (overwrites v)
!
	real x(n),u(m),v(n),w(n)
	real aaa(nz)
	integer ncolrow(nz),ncol(m)
	x=0.
	v=0.
	call normlz(m,u,beta)
	b1=beta 
	!write(*,*)' b1=',b1
	call atupv(m,n,u,v,nz,aaa,ncolrow,ncol)
	!write(*,*)' b1=',b1
	call normlz(n,v,alfa)
	!write(*,*)' alfa=',alfa,' beta=',beta
	rhobar=alfa
	phibar=beta 
	w=v
	do iter =1,itmax 		! repeat for fixed nr of iterations
		!if(mod(iter,20).eq.0) write(*,*)' iter=',iter
   	   a=-alfa
   	   u=a*u			! bidiagonalization
   	   call avpu(m,n,u,v,nz,aaa,ncolrow,ncol)
   	   call normlz(m,u,beta) 
   	   b=-beta
   	   v=b*v
   	   call atupv(m,n,u,v,nz,aaa,ncolrow,ncol) 
   	   call normlz(n,v,alfa)
!	write(*,*)' rhobar*rhobar+beta*beta=',rhobar*rhobar+beta*beta
	if(rhobar*rhobar+beta*beta.lt.0.) then
		!write(*,*)' rhobar*rhobar+beta*beta=',rhobar*rhobar+beta*beta
	end if
   	   rho=sqrt(rhobar*rhobar+beta*beta)	! modified QR factorization
   	   c=rhobar/rho
   	   s=beta/rho 
   	   teta=s*alfa
   	   rhobar=-c*alfa
   	   phi=c*phibar
   	   phibar=s*phibar 
   	   t1=phi/rho
   	   t2=-teta/rho
   	   xnorm=0.0
   	   do i=1,n				! update solution x and storage vector w
      	      x(i)=t1*w(i)+x(i)
      	      xnorm=xnorm+abs(x(i))
      	      w(i)=t2*w(i)+v(i)
           end do
   	   r=phibar/b1 
   	   xnorm=xnorm/n
	end do
	return
	end


	subroutine avpu(m,n,u,v,nz,aaa,ncolrow,ncol)		! computes u=u+a*v
	real u(m),v(n)
	real aaa(nz)
	integer ncolrow(nz),ncol(m)
	kount=0
	do i=1,m 				! work row by row
  	   nc=ncol(i)
  	   if (nc.eq.0) cycle
  	   do j=1,nc 
    	      kount=kount+1
              u(i)=u(i)+aaa(kount)*v(ncolrow(kount))
 	   end do
	end do 	   
	return
	end
	
	
	subroutine atupv(m,n,u,v,nz,aaa,ncolrow,ncol)		! computes v=v+a(transpose)*u
	real u(m),v(n)
	real aaa(nz)
	integer ncolrow(nz),ncol(m)
	kount=0
	!write(*,*)' m=',m,' nz=',nz,' n=',n
	do i=1,m 				! work row by row (here too!)
		!if(mod(i,100).eq.0)write(*,*)' i=',i
  		nc=ncol(i)
!if(i.ge.10920)write(*,*)' i=',i,' nc=',nc
!		if(i.ge.24662)write(*,*)' i=',i
  	   if (nc.eq.0) cycle
  	   do j=1,nc 
    		kount=kount+1
!if(i.eq.10920)write(*,*)' kount=',kount,' j=',j
    		jj=ncolrow(kount)
!if(i.eq.10920) write(*,*)' i=',i,kount,' aaa=',aaa(kount),' u=',u(i)
!	if(i.ge.24662)write(*,*)' jj=',jj
    		v(jj)=v(jj)+aaa(kount)*u(i)	! add to the right vector element
!	if(i.ge.24662)write(*,*)' v(jj)=',v(jj)
	   end do
	end do
	!write(*,'(15f5.1)')(v(i),i=300,500)
	!do i=300,500
		!write(*,*)i,v(i)
	!end do
	!write(*,*)' n=',n
	return
	end

	subroutine normlz(n,x,s)
	real x(n)
	s=0.
	!write(*,'(10f5.1)')(x(i),i=1,n)
	do i=1,n
	   s=s+x(i)**2
	end do
	!write(*,*)' s=',s
	s=sqrt(s)
	ss=1./s
	x=x*ss
	return
	end


	function vref(z, nma,hma,vma)
	real hma(300),vma(300)
	do ilay=2,nma
		z1=hma(ilay-1)
		z2=hma(ilay)
		if(z.ge.z1.and.z.le.z2) then
			v1=vma(ilay-1)
			v2=vma(ilay)
			vref=v1+(v2-v1)/(z2-z1)*(z-z1)
		end if
	end do
	return
	end
