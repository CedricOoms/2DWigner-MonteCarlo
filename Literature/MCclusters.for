C Program made for different confining potential (power n called 'xmacht') 
C and different power of interparticle interaction (power m called 'power'), 
C with possible Yukawa screening. 
C Program calculates the groundstate energies and configurations (in reduced 
C units) of charged particles in a confinement potential, which can be 
C different along (x,y) and z directions.   
C The accompanying datafile 'input.dat' contains all the necessary parameters 
C to run the program.
C -------------------------------------------------------------------

      PROGRAM ARTIFICIAL_ATOMS_2D

      IMPLICIT REAL*8 (a-h,o-z)
      COMMON /grid/x(1000),y(1000)
	common /bart/power,xmacht,yuk,dkapp,dpar
      DIMENSION energy_t(20000),nrt(20000),npsh(1000),xx(1000),yy(1000),
     :step(1000),stepz(1000),r1(1000),rr(1000),ic(1000),diffrr(1000)
      real time_begin, time_end, time_beginn, time_endn
	REAL*4 rndm 
      
	CHARACTER name*8
	
	
      OPEN (1,file='input.dat')      
      read (1,*) number         ! Number of particles
	read (1,*) ntests         ! Number of initial tests
      read (1,*) mont	        ! Number of Monte Carlo steps
      read (1,*) icc	        ! Parameter for random generator
      read (1,*) xmacht         ! u= xmacht potential' xmacht 
	read (1,*) power          ! Coulomb ' power 
      read (1,*) yuk            ! if yuk=1 Yukawa interaction with parameter dkapp
	read (1,*) dkapp		    ! kappa
	read (1,*) dpar			! confinement - yes or no
	CLOSE (1)                
 
    	write (name(1:8),'(a2,i4)') 'd0',number

      OPEN (3,file=name//'.coo')
	OPEN (8,file='shell.dat')  ! this file is exactly same position .coo
      open (11,file='cputime.dat')

      pi=4.d0*atan(1.d0)
      jtest=0
      energy_opt=1.e10      
             
C Throw particles *ntests* times
      DO itest=1,ntests
 
      ! Give every particle random coordinates
      DO n=1,number
         CALL RANDOM(rndm,icc)
         a=sqrt(0.25*number*rndm)
	   CALL RANDOM(rndm,icc)
         x(n)=a*cos(2*pi*rndm)
         y(n)=a*sin(2*pi*rndm)
	   step(n)=0.2d0
      END DO 

      ! calculate and write energy of random configuration
      
	CALL ENERGY_Q(number,energy_1,x,y)
      write (*,'(a5,i4,a7,i3)') 'Test=',itest,'Numb=',number
      write (*,'(a16,1pe21.14)') 'Starting energy=',energy_1
      nreplace=0
 
           ! MONTE CARLO
      DO im=1,mont      
         ! Vary coordinates of every particle
      	 DO n=1,number   

            CALL RANDOM(rndm,icc)
            xn=x(n)+step(n)*(2*rndm-1)
            CALL RANDOM(rndm,icc)
            yn=y(n)+step(n)*(2*rndm-1)

            ! DELTA_Q calculates the difference in energy
            
		  CALL DELTA_Q(number,delta,x,y,xn,yn,n)
	
	  ! If new coordinates represent a lower energy
            IF (delta.lt.0.d0) THEN
               ! Replace old coordinates by new ones               
               x(n)=xn
               y(n)=yn
	       ! nreplace = # M.C. steps --> better coordinates
               nreplace=nreplace+1               
	       ! Increase variation of coordinates
               step(n)=dmin1(1.d0,1.2d0*step(n))
	      ELSE
	         step(n)=dmax1(1.d-1,step(n)/1.2d0)
            END IF

         END DO   ! Vary coordinates
        
         ! After every 100 Monte Carlo steps
         IF (mod(im,100).eq.0) THEN    

            ! ENERGY_Q calculates the energy per particle
            CALL ENERGY_Q(number,energy,x,y)
            ! write to screen       
3           FORMAT (a5,i4,a6,1pe21.14)
	    write (*,3) 'Mont=',im,'U/N=',energy
            write (*,'(a10,i8)') 'nreplace=',nreplace
            nreplace=0

         END IF
      END DO   ! MONTE CARLO

	nn=2*number
      nnn=nn+1       
      
      ! Calculate the difference of new energy with those previously found
      delte=1.d0
      DO n=1,jtest
         delte=dmin1(delte,abs(energy-energy_t(n)))
      END DO      
	
      ! If a new energy value is found
      IF (delte.gt.1.d-8) THEN   

         jtest=jtest+1
         IF (jtest.gt.10000) goto 4
         nrt(jtest)=itest
         energy_t(jtest)=energy
 
         ! write to screen
         write (*,'(a5,i3,a7,i4,a6,1pe13.6)') 
     :         'Numb=',number,'jtes=',jtest,'E/N=',energy        
        
               
         ! Name.coo file contains coordinates of all energies found
         DO n=1,number
	        r1(n)=sqrt(x(n)**2+y(n)**2)
	    write (3,'(i4,4f15.10)')jtest,energy,x(n),y(n)
         END DO

        ! add this part to calculate how many particles in each shell for every state
	   call shell(number,x,y,ic,imax)	  
	
         do i=1, number
	      npsh(i)=0   ! npsh means number of particles in one shell
	      rr(i)=0
	   end do
             
	       do i=1,imax
	       rrr=0.d0
		      do j=1, number
	          if(ic(j).eq.i) then 
			  npsh(i)=npsh(i)+1
	          rrr=rrr+r1(j)
			  endif   
			  end do
	          rr(i)=rrr/npsh(i)
	       end do  
	do i=1,imax
	    if (i.ge.2)  diffrr(i)=rr(i)-rr(i-1)
		write(8,'(f4.2,1x,i3,1x,a7,e13.6,1x,a6,i3,1x,a9,i3,1x,f10.5,
     :1x,f10.5)')
     :xmacht,jtest,'energy=',energy,'shell=',i,'particle=',npsh(i),
     :rr(i),diffrr(i)
	enddo
	  
	   
	   ! If new energy is lower than those previously found
         IF (energy.le.energy_opt) THEN

            energy_opt=energy
          
            ! Name file contains groundstate energy and coordinates
            OPEN (1,file=name)	
            write (1,'(a20,i3)') 'Number of particles=',number
            write (1,'(a7,1pe21.14,2(a8,1pe9.2))')
     :            'Energy=',energy
            write (1,*)'  x                     y  '  
            DO n=1,number
	       write (1,'(2(1x,1pe21.14))')x(n),y(n)
            END DO 
		  CLOSE (1)
				      		
         END IF   ! Lower energy found
 
      END IF   ! New energy found

      END DO   ! Test 
       
      ! Create ALL-file to store all energies found
4     OPEN (1,file=name//'.all')
      write (1,'(a25,i4)') 'Number of energies found=',jtest
      DO n=1,jtest
         write (1,'(i4,a7,i4,2(a7,1pe21.14))')
     :         n,'Test=',nrt(n),'U/N=',energy_t(n),
     :         'E-E°=',energy_t(n)-energy_opt
      END DO
      CLOSE (1)
      
      CLOSE (3)         ! name.coo
	close (8)         ! shell.dat
      close (11)        ! cputime.dat
      

       END 

****************************   RANDOM   ***************************

      SUBROUTINE RANDOM(r,i)

  1   r=RAN1(i)
      IF (r.eq.0.0.or.r.eq.1.0) goto 1
      RETURN
      END
      
      FUNCTION RAN1(IDUM)
    
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=242000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END

**************************   ENERGY_Q   ***************************

      SUBROUTINE ENERGY_Q(number,energy,x,y)
     
	 
      IMPLICIT REAL*8 (a-h,o-z)
	common /bart/power,xmacht,yuk,dkapp,dpar
      DIMENSION x(number),y(number)

      energy=0.d0

	if (yuk.eq.1)then

      do i=1,number
         energy=energy+dpar*(dsqrt(x(i)**2+y(i)**2))**xmacht
         do j=i+1,number
         energy=energy+dexp(-dkapp*dsqrt((x(i)-x(j))**2+(y(i)-y(j))**2))
     :/dsqrt((x(i)-x(j))**2+(y(i)-y(j))**2)
         end do
      end do

      else
      DO i=1,number
         energy=energy+dpar*(sqrt(x(i)**2+y(i)**2))**xmacht
         DO j=i+1,number
          energy=energy+1/(sqrt((x(i)-x(j))**2+(y(i)-y(j))**2)**power) 
         END DO
      END DO
      endif

      energy=energy/number

      END 

***************************   DELTA_Q   ***************************

      SUBROUTINE DELTA_Q(number,delta,x,y,xx,yy,i)
      
      IMPLICIT REAL*8 (a-h,o-z)
	common /bart/power,xmacht,yuk,dkapp,dpar
      DIMENSION x(number),y(number)

      delta=dpar*(-(x(i)**2+y(i)**2)**(xmacht/2)
     :+(xx**2+yy**2)**(xmacht/2))  

	if(yuk.eq.1)then
      do j=1,number
      IF (i.ne.j) THEN   ! don't forget!!! only here need to put 'if'. 
	delta=delta
     :-dexp(-dkapp*dsqrt((x(i)-x(j))**2+(y(i)-y(j))**2))
     :/dsqrt((x(i)-x(j))**2+(y(i)-y(j))**2)
     :+dexp(-dkapp*dsqrt((x(j)-xx  )**2+(y(j)-yy  )**2))
     :/dsqrt((x(j)-xx  )**2+(y(j)-yy  )**2)
      endif
	end do
	
	else
  
      DO j=1,number
         IF (i.ne.j) THEN
        delta=delta
     :  -(1/sqrt((x(i)-x(j))**2+(y(i)-y(j))**2))**power
     :  +(1/sqrt((x(j)-xx  )**2+(y(j)-yy  )**2))**power

         END IF
      END DO
      endif
      END 

**********************    SHELL    ************************
      subroutine shell(number,x,y,ic,imax)
      implicit real*8(a-h,o-z)
	dimension   x(1000),y(1000),ic(1000)
      rmin=-rmax
      do n=1,number
      ic(n)=0
      rmax=dmax1(rmax,sqrt(x(n)**2+y(n)**2)) 
      rmin=dmin1(rmin,sqrt(x(n)**2+y(n)**2)) 
      end do
      delta=1.3*sqrt((rmax**2-rmin**2)/number)
      do i=1,number
      a=1.d10
      icmin=number
      do j=1,number
      rr=sqrt(x(j)**2+y(j)**2)
      if(rr.le.a.and.ic(j).eq.0)a=rr
      icmin=min(icmin,ic(j))
      end do
      if(icmin.ne.0)then
      imax=i-1
      return
      end if
      a=0.95*a
      b=1.1*a+delta
      do j=1,number
      rr=sqrt(x(j)**2+y(j)**2)
      if(rr.ge.a.and.rr.le.b.and.ic(j).eq.0)ic(j)=i
      end do
      end do
	end	
*******************************************************************


