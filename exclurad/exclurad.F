      program exclurad
c	latest update: Jan.7, 2003
c      version 1.1  December 31, 2000
c
c      Andrei Afanasev	 afanas@jlab.org
c      Igor Akushevich	    aku@jlab.org
c      Volker Burkert	burkert@jlab.org
c      Kyungseon Joo	   kjoo@jlab.org
c
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      common/keys/ivec,iphy
      common/born/w,w2,q2,s,x,sx,sxp,sqlw0,ssh0,xxh0,sxtm0,sqll
     .	   ,cosh,phih,sinh,vmax,vv10,vv20,vvp0,vvm0,t0	!  /born/
     .	   ,cos1,sin1,cos2,sin2,ee1,ee2,pp1,pp2,eeq,ppq !  /born/
      common/isfisf/isf1,isf2

      parameter(npoimax=10000)
      dimension q2_input(npoimax),w_input(npoimax)
      dimension cscm_input(npoimax),phicm_input(npoimax)
      dimension sibar(2),sigar(2),silar(2)
      parameter (netx=1)
      parameter (nety=6)
      parameter (nett=2)
      dimension wnet(netx),q2net(nety),ctnet(nett)
      data wnet/1.5d0/
     .	   q2net/0.7d0,1d0,1.3d0,1.5d0,1.8d0,2.0d0/
     .	   ctnet/-.1d0,-.3d0/
      data ipolmax/2/

	 open(unit=22,file='radcor.dat',status='unknown')
	 open(unit=41,file='radsigpl.dat',status='unknown')
	 open(unit=42,file='radsigmi.dat',status='unknown')
	 open(unit=27,file='radasm.dat',status='unknown')
	 open(unit=28,file='radtot.dat',status='unknown')
	 open(unit=5,file='input.dat',status='old')

	  read(5,*)iphy
	  read(5,*)ll
	  read(5,*)bmom
	  read(5,*)tmom
	  read(5,*)lepton
	  read(5,*)ivec
	  read(5,*)cutv
	  read(5,*)npoi

	  if (npoi.ge.1)then
	    read(5,*)(w_input(i),i=1,npoi)
	    read(5,*)(Q2_input(i),i=1,npoi)
	    read(5,*)(cscm_input(i),i=1,npoi)
	    read(5,*)(phicm_input(i),i=1,npoi)
	  elseif (npoi.eq.-100) then
	    ii = 0
	    read(5,*)w_dummy
	    read(5,*)q2_dummy
	    do ic = 1, 10
	       csthcm_dummy = -0.9+(ic-1)*0.2
	       do ip = 1, 1
		  ii = ii + 1
		  phicm_dummy = 15.+(ip-1)*30
		  q2_input(ii) = q2_dummy
		  w_input(ii)  = w_dummy
		  cscm_input(ii) = csthcm_dummy
		  phicm_input(ii) = phicm_dummy
	       enddo
	    enddo
	    npoi = ii
	  else
c	     if(intphi.ne.1)stop 'intphi must be =1 for npoi=0'
	    do it=1,nett
	      do ix=1,netx
		do iy=1,nety
		  npoi=npoi+1
		  w_input(npoi)=wnet(ix)**2
		  q2_input(npoi)=q2net(iy)
		  cscm_input(npoi)=ctnet(it)
		enddo
	      enddo
	    enddo
	  endif

	 call titout('03.31.2000')
	 write(9,'(a8,f8.3)')' bmom  =',bmom
	 write(9,'(a8,f8.3)')' tmom  =',tmom
	 write(9,'(a8,i2)')' lepton=',lepton
	 write(9,'(a8,i2)')' ivec  =',ivec
	 write(9,'(a8,f6.3)')' cutv  =',cutv
	 write(9,'(a8,i3)')' npoi  =',npoi

	 call setcon(ivec,lepton)

      s=2.*(sqrt(tmom**2+amp*amp)*sqrt(bmom**2+aml2)+bmom*tmom)

      print *, 'npoi=', npoi, ' data points'

      do  i=1,npoi
	w=w_input(i)
	q2=q2_input(i)

	if(ivec.eq.1)then
	  cosh= -cscm_input(i)
	  phih=pi+phicm_input(i)*pi/180d0
	elseif(ivec.eq.2)then
	  cosh=  cscm_input(i)
	  phih=phicm_input(i)*pi/180d0
	endif


       call borntest(bmom,w,q2,cscm_input(i),phicm_input(i)*pi/180.
     .			 ,sibtu,sibtp,eps_d)

       call bornkin(cutv)

       do ipol=1,ipolmax
	 if(ipol.eq.1)then
	   isf1=1
	   isf2=4
	   sibt=sibtu
	 elseif(ipol.eq.2)then
	   isf1=5
	   isf2=5
	   sibt=sibtp
	 endif

       call bornin(sib)

       sum=vacpol(q2)

       dlm=log(q2/aml2)

       deltavr=(1.5d0*dlm-2.d0-.5d0*log(xxh0/ssh0)**2
     .		    +fspen(1d0-amp2*q2/ssh0/xxh0)-pi**2/6.d0)


       delinf=(dlm-1.d0)*log(vmax**2/ssh0/xxh0)
       extai1=exp(alpha/pi*delinf)

	tai=0.
c	 taill=0.
	if(ll.ne.1)then 
	  call qqt(tai)
	endif
	
	call qqtll(taill)

       sig=sib*extai1*(1.d0+alpha/pi*(deltavr+sum))+tai
       sigll=sib*extai1*(1.d0+alpha/pi*(deltavr+sum))+taill
c
       write(*,'(a8,8g11.3)') ' test : ',sib,sibt,sib/sibt

       if(ipol.eq.1)ddcr=sig/sib

       write(22,'(8f8.3)')w,q2,eps_d,cscm_input(i),phicm_input(i)
     . ,sig/sib,sib/sibt,sigll/sib

       write(*,'(8f8.3)') w,q2,eps_d,cscm_input(i),phicm_input(i)
     . ,sig/sib,sib/sibt,sigll/sib

c      write(*,'(8f8.5)') extai1-1.,alpha/pi*(deltavr+sum),tai/sib

       silar(ipol)=sigll
       sigar(ipol)=sig
       sibar(ipol)=sib

       enddo

       dd=1.
       if(ipolmax.eq.2)then
	 asb=sibar(2)/sibar(1)*1d2
	 as =sigar(2)/sigar(1)*1d2
	 asl=silar(2)/silar(1)*1d2
	 dd=(as-asb)/asb*1d2
	 dl=(asl-asb)/asb*1d2


	 write(27,'(5f6.2,5F8.3)')
     .	 w,q2,eps_d,cscm_input(i),phicm_input(i),as,asb,dd
	 write(41,'(5f6.2,2g11.3)')
     . 	 w,q2,eps_d,cscm_input(i),phicm_input(i)
     .	,sibar(1)+sibar(2),sigar(1)+sigar(2)
	 write(42,'(5f6.2,2g11.3)')
     .	 w,q2,eps_d,cscm_input(i),phicm_input(i)
     .	,sibar(1)-sibar(2),sigar(1)-sigar(2)
       endif

       write(28,'(5f6.2,5F8.3)')
     .	   e1,w,q2,cscm_input(i),phicm_input(i),ddcr,dd

      enddo
      close(22)
      end



****************** setcon *************************************

      subroutine setcon(ivec,lepton)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      parameter(ncases=2)
      dimension amhad(ncases),amunm(ncases)
      data amhad/0.938272d0,0.13957d0/
      data amunm/0.13957d0,0.938272d0/

      if(lepton.eq.1)aml2=.261112d-6
      if(lepton.eq.2)aml2=.111637d-1
      pi=3.1415926d0
      alpha=.729735d-2
      barn=.389379d6
      amh2=amhad(ivec)**2
      amu2=amunm(ivec)**2

      amp=.938272d0
      amp2=amp**2
      ap=2d0*amp
      ap2=2d0*amp2


      end



****************** bornkin *************************************

      subroutine bornkin(cutv)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      common/born/w,w2,q2,s,x,sx,sxp,sqlw0,ssh0,xxh0,sxtm0,sqll
     .	   ,cosh,phih,sinh,vmax,vv10,vv20,vvp0,vvm0,t0 !  /born/
     .	   ,cos1,sin1,cos2,sin2,ee1,ee2,pp1,pp2,eeq,ppq !  /born/

      w2=w**2
      sx=w2+q2-amp2
      x=s-sx
      sxp=s+x
      aly=sx**2+4.*amp2*q2
      sqly=dsqrt(aly)
      eeq=(sx-2.*q2)/2./w
      ppq=sqly/2./w

      vmax=w2-amu2+amh2-2.*sqrt(w2)*sqrt(amh2)	- 1d-8
      if(cutv.gt.1d-12)vmax=min(vmax,cutv)


      u1=s-q2
      u2=x+q2
      ee1=u1/2./w
      ee2=u2/2./w
      sql1=sqrt(u1**2-4.*w2*aml2)
      sql2=sqrt(u2**2-4.*w2*aml2)
      pp1=sql1/2./w
      pp2=sql2/2./w
      cos1=(u1*(sx-2.*q2)+2.*q2*w2)/sqly/sql1
      cos2=(u2*(sx-2.*q2)-2.*q2*w2)/sqly/sql2
      sql12=sqrt(q2*u1*u2-q2**2*w2-aml2*aly)
      sin1=2.*w*sql12/sql1/sqly
      sin2=2.*w*sql12/sql2/sqly
      sqll=sqrt(q2*u1*u2-q2**2*w2-aml2*aly)


      sinh=sqrt(1d0-cosh**2)
      eeh0=(w2+amh2-amu2)/2./w
      sqlw0=sqrt((w2+amh2-amu2)**2-4.*w2*amh2)
      pph0=sqlw0/2./w

      vv10=2d0*(ee1*eeh0-pp1*pph0*(cosh*cos1+sinh*sin1*cos(phih)))
      vv20=2d0*(ee2*eeh0-pp2*pph0*(cosh*cos2+sinh*sin2*cos(phih)))

      ssh0=x+q2-vv20
      xxh0=s-q2-vv10
      t0=-q2-vv10+vv20+amh2

      sxtm0=sx+t0+amp2-amu2

c     write(*,'(a8,2g11.3)')' test1: ',sin1**2+cos1**2
c     write(*,'(a8,2g11.3)')' test2: ',sin2**2+cos2**2
c     write(*,'(a8,2g11.3)')' test3: ',sinh**2+cosh**2
c     write(*,'(a8,2g11.3)')' test4: ',pph0**2+amh2,eeh0**2
c     write(*,'(a8,2g11.3)')' test5: ',pp1**2+aml2,ee1**2
c     write(*,'(a8,2g11.3)')' test6: ',pp2**2+aml2,ee2**2
c     write(*,'(a8,2g11.3)')' test7: ',ppq**2-eeq**2,q2
c     write(*,'(a8,2g11.3)')' test8: ',w2+amh2-2.*w*eeh0,amu2

      end

****************** borntest *************************************
      subroutine borntest(bmom,w,q2,csin,phih,sig0bu,sig0bp,eps_d)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      common/keys/ivec,iphy

	enu = 0.5d0*(w**2+q2-amp2)/amp
	bmom_prime = bmom - enu
	cthe = 1.d0-q2/(2.0*bmom*bmom_prime)
	the_e = acos(cthe)
	eps_d= 1d0/(1d0+2d0*(1d0+enu*enu/q2)*tan(0.5d0*the_e)**2)

	call phy_model(iphy,q2,w,csin
     .	 ,2*ivec-1,sig_t,sig_l,sig_tt,sig_lt,sig_ltp)

c	write(*,'(5g12.4)')q2,w,csin
c	 write(*,'(5g12.4)')sig_t,sig_l,sig_tt,sig_lt,sig_ltp


c	 sig_t=0.
c	 sig_l=0.
c	 sig_tt=0.
c	 sig_lt=0.

	sig0b = sig_t + eps_d*sig_l + eps_d*sig_tt*cos(2d0*phih)
     .	 + sqrt(eps_d*(eps_d+1d0)/2d0)*sig_lt*cos(phih)
	sig0b=sig0b * alpha/2d0/pi*bmom_prime/bmom*(w**2-amp2)
     .	 /2./amp/q2/(1d0-eps_d)

	sig0bu=sig0b*pi/2./amp/bmom/bmom_prime

	sig0b =
     .	 + sqrt(eps_d*(-eps_d+1d0)/2d0)*sig_ltp*sin(phih)
	sig0b=sig0b * alpha/2d0/pi*bmom_prime/bmom*(w**2-amp2)
     .	 /2./amp/q2/(1d0-eps_d)

	sig0bp=sig0b*pi/2./amp/bmom/bmom_prime

      end


****************** bornin *************************************

      subroutine bornin(sibor)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      common/born/w,w2,q2,s,x,sx,sxp,sqlw0,ssh0,xxh0,sxtm0,sqll
     .	   ,cosh,phih,sinh,vmax,vv10,vv20,vvp0,vvm0,t0	!  /born/
     .	   ,cos1,sin1,cos2,sin2,ee1,ee2,pp1,pp2,eeq,ppq !  /born/
      dimension tb(5),sfm0(5)
      common/isfisf/isf1,isf2

       call subtm0(5,tb)
       call sffun(sfm0,q2,w2,t0)

       ssu=0d0
       do isf=isf1,isf2
	  ssu=ssu+tb(isf)*sfm0(isf)
       enddo

       sibor=alpha**2*sqlw0/(32.*pi**2*w2*s**2*q2**2) * ssu

       end



      subroutine subtm0(ntm,tb)
      implicit real*8(a-h,o-z)
      dimension tb(ntm)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      common/born/w,w2,q2,s,x,sx,sxp,sqlw0,ssh0,xxh0,sxtm0,sqll
     .	   ,cosh,phih,sinh,vmax,vv10,vv20,vvp0,vvm0,t0	!  /born/
     .	   ,cos1,sin1,cos2,sin2,ee1,ee2,pp1,pp2,eeq,ppq !  /born/

      TB(1)=Q2
      TB(2)=(-AMP2*Q2+S*X)/2.
      TB(3)=(-AMH2*Q2+Vv10*Vv20)/2.
      TB(4)=(-Q2*SXTM0+S*Vv20+X*Vv10)/2.

      sqla=sqrt(q2*(s-q2)*(x+q2)-q2**2*W2)
      eps012=sqlw0*sqla/4./w*sinh*sin(phih)

c      tt=(-q2+amh2-vv10+vv20)
c      sxtm=s-x+tt+amp2-amu2
c      eee12=SQRT(-((vv10*x-q2*sxtm+s*vv20)**2
c     .        +4.*(s*x-amp2*q2)*(q2*amh2-vv10*vv20))/16.)
c
c	WRITE(*,*)EPS012,eee12

      TB(5)=-2.*eps012

      end

*************************** qqt *****************************
      subroutine qqt(tai)
      implicit real*8(a-h,o-z)
      external rv2tr
      real*8 rv2tr
      dimension am(3),bm(3),wrk(500)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      common/born/w,w2,q2,s,x,sx,sxp,sqlw0,ssh0,xxh0,sxtm0,sqll
     .     ,cosh,phih,sinh,vmax,vv10,vv20,vvp0,vvm0,t0  !  /born/
     .     ,cos1,sin1,cos2,sin2,ee1,ee2,pp1,pp2,eeq,ppq !  /born/
      dimension phiar(4),coar(6)

      am(1)=0d0
      bm(1)=vmax

      phiar(1)=0.
      phiar(2)=0.01*pi
      phiar(3)=2d0*pi-0.01*pi
      phiar(4)=2d0*pi

      coar(1)=-1d0
      coar(2)=cos2-0.05d0
      coar(3)=cos2+0.05d0
      coar(4)=cos1-0.05d0
      coar(5)=cos1+0.05d0
      coar(6)=1d0


      ma=100000
      ot=1d-2

      rere=0.
      do iph=1,3
      do ico=1,5

      am(2)=coar(ico)
      bm(2)=coar(ico+1)
      am(3)=phiar(iph)
      bm(3)=phiar(iph+1)

       if(am(2).gt.bm(2))write(*,*)' am(2)<bm(2)'

      id=1
      mir=10000
      call d01fce(3,am,bm,mir,ma,rv2tr,ot,otr,500,wrk,re,id)
      write(9,'(1x,''tai:'',2i3,g13.4,f8.4,i9,i3)')ico,iph,re,otr,mir,id
      write(*,'(1x,''tai:'',2i3,g13.4,f8.4,i9,i3)')ico,iph,re,otr,mir,id
      rere=rere+re
      enddo
      enddo

c     

      tai=-alpha**3/(512.*pi**4*s**2*w2**2) * rere

      write(9,'(1x,''tai:'',2g13.4,2i5)')tai,otr,mir,id
      write(*,'(1x,''tai:'',2g13.4,2i5)')tai,otr,mir,id

      end


****************** qqtll **************************************

      subroutine qqtll(taill)
      implicit real*8(a-h,o-z)
      external fvll
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      common/born/w,w2,q2,s,x,sx,sxp,sqlw0,ssh0,xxh0,sxtm0,sqll
     .	   ,cosh,phih,sinh,vmax,vv10,vv20,vvp0,vvm0,t0	!  /born/
     .	   ,cos1,sin1,cos2,sin2,ee1,ee2,pp1,pp2,eeq,ppq !  /born/

	call simpsx(1d-6,vmax,150,1d-3,fvll,taill)

      end

****************** fvll **************************************

      real*8 function fvll(v)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      common/born/w,w2,q2,s,x,sx,sxp,sqlw0,ssh0,xxh0,sxtm0,sqll
     .	   ,cosh,phih,sinh,vmax,vv10,vv20,vvp0,vvm0,t0 !  /born/
     .	   ,cos1,sin1,cos2,sin2,ee1,ee2,pp1,pp2,eeq,ppq !  /born/

      eeh=(w2+amh2-amu2-v)/2./w
      sqlw=sqrt((w2+amh2-amu2-v)**2-4.*w2*amh2)
      pph=sqlw/2./w
      vv1=2d0*(ee1*eeh-pp1*pph*(cosh*cos1+sinh*sin1*cos(phih)))
      vv2=2d0*(ee2*eeh-pp2*pph*(cosh*cos2+sinh*sin2*cos(phih)))

      tt=(-q2+amh2-vv1+vv2)

      sxtm=sx+tt-v+amp2-amu2


      xx=x+q2-vv2
      ss=s-q2-vv1

      z1=1d0-v/ss
      z2=1d0/(1d0+v/xx)

       sigmas=sigm5(z1*s,x,z1*q2,z1*vv1,vv2,sinh,phih)
       sigmax=sigm5(s,x/z2,q2/z2,vv1,vv2/z2,sinh,phih)
       sigma0=sigm5(s,x,q2,vv1,vv2,sinh,phih)

       w2s=z1*s-x-z1*q2+amp2
       w2x=s-x/z2-q2/z2+amp2


      sqlws0=sqrt((w2s+amh2-amu2)**2-4.*w2s*amh2)
      sqlwx0=sqrt((w2x+amh2-amu2)**2-4.*w2x*amh2)

       ccs=sqlw/sqlws0*w2s/w2
       ccx=sqlw/sqlwx0*w2x/w2

      sis=((1d0+z1**2)*ccs*sigmas-2d0*sigma0)/ss/(1d0-z1)
      six=((1d0+z2**2)*ccx*sigmax-2d0*sigma0)*z2/xx/(1d0-z2)

      fvll=alpha/2d0/pi*(log(Q2/aml2)-1d0)*(sis+six)

c     write(33,'(5g12.3)')z1,z2,v
c      pause

      end



       real*8 function sigm5(sz,xz,q2z,v1z,v2z,sinh,phih)
       implicit real*8(a-h,o-z)
       common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
       dimension tb(5),sfm0(5)
      common/isfisf/isf1,isf2

      s=sz
      x=xz
      q2=q2z
      vv1=v1z
      vv2=v2z


      w2=s-x-q2+amp2
      tt=(-q2+amh2-vv1+vv2)
      sqlw0=sqrt((w2+amh2-amu2)**2-4.*w2*amh2)

      sxtm=s-x+tt+amp2-amu2

       call sffun(sfm0,q2,w2,tt)

      TB(1)=Q2
      TB(2)=(-AMP2*Q2+S*X)/2.
      TB(3)=(-AMH2*Q2+Vv1*Vv2)/2.
      TB(4)=(-Q2*SXTM+S*Vv2+X*Vv1)/2.


c      eee12=sqrt(-((vv1*x-q2*sxtm+s*vv2)**2
c     .        +4.*(s*x-amp2*q2)*(q2*amh2-vv1*vv2))/16.)
c	rd12:=moph*sqll/2*sih*sin(phih)$


	  sqll=sqrt(q2*s*x-q2**2*amp2)
	  eee12=sqlw0*sqll/4./sqrt(w2)*sinh*sin(phih)


       tb(5)=-2*eee12
c	write(33,'(5g12.3)')eee12

       ssu=0d0
       do isf=isf1,isf2
	  ssu=ssu+tb(isf)*sfm0(isf)
c	  write(*,'(i5,5g12.4)')isf,ssu,tb(isf),sfm0(isf)
       enddo

       sigm5=alpha**2*sqlw0/(32.*pi**2*w2*s**2*q2**2) * ssu

c	 write(*,'(5g11.3)')w2,q2,tt,sigm5

       end






****************** rv2tr **************************************

      double precision function rv2tr(ndim,arg)
      implicit real*8(a-h,o-z)
      dimension arg(ndim),tm(5,3),tb(5),sfm(5),sfm0(5)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      common/born/w,w2,q2,s,x,sx,sxp,sqlw0,ssh0,xxh0,sxtm0,sqll
     .	   ,cosh,phih,sinh,vmax,vv10,vv20,vvp0,vvm0,t0	!  /born/
     .	   ,cos1,sin1,cos2,sin2,ee1,ee2,pp1,pp2,eeq,ppq !  /born/
      common/isfisf/isf1,isf2


      phik=arg(3)
      cosk=arg(2)
      vv=arg(1)

      sink2=1.d0-cosk**2
      sink=sqrt(max(0d0,sink2))

      z1=(ee1-pp1*(cos1*cosk+sin1*sink*cos(phik)))/w
      z2=(ee2-pp2*(cos2*cosk+sin2*sink*cos(phik)))/w
       bb=1.
       bi12=bb/(z1*z2)
       bi1pi2=bb/z2+bb/z1
       bis=(bb/z2**2+bb/z1**2)*aml2
       bir=(bb/z2**2-bb/z1**2)*aml2
       bi22m=bb/z2**2*aml2
       bi11m=bb/z1**2*aml2
       hi2=bis-q2*bi12

      eeh=(w2+amh2-amu2-vv)/2./w
      sqlw=sqrt((w2+amh2-amu2-vv)**2-4.*w2*amh2)
      pph=sqlw/2./w
      vv1=2d0*(ee1*eeh-pp1*pph*(cosh*cos1+sinh*sin1*cos(phih)))
      vv2=2d0*(ee2*eeh-pp2*pph*(cosh*cos2+sinh*sin2*cos(phih)))
      fwiw=1d0-(eeh-pph*(cosh*cosk+sinh*sink*cos(phih-phik)))/w

      ta=(eeq-ppq*cosk)/w
      dmu=1d0-fwiw
      r=vv/fwiw

      tt=-q2+amh2-vv1+vv2

      tldq2=q2+r*ta
      tldw2=w2-r
      tldtm=(-q2+amh2-vv1+vv2)-r*(ta-dmu)
      tldtm=tt-r*(ta-dmu)

c	write(22,'(a5,7g11.3)')' tri ',eeh/w,pph/w
c     .  ,w-eeh-r/2./w,sqrt(amu2)
c     .  ,fwiw

       call sffun(sfm,tldq2,tldw2,tldtm)
       call sffun(sfm0,q2,w2,t0)


      sxtm=sx+tt-vv+amp2-amu2

      vpl=vv1+vv2
      vmi=vv1-vv2

       call subtm0(5,tb)

c     TM(1,1)=4.*HI2*Q2
c     TM(1,2)=4.*HI2*TA
c     TM(1,3)=2.*(-2.*BB-BI12*TA**2)
c     TM(2,1)=(HI2*(-4.*AMP2*Q2+SXP**2-SX**2))/2.
c     TM(2,2)=(-4.*AMP2*HI2*TA-BI12*SXP**2*TA+BI1PI2*SXP*SX+2.*
c    . BIR*SXP+2.*HI2*SX)/2.
c     TM(2,3)=(4.*AMP2*BB+2.*AMP2*BI12*TA**2-BI12*SX*TA-BI1PI2*
c    . SXP)/2.
c     TM(3,1)=(HI2*(-4.*AMH2*Q2-VMI**2+VPL**2))/2.
c     TM(3,2)=(-4.*AMH2*HI2*TA-BI12*TA*VPL**2+BI1PI2*VMI*VPL+2.*
c    . BIR*DMU*VPL+2.*DMU*HI2*VMI)/2.
c     TM(3,3)=(4.*AMH2*BB+2.*AMH2*BI12*TA**2-BI12*DMU*TA*VMI-
c    . BI1PI2*DMU*VPL)/2.
c     TM(4,1)=HI2*(-2.*Q2*SXTM+SXP*VPL-SX*VMI)
c     TM(4,2)=(-2.*BI12*SXP*TA*VPL+BI1PI2*SXP*VMI+BI1PI2*SX*VPL+
c    . 2.*BIR*DMU*SXP+2.*BIR*VPL+2.*DMU*HI2*SX-4.*HI2*SXTM*TA+2.*
c    . HI2*VMI)/2.
c     TM(4,3)=(4.*BB*SXTM-BI12*DMU*SX*TA+2.*BI12*SXTM*TA**2-BI12
c    . *TA*VMI-BI1PI2*DMU*SXP-BI1PI2*VPL)/2.


      v1=vv1
      v2=vv2
      tm(1,1)=4.*hi2*q2
      tm(1,2)=4.*hi2*ta
      tm(1,3)=2.*(-2.*bb-bi12*ta**2)
      tm(2,1)=2.*hi2*(-amp2*q2+s*x)
      tm(2,2)=(-4.*amp2*hi2*ta-bi12*s**2*ta-2.*bi12*s*ta*x
     . -bi12*ta*x**2+bi1pi2*s**2-bi1pi2*x**2-2.*bir*s*ta+
     . 2.*bir*s-2.*bir*ta*x+2.*bir*x-2.*hi2*s*ta+2.*hi2*s+
     . 2.*hi2*ta*x-2.*hi2*x)/2.
      tm(2,3)=(4.*amp2*bb+2.*amp2*bi12*ta**2+bi12*s*ta**2-
     . bi12*s*ta-bi12*ta**2*x+bi12*ta*x+bi1pi2*s*ta-bi1pi2
     . *s+bi1pi2*ta*x-bi1pi2*x)/2.
      tm(3,1)=2.*hi2*(-amh2*q2+v1*v2)
      tm(3,2)=(-4.*amh2*hi2*ta-bi12*ta*v1**2-2.*bi12*ta*v1
     . *v2-bi12*ta*v2**2+bi1pi2*v1**2-bi1pi2*v2**2+2.*bir*
     . dmu*v1+2.*bir*dmu*v2+2.*dmu*hi2*v1-2.*dmu*hi2*v2)/
     . 2.
      tm(3,3)=(4.*amh2*bb+2.*amh2*bi12*ta**2-bi12*dmu*ta*
     . v1+bi12*dmu*ta*v2-bi1pi2*dmu*v1-bi1pi2*dmu*v2)/2.
      tm(4,1)=2.*hi2*(-q2*sxtm+s*v2+v1*x)
      tm(4,2)=-bi12*s*ta*v1-bi12*s*ta*v2-bi12*ta*v1*x-bi12
     . *ta*v2*x+bi1pi2*s*v1-bi1pi2*v2*x+bir*dmu*s+bir*dmu*
     . x-bir*ta*v1-bir*ta*v2+bir*v1+bir*v2+dmu*hi2*s-dmu*
     . hi2*x-2.*hi2*sxtm*ta-hi2*ta*v1+hi2*ta*v2+hi2*v1-hi2
     . *v2
      tm(4,3)=(4.*bb*sxtm-bi12*dmu*s*ta+bi12*dmu*ta*x+2.*
     . bi12*sxtm*ta**2+bi12*ta**2*v1-bi12*ta**2*v2-bi12*ta
     . *v1+bi12*ta*v2-bi1pi2*dmu*s-bi1pi2*dmu*x+bi1pi2*ta*
     . v1+bi1pi2*ta*v2-bi1pi2*v1-bi1pi2*v2)/2.



      eee12i=sqrt(-((vv1*x-q2*sxtm+s*vv2)**2
     .	      +4.*(s*x-amp2*q2)*(q2*amh2-vv1*vv2))/16.)

      eee2i=sqrt(- ((vv2*(ta-1.)-dmu*x+sxtm*z2)**2
     .	      +4.*(z2*amp2+(ta-1.)*x)*(-z2*amh2+vv2*dmu))/16.)

      eee1i=sqrt(- ((vv1*(ta-1.)-dmu*s+sxtm*z1)**2
     .	      +4.*(z1*amp2+(ta-1.)*s)*(-z1*amh2+vv1*dmu))/16.)


       eee12=pph*sqll*sinh*sin(phih)/2.
       om=r/2./w
       sqlq=2.*w*ppq
       bee1=(sx*s+ap2*q2)*sin(phik-phih)*pph*sinh*sink/(4.*w*sqlq)
       bee2=(sx*x-ap2*q2)*sin(phik-phih)*pph*sinh*sink/(4.*w*sqlq)
       beee=(sx+ap2)*(sin(phih)*cosk*sinh - sin(phik)*cosh*sink)
     .	   *pph*sqll/(4.*w2*sqlq)
       bbbb = sqll*(sin(phih)*pph*sinh - sin(phik)*eeh*sink)/(4.*w**2)
       eee1=beee+bbbb+bee1
       eee2=beee+bbbb+bee2

c
c	write(*,'(6g12.3)')eee12,eee12i,eee1/eee1i,eee2/eee2i
c	write(*,'(6g12.3)')eee1i/eee2i,eee1/eee2
c	write(*,'(6g12.3)')eee1i,eee2i,eee1,eee2

c bbbb := sqll*(sin(phih)*moph*sih - sin(phik)*eeh*sik)/(4*w**2)$
c rd12:=moph*sqll/2*sih*sin(phih)$
c sin(phih-phik):=sin(phih)*cos(phik)-sin(phik)*cos(phih);
c bee := (sx + 2*m1**2)*
c (sin(phih)*cok*sih - sin(phik)*coh*sik)
c *mok*moph*sqll/(2*w*sqlq)$
c be1 := (sx*s + 2*m1**2*q2)*sin(phik-phih)*mok*moph*sih*sik/(2*sqlq)$
c be2 := (sx*x - 2*m1**2*q2)*sin(phik-phih)*mok*moph*sih*sik/(2*sqlq)$
c rd1:=bee+be1;
c rd2:=bee+be2;

c      bi11m=0d0

       tm(5,1)=-8.*hi2*eee12
       tm(5,2)=4.*(2.*bi11m*eee2 + bi12*(eee12*ta-q2*(eee1+eee2)) +
     .	2.*bi22m*eee1+2.*bi11m*(eee1*s - eee12*ta + eee12 - eee2*s)/s)
       tm(5,3)=2.*( bi1pi2*(eee2-eee1)-bi12*(eee1+eee2)*ta
     .	+4.*bi11m*eee2*(ta-1.)/s)

      podinl=0.
      do  isf=isf1,isf2
      do  irr=1,3
	 pp=sfm(isf)*tm(isf,irr)*sqlw/tldq2**2
	 if(irr.eq.1)
c    .	 pp=pp-sqlw0*tm(isf,1)*sfm0(isf)/q2**2
     .	 pp=pp-4d0*sqlw0*hi2*tb(isf)*sfm0(isf)/q2**2
	 pres=pp*r**(irr-2)
	 podinl=podinl+pres
      if(isf.eq.5.and.r.lt.1d-8.and.irr.eq.1)
     . write(*,'(2i4,7g11.3)')isf,irr,pres,r,sfm(isf),tm(isf,irr)
     . ,tb(isf)*4*hi2,sfm0(isf)
      enddo
      enddo

      rv2tr=podinl/fwiw

c      write(*,'(7g11.3)')podinl,r,tldq2,sqlw
c      if(abs(rv2tr).gt.1d13)then
c      write(22,'(7g11.3)')rv2tr,r,cosk,phik,tldq2,tldw2,tldtm
c      endif
c     pause



      end




      subroutine sffun(sfm,q2,w2,t)
      implicit real*8(a-h,o-z)
      dimension sfm(5)
      common/keys/ivec,iphy
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn

      sqw2=sqrt(w2)
      iepi=2*ivec-1

      sx=w2+q2-amp2
      aly=sx**2+4.*amp2*q2
      sxt=sx+t+amp2-amu2
      tq=t+q2-amh2

      sqlw=sqrt((W2-amu2-amh2)**2-4.*amu2*amh2)
      sqll=sqrt(q2*sxt**2-sxt*sx*tq-amp2*tq**2-amh2*aly)

      csthcm=(2.*tq*w2+(sx-2.*q2)*(w2+amh2-amu2))/sqlw/sqrt(aly)

c      write(22,'(7g11.3)')q2,w2,t

      if(ivec.eq.1)then
	cspion=-csthcm
      elseif(ivec.eq.2)then
	cspion= csthcm
      endif

      call phy_model(iphy,q2,sqw2,cspion,iepi,st,sl,stt,slt,sltp)

cc	write(*,'(5g12.4)')q2,sqw2,cspion
c	write(*,'(5g12.4)')st,sl,stt,slt,sltp
c      write(22,'(7g11.3)')q2,sqw2,cspion,st,sl,stt,slt

c	stt=0.
c      slt=0.
c	st=0.
c	sl=0.


      if(ivec.eq.1)then
	 slt=-slt
	 sltp=-sltp
      endif

      sfm10=(st-stt)
      sfm20=4.*(st+sl)*q2/aly
      sfm2tl=2.*slt*sqrt(q2)*(-sx*tq+2*q2*sxt)/(aly*sqll)
      sfm4tl=-slt*sqrt(q2)/sqll
      sfm4tt=-2.*stt*(-sx*tq+2*q2*sxt)/sqll**2
      sfm3tt=2.*stt*aly/sqll**2
      sfm2tt=2.*stt*((-sx*tq+2*q2*sxt)**2-2*q2*sqll**2)/(aly*sqll**2)
      sfm5tl=-sltp*sqrt(q2)/sqll


       coetr=16.*pi**2*(w2-amp2)*w2/(alpha*sqlw)

       sfm(1)=coetr*sfm10
       sfm(2)=coetr*(sfm20+sfm2tl+sfm2tt)
       sfm(3)=coetr*sfm3tt
       sfm(4)=coetr*(sfm4tl+sfm4tt)
       sfm(5)=coetr*sfm5tl

c	 write(*,*)q2,w2,t,sfm
c	pause

      end



      subroutine phy_model(iphy,q2,w,csthcm,iepi,st,sl,stt,slt,sltp)
c
      real*8 q2,w,csthcm,st,sl,stt,slt,sltp
      integer iphy
c
      hq2=sngl(q2)
      hw =sngl(w)
      hcs=sngl(csthcm)

      if (iphy.eq.1) then
	call ao_model(hQ2,hW,0.5,hcs,0.0,iepi,
     *	     sigma0,sig_t,sig_tt,sig_l,sig_lt,sig_ltp)
      else
	call maid_lee(iphy,hQ2,hW,0.5,hcs,0.0,iepi,
     *	     sigma0,sig_t,sig_tt,sig_l,sig_lt,sig_ltp)
c	if (abs(sig_t).gt.50 .or. abs(sig_tt).gt.50 .or.
c    *	    abs(sig_l).gt.10 .or. abs(sig_lt).gt.10) then
c	  print *, hW, sig_t,sig_tt,sig_l,sig_lt
c	  sig_t = 0.0
c	  sig_tt = 0.0
c	  sig_l = 0.0
c	  sig_lt = 0.0
c	endif
      endif
      st=dble(sig_t)
      sl=dble(sig_l)
      stt=dble(sig_tt)
      slt=dble(sig_lt)
      sltp=dble(sig_ltp)

      end


      subroutine ao_model(Q2,W,epsilon,csthcm,phicm,epirea,
     *sigma0,sig_t,sig_tt,sig_l,sig_lt,sig_ltp)

      implicit none

      real Q2, W, epsilon, csthcm, phicm
      real sig_t, sig_l, sig_lt, sig_tt, sig_ltp
      real pi
      real epq2,epw,epeps,epcos,epphi,sigma0,sigu,sigt,sigl,sigi,sigip
      real tandd,fkt,sig0,sigmao
      real xmprot

***************************************************************************
*      AUTHOR:	      V. BURKERT AND Z. LI
*      FIRST VERSION: SUMMER, 1991.  RECENT UPDATE: MAY.1993
* AO IS WRITTEN BASED ON THE ORIGINAL PROGRAM A_AND_O FROM V.BURKERT
* THIS PROGRAM SHOULD BE LINKED TO EITHER QKMC OR QKXM FOR THE CALCULATION
* OF QUARK MODEL.
* QKMC AND QKXM USE DIFFERENT FORM FACTORS:
* QKMC USES THE TREATMENT OF FOSTER ET. AL.
* QKXM USES THE DIPOLE FORM FACTOR
* QKMC IS THE DEFAULT CHOICE IN THIS PACKAGE.
* AO CAN BE USED TO EXAM: Q2-DEPENDENCE OF HELICITY AMP. AT RES. POSITION (GO1)
*			  OUTPUT:GDH.TOP: A B CA(A_1/2) CB(A_3/2)
*			  GDH SUM RULE (GO2):OUTPUT GDH.TOP
*			  OBERVABLES  (RETURN):OUTPUT TEST.TOP
**************************************************************************
*	    AO.FOR
* AO CONTAINS THE FOLLOWING SUBROUTINES/FUNCTIONS:
* SIGMA--CALCULATES OBSERVABLES
* EPRES--CALCULATES BREIT-WIGNER RESONANCE AMPLITUDES
* EXPA --CALCULATES THE HELICITY AMPLITUTES FROM EXPT (V.BURKERT ORIGINAL)
* BORNT--CALCULATES THE BORN TERM CONTRIBUTIONS
* BACK --CALCULATES BORN AND NON-BORN BACK GROUND TERMSS
* QKMA --CALCULATES THE HELICITY AMPLITUTES FROM QUARK MODEL
* RAMPF --CALCULATES THE Q2-DEPENDENCE OF THE HELICITY AMP. AT RES. POSITIONS
* HAMP --CALCULATES THE ENERGY-DEPENDECE OF THE RESONANCE HELICITY AMPLITUDES
* QKMC --CALCULATES THE COUPLING CONSTANTS FROM QUARK MODEL
************************************************************************
**************************************************************************
* UPDTATED: NOV., 1992
*	  * WITH NEW OPTIONS TO TURN THE BORN TERM ON AND OFF
*	  * BORN TERM ARE MODIFIED WITH A CUT OFF FACTOR AT HIGHER
*	    ENERGIES (WCM>1.3 GEV)
*	  * SOME CORRECTIONS HAVE BEEN MADE TO THE EXPA SUBROUTINE
***************************************************************************

      COMMON /IP/IT,IB,NF,IBORN,CUT,IP11

c     COMMON /CA/CS11P,CP11P,CP333,CP331
C* CA HERE USED FOR ADJUSTING THE COUPLING CONSTANT CALCULATED
C* FROM QUARK MODEL

C* IT:		DETERMINE WHICH MODEL TO USE,
C* IB=1,0:	TURNS THE NON-BORN BACKGROUND ON AND OFF
C* IBORN=1,0:	TURNS THE BORN	TERM ON AND OFF
C* NF=1,3:	DETERMINE THE FORM OF P11_1440
C* IP11=1,4:	DETERMINES THE USE OF P11_1440 AMPLITUDE IN EXPT

      real cut
      integer it,ib,nf,iborn,ip11,epirea

       real
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710,
     * PIBP13_1720, PIBD13_1520, PIBD13_1700,
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600,
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975

       COMMON/RCONST/
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710,
     * PIBP13_1720, PIBD13_1520, PIBD13_1700,
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600,
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975
C
C      Branching ratios into N-pion channel
C
      DATA PIBS11_1535	    /0.40    /
      DATA PIBS11_1650	    /0.60    /
      DATA PIBP11_1440	    /0.60    /
      DATA PIBP11_1710	    /0.15    /
      DATA PIBP13_1720	    /0.15    /
      DATA PIBD13_1520	    /0.55    /
      DATA PIBD13_1700	    /0.10    /
      DATA PIBD15_1675	    /0.40    /
      DATA PIBF15_1680	    /0.60    /
      DATA PIBG17_2190	    /0.14    /
      DATA PIBG19_2250	    /0.10    /
      DATA PIBH19_2220	    /0.18    /
      DATA PIBI111_2600     /0.05    /

      DATA PIBS31_1620	    /0.30    /
      DATA PIBS31_1900	    /0.10    /
      DATA PIBP31_1910	    /0.20    /
      DATA PIBP33_1232	    /0.994   /
      DATA PIBP33_1920	    /0.20    /
      DATA PIBD33_1700	    /0.15    /
      DATA PIBD35_1930	    /0.10    /
      DATA PIBF35_1905	    /0.10    /
      DATA PIBF37_1950	    /0.40    /
      DATA PIBH311_2420     /0.10    /

      DATA PIBP33_1600	    /0.19    /
      DATA PIBF17_1990	    /0.05    /
      DATA PIBF15_2000	    /0.12    /
      DATA PIBP11_2100	    /0.11    /
      DATA PIBF35_2000	    /0.02    /

      DATA PIBP13_1870	    /0.02   /
      DATA PIBP31_1925	    /0.02   /
      DATA PIBP13_1980	    /0.02   /
      DATA PIBF15_1955	    /0.02   /
      DATA PIBP13_1955	    /0.02   /
      DATA PIBP33_1975	    /0.02   /

      DATA LS11_1535	  /1.532    /
      DATA LS11_1650	  /1.650    /
      DATA LP11_1440	  /1.440    /
      DATA LP11_1710	  /1.710    /
      DATA LP13_1720	  /1.740    /
      DATA LD13_1520	  /1.520    /
      DATA LD13_1700	  /1.700    /
      DATA LD15_1675	  /1.675    /
      DATA LF15_1680	  /1.688    /
      DATA LG17_2190	  /2.175    /
      DATA LG19_2250	  /2.250    /
      DATA LH19_2220	  /2.220    /
      DATA LI111_2600	  /2.640    /
      DATA LS31_1620	  /1.620    /
      DATA LS31_1900	  /1.925    /
      DATA LP31_1910	  /1.905    /
      DATA LP33_1232	  /1.232    /
      DATA LP33_1920	  /1.920    /
      DATA LD33_1700	  /1.670    /
      DATA LD35_1930	  /1.930    /
      DATA LF35_1905	  /1.905    /
      DATA LF37_1950	  /1.950    /
      DATA LH311_2420	  /2.420    /
      DATA LP33_1600	  /1.600   /
      DATA LF17_1990	  /1.990   /
      DATA LF15_2000	  /2.000   /
      DATA LP11_2100	  /2.100   /
      DATA LF35_2000	  /2.000   /
      DATA LP13_1870	  /1.870   /
      DATA LP31_1925	  /1.925   /
      DATA LP13_1980	  /1.980   /
      DATA LF15_1955	  /1.955   /
      DATA LP13_1955	  /1.955   /
      DATA LP33_1975	  /1.975   /

      DATA WS11_1535	  /.150     /
      DATA WS11_1650	  /.150     /
      DATA WP11_1440	  /.300     /
      DATA WP11_1710	  /.110     /
      DATA WP13_1720	  /.200     /
      DATA WD13_1520	  /.125     /
      DATA WD13_1700	  /.100     /
      DATA WD15_1675	  /.155     /
      DATA WF15_1680	  /.125     /
      DATA WG17_2190	  /.350     /
      DATA WG19_2250	  /.300     /
      DATA WH19_2220	  /.400     /
      DATA WI111_2600	  /.400     /
      DATA WS31_1620	  /.140     /
      DATA WS31_1900	  /.150     /
      DATA WP31_1910	  /.220     /
      DATA WP33_1232	  /.115     /
      DATA WP33_1920	  /.250     /
      DATA WD33_1700	  /.250     /
      DATA WD35_1930	  /.250     /
      DATA WF35_1905	  /.300     /
      DATA WF37_1950	  /.240     /
      DATA WH311_2420	  /.300     /
      DATA WP33_1600	  /.230     /
      DATA WF17_1990	  /.300     /
      DATA WF15_2000	  /.135     /
      DATA WP11_2100	  /.230     /
      DATA WF35_2000	  /.300     /

      DATA WP13_1870	  /.300     /
      DATA WP31_1925	  /.300     /
      DATA WP13_1980	  /.300     /
      DATA WF15_1955	  /.300     /
      DATA WP13_1955	  /.300     /
      DATA WP33_1975	  /.300     /
      DATA XMPROT     /.938	/
c      DATA XMPIP      /.1395	 /
c      DATA XMPI0      /.1349	 /
c      DATA XMETA      /.5488	 /

      pi	= 4.*atan(1.)
      nf	= 1
      it	= 1
      ip11	= 1
      ib	= 0
      iborn	= 1
      cut	= 0.4
      sigma0	= 0.
c
      EPW   = W
      EPQ2  = Q2
      EPEPS = epsilon
      EPCOS = csthcm
      EPPHI = phicm
c
      sigma0 = SIGMAO(EPW,EPQ2,EPEPS,EPCOS,EPPHI,EPIREA,
     *	 0.,0.,0.,0.,sig0,sigu,sigt,sigl,sigi,sigip,fkt)
c
      sig_L  = sigl
      sig_T  = sigu
      sig_LT = sigi
      sig_TT = sigt
      sig_LTP = sigip
c
      end
C**************************************************************************
	 COMPLEX FUNCTION EPRES(IREA,XMP,XMPI,W,W0,GAM,X,L,J,ISO,
     *	 PIMP,GENER,TSTS11)
	 COMPLEX QETA,GAMMAT,Q0ETA
	 INTEGER L,J,ISO,IREA
	 REAL XMP,XMPI,W,W0,GAM,X,PIMP0
	 REAL PENER0,GENER0,XMETA,GAMMA
	 LOGICAL TSTS11
	 COMPLEX IM
	 DATA IM/(0.,1.)/
	 DATA XMETA/.5488/

	 PENER0=(W0**2-XMP**2+XMPI**2)/(2.*W0)
	 PIMP0 = SQRT ( PENER0**2 - XMPI**2 )
	 GENER0 = (W0**2-XMP**2)/(2.*W0)
	 if (PIMP0 .le. 0.)then
	    write(6,*)' epres: pimp0 =',PIMP0
	    PIMP0=abs(PIMP0)
	 endif

	 GAMMA=GAM * ((PIMP/PIMP0)**(2.*L+1))*
     *	 (((PIMP0**2+X**2)/(PIMP**2+X**2))**L)

	 GAMMAG=GAM * ((GENER/GENER0)**(2.*J))*
     *	 (((GENER0**2+X**2)/(GENER**2+X**2))**J)

	 IF(TSTS11)THEN
	   XETA=(W**2-XMP**2+XMETA**2)/(2.*W)
	   XETA=XETA**2-XMETA**2
	   IF(XETA.GE.0.)THEN
	  QETA=CMPLX(SQRT(XETA),0.)
	   ELSE
	     QETA=CMPLX(0.,SQRT(-XETA))
	   END IF
	  Q0ETA=(W0**2-XMP**2+XMETA**2)/(2.*W0)
	  Q0ETA=SQRT(Q0ETA**2-XMETA**2)
	  GAMMAT=GAM*(.55*(QETA/Q0ETA)+.4*(PIMP/PIMP0)+.05)
C	     GAMMAT=GAM*(PIMP/PIMP0)
	   ELSE
	    GAMMAT=CMPLX(GAMMA,0.)
	  END IF
	 EPRES=SQRT((GENER0*PIMP0)/(GENER*PIMP))
	 EPRES=EPRES*W0*SQRT(GAMMA*GAMMAG)
	 EPRES=EPRES/(W0**2-W**2-IM*W0*GAMMAT)
	 IF(IREA.EQ.1.OR.IREA.EQ.2)THEN
	   IF(ISO.EQ.1)THEN
	     EPRES=-EPRES/SQRT(2.)
	   ELSE
	     EPRES=EPRES*SQRT(2.)
	   END IF
	 END IF
99	  RETURN
	 END

***************************************************************************
      SUBROUTINE EXPA(EPIREA,EPQ2)
      implicit none
c
c     Amplitudes from Experiment
c

      COMMON /IP/IT,IB,NF,IBORN,CUT,IP11

      real cut
      real pi,xmprot,xmpip,xmpi0,xmpi,epw
      real elchar,epq2,qme11,qmm11,qmm12,qme22,qmm21,qmm22,qmm23
      real fact,cpin,pvec_gam,pvec_pi,pvec_gam0,q2evf,q2evf0,dip_evf
      real q0_gam,ffac
      real cosdd,sindd

      integer it,ib,nf,iborn,ip11
      INTEGER EPIREA


       real
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710,
     * PIBP13_1720, PIBD13_1520, PIBD13_1700,
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600,
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975

       COMMON/RCONST/
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710,
     * PIBP13_1720, PIBD13_1520, PIBD13_1700,
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600,
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975

	 COMMON/RAMP/
     * RAS11_1535,	       RCS11_1535,
     * RAS11_1650,	       RCS11_1650,
     * RAP11_1440,	       RCP11_1440,
     * RAP11_1710,	       RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,
     * RAP13_1910, RBP13_1910, RCP13_1910,
     * RAD13_1520, RBD13_1520, RCD13_1520,
     * RAD13_1700, RBD13_1700, RCD13_1700,
     * RAD15_1675, RBD15_1675, RCD15_1675,
     * RAF15_1680, RBF15_1680, RCF15_1680,
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,	       RCS31_1620,
     * RAS31_1900,	       RCS31_1900,
     * RAP31_1910,	       RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232,
     * RAP33_1920, RBP33_1920, RCP33_1920,
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1975, RCP31_1975,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975

	 real
     * RAS11_1535,	       RCS11_1535,
     * RAS11_1650,	       RCS11_1650,
     * RAP11_1440,	       RCP11_1440,
     * RAP11_1710,	       RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,
     * RAP13_1910, RBP13_1910, RCP13_1910,
     * RAD13_1520, RBD13_1520, RCD13_1520,
     * RAD13_1700, RBD13_1700, RCD13_1700,
     * RAD15_1675, RBD15_1675, RCD15_1675,
     * RAF15_1680, RBF15_1680, RCF15_1680,
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,	       RCS31_1620,
     * RAS31_1900,	       RCS31_1900,
     * RAP31_1910,	       RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232,
     * RAP33_1920, RBP33_1920, RCP33_1920,
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1975, RCP31_1975,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975


      REAL A12S11_1535, 	     S12S11_1535
      REAL A12S11_1650, 	     S12S11_1650
      REAL A12P11_1440, 	     S12P11_1440
      REAL A12P11_1710, 	     S12P11_1710
      REAL A12P13_1720, A32P13_1720, S12P13_1720
      REAL A12P13_1910, A32P13_1910, S12P13_1910
      REAL A12D13_1520, A32D13_1520, S12D13_1520
      REAL A12D13_1700, A32D13_1700, S12D13_1700
      REAL A12D15_1675, A32D15_1675, S12D15_1675
      REAL A12F15_1680, A32F15_1680, S12F15_1680
      REAL A12G17_2190, A32G17_2190, S12G17_2190
      REAL A12G19_2250, A32G19_2250, S12G19_2250
      REAL A12H19_2220, A32H19_2220, S12H19_2220
      REAL A12I111_2600,A32I111_2600,S12I111_2600
      REAL A12S31_1620, 	     S12S31_1620
      REAL A12S31_1900, 	     S12S31_1900
      REAL A12P31_1910, 	     S12P31_1910
      REAL A12P33_1232, A32P33_1232, S12P33_1232
      REAL A12P33_1920, A32P33_1920, S12P33_1920
      REAL A12D33_1700, A32D33_1700, S12D33_1700
      REAL A12D35_1930, A32D35_1930, S12D35_1930
      REAL A12F35_1905, A32F35_1905, S12F35_1905
      REAL A12F37_1950, A32F37_1950, S12F37_1950
      REAL A12H311_2420,A32H311_2420,S12H311_2420
      REAL THETA_S11
c      REAL THETA_D13
c      REAL THETA_D15

      DATA A12S11_1535, 	     S12S11_1535  /0.051  ,0.	/
      DATA A12S11_1650, 	     S12S11_1650  /0.061  ,0.	/
      DATA A12P11_1440, 	     S12P11_1440  /-0.071 ,0. /
      DATA A12P11_1710, 	     S12P11_1710  /-0.039  ,0.	 /
      DATA A12P13_1720, A32P13_1720, S12P13_1720  /0.029  ,-0.034 ,0./
      DATA A12P13_1910, A32P13_1910, S12P13_1910  /0.  ,0.   ,0./
      DATA A12D13_1520, A32D13_1520, S12D13_1520  /-0.021  ,0.164 ,0./
      DATA A12D13_1700, A32D13_1700, S12D13_1700  /0.  ,0.   ,0./
      DATA A12D15_1675, A32D15_1675, S12D15_1675  /0.004  ,0.025  ,0./
      DATA A12F15_1680, A32F15_1680, S12F15_1680  /-0.016  ,0.146 ,0./
      DATA A12G17_2190, A32G17_2190, S12G17_2190  /0.  ,0.   ,0./
      DATA A12G19_2250, A32G19_2250, S12G19_2250  /0.  ,0.   ,0./
      DATA A12H19_2220, A32H19_2220, S12H19_2220  /0.  ,0.   ,0./
      DATA A12I111_2600,A32I111_2600,S12I111_2600 /0.  ,0.   ,0./
      DATA A12S31_1620, 	     S12S31_1620  /0.040  ,0.	/
      DATA A12S31_1900, 	     S12S31_1900  /0.  ,0.   /
      DATA A12P31_1910, 	     S12P31_1910  /0.033  ,0.	/
      DATA A12P33_1232, A32P33_1232, S12P33_1232  /-0.143,-0.259,0./
      DATA A12P33_1920, A32P33_1920, S12P33_1920  /0.  ,0.   ,0.   /
      DATA A12D33_1700, A32D33_1700, S12D33_1700  /0.110,0.093 ,0. /
      DATA A12D35_1930, A32D35_1930, S12D35_1930  /-0.021,0.004,0. /
      DATA A12F35_1905, A32F35_1905, S12F35_1905  /0.049 ,-0.011,0./
      DATA A12F37_1950, A32F37_1950, S12F37_1950  /-0.090,-0.113,0./
      DATA A12H311_2420,A32H311_2420,S12H311_2420 /0.  ,0.   ,0.   /

      DATA THETA_S11 / 38. /
c      DATA THETA_D13 /  0. /
c      DATA THETA_D15 /  0. /


c      DATA SQ2        /1.41421  /
      DATA PI	      /3.14159	/
      DATA XMPROT     /.938	/
      DATA XMPIP      /.1395	/
      DATA XMPI0      /.1349	/
c      DATA XMETA      /.5488	 /




C	Now calculate the resonant helicity amplitudes
C
	 IF(EPIREA.EQ.1.OR.EPIREA.EQ.2)XMPI=XMPI0
	 IF(EPIREA.EQ.3.OR.EPIREA.EQ.4)XMPI=XMPIP
	 ELCHAR = SQRT(4.*PI/137.)

	IF (EPIREA.EQ.1.OR.EPIREA.EQ.3) THEN
	CALL QKEM(LS11_1535,EPQ2,QME11,QMM11,QMM12)
	A12S11_1535 = ELCHAR/2./SQRT(XMPROT*(LS11_1535**2-XMPROT**2))*
     *	     (SQRT(1./3.)*QME11 - SQRT(2./3.)*QMM11)*COSDD(THETA_S11)
	CALL QKEM(LS11_1650,EPQ2,QME11,QMM11,QMM12)
	A12S11_1650 = ELCHAR/2./SQRT(XMPROT*(LS11_1650**2-XMPROT**2)) *
     *	      (SQRT(1./3.)*QME11 - SQRT(2./3.)*QMM11)*SINDD(THETA_S11)

C
C	 Assume no mixing for D13 and D15
C
	CALL QKEM(LD13_1520,EPQ2,QME11,QMM11,QMM12)
	A12D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *	 (SQRT(1./6.)*QME11 + SQRT(1./12.)*QMM11 - SQRT(3./4.)*QMM12)
	A32D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *		 (SQRT(1./2.)*QME11 + 0.5*QMM11 + 0.5*QMM12)

	 A12D13_1700 = 0.
	 A32D13_1700 = 0.

	 A12D15_1675 = 0.
	 A32D15_1675 = 0.

	CALL QKEM(LD33_1700,EPQ2,QME11,QMM11,QMM12)
	A12D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *	 (-SQRT(1./12.) * QMM12 - SQRT(1./6.) * QME11 +
     *	  0.5*SQRT(1./27.)*QMM11)
	A32D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *	      (1./6.*QMM12 - SQRT(1./2.)*QME11 + 1./6.*QMM11)

	CALL QKEM(LS31_1620,EPQ2,QME11,QMM11,QMM12)
	A12S31_1620 = ELCHAR/2./SQRT(XMPROT*(LS31_1620**2-XMPROT**2))*
     *	      ( SQRT(1./3.)*QME11 + SQRT(2./27.)*QMM11)

C*  SET ALL TO ZERO EXCEPT S11(1535) AND D13(1520)
C*	   A12S11_1650=0.
C*	   A12S31_1620=0.
C*	   A12D33_1700=0.
C*	   A32D33_1700=0.


C
C	  Now do the {56,0+} ---> {56,2+} transition.
C	  Here we use the quark multipole moments of Cottingham & Dunbar
C	  directly. Note thet quark multipole QMM21 is assumed = 0.
C	  This assumption has no justification (data do not allow its
C	  determination).
C
	  CALL QKEM2(LF15_1680,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LF15_1680**2 - XMPROT**2))
	  A12F15_1680 = FACT*(-SQRT(2./3.)*QMM23 + SQRT(1./5.)*QME22+
     *			SQRT(2./15.)*QMM22)
	  A32F15_1680 = FACT*(SQRT(1./3.)*QMM23 + SQRT(2./5.)*QME22+
     *			SQRT(4./15.)*QMM22)

	  CALL QKEM2(LP13_1720,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LP13_1720**2 - XMPROT**2))
	  A12P13_1720 = -FACT*(SQRT(3./4.)*QMM21 + SQRT(1./10.)*QME22-
     *			 SQRT(3./20.)*QMM22)
	  A32P13_1720 =  FACT*(0.5*QMM21 - SQRT(3./10.)*QME22 +
     *			 SQRT(9./20.)*QMM22)

	  CALL QKEM2(LF35_1905,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LF35_1905**2 - XMPROT**2))
	  A12F35_1905 =  FACT*(SQRT(16./189.)*QMM23 +
     *			SQRT(28./135.)*QMM22)
	  A32F35_1905 = -FACT*(SQRT(8./189.)*QMM23 -
     *			SQRT(56./135.)*QMM22)

	  CALL QKEM2(LP33_1920,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LP33_1920**2 - XMPROT**2))
	  A12P33_1920 = -FACT*(SQRT(3./15.)*QMM22 + 1./3.*QMM21)
	  A32P33_1920 =  FACT*(SQRT(1./15.)*QMM22 - SQRT(1./3.)*QMM21)

	  CALL QKEM2(LF37_1950,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LF37_1950**2 - XMPROT**2))
	  A12F37_1950 = -FACT*SQRT(2./7.)*QMM23
	  A32F37_1950 = -FACT*SQRT(10./21.)*QMM23
cccc
cccc	    multiply A12F37_1950 and A32F37_1950 by factor 2 (temporarily).
cccc
	  A12F37_1950 = A12F37_1950 * 2.
	  A32F37_1950 = A32F37_1950 * 2.
cccc
cccc
cccc
	  CALL QKEM2(LP31_1910,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LP31_1910**2 - XMPROT**2))
	  A12P31_1910 = -FACT*2./3.*QMM21

C*   SET ALL STATES ZERO EXCEPT F15(1680)
C	   A12P13_1720=0.
C	   A32P13_1720=0.
C	   A12F37_1950=0.
C	   A32F37_1950=0.
C	   A12F35_1905=0.
C	   A32F35_1905=0.
C	   A12P33_1920=0.
C	   A32P33_1920=0.
C	   A12P31_1910=0.
C
	  ELSE
C* THE FOLLOWING IS FOR NEUTRON TARGET
	A12P11_1440=0.056

	CALL QKEM(LS11_1535,EPQ2,QME11,QMM11,QMM12)
	A12S11_1535 = ELCHAR/2./SQRT(XMPROT*(LS11_1535**2-XMPROT**2))*
     *	     ((-SQRT(1./3.)*QME11 +SQRT(2./27.)*QMM11)*COSDD(THETA_S11)
     *	    +SINDD(THETA_S11)*(-SQRT(2./27.)*QMM11))
	CALL QKEM(LS11_1650,EPQ2,QME11,QMM11,QMM12)
	A12S11_1650 = ELCHAR/2./SQRT(XMPROT*(LS11_1650**2-XMPROT**2)) *
     *	      ((-SQRT(1./3.)*QME11 +SQRT(2./27.)*QMM11)*
     *	      (SINDD(THETA_S11))
     *	    +COSDD(THETA_S11)*(-SQRT(2./27.)*QMM11))
C
C	 Assume no mixing for D13 and D15
C
	CALL QKEM(LD13_1520,EPQ2,QME11,QMM11,QMM12)
	A12D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *	 (-SQRT(1./6.)*QME11-SQRT(1./108.)*QMM11 +SQRT(1./12.)*QMM12)
	A32D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *		 (-SQRT(1./2.)*QME11 - 1./6.*QMM11 -1./6.*QMM12)


	CALL QKEM(LD13_1700,EPQ2,QME11,QMM11,QMM12)
	 A12D13_1700 =ELCHAR/2./SQRT(XMPROT*(LD13_1700**2-XMPROT**2))*
     *	  (-SQRT(5./54.)*QMM11-SQRT(1./120.)*QMM12)
	 A32D13_1700 =ELCHAR/2./SQRT(XMPROT*(LD13_1700**2-XMPROT**2))*
     *	  (-SQRT(5./18.)*QMM11+SQRT(1./360.)*QMM12)

	CALL QKEM(LD15_1675,EPQ2,QME11,QMM11,QMM12)
	 A12D15_1675 =ELCHAR/2./SQRT(XMPROT*(LD15_1675**2-XMPROT**2))*
     *	  (-SQRT(2./15.)*QMM12)
	 A32D15_1675 =ELCHAR/2./SQRT(XMPROT*(LD15_1675**2-XMPROT**2))*
     *	   (-SQRT(4./15.)*QMM12)
	CALL QKEM(LD33_1700,EPQ2,QME11,QMM11,QMM12)
	A12D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *	 (-SQRT(1./12.) * QMM12 - SQRT(1./6.) * QME11 +
     *	  0.5*SQRT(1./27.)*QMM11)
	A32D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *	      (1./6.*QMM12 - SQRT(1./2.)*QME11 + 1./6.*QMM11)

	CALL QKEM(LS31_1620,EPQ2,QME11,QMM11,QMM12)
	A12S31_1620 = ELCHAR/2./SQRT(XMPROT*(LS31_1620**2-XMPROT**2))*
     *	      ( SQRT(1./3.)*QME11 + SQRT(2./27.)*QMM11)

C*  SET ALL TO ZERO EXCEPT S11(1535) AND D13(1520)
C*	   A12S11_1650=0.
C*	   A12S31_1620=0.
C*	   A12D33_1700=0.
C*	   A32D33_1700=0.


C
C	  Now do the {56,0+} ---> {56,2+} transition.
C	  Here we use the quark multipole moments of Cottingham & Dunbar
C	  directly. Note thet quark multipole QMM21 is assumed = 0.
C	  This assumption has no justification (data do not allow its
C	  determination).
C

	  CALL QKEM2(LF15_1680,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LF15_1680**2 - XMPROT**2))
	  A12F15_1680 = FACT*(SQRT(8./27.)*QMM23-SQRT(8./135.)*QMM22)
	  A32F15_1680 = FACT*(-SQRT(4./27.)*QMM23-SQRT(16./135.)*QMM22)

	  CALL QKEM2(LP13_1720,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LP13_1720**2 - XMPROT**2))
	  A12P13_1720 = -FACT*(-SQRT(1./3.)*QMM21 + SQRT(1./15.)*QMM22)
	  A32P13_1720 = -FACT*(1./3.*QMM21+SQRT(1./5.)*QMM22)

	  CALL QKEM2(LF35_1905,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LF35_1905**2 - XMPROT**2))
	  A12F35_1905 =  FACT*(SQRT(16./189.)*QMM23 +
     *			SQRT(28./135.)*QMM22)
	  A32F35_1905 = -FACT*(SQRT(8./189.)*QMM23 -
     *			SQRT(56./135.)*QMM22)

	  CALL QKEM2(LP33_1920,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LP33_1920**2 - XMPROT**2))
	  A12P33_1920 = -FACT*(SQRT(3./15.)*QMM22 + 1./3.*QMM21)
	  A32P33_1920 =  FACT*(SQRT(1./15.)*QMM22 - SQRT(1./3.)*QMM21)

	  CALL QKEM2(LF37_1950,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LF37_1950**2 - XMPROT**2))
	  A12F37_1950 = -FACT*SQRT(2./7.)*QMM23
	  A32F37_1950 = -FACT*SQRT(10./21.)*QMM23
cccc
cccc	    multiply A12F37_1950 and A32F37_1950 by factor 2 (temporarily !!!)
cccc	    to bring them in quantitative agreement with photoproduction
cccc	    analysis
cccc
	  A12F37_1950 = A12F37_1950 * 2.
	  A32F37_1950 = A32F37_1950 * 2.
cccc
cccc
cccc
	  CALL QKEM2(LP31_1910,EPQ2,QME22,QMM21,QMM22,QMM23)
	  FACT = ELCHAR/2./SQRT(XMPROT*(LP31_1910**2 - XMPROT**2))
	  A12P31_1910 = -FACT*2./3.*QMM21

	  END IF

C
C	Now calculate the resonant partial wave helicity elements
C

	  CPIN=SQRT(2./3.)
	  IF(EPIREA.EQ.1.OR.EPIREA.EQ.3)CPIN=-CPIN
98	   PVEC_GAM = SQRT((EPQ2 + (LS11_1535+XMPROT)**2)*
     *	 (EPQ2+(LS11_1535-XMPROT)**2))/(2.*LS11_1535)
	 PVEC_PI= SQRT((LS11_1535**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LS11_1535**2) - XMPI**2)
	 FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LS11_1535*PIBS11_1535/WS11_1535)
	 RAS11_1535 =-FACT*A12S11_1535*19.73*CPIN
	 RCS11_1535 = 0.

	 PVEC_GAM = SQRT((EPQ2 + (LS11_1650+XMPROT)**2)*
     *	 (EPQ2+(LS11_1650-XMPROT)**2))/(2.*LS11_1650)
	 PVEC_PI= SQRT((LS11_1650**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LS11_1650**2) - XMPI**2)
	 FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LS11_1650*PIBS11_1650/WS11_1650)
	 RAS11_1650 = -FACT*A12S11_1650*19.73*CPIN
	 RCS11_1650 = 0.

	 PVEC_GAM = SQRT((EPQ2 + (LP11_1440+XMPROT)**2)*
     *	 (EPQ2+(LP11_1440-XMPROT)**2))/(2.*LP11_1440)
	 PVEC_GAM0 = SQRT(((LP11_1440+XMPROT)**2)*
     *	 ((LP11_1440-XMPROT)**2))/(2.*LP11_1440)
	 PVEC_PI= SQRT((LP11_1440**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LP11_1440**2) - XMPI**2)
	 FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LP11_1440*PIBP11_1440/WP11_1440)
	 EPW=LP11_1440
	 Q2EVF=((EPW**2-XMPROT**2)**2+EPQ2*(EPW+XMPROT)**2)/
     *	 (4.*EPW*XMPROT)
	 Q2EVF0=((EPW**2-XMPROT**2)**2)/(4.*EPW*XMPROT)
	 DIP_EVF = (1.+ Q2EVF/0.71)**(-2.)

	 RAP11_1440 = FACT*A12P11_1440*(Q2EVF/Q2EVF0)*
     *	 DIP_EVF*19.73*CPIN
	 RCP11_1440 =		    0.0000D+00

	 RAP11_1710 =		    0.0000D+00
	 RCP11_1710 =		    0.0000D+00

	 PVEC_GAM = SQRT((EPQ2 + (LP13_1720+XMPROT)**2)*
     *	 (EPQ2+(LP13_1720-XMPROT)**2))/(2.*LP13_1720)
	 PVEC_PI= SQRT((LP13_1720**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LP13_1720**2) - XMPI**2)
	 FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LP13_1720*PIBP13_1720/WP13_1720)
	 RAP13_1720 = - FACT*A12P13_1720*19.73*CPIN
	 RBP13_1720 = FACT*SQRT(16./12.)*A32P13_1720*19.73*
     *		       CPIN
	 RCP13_1720 =		    0.0000D+00

	 PVEC_GAM = SQRT((EPQ2 + (LD13_1520+XMPROT)**2)*
     *	 (EPQ2+(LD13_1520-XMPROT)**2))/(2.*LD13_1520)
	 PVEC_PI= SQRT((LD13_1520**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LD13_1520**2) - XMPI**2)
	 FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LD13_1520*PIBD13_1520/WD13_1520)
	 RAD13_1520 = FACT*A12D13_1520*19.73*CPIN
	 RBD13_1520 =  -FACT*SQRT(16./12.)*A32D13_1520*19.73*
     *		       CPIN
	 RCD13_1520 =		    0.0000D+00

	 PVEC_GAM = SQRT((EPQ2 + (LD13_1700+XMPROT)**2)*
     *	 (EPQ2+(LD13_1700-XMPROT)**2))/(2.*LD13_1700)
	 PVEC_PI= SQRT((LD13_1700**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LD13_1700**2) - XMPI**2)
	 FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LD13_1700*PIBD13_1700/WD13_1700)
	 RAD13_1700 = FACT*A12D13_1700*19.73*CPIN
	 RBD13_1700 = -FACT*SQRT(16./12.)*A32D13_1700*19.73*CPIN
	 RCD13_1700 =		    0.0000D+00

	 PVEC_GAM = SQRT((EPQ2 + (LD15_1675+XMPROT)**2)*
     *	 (EPQ2+(LD15_1675-XMPROT)**2))/(2.*LD15_1675)
	 PVEC_PI= SQRT((LD15_1675**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LD15_1675**2) - XMPI**2)
	 FACT = SQRT(1./6./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LD15_1675*PIBD15_1675/WD15_1675)
	 RAD15_1675 =-FACT*A12D15_1675*19.73*CPIN
	 RBD15_1675 = FACT*SQRT(16./32.)*A32D15_1675*19.73*CPIN
	 RCD15_1675 =		    0.0000D+00

	 PVEC_GAM = SQRT((EPQ2 + (LF15_1680+XMPROT)**2)*
     *	 (EPQ2+(LF15_1680-XMPROT)**2))/(2.*LF15_1680)
	 PVEC_PI= SQRT((LF15_1680**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LF15_1680**2) - XMPI**2)
	 FACT = SQRT(1./6./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LF15_1680*PIBF15_1680/WF15_1680)
	 RAF15_1680 = FACT*A12F15_1680*19.73*CPIN
	 RBF15_1680 = -FACT*SQRT(16./32.)*A32F15_1680*19.73*
     *		       CPIN

	 RCF15_1680 =		   0.0000D+00

	 RAG17_2190 =		   0.0000D+00
	 RBG17_2190 =		   0.0000D+00
	 RCG17_2190 =		   0.0000D+00

	 RAG19_2250 =		   0.0000D+00
	 RBG19_2250 =		   0.0000D+00
	 RCG19_2250 =		   0.0000D+00

	 RAH19_2220 =		   0.0000D+00
	 RBH19_2220 =		   0.0000D+00
	 RCH19_2220 =		   0.0000D+00

	 RAI111_2600 =		   0.0000D+00
	 RBI111_2600 =		   0.0000D+00
	 RCI111_2600 =		   0.0000D+00
C
	 PVEC_GAM = SQRT((EPQ2 + (LS31_1620+XMPROT)**2)*
     *	 (EPQ2+(LS31_1620-XMPROT)**2))/(2.*LS31_1620)
	 PVEC_PI= SQRT((LS31_1620**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LS31_1620**2) - XMPI**2)
	 FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LS31_1620*PIBS31_1620/WS31_1620)
	 RAS31_1620 =  FACT*A12S31_1620*19.73*SQRT(1./3.)
	 RCS31_1620 =		   0.0000D+00

	 RAS31_1900 =		   0.0000D+00
	 RCS31_1900 =		   0.0000D+00

	 PVEC_GAM = SQRT((EPQ2 + (LP31_1910+XMPROT)**2)*
     *	 (EPQ2+(LP31_1910-XMPROT)**2))/(2.*LP31_1910)
	 PVEC_PI= SQRT((LP31_1910**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LP31_1910**2) - XMPI**2)
	 FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LP31_1910*PIBP31_1910/WP31_1910)
	 RAP31_1910 = -FACT*A12P31_1910*19.73*SQRT(1./3.)
	 RCP31_1910 =		   0.0000D+00


	 PVEC_GAM = SQRT((EPQ2 + (LP33_1232+XMPROT)**2)*
     *	 (EPQ2+(LP33_1232-XMPROT)**2))/(2.*LP33_1232)
	 PVEC_GAM0 = SQRT(((LP33_1232+XMPROT)**2)*
     *	 ((LP33_1232-XMPROT)**2))/(2.*LP33_1232)
	 Q0_GAM=(LP33_1232**2 - XMPROT**2+EPQ2)/(2.*LP33_1232)
	 PVEC_PI= SQRT((LP33_1232**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LP33_1232**2) - XMPI**2)
	   FACT = SQRT(1./137.*LP33_1232*PIBP33_1232/6./XMPROT/PVEC_PI/
     *		    XMPROT**2/WP33_1232)*PVEC_GAM
	 FFAC = 1./(1+EPQ2/0.71)**2*EXP(-0.21*EPQ2)
c	  RAP33_1232 =	-0.44*3.1*FFAC*FACT*SQRT(1./3.)*19.73
	 RAP33_1232 =  -0.5*3.1*FFAC*FACT*SQRT(1./3.)*19.73
	 RBP33_1232 =	    3.1*FFAC*FACT*SQRT(1./3.)*19.73
c	  RCP33_1232 =	-0.03*2.*SQRT(EPQ2)/Q0_GAM*3.*FFAC*
c     * 		 FACT*SQRT(1./3.)*19.73

C*  SET THE AMPLITUDE OF P33 TO ZERO
ccc	     RAP33_1232=0.
ccc	     RBP33_1232=0.
cccc	     RCP33_1232 =  -0.07*2.*SQRT(EPQ2)/Q0_GAM*3.*FFAC*
cccc	 *		    FACT*SQRT(1./3.)*19.73


	 PVEC_GAM = SQRT(EPQ2+ (LP33_1920+XMPROT)**2)*
     *	 (EPQ2+(LP33_1920-XMPROT)**2)/(2.*LP33_1920)
	 PVEC_PI= SQRT(((LP33_1920**2+XMPI**2-XMPROT**2)/
     *	 (4*LP33_1920**2) - XMPI**2))
	 FACT = SQRT((1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LP33_1920*PIBP33_1920/WP33_1920))
	 RAP33_1920 =  FACT*A12P33_1920*SQRT(1./3.)
	 RBP33_1920 = -FACT*SQRT(16./12.)*A32P33_1920*SQRT(1./3.)
	 RCP33_1920 =		    0.0000D+00

	 PVEC_GAM = SQRT((EPQ2 + (LD33_1700+XMPROT)**2)*
     *	 (EPQ2+(LD33_1700-XMPROT)**2))/(2.*LD33_1700)
	 PVEC_PI= SQRT((LD33_1700**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LD33_1700**2) - XMPI**2)
	 FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LD33_1700*PIBD33_1700/WD33_1700)
	 RAD33_1700 = -FACT*A12D33_1700*19.73*SQRT(1./3.)
	 RBD33_1700 =  FACT*SQRT(16./12.)*A32D33_1700*19.73*
     *		      SQRT(1./3.)
	 RCD33_1700 =		 0.0000D+00

	 RAD35_1930 =		    0.0000D+00
	 RBD35_1930 =		    0.0000D+00
	 RCD35_1930 =		    0.0000D+00

	 PVEC_GAM = SQRT((EPQ2 + (LF35_1905+XMPROT)**2)*
     *	 (EPQ2+(LF35_1905-XMPROT)**2))/(2.*LF35_1905)
	 PVEC_PI= SQRT((LF35_1905**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LF35_1905**2) - XMPI**2)
	 FACT = SQRT(1./6./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LF35_1905*PIBF35_1905/WF35_1905)
	 RAF35_1905 = -FACT*A12F35_1905*19.73*SQRT(1./3.)
	 RBF35_1905 =  FACT*SQRT(16./32.)*A32F35_1905*19.73*
     *		       SQRT(1./3.)
	 RCF35_1905 =		    0.0000D+00

	 PVEC_GAM = SQRT((EPQ2 + (LF37_1950+XMPROT)**2)*
     *	 (EPQ2+(LF37_1950-XMPROT)**2))/(2.*LF37_1950)
	 PVEC_PI= SQRT((LF37_1950**2+XMPI**2-XMPROT**2)**2/
     *	 (4.*LF37_1950**2) - XMPI**2)
	 FACT = SQRT(1./8./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *	 LF37_1950*PIBF37_1950/WF37_1950)
	 RAF37_1950 =  FACT*A12F37_1950*19.73*SQRT(1./3.)
	 RBF37_1950 = -FACT*SQRT(16./60.)*A32F37_1950*19.73*
     *		       SQRT(1./3.)
	 RCF37_1950 =		    0.0000D+00

	 RAH311_2420 =		    0.0000D+00
	 RBH311_2420 =		    0.0000D+00
	 RCH311_2420 =		    0.0000D+00
c      rm1m=-0.19
c      RAP11_1440=RM1M/sqrt(0.60)
c      RE0P=-0.11
c      RAS11_1535=RE0P/sqrt(0.40)
c      RE0P=0.090
c      RAS11_1650=RE0P/sqrt(0.60)
c      RE2M=-0.55
c      RM2M=-0.16
c      RAD13_1520=(3./2.*RM2M-1./2.*RE2M)/sqrt(0.55)
c      RBD13_1520=(RE2M+RM2M)/sqrt(0.55)
c      RE2M=-0.034
c      RM2M=-0.016
c      RAD13_1700=(3./2.*RM2M-1./2.*RE2M)/sqrt(0.10)
c      RBD13_1700=(RE2M+RM2M)/sqrt(0.10)
c      RE3M=0.039
c      RM3M=0.0
c      RAF15_1680=(2.*RM3M-RE3M)/sqrt(0.60)
c      RBF15_1680=(RE3M+RM3M)/sqrt(0.60)
c      RE3M=0.23
c      RM3M=-0.17
c      RAD33_1700=(3./2.*RM2M-1./2.*RE2M)/sqrt(0.15)
c      RBD33_1700=(RE2M+RM2M)/sqrt(0.15)

	 RETURN
	 END
***************************************************************************
***************************************************************************
* DEVELOPED IN END OF SEPT. 1992 Z.LI
* With IFLAG and VFLAG both equal to 0, we calculate electric born only
* Both neutron and proton form factor took from original program hborn.
* The axial form factor may not correct.
* A cut off factor is added for the born terms when the center of
* mass energy is larger than 1.3 GEV:
* Fcut=cut**2/(cut**2-(w-1.3)**2), cut is set to 0.4
**************************************************************************
       SUBROUTINE BORNT(Q2,W,SINX,COSX,SINX2,COSX2,BORN,IR,CUT)
       implicit none
       real q2,w,sinx,cosx,sinx2,cosx2,cut
       integer ir
       real xmpip,xmpi0,xmpi,xmp,up,un,ec,gn
       real BORN(6),GAMMA(6,3)
       real A(6,3),AR(6),FR(6),F1(3),F2(3),FP(3),F2A(3)
       real C(6),B(6,6),ETA(6),EPS(3)
       real pi
       real s,q0,e1,e2,p,p0,t,u,sq2
       real gg,gep,cap,f1p,f1n,f2p,f2n,f1v,f1s,f2v,f2s
       real q,z1,z2,fc,pq,fcut
       integer iflag
       real vflag
       integer i,j

       DATA XMPIP/0.139563/,XMPI0/0.1349630/,XMP/0.93827/
       DATA UP,UN/1.793,-1.913/
       DATA EC/.302862/,GN/13.5/
       DATA ETA/1.,1.,-1.,1.,-1.,-1./,EPS/1.,1.,-1./
       DATA IFLAG/0/,VFLAG/0./
C* IFLAG=0, VFLAG=0. WOULD GET RIDE OF MAGNETIC AND PV EXTRA TERM
C* RESPECTIVELY.
       IF(IR.EQ.1.OR.IR.EQ.2)XMPI=XMPI0
       IF(IR.EQ.3.OR.IR.EQ.4)XMPI=XMPIP
       PI=4.0*ATAN(1.0)
       GG=SQRT(4.*PI*GN)
C*  THE PI-N COUPLING CONSTANT
       GEP=1./(1.+3.04*Q2+1.54*Q2*Q2+0.068*Q2*Q2*Q2)
       CAP=Q2/(4.*XMP*XMP)
       F1P=EC*GEP*(1.+(1.+UP)*CAP)/(1.+CAP)
       F1N=EC*GEP*UN*CAP/(1.+CAP)
       IF(IFLAG.EQ.0) THEN
	 F2P=0.
	 F2N=0.
       ELSE
	 F2P=EC/(2.*XMP)*GEP*UP/(1.+CAP)
	 F2N=EC/(2.*XMP)*GEP*UN/(1.+CAP)
       END IF
       F1V=F1P-F1N
       F1S=F1P+F1N
       F2V=F2P-F2N
       F2S=F2P+F2N
       F1(1)=F1V
       F1(2)=F1S
       F1(3)=F1V
       F2(1)=F2V
       F2(2)=F2S
       F2(3)=F2V
C* 1,3,2, CORESSPONDING TO THE ISOVECTOR (+,-) AND ISOSCALAR (0).
       FP(3)=EC/(1.+Q2/0.5)
       FP(1)=0.
       FP(2)=0.
C* PION FORM FACTOR.
       S=W**2
       Q0=(S-XMP**2-Q2)/(2.*W)
       Q=SQRT(Q2+Q0**2)
       E1=(S+XMP**2+Q2)/(2.*W)
       E2=(S+XMP**2-XMPI**2)/(2.*W)
       P=SQRT(E2**2-XMP**2)
       P0=SQRT(P**2+XMPI**2)
       T=2.*Q*P*COSX-2.*Q0*P0+XMPI**2-Q2
       U=-2.*Q*P*COSX-2.*Q0*E2+XMP**2-Q2
       DO 15,J=1,3
	 GAMMA(1,J)=0.5*GG*F1(J)
	 GAMMA(2,J)=-GG*F1(J)/(T-XMPI**2)
	 GAMMA(3,J)=-0.5*GG*F2(J)
	 GAMMA(4,J)=GAMMA(3,J)
	 GAMMA(5,J)=0.5*GAMMA(2,J)
	 GAMMA(6,J)=0.
15     CONTINUE
       DO 25,I=1,4
       DO 35,J=1,3
	 A(I,J)=(1./(S-XMP**2)+EPS(J)*ETA(I)/(U-XMP**2))*GAMMA(I,J)
35     CONTINUE
25     CONTINUE
       F2A(1)=F2(1)
       F2A(2)=F2(2)
       A(1,1)=VFLAG*0.5*GG/XMP*F2A(1)+A(1,1)
       A(1,2)=VFLAG*0.5*GG/XMP*F2A(2)+A(1,2)
c* Add in the pseudovector extra term, here we simply used F2(1,2) same
c* as in non-axial term.
       DO 45,J=1,3
       IF (Q2.EQ.0.) THEN
	 A(5,J)=0.
	 A(6,J)=0.
       ELSE
	 A(5,J)=(1./(S-XMP**2)+ETA(5)/(U-XMP**2))*GAMMA(5,J)
     C	 +0.5*(1.-EPS(J))*2.*GG/Q2*(FP(J)-F1(J))/(T-XMPI**2)
	 A(6,J)=(1./(S-XMP**2)+ETA(6)/(U-XMP**2))*GAMMA(6,J)
       END IF
45     CONTINUE
       SQ2=SQRT(2.0)
       DO 50, I=1,6
	 IF(IR.EQ.1) AR(I)=(A(I,2)+A(I,1))
	 IF(IR.EQ.2) AR(I)=(-A(I,2)+A(I,1))
	 IF(IR.EQ.3) AR(I)=-SQ2*(A(I,2)+A(I,3))
	 IF(IR.EQ.4) AR(I)=SQ2*(A(I,2)-A(I,3))
50     CONTINUE
       Z1=SQRT(E1+XMP)
       Z2=SQRT(E2+XMP)
       FC=1./(8.*PI*W)
       C(1)=FC*(W-XMP)*Z1*Z2
       C(2)=FC*(W+XMP)*P*Q/Z1/Z2
       C(3)=FC*(W+XMP)*P*Q*Z2/Z1
       C(4)=FC*(W-XMP)*P**2*Z1/Z2
       C(5)=FC*Z1/Z2*P
       C(6)=FC*Q*Z2/Z1

c     The following is not used:
c	XX=C(1)/C(2)

       PQ=(T-XMPI**2+Q2)/2.
       B(1,1)=1.
       B(1,2)=0.
       B(1,3)=-PQ/(W-XMP)
       B(1,4)=W-XMP+PQ/(W-XMP)
       B(1,5)=0.
       B(1,6)=Q2/(W-XMP)
       B(2,1)=-1.
       B(2,2)=0.
       B(2,3)=-PQ/(W+XMP)
       B(2,4)=W+XMP+PQ/(W+XMP)
       B(2,5)=0.
       B(2,6)=Q2/(W+XMP)
       B(3,1)=0.
       B(3,2)=(S-XMP**2+Q2/2.)/(W+XMP)
       B(3,3)=1.
       B(3,4)=-1.
       B(3,5)=-Q2/(W+XMP)
       B(3,6)=0.
       B(4,1)=0.
       B(4,2)=-(S-XMP**2+Q2/2.)/(W-XMP)
       B(4,3)=1.
       B(4,4)=-1.
       B(4,5)=Q2/(W-XMP)
       B(4,6)=0.
       B(5,1)=-E1+XMP
       B(5,2)=1./(2.*Q0)*(Q**2*(3.*PQ+2.*Q0*W)-(PQ+P0*Q0)*
     C	 (2.*(S-XMP**2)+Q2))
       B(5,3)=PQ+P0*(W-XMP)
       B(5,4)=-PQ-P0*(W+XMP)+(E1-XMP)*(W+XMP)
       B(5,5)=P0*Q2-Q0*PQ
       B(5,6)=-(E1-XMP)*(W+XMP)
       B(6,1)=E1+XMP
       B(6,2)=-B(5,2)
       B(6,3)=PQ+P0*(W+XMP)
       B(6,4)=-PQ-P0*(W+XMP)+(E1+XMP)*(W-XMP)
       B(6,5)=-B(5,5)
       B(6,6)=-(E1+XMP)*(W-XMP)
       DO 55, I=1,6
       FR(I)=0.
       DO 65, J=1,6
	 FR(I)=FR(I)+C(I)*B(I,J)*AR(J)*19.732857
65     CONTINUE
55     CONTINUE
       IF(W.LE.1.3) FCUT=1.
       IF(W.GT.1.3) FCUT=CUT**2/(CUT**2+(W-1.3)**2)
c	fcut=1
       BORN(1)=-SINX*COSX2/SQ2*(FR(3)+FR(4))
       BORN(2)=-SQ2*COSX2*(FR(1)-FR(2))+SINX2*SINX*(FR(3)-FR(4))/SQ2
       BORN(3)=SINX2*SINX/SQ2*(FR(3)-FR(4))
       BORN(4)=SQ2*SINX2*(FR(1)+FR(2))+COSX2*SINX*(FR(3)+FR(4))/SQ2
       BORN(5)=COSX2*(FR(5)+FR(6))*Q0/Q
       BORN(6)=SINX2*(FR(5)-FR(6))*Q0/Q
       DO 75, I=1,6
       BORN(I)=BORN(I)*FCUT
75     CONTINUE
       RETURN
       END
c*********************************************************************
	 SUBROUTINE QKEM(EPW,EPQ2,QME11,QMM11,QMM12)
	 implicit none
	 real epw,epq2,qme11,qmm11,qmm12,xmprot
	 real q2evf,q2evf0,dip_evf,dip_evf0
	 DATA XMPROT/0.938/
	 Q2EVF=((EPW**2-XMPROT**2)**2+EPQ2*(EPW+XMPROT)**2)/
     *	 (4.*EPW*XMPROT)
	 Q2EVF0=((EPW**2-XMPROT**2)**2)/(4.*EPW*XMPROT)
	 DIP_EVF = (1.+ Q2EVF/0.71)**(-2.)
	 dip_evf0 = (1.+ Q2EVF0/0.71)**(-2.)
C* PUT THE FOLLOWING TO USE THE COUPLING CONSTANT GIVEN IN DATA STATEMENT
C	   GO TO 98
		  QME11 = 3.  * DIP_EVF
		  QMM11 =  3.8 * (0.4 - Q2EVF) * DIP_EVF
	 IF(Q2EVF.LE.1.2)THEN
		  QMM12 = 5.* Q2EVF * DIP_EVF
	  ELSE
		  QMM12 = ( 7.2 - 1.0*Q2EVF ) * DIP_EVF
	 ENDIF
	 RETURN
	 END

	 SUBROUTINE QKEM2(EPW,EPQ2,QME22,QMM21,QMM22,QMM23)
	 implicit none
	 real epw,epq2,qme22,qmm21,qmm22,qmm23
	 real xmprot,q2evf,q2evf0,dip_evf,dip_evf0

	 DATA XMPROT/0.938/
	 Q2EVF=((EPW**2-XMPROT**2)**2+EPQ2*(EPW+XMPROT)**2)/
     *	 (4.*EPW*XMPROT)
	 Q2EVF0=((EPW**2-XMPROT**2)**2)/(4.*EPW*XMPROT)
	 DIP_EVF = (1.+ Q2EVF/0.71)**(-2.)
	 DIP_EVF0 = (1.+ Q2EVF0/0.71)**(-2)
CCC
CCC	  NOTE THAT QME22 HAS TO BE RENORMALIZED TO REPRODUCE DATA
CCC	      REASON IS UNCLEAR !!!!
CCC
	 QME22 = 0.99*SQRT(Q2EVF)*DIP_EVF
     *		  /SQRT(Q2EVF0)/DIP_EVF0
C
	  QMM21 = 0.* SQRT(Q2EVF)*DIP_EVF
	  QMM22 = (0.75*SQRT(5.)-1.5*SQRT(5.)*Q2EVF)*
     *		  SQRT(Q2EVF)*DIP_EVF
	  IF(Q2EVF.LE.1.2)THEN
	  QMM23 = 5.*Q2EVF*SQRT(Q2EVF)*DIP_EVF
	  ELSE
	  QMM23 = (7.2 - Q2EVF)*SQRT(Q2EVF)*DIP_EVF
	  ENDIF
	  RETURN
	  END
*************************************************************************
	 REAL FUNCTION SIGMAO(EPW,EPQ2,EPEPS,EPCOS,EPPHI,EPIREA,
     *	 POL_ELEC,POL_TARG,POL_TARG_THETA,POL_TARG_PHI,sig0
     *	 ,sigu,sigt,sigl,sigi,sigip,fkt)

	 implicit none
	 real epw,epq2,epeps,epcos,EPPHI,pol_elec,pol_targ
	 real pol_targ_theta,pol_targ_phi,sig0,sigu,sigt,sigl,sigi,sigip
	 real fkt
	 real sindd,cosdd,pol_x,pol_y,pol_z,sigpt,siget
	 real xila,sige,theta_gam,enue,theta_elec,e_elec_scatt
	 real e_elec_prim
	 real c2,c3
	 real q2,sinx,cosx,sinx2,cosx2
	 real eph1b,eph2b,eph3b,eph4b,eph5b,eph6b
	 real xmpot
	 real sq2,pi,xmprot,xmpip,xmpi0,xmpi,dum
	 real epsin,epsin2,epcos2,pener,pimp,gener,gimp,gamk
	 real p1s,p2s,p3s,p4s,p5s,p6s,p7s,p2s2,p3s2,p4s2,p5s2
	 real p6s2,p7s2
	 real x,x0,f
	 real bpipy0,bpipy1,bpipy2,bpipy3
	 integer irea
c	The folowing were not declared in the original AO (RCM)
c	I think they are supposed to be complex
	 complex xa4m,xb4m,xc4m,xa4p,xb4p,xc4p,xa5m,xb5m,xc5m,xa5p
	 complex xb5p,xc5p,xa6m,xb6m,xc6m

c      CHARACTER*3 GON
c      LOGICAL TEST
      COMMON /IP/IT,IB,NF,IBORN,CUT,IP11
c      COMMON /GOCOM/TEST,GON
      real cut
      integer it,ib,nf,iborn,ip11

	 COMMON /BACKG/
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0,
     * BGB2M0,BGB3M0,
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP
       COMMON/RCONST/
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710,
     * PIBP13_1720, PIBD13_1520, PIBD13_1700,
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600,
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975
	 COMMON/RAMP/
     * RAS11_1535,	       RCS11_1535,
     * RAS11_1650,	       RCS11_1650,
     * RAP11_1440,	       RCP11_1440,
     * RAP11_1710,	       RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,
     * RAP13_1910, RBP13_1910, RCP13_1910,
     * RAD13_1520, RBD13_1520, RCD13_1520,
     * RAD13_1700, RBD13_1700, RCD13_1700,
     * RAD15_1675, RBD15_1675, RCD15_1675,
     * RAF15_1680, RBF15_1680, RCF15_1680,
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,	       RCS31_1620,
     * RAS31_1900,	       RCS31_1900,
     * RAP31_1910,	       RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232,
     * RAP33_1920, RBP33_1920, RCP33_1920,
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1925, RCP31_1925,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975

	 real
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0,
     * BGB2M0,BGB3M0,
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP

       real
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710,
     * PIBP13_1720, PIBD13_1520, PIBD13_1700,
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600,
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975

       real
     * RAS11_1535,	       RCS11_1535,
     * RAS11_1650,	       RCS11_1650,
     * RAP11_1440,	       RCP11_1440,
     * RAP11_1710,	       RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,
     * RAP13_1910, RBP13_1910, RCP13_1910,
     * RAD13_1520, RBD13_1520, RCD13_1520,
     * RAD13_1700, RBD13_1700, RCD13_1700,
     * RAD15_1675, RBD15_1675, RCD15_1675,
     * RAF15_1680, RBF15_1680, RCF15_1680,
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,	       RCS31_1620,
     * RAS31_1900,	       RCS31_1900,
     * RAP31_1910,	       RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232,
     * RAP33_1920, RBP33_1920, RCP33_1920,
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1925, RCP31_1925,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975


      COMPLEX
     * QAS11_1535, QCS11_1535,
     * QAS11_1650, QCS11_1650,
     * QAP11_1440, QCP11_1440,
     * QAP11_1710, QCP11_1710,
     * QAP13_1720, qbp13_1720, QCP13_1720,
     * QAD13_1520, QBD13_1520, QCD13_1520,
     * QAD13_1700, QBD13_1700, QCD13_1700,
     * QAD15_1675, QBD15_1675, QCD15_1675,
     * QAF15_1680, QBF15_1680, QCF15_1680,
     * QAG17_2190, QBG17_2190, QCG17_2190,
     * QAG19_2250, QBG19_2250, QCG19_2250,
     * QAH19_2220, QBH19_2220, QCH19_2220,
     * QAI111_2600,QBI111_2600,QCI111_2600,
     * QAS31_1620, QCS31_1620,
     * QAP31_1910, QCP31_1910,
     * QAP33_1232, QBP33_1232, QCP33_1232,
     * QAP33_1920, QBP33_1920, QCP33_1920,
     * QAD33_1700, QBD33_1700, QCD33_1700,
     * QAD35_1930, QBD35_1930, QCD35_1930,
     * QAF35_1905, QBF35_1905, QCF35_1905,
     * QAF37_1950, QBF37_1950, QCF37_1950,
     * QAH311_2420,QBH311_2420,QCH311_2420,
     * QAP33_1600, QBP33_1600, QCP33_1600,
     * QAF17_1990, QBF17_1990, QCF17_1990,
     * QAF15_2000, QBF15_2000, QCF15_2000,
     * QAP11_2100, QBP11_2100, QCP11_2100,
     * QAF35_2000, QBF35_2000, QCF35_2000,
     * QAP13_1870, QBP13_1870, QCP13_1870,
     * QAP31_1925, QBP31_1925, QCP31_1925,
     * QAP13_1980, QBP13_1980, QCP13_1980,
     * QAF15_1955, QBF15_1955, QCF15_1955,
     * QAP13_1955, QBP13_1955, QCP13_1955,
     * QAP33_1975, QBP33_1975, QCP33_1975,
     * qas31_1900, qcs31_1900,
     * XA0P  ,XA1P  ,XA2P  ,XA3P  ,
     * XA1M  ,XA2M  ,XA3M  ,
     * XB1P  ,XB2P  ,XB3P  ,
     * XB2M  ,XB3M  ,
     * XC0P  ,XC1P  ,XC2P  ,XC3P  ,
     * XC1M  ,XC2M  ,XC3M
      COMPLEX
     * RS11_1535, RS11_1650, RP11_1440, RP11_1710,
     * RP13_1720, RD13_1520, RD13_1700,
     * RD15_1675, RF15_1680,
     * RG17_2190,
     * RG19_2250, RH19_2220,
     * RI111_2600,
     * RS31_1620, RS31_1900, RP31_1910,
     * RP33_1232, RP33_1920, RD33_1700,
     * RD35_1930, RF35_1905,
     * RF37_1950,
     * RH311_2420,
     * RP33_1600, RF17_1990, RF15_2000, RP11_2100, RF35_2000,
     * RP13_1870,RP31_1925,RP13_1980,RF15_1955,RP13_1955,RP33_1975

      COMPLEX EPH1,EPH2,EPH3,EPH4,EPH5,EPH6
      COMPLEX ADDA0,ADDA1,ADDA2,ADDA3,ADDA4,ADDA5
      COMPLEX ADDB1,ADDB2,ADDB3,ADDB4,ADDB5
      COMPLEX ADDC0,ADDC1,ADDC2,ADDC3,ADDC4,ADDC5
      COMPLEX SUBA0,SUBA1,SUBA2,SUBA3,SUBA4,SUBA5
      COMPLEX SUBB1,SUBB2,SUBB3,SUBB4,SUBB5
      COMPLEX SUBC0,SUBC1,SUBC2,SUBC3,SUBC4,SUBC5
      COMPLEX HNP,HNM,HFP,HFM,HN0,HF0
      COMPLEX X1,X2,Y1,Y2,Y3,Y4,Z1,Z2
      COMPLEX R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11
      REAL HB(6)

	 COMPLEX ZERO
	 INTEGER EPIREA
	 INTEGER IR
	 COMPLEX EPRES
	 REAL YMPI,W,PIM,GK

	 DATA ZERO/(0.,0.)/
      DATA SQ2	      /1.41421	/
      DATA PI	      /3.14159	/
      DATA XMPROT     /.938	/
      DATA XMPIP      /.1395	/
      DATA XMPI0      /.1349	/
c      DATA XMETA      /.5488	 /


	sig0=0.

	 DUM=ACOS(EPCOS)
	 EPSIN=SIN(DUM)
	 EPSIN2=SIN(DUM/2.)
	 EPCOS2=COS(DUM/2.)
	 IF(EPIREA.EQ.1.OR.EPIREA.EQ.2)XMPI=XMPI0
	 IF(EPIREA.EQ.3.OR.EPIREA.EQ.4)XMPI=XMPIP

	   PENER  =(EPW**2-XMPROT**2+XMPI   **2)/(2.*EPW   )
	   PIMP   =SQRT(PENER	**2 - XMPI   **2)
	   GENER  =(EPW   **2-XMPROT**2-EPQ2   )/(2.*EPW   )
	   GIMP   =SQRT(GENER**2 + EPQ2   )
	   GAMK   =(EPW**2-XMPROT**2)/(2.*EPW	)
	   X=EPCOS
	   P1S	=1.
	   P2S	=3.*X
	   P3S	=(15.*X*X-3.)/2.
	   P4S	=(35.*X*X*X-15.*X)/2.
	   P5S	=(315.*X*X*X*X-210.*X*X+15.)/8.
	   P6S	=(693.*X*X*X*X*X-630.*X*X*X)/8.
	   P7S	=(3003.*X*X*X*X*X*X-3465.*X*X*X*X)/16.
	   P2S2 =3.
	   P3S2 =15.*X
	   P4S2 =(105.*X*X-15.)/2.
	   P5S2 =(315.*X*X*X-105.*X)/2.
	   P6S2 =(3465.*X*X*X*X-1890.*X*X)/8.
	   P7S2 =(18018.*X*X*X*X*X-13860.*X*X*X)/16.
	   X0=1./SQRT((PIMP   **2)+(.35**2))
	   F=(1.+(PIMP	 **2)/.71)
	   BPIPY0   =1./F
	   BPIPY1   =X0*PIMP/F
	   BPIPY2   =X0*X0*PIMP*PIMP/F
	   BPIPY3   =X0*X0*X0*PIMP*PIMP*PIMP/F

	   if (gamk .gt. 0.)then
	      FKT=PIMP/GAMK
	   else
	      FKT=0.
	   endif


	   Q2=EPQ2
	   W =EPW
	   SINX =EPSIN
	   COSX =EPCOS
	   SINX2=EPSIN2
	   COSX2=EPCOS2
	   IREA=EPIREA
	   IBORN = 1
	   IF(IBORN.EQ.1) THEN
	    CALL BORNT(Q2,W,SINX,COSX,SINX2,COSX2,HB,IREA,CUT)
	    EPH1B   =HB(1)
	    EPH2B   =HB(2)
	    EPH3B   =HB(3)
	    EPH4B   =HB(4)
	    EPH5B   =HB(5)
	    EPH6B   =HB(6)
	   ELSE
	    EPH1B=0.
	    EPH2B=0.
	    EPH3B=0.
	    EPH4B=0.
	    EPH5B=0.
	    EPH6B=0.
	   END IF

	   IR=EPIREA
	   YMPI  =  XMPI
	   W	 =  EPW
	   PIM	 =  PIMP
	   GK	 =  GAMK

      RS11_1535 = EPRES(IR,XMPROT,YMPI,W,LS11_1535,WS11_1535,
     * .35,0,1,1,PIM,GK,.TRUE.)
       RS11_1650  = EPRES(IR,XMPROT,YMPI,W,LS11_1650,WS11_1650,
     * .35,0,1,1,PIM,GK,.FALSE.)
       RP11_1440 = EPRES(IR,XMPROT,YMPI,W,LP11_1440,WP11_1440,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP11_1710 = EPRES(IR,XMPROT,YMPI,W,LP11_1710,WP11_1710,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP13_1720 = EPRES(IR,XMPROT,YMPI,W,LP13_1720,WP13_1720,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RD13_1520 = EPRES(IR,XMPROT,YMPI,W,LD13_1520,WD13_1520,
     * .35,2,1,1,PIM,GK,.FALSE.)
       RD13_1700 = EPRES(IR,XMPROT,YMPI,W,LD13_1700,WD13_1700,
     * .35,2,1,1,PIM,GK,.FALSE.)
       RD15_1675 = EPRES(IR,XMPROT,YMPI,W,LD15_1675,WD15_1675,
     * .35,2,2,1,PIM,GK,.FALSE.)
       RF15_1680 = EPRES(IR,XMPROT,YMPI,W,LF15_1680,WF15_1680,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RG17_2190 = EPRES(IR,XMPROT,YMPI,W,LG17_2190,WG17_2190,
     * .35,4,3,1,PIM,GK,.FALSE.)
       RG19_2250 = EPRES(IR,XMPROT,YMPI,W,LG19_2250,WG19_2250,
     * .35,4,4,1,PIM,GK,.FALSE.)
       RH19_2220 = EPRES(IR,XMPROT,YMPI,W,LH19_2220,WH19_2220,
     * .35,5,4,1,PIM,GK,.FALSE.)
       RI111_2600 = EPRES(IR,XMPROT,YMPI,W,LI111_2600,WI111_2600,
     * .35,6,5,1,PIM,GK,.FALSE.)


       RS31_1620 = EPRES(IR,XMPROT,YMPI,W,LS31_1620,WS31_1620,
     * .35,0,1,3,PIM,GK,.FALSE.)
       RS31_1900 = EPRES(IR,XMPROT,YMPI,W,LS31_1900,WS31_1900,
     * .35,0,1,3,PIM,GK,.FALSE.)
       RP31_1910 = EPRES(IR,XMPROT,YMPI,W,LP31_1910,WP31_1910,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RP33_1232 = EPRES(IR,XMPROT,YMPI,W,LP33_1232,WP33_1232,
     * .185,1,1,3,PIM,GK,.FALSE.)
       RP33_1920 = EPRES(IR,XMPROT,YMPI,W,LP33_1920,WP33_1920,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RD33_1700 = EPRES(IR,XMPROT,YMPI,W,LD33_1700,WD33_1700,
     * .35,2,1,3,PIM,GK,.FALSE.)
       RD35_1930 = EPRES(IR,XMPROT,YMPI,W,LD35_1930,WD35_1930,
     * .35,2,2,3,PIM,GK,.FALSE.)
       RF35_1905 = EPRES(IR,XMPROT,YMPI,W,LF35_1905,WF35_1905,
     * .35,3,2,3,PIM,GK,.FALSE.)
       RF37_1950 = EPRES(IR,XMPROT,YMPI,W,LF37_1950,WF37_1950,
     * .35,3,3,3,PIM,GK,.FALSE.)
       RH311_2420 = EPRES(IR,XMPROT,YMPI,W,LH311_2420,WH311_2420,
     * .35,5,5,3,PIM,GK,.FALSE.)

C* THE FOLOWING ARE ONE OR TWO STAR RESONANCES NOT INCLUDED IN THE ABOVE
       RP33_1600 = EPRES(IR,XMPROT,YMPI,W,LP33_1600,WP33_1600,
     * .35,1,1,3,PIM,GK,.FALSE.)
! aku
c	 RF17_1990=0d0
c	RF15_2000=0d0
       RF17_1990 = EPRES(IR,XMPOT,YMPI,W,LF17_1990,WF17_1990,
     * .35,3,3,1,PIM,GK,.FALSE.)
       RF15_2000 = EPRES(IR,XMPOT,YMPI,W,LF15_2000,WF15_2000,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RP11_2100 = EPRES(IR,XMPROT,YMPI,W,LP11_2100,WP11_2100,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RF35_2000 = EPRES(IR,XMPROT,YMPI,W,LF35_2000,WF35_2000,
     * .35,3,2,3,PIM,GK,.FALSE.)

C*  THE FOLOWING ARE RESONANCES PREDICTED BY THE QUARK MODEL ONLY
       RP13_1870 = EPRES(IR,XMPROT,YMPI,W,LP13_1870,WP13_1870,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP31_1925 = EPRES(IR,XMPROT,YMPI,W,LP31_1925,WP31_1925,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RP13_1980 = EPRES(IR,XMPROT,YMPI,W,LP13_1980,WP13_1980,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RF15_1955 = EPRES(IR,XMPROT,YMPI,W,LF15_1955,WF15_1955,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RP13_1955 = EPRES(IR,XMPROT,YMPI,W,LP13_1955,WP13_1955,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP33_1975 = EPRES(IR,XMPROT,YMPI,W,LP33_1975,WP33_1975,
     * .35,1,1,3,PIM,GK,.FALSE.)

c	    CALL BACK(IB)

	   IF(EPIREA   .EQ.1)THEN
	     XA0P    = BGA0P0*BPIPY0   +ZERO
	     XA1P    = BGA1P0*BPIPY1   +ZERO
	     XA2P    = BGA2P0*BPIPY2   +ZERO
	     XA3P    = BGA3P0*BPIPY3   +ZERO
	     xa4p    = zero
	     xa5p    = zero

	     XB1P    = BGB1P0*BPIPY1   +ZERO
	     XB2P    = BGB2P0*BPIPY2   +ZERO
	     XB3P    = BGB3P0*BPIPY3   +ZERO
	     xb4p    = zero
	     xb5p    = zero

	     XC0P    = BGC0P0*BPIPY0   +ZERO
	     XC1P    = BGC1P0*BPIPY1   +ZERO
	     XC2P    = BGC2P0*BPIPY2   +ZERO
	     XC3P    = BGC3P0*BPIPY3   +ZERO
	     xc4p     =  zero
	     xc5p     = zero

	     XA1M    = BGA1M0*BPIPY1   +ZERO
	     XA2M    = BGA2M0*BPIPY2   +ZERO
	     XA3M    = BGA3M0*BPIPY3   +ZERO
	     xa4m    = zero
	     xa5m     = zero
	     xa6m     = zero

	     XB2M    = BGB2M0*BPIPY2   +ZERO
	     XB3M    = BGB3M0*BPIPY3   +ZERO
	     xb4m     = zero
	     xb5m    = zero
	     xb6m     = zero

	     XC1M    = BGC1M0*BPIPY1   +ZERO
	     XC2M    = BGC2M0*BPIPY2   +ZERO
	     XC3M    = BGC3M0*BPIPY3   +ZERO
	     xc4m    =	zero
	     xc5m     = zero
	     xc6m    = zero

	   ELSE

	     XA0P    = BGA0PP*BPIPY0   +ZERO
	     XA1P    = BGA1PP*BPIPY1   +ZERO
	     XA2P    = BGA2PP*BPIPY2   +ZERO
	     XA3P    = BGA3PP*BPIPY3   +ZERO

	     XB1P    = BGB1PP*BPIPY1   +ZERO
	     XB2P    = BGB2PP*BPIPY2   +ZERO
	     XB3P    = BGB3PP*BPIPY3   +ZERO

	     XC0P    = BGC0PP*BPIPY0   +ZERO
	     XC1P    = BGC1PP*BPIPY1   +ZERO
	     XC2P    = BGC2PP*BPIPY2   +ZERO
	     XC3P    = BGC3PP*BPIPY3   +ZERO

	     XA1M    = BGA1MP*BPIPY1   +ZERO
	     XA2M    = BGA2MP*BPIPY2   +ZERO
	     XA3M    = BGA3MP*BPIPY3   +ZERO

	     XB2M    = BGB2MP*BPIPY2   +ZERO
	     XB3M    = BGB3MP*BPIPY3   +ZERO

	     XC1M    = BGC1MP*BPIPY1   +ZERO
	     XC2M    = BGC2MP*BPIPY2   +ZERO
	     XC3M    = BGC3MP*BPIPY3   +ZERO

	   END IF

C* IT=1 FOR EXPT ONLY, IT=2 FOR EXPT+ONE OR TWO STAR FROM QKM
C* IT=2 FOR EXPT+ALL THE OTHER FROM QKM
C* IT=4 FOR ALL FROM QKM
C* C2,C3 ARE USED TO TURN THE QKM PREDICTION OFF AND ON

	 IF(IT.EQ.1) THEN
	   CALL EXPA(EPIREA,EPQ2)
	 IF(IP11.EQ.2) THEN
	   RAP11_1440=0.0
	 ENDIF
	 IF(IP11.EQ.3.OR.IP11.EQ.4) THEN
c	    CALL QKMA(EPIREA,EPQ2)
	 ENDIF
	 C2=0.
	 C3=0.
	 GO TO 35
	 ENDIF

	 IF(IT.EQ.2) THEN
	 C2=1.
	 C3=0.
c	  CALL QKMA(EPIREA,EPQ2)
	 CALL EXPA(EPIREA,EPQ2)
	 GO TO 25
	 ENDIF

	 IF (IT.EQ.3) THEN
c	  CALL QKMA(EPIREA,EPQ2)
	 CALL EXPA(EPIREA,EPQ2)
	 C2=1.
	 C3=1.
	 ENDIF

	 IF(IT.EQ.4) THEN
	 C2=0.
	 C3=0.
C* here I get ride of the terms which is not used in
c* the EXPT
c	  CALL QKMA(EPIREA,EPQ2)
C*	   GO TO 35
	 ENDIF

C* THIS FOLLOWING RESONANCES ARE PREDICTED BY QK MODEL ONLY
15	 QAP13_1870   =RP13_1870   * RAP13_1870
	 QBP13_1870   =RP13_1870   * RBP13_1870
	 QCP13_1870   =RP13_1870   * 0.0

	 QAP31_1925   =RP31_1925   * RAP31_1925
	 QBP31_1925   =RP31_1925   * RBP31_1925
	 QCP31_1925   =RP31_1925   * 0.0

	 QAP13_1980   =RP13_1980   * RAP13_1980
	 QBP13_1980   =RP13_1980   * RBP13_1980
	 QCP13_1980   =RP13_1980   * 0.0

	 QAF15_1955   =RF15_1955   * RAF15_1955
	 QBF15_1955   =RF15_1955   * RBF15_1955
	 QCF15_1955   =RF15_1955   * 0.0

	 QAP13_1955   =RP13_1955   * RAP13_1955
	 QBP13_1955   =RP13_1955   * RBP13_1955
	 QCP13_1955   =RP13_1955   * 0.0

	 QAP33_1975   =RP33_1975   * RAP33_1975
	 QBP33_1975   =RP33_1975   * RBP33_1975
	 QCP33_1975   =RP33_1975   * 0.0


C*  THE FOLLOWING RESONANCES ARE SEEN AS ONE OR TWO STAR
25	 QAP33_1600   =RP33_1600   * RAP33_1600
	 QBP33_1600   =RP33_1600   * RBP33_1600
	 QCP33_1600   =RP33_1600   * 0.0

	 QAF17_1990   =RF17_1990   * RAF17_1990
	 QBF17_1990   =RF17_1990   * RBF17_1990
	 QCF17_1990   =RF17_1990   * 0.0

	 QAF15_2000   =RF15_2000   * RAF15_2000
	 QBF15_2000   =RF15_2000   * RBF15_2000
	 QCF15_2000   =RF15_2000   * 0.0

	 QAP11_2100   =RP11_2100   * RAP11_2100
	 QBP11_2100   =RP11_2100   * RBP11_2100
	 QCP11_2100   =RP11_2100   * 0.0

	 QAF35_2000   =RF35_2000   * RAF35_2000
	 QBF35_2000   =RF35_2000   * RBF35_2000
	 QCF35_2000   =RF35_2000   * RCF35_2000


C* THE FOLLWING RESONANCES ARE SEEN BY EXPERIMENT
35	 QAS11_1535   =RS11_1535   * RAS11_1535
	 QCS11_1535   =RS11_1535   * RCS11_1535

	 QAS11_1650   =RS11_1650   * RAS11_1650
	 QCS11_1650   =RS11_1650   * RCS11_1650

	 QAP11_1440   =RP11_1440   * RAP11_1440
	 QCP11_1440   =RP11_1440   * RCP11_1440

	 QAP13_1720   =RP13_1720   * RAP13_1720
	 QBP13_1720   =RP13_1720   * RBP13_1720
	 QCP13_1720   =RP13_1720   * RCP13_1720

	 QAD13_1520   =RD13_1520   * RAD13_1520
	 QBD13_1520   =RD13_1520   * RBD13_1520
	 QCD13_1520   =RD13_1520   * RCD13_1520

	 QAD13_1700   =RD13_1700   * RAD13_1700
	 QBD13_1700   =RD13_1700   * RBD13_1700
	 QCD13_1700   =RD13_1700   * RCD13_1700

	 QAD15_1675   =RD15_1675   * RAD15_1675
	 QBD15_1675   =RD15_1675   * RBD15_1675
	 QCD15_1675   =RD15_1675   * RCD15_1675

	 QAF15_1680   =RF15_1680   * RAF15_1680
	 QBF15_1680   =RF15_1680   * RBF15_1680
	 QCF15_1680   =RF15_1680   * RCF15_1680

	 QAG17_2190   =RG17_2190   * RAG17_2190
	 QBG17_2190   =RG17_2190   * RBG17_2190
	 QCG17_2190   =RG17_2190   * RCG17_2190

	 QAG19_2250   =RG19_2250   * RAG19_2250
	 QBG19_2250   =RG19_2250   * RBG19_2250
	 QCG19_2250   =RG19_2250   * RCG19_2250

	 QAH19_2220   =RH19_2220   * RAH19_2220
	 QBH19_2220   =RH19_2220   * RBH19_2220
	 QCH19_2220   =RH19_2220   * RCH19_2220

	 QAI111_2600  =RI111_2600  * RAI111_2600
	 QBI111_2600  =RI111_2600  * RBI111_2600
	 QCI111_2600  =RI111_2600  * RCI111_2600

	 QAS31_1620   =RS31_1620   * RAS31_1620
	 QCS31_1620   =RS31_1620   * RCS31_1620

	 QAS31_1900   =RS31_1900   * RAS31_1900
	 QCS31_1900   =RS31_1900   * RCS31_1900

	 QAP31_1910   =RP31_1910   * RAP31_1910
	 QCP31_1910   =RP31_1910   * RCP31_1910

	 QAP33_1232   =RP33_1232   * RAP33_1232
	 QBP33_1232   =RP33_1232   * RBP33_1232
	 QCP33_1232   =RP33_1232   * RCP33_1232

	 QAP33_1920   =RP33_1920   * RAP33_1920
	 QBP33_1920   =RP33_1920   * RBP33_1920
	 QCP33_1920   =RP33_1920   * RCP33_1920


	 QAD33_1700   =RD33_1700   * RAD33_1700
	 QBD33_1700   =RD33_1700   * RBD33_1700
	 QCD33_1700   =RD33_1700   * RCD33_1700

	 QAD35_1930   =RD35_1930   * RAD35_1930
	 QBD35_1930   =RD35_1930   * RBD35_1930
	 QCD35_1930   =RD35_1930   * RCD35_1930

	 QAF35_1905   =RF35_1905   * RAF35_1905
	 QBF35_1905   =RF35_1905   * RBF35_1905
	 QCF35_1905   =RF35_1905   * RCF35_1905

	 QAF37_1950   =RF37_1950   * RAF37_1950
	 QBF37_1950   =RF37_1950   * RBF37_1950
	 QCF37_1950   =RF37_1950   * RCF37_1950

	 QAH311_2420  =RH311_2420  * RAH311_2420
	 QBH311_2420  =RH311_2420  * RBH311_2420
	 QCH311_2420  =RH311_2420  * RCH311_2420

c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c	The following were not initialized in my copy of AO (RCM)
	qap11_1710=zero
	qcp11_1710=zero
c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c



!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 if(abs(c3).lt.1d-12.and.abs(c2).lt.1d-12)then

	 XA0P=XA0P + QAS11_1535 + QAS11_1650 + QAS31_1620 + QAS31_1900
	 XC0P=XC0P + QCS11_1535 + QCS11_1650 + QCS31_1620 + QCS31_1900

	 XA1M=XA1M + QAP11_1440 + QAP11_1710 + QAP31_1910
c     *      +C2*(QAP11_2100)+C3*(QAP31_1925)
	 XC1M=XC1M + QCP11_1440 + QCP11_1710 + QCP31_1910
c     *      +C2*(QCP11_2100)+C3*(QCP31_1925)

	 XA1P=XA1P + QAP13_1720 + QAP33_1232 + QAP33_1920
c     *+C2*(QAP33_1600)+C3*(QAP13_1870+QAP13_1980+QAP13_1955+QAP33_1975)
	 XB1P=XB1P + QBP13_1720 + QBP33_1232 + QBP33_1920
c     *+C2*(QBP33_1600)+C3*(QBP13_1870+QBP13_1980+QBP13_1955+QBP33_1975)
	 XC1P=XC1P + QCP13_1720 + QCP33_1232 + QCP33_1920
c     *+C2*(QCP33_1600)+C3*(QCP13_1870+QCP13_1980+QCP13_1955+QCP33_1975)

	 XA2M = XA2M + QAD13_1520 + QAD13_1700 + QAD33_1700
	 XB2M = XB2M + QBD13_1520 + QBD13_1700 + QBD33_1700
	 XC2M = XC2M + QCD13_1520 + QCD13_1700 + QCD33_1700

	 XA2P = XA2P + QAD15_1675 + QAD35_1930
	 XB2P = XB2P + QBD15_1675 + QBD35_1930
	 XC2P = XC2P + QCD15_1675 + QCD35_1930

	 XA3M = XA3M + QAF15_1680 + QAF35_1905
c     *       +C2*(QAF15_2000+QAF35_2000)+C3*(QAF15_1955)
	 XB3M = XB3M + QBF15_1680 + QBF35_1905
c     *       +C2*(QBF15_2000+QBF35_2000)+C3*(QBF15_1955)

	 XC3M = XC3M + QCF15_1680 + QCF35_1905
c     *       +C2*(QCF15_2000+QCF35_2000)+C3*(QCF15_1955)

	 XA3P = XA3P + QAF37_1950
c     *       +C2*(QAF17_1990)
	 XB3P = XB3P + QBF37_1950
c     *       +C2*(QBF17_1990)
	 XC3P = XC3P + QCF37_1950
c     *       +C2*(QCF17_1990)

	 else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 XA0P=XA0P + QAS11_1535 + QAS11_1650 + QAS31_1620 + QAS31_1900
	 XC0P=XC0P + QCS11_1535 + QCS11_1650 + QCS31_1620 + QCS31_1900

	 XA1M=XA1M + QAP11_1440 + QAP11_1710 + QAP31_1910
     *	     +C2*(QAP11_2100)+C3*(QAP31_1925)
	 XC1M=XC1M + QCP11_1440 + QCP11_1710 + QCP31_1910
     *	     +C2*(QCP11_2100)+C3*(QCP31_1925)

	 XA1P=XA1P + QAP13_1720 + QAP33_1232 + QAP33_1920
     *+C2*(QAP33_1600)+C3*(QAP13_1870+QAP13_1980+QAP13_1955+QAP33_1975)
	 XB1P=XB1P + QBP13_1720 + QBP33_1232 + QBP33_1920
     *+C2*(QBP33_1600)+C3*(QBP13_1870+QBP13_1980+QBP13_1955+QBP33_1975)
	 XC1P=XC1P + QCP13_1720 + QCP33_1232 + QCP33_1920
     *+C2*(QCP33_1600)+C3*(QCP13_1870+QCP13_1980+QCP13_1955+QCP33_1975)

	 XA2M = XA2M + QAD13_1520 + QAD13_1700 + QAD33_1700
	 XB2M = XB2M + QBD13_1520 + QBD13_1700 + QBD33_1700
	 XC2M = XC2M + QCD13_1520 + QCD13_1700 + QCD33_1700

	 XA2P = XA2P + QAD15_1675 + QAD35_1930
	 XB2P = XB2P + QBD15_1675 + QBD35_1930
	 XC2P = XC2P + QCD15_1675 + QCD35_1930

	 XA3M = XA3M + QAF15_1680 + QAF35_1905
     *	      +C2*(QAF15_2000+QAF35_2000)+C3*(QAF15_1955)
	 XB3M = XB3M + QBF15_1680 + QBF35_1905
c     *       +C2*(QBF15_2000+QBF35_2000)+C3*(QBF15_1955)

	 XC3M = XC3M + QCF15_1680 + QCF35_1905
     *	      +C2*(QCF15_2000+QCF35_2000)+C3*(QCF15_1955)

	 XA3P = XA3P + QAF37_1950
     *	      +C2*(QAF17_1990)
	 XB3P = XB3P + QBF37_1950
     *	      +C2*(QBF17_1990)
	 XC3P = XC3P + QCF37_1950
     *	      +C2*(QCF17_1990)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	 XA4M = XA4M + QAG17_2190
	 XB4M = XB4M + QBG17_2190
	 XC4M = XC4M + QCG17_2190

	 XA4P = XA4P + QAG19_2250
	 XB4P = XB4P + QBG19_2250
	 XC4P = XC4P + QCG19_2250

	 XA5M = XA5M + QAH19_2220
	 XB5M = XB5M + QBH19_2220
	 XC5M = XC5M + QCH19_2220

	 XA5P = XA5P + QAH311_2420
	 XB5P = XB5P + QBH311_2420
	 XC5P = XC5P + QCH311_2420

	 XA6M = XA6M + QAI111_2600
	 XB6M = XB6M + QBI111_2600
	 XC6M = XC6M + QCI111_2600

	 SUBA0=XA0P-XA1M
	 SUBA1=XA1P-XA2M
	 SUBA2=XA2P-XA3M
	 SUBA3=XA3P-XA4M
	 SUBA4=XA4P-XA5M
	 SUBA5=XA5P-XA6M

	 SUBC0=XC0P-XC1M
	 SUBC1=XC1P-XC2M
	 SUBC2=XC2P-XC3M
	 SUBC3=XC3P-XC4M
	 SUBC4=XC4P-XC5M
	 SUBC5=XC5P-XC6M

	 SUBB1=XB1P-XB2M
	 SUBB2=XB2P-XB3M
	 SUBB3=XB3P-XB4M
	 SUBB4=XB4P-XB5M
	 SUBB5=XB5P-XB6M

	 ADDA0=XA0P+XA1M
	 ADDA1=XA1P+XA2M
	 ADDA2=XA2P+XA3M
	 ADDA3=XA3P+XA4M
	 ADDA4=XA4P+XA5M
	 ADDA5=XA5P+XA6M

	 ADDC0=XC0P+XC1M
	 ADDC1=XC1P+XC2M
	 ADDC2=XC2P+XC3M
	 ADDC3=XC3P+XC4M
	 ADDC4=XC4P+XC5M
	 ADDC5=XC5P+XC6M

	 ADDB1=XB1P+XB2M
	 ADDB2=XB2P+XB3M
	 ADDB3=XB3P+XB4M
	 ADDB4=XB4P+XB5M
	 ADDB5=XB5P+XB6M
c
C   Now do the partial wave expansion for the Walker amplitudes
C
       EPH1   =SUBB1*(	     -P2S2   )
     *	     + SUBB2*(P2S2   -P3S2   )
     *	     + SUBB3*(P3S2   -P4S2   )
     *	     + SUBB4*(P4S2   -P5S2   )
     *	     + SUBB5*(P5S2   -P6S2   )
       EPH1   = EPH1   * EPSIN	 * EPCOS2   /SQ2


       EPH2  = SUBA0   *(	-P1S   )
     *	     + SUBA1   *(P1S	-P2S   )
     *	     + SUBA2   *(P2S	-P3S   )
     *	     + SUBA3   *(P3S	-P4S   )
     *	     + SUBA4   *(P4S	-P5S   )
     *	     + SUBA5   *(P5S	-P6S   )
       EPH2  = EPH2    * EPCOS2   * SQ2


       EPH3   = ADDB1	 *(	  +P2S2   )
     *	     +	ADDB2	 *(P2S2   +P3S2   )
     *	     +	ADDB3	 *(P3S2   +P4S2   )
     *	     +	ADDB4	 *(P4S2   +P5S2   )
     *	     +	ADDB5	 *(P5S2   +P6S2   )
       EPH3   = EPH3   *  EPSIN   * EPSIN2   /SQ2


       EPH4  = ADDA0	 * (	   +P1S   )
     *	     + ADDA1	 * (P1S    +P2S   )
     *	     + ADDA2	 * (P2S    +P3S   )
     *	     + ADDA3	 * (P3S    +P4S   )
     *	     + ADDA4	 * (P4S    +P5S   )
     *	     + ADDA5	 * (P5S    +P6S   )
       EPH4   = EPH4   * EPSIN2   * SQ2


       EPH5   =SUBC0	 *(	  -P1S	 )
     *	     + SUBC1	 *(P1S	  -P2S	 )
     *	     + SUBC2	 *(P2S	  -P3S	 )
     *	     + SUBC3	 *(P3S	  -P4S	 )
     *	     + SUBC4	 *(P4S	  -P5S	 )
     *	     + SUBC5	 *(P5S	  -P6S	 )
       EPH5   = EPH5   * EPCOS2   * SQ2


       EPH6  = ADDC0	 * (	   +P1S   )
     *	     + ADDC1	 * (P1S    +P2S   )
     *	     + ADDC2	 * (P2S    +P3S   )
     *	     + ADDC3	 * (P3S    +P4S   )
     *	     + ADDC4	 * (P4S    +P5S   )
     *	     + ADDC5	 * (P5S    +P6S   )
	EPH6   =  EPH6	 * EPSIN2   * SQ2

C
C      Now add the Born term contributions. Remember, these are real functions!
C
	 EPH1	=EPH1	+CMPLX(EPH1B   ,0.)
	 EPH2	=EPH2	+CMPLX(EPH2B   ,0.)
	 EPH3	=EPH3	+CMPLX(EPH3B   ,0.)
	 EPH4	=EPH4	+CMPLX(EPH4B   ,0.)
	 EPH5	=EPH5	+CMPLX(EPH5B   ,0.)
	 EPH6	=EPH6	+CMPLX(EPH6B   ,0.)


C
C	Now we are ready to calculate all the observables
C
C	   First the unpolarized transverse term
C

	   SIGU   =(CABS(EPH1	)**2
     *		  +CABS(EPH2   )**2
     *		  +CABS(EPH3   )**2
     *		  +CABS(EPH4   )**2)/2
C
C	     Next the transverse polarized term
C
	   SIGT   =REAL(EPH2   *CONJG(EPH3   )
     *		  -	EPH1   *CONJG(EPH4   ))
C
C	      And the longitudinal term
C
	   SIGL   =CABS(EPH5   )**2
     *		  +CABS(EPH6   )**2
C
C	      Finally, the longitudinal-transverse interference term
C
	   SIGI   =SQ2*REAL(EPH5   *CONJG(EPH4	 -EPH1	 )
     *		  +	    EPH6   *CONJG(EPH3	 +EPH2	 ))
C
       SIGIP  =SQ2*AIMAG(EPH5	*CONJG(EPH4   -EPH1   )
     *	      + 	EPH6   *CONJG(EPH3   +EPH2   ))
C
C
C	       Sum up all the unpolarized cross section terms
C

	     sig0  = SIGU + EPEPS*SIGL
c	     print *,SIGT,SIGI,SIGU,SIGL
	   SIGMAO   = sig0
     *		   +EPEPS*SIGT*COS(EPPHI*PI/90.)
     *		   +SQRT(EPEPS*(1+EPEPS)/2.)*SIGI
     *		   *COS(EPPHI*PI/180.)



C
	 if(POL_ELEC .eq. 0. .and. POL_TARG .eq. 0.)then
	    sigmao=sigmao*FKT

	    return
	 endif



C
C	       Convert Walker amplitudes to Bartl & Majerotto spin-flip and non spin-flip amplitudes
C
	   HNP=(EPH4+EPH1)/SQ2
	   HNM=(EPH4-EPH1)/SQ2
	   HFP=(EPH3-EPH2)/SQ2
	   HFM=(EPH3+EPH2)/SQ2
	   HN0=EPH5
	   HF0=EPH6
C
C		Define the appreviations for the various polarized cross section terms
C		Polarized beam, polarized target, polarized beam - polarized target
C
	   X1=HF0*CONJG(HNP) + HN0*CONJG(HFP)
	   X2=HFM*CONJG(HNP) + HNM*CONJG(HFP)
	   Y1=HNP*CONJG(HFP) + HNM*CONJG(HFM)
	   Y2=HNM*CONJG(HFM) - HNP*CONJG(HFP)
	   Y3=HN0*CONJG(HF0)
	   Y4=HN0*CONJG(HFM) - HF0*CONJG(HNM)
	   Z1=HN0*CONJG(HNP) - HF0*CONJG(HFP)
	   Z2=HNM*CONJG(HNP) - HFM*CONJG(HFP)
C
C		These are the terms for the recoil polarizations
C
	   R1=CONJG(HN0)*HF0
	   R2=CONJG(HNP)*HFP - CONJG(HNM)*HFM
	   R3=CONJG(HN0)*HFM + CONJG(HF0)*HNP
	   R4=CONJG(HF0)*HNM - CONJG(HN0)*HFM
	   R5=CONJG(HN0)*HNP + CONJG(HF0)*HFP
	   R6=CONJG(HNP)*HFP + CONJG(HNM)*HFM
	   R7=CONJG(HFM)*HNP + CONJG(HFP)*HNM
	   R8=CONJG(HFP)*HFM + CONJG(HNP)*HNM
	   R9=CONJG(HNP)*HFM - CONJG(HNM)*HFP
	   R10=CONJG(HNP)*HNM + CONJG(HFP)*HFM
	   R11=CONJG(HN0)*HNM + CONJG(HF0)*HFM


	   XILA=SQRT(2.*EPEPS*(EPEPS+1))
	   SIGE=-POL_ELEC*SQRT(2*EPEPS*(1-EPEPS))*SIN(EPPHI*PI/180.)*
     *	   AIMAG(HN0*CONJG(HNM)+HF0*CONJG(HFM))


C	    SIGPT=PX*(-XILA*SIN(EPPHI*PI/180)*AIMAG(X1)-EPEPS*
C     * 	SIN(EPPHI*PI/90.)*AIMAG(X2))
C     *       -PY*(AIMAG(Y1)+EPEPS*COS(EPPHI*PI/90.)*AIMAG(Y2)
C     * 	 +2.*EPEPS*AIMAG(Y3)+XILA*COS(EPPHI*PI/180)*AIMAG(Y4))
C     * 	+PZ*(EPEPS*SIN(EPPHI*PI/90.)*AIMAG(Z2)+
C     * 	     XILA*SIN(EPPHI*PI/180.)*AIMAG(Z1))
C	SIGET=-POL_ELEC*(-PX*(SQRT(2.*EPEPS*(1-EPEPS))*COS(EPPHI*PI/180)
C     * 	      *REAL(X1)+SQRT(1-EPEPS**2)*REAL(X2))
C     *       +PY*SQRT(2.*EPEPS*(1-EPEPS))*SIN(EPPHI*PI/180.)*REAL(Y4)
C     *       +PZ*(SQRT(1-EPEPS**2)*REAL(Z2)+SQRT(2*EPEPS*(1-EPEPS))*
C     *        COS(EPPHI*PI/180.)*REAL(Z1)))
C	    SIGMA = SIGMA + SIGE + SIGPT + SIGET
C	    SIGMA=SIGMA*FKT
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C	  The target is polarized along the beam axis.
c					PS = -1 anti-parallel to beam
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	   IF(EPQ2.LT.1.E-06.OR.EPEPS.LT.1.E-06)THEN
		THETA_GAM = 0.
	   ELSE
	   ENUE = (EPW**2 - XMPROT**2 + EPQ2)/2./XMPROT
	   THETA_ELEC = 2.*ATAN(SQRT(EPQ2/2./(EPQ2+ENUE**2)*
     *			(1.-EPEPS)/EPEPS))
	   E_ELEC_SCATT = -ENUE/2 + SQRT((ENUE/2.)**2 +
     *			 EPQ2/4./SIN(THETA_ELEC/2.)**2)
	   E_ELEC_PRIM = E_ELEC_SCATT + ENUE
	   THETA_GAM = ATAN(SIN(THETA_ELEC)/
     *		       (E_ELEC_PRIM/E_ELEC_SCATT - COS(THETA_ELEC)))
	   THETA_GAM = THETA_GAM*180./PI
	   ENDIF

	   POL_X = POL_TARG*SINDD(THETA_GAM)*COSDD(EPPHI)
	   POL_Y = POL_TARG*SINDD(THETA_GAM)*SINDD(EPPHI)
	   POL_Z = POL_TARG*COSDD(THETA_GAM)

	   SIGPT = - POL_X*(SINDD(EPPHI)*
     1	SQRT(2.*EPEPS*(1. + EPEPS))*AIMAG(X1) +
     2	EPEPS*SINDD(2.*EPPHI)*AIMAG(X2)) -
     3	POL_Y*(AIMAG(Y1) +
     4	EPEPS*COSDD(2.*EPPHI)*AIMAG(Y2) +
     5	2.*EPEPS*AIMAG(Y3) +
     6	SQRT(2.*EPEPS*(1.+EPEPS))*COSDD(EPPHI)*AIMAG(Y4)) +
     7	POL_Z*(EPEPS*SINDD(2.*EPPHI)*AIMAG(Z2) +
     8	SQRT(2.*EPEPS*(1. + EPEPS))*SINDD(EPPHI)*AIMAG(Z1))

	SIGET = POL_ELEC*POL_TARG*(SINDD(THETA_GAM)*(COSDD(EPPHI))**2*
     1	(REAL(Y4) + SQRT(2.*EPEPS*(1.-EPEPS))*REAL(X1))+
     2	COSDD(EPPHI)*(SINDD(THETA_GAM)*SQRT(1.-EPEPS**2)*REAL(X2) -
     3	COSDD(THETA_GAM)*SQRT(2.*EPEPS*(1.-EPEPS))*REAL(Z1)) -
     4	SINDD(THETA_GAM)*REAL(Y4) - COSDD(THETA_GAM)*
     5	SQRT(1.-EPEPS**2)*REAL(Z2))

	   SIGMAO = SIGMAO + SIGE + SIGPT + SIGET

c	    print *,SIGMAO , SIGE , SIGPT , SIGET


 99	   SIGMAO = SIGMAO*FKT

	   RETURN
	   END

*****************************************************************************

      real function sindd(x)
      implicit none

      real pi,x
      pi=4.*atan(1.)
      sindd=sin(pi*x/180.)
      return
      end
*****************************************************************************
*****************************************************************************


*****************************************************************************
*****************************************************************************

      real function cosdd(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      cosdd=cos(pi*x/180.)
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************

      real function tandd(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      tandd=tan(pi*x/180.)
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************

      real function asind(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      asind=180.*asin(x)/pi
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************

      real function acosd(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      acosd=180.*acos(x)/pi
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************

      real function atand(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      atand=180.*atan(x)/pi
      return
      end
*****************************************************************************





****************** vacpol *************************************

      double precision function vacpol(t)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amu2,amh2,barn
      dimension am2(3)
c
c    am2 : squared masses of charge leptons
c
      data am2/.26110d-6,.111637d-1,3.18301d0/

      suml=0.
      do 10 i=1,3
	 a2=2.*am2(i)
	 sqlmi=dsqrt(t*t+2.*a2*t)
	 allmi=dlog((sqlmi+t)/(sqlmi-t))/sqlmi
  10  suml=suml+2.*(t+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./t
      if(t.lt.1.d0)then
	aaa = -1.345d-9
	bbb = -2.302d-3
	ccc = 4.091
      elseif(t.lt.64d0)then
	aaa = -1.512d-3
	bbb =  -2.822d-3
	ccc = 1.218
      else
	aaa = -1.1344d-3
	bbb = -3.0680d-3
	ccc = 9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.+ccc*t)) *2*pi/alpha

      vacpol=suml+sumh

      end

****************** fspens *************************************

      double precision function fspens(x)
c
c    spence function
c
      implicit real*8(a-h,o-z)
      f=0.d0
      a=1.d0
      an=0.d0
      tch=1.d-16
  1   an=an+1.d0
      a=a*x
      b=a/an**2
      f=f+b
      if(b-tch)2,2,1
  2   fspens=f
      return
      end

      double precision function fspen(x)
      implicit real*8(a-h,o-z)
      data f1/1.644934d0/
      if(x)8,1,1
  1   if(x-.5d0)2,2,3
    2 fspen=fspens(x)
      return
    3 if(x-1d0)4,4,5
    4 fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
      return
    5 if(x-2d0)6,6,7
    6 fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
    7 fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
      return
    8 if(x+1d0)10,9,9
   9  fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
      return
  10  fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
      return
      end


************ titout **************************
      subroutine titout(date)
      implicit real*8(a-h,o-z)
      character*10 date

      open(unit=9,file='all.dat')
      open(unit=21,file='allu.dat')
	 write(9,1)date
	 write(9,3)
3     format(/' the file contains work information about'/
     .' contribution of tails')
1     format(1x,'program DIFFRAD 1.0 version from ',a10)
      end


      subroutine simpsx(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end
      subroutine simptx(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simpt(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end
      subroutine simpux(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simpu(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

      subroutine dqg32(xl,xu,fct,y)
c
c  computation of integrals by means of 32-point gauss quadrature
c  formula, which integrates polynomials up to degree 63.
c
c
      double precision xl,xu,y,a,b,c,fct
c      external fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      return
      end
      subroutine dqm32(xl,xu,fct,y)
c
c  computation of integrals by means of 32-point gauss quadrature
c  formula, which integrates polynomials up to degree 63.
c
c
      double precision xl,xu,y,a,b,c,fct
c      external fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      return
      end
      subroutine dql32(xl,xu,fct,y)
c
c  computation of integrals by means of 32-point gauss quadrature
c  formula, which integrates polynomials up to degree 63.
c
c
      double precision xl,xu,y,a,b,c,fct
c      external fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      return
      end
      subroutine dqn32(xl,xu,fct,y)
c
c  computation of integrals by means of 32-point gauss quadrature
c  formula, which integrates polynomials up to degree 63.
c
c
      double precision xl,xu,y,a,b,c,fct
c      external fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      return
      end

      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
c simps
c a1,b1 -the limits of integration
c h1 -an initial step of integration
c reps1,aeps1 - relative and absolute precision of integration
c funct -a name of function subprogram for calculation of integrand +
c x - an argument of the integrand
c ai - the value of integral
c aih- the value of integral with the step of integration
c aiabs- the value of integral for module of the integrand
c this subrogram calculates the definite integral with the relative or
c absolute precision by simpson+s method with the automatical choice
c of the step of integration
c if aeps1    is very small(like 1.e-17),then calculation of integral
c with reps1,and if reps1 is very small (like 1.e-10),then calculation
c of integral with aeps1
c when aeps1=reps1=0. then calculation with the constant step h1
c
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.d0
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end


      subroutine simpt(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpu(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

	subroutine cgln_amps

	implicit none

	include 'mpintp.inc'
	include 'spp.inc'

	integer l

	call multipole_amps

	call legendre

	if (method_helicity.eq.1) then

	ff1 = 0.0
	ff2 = 0.0
	ff3 = 0.0
	ff4 = 0.0
	ff5 = 0.0
	ff6 = 0.0

	do l = 0, wave_L

	  if(l.lt.2) then
	     ff1 = ff1 + (    l*mp(l) +  ep(l))*pol(l+1,1)
	  else
	     ff1 = ff1 + (    l*mp(l) +  ep(l))*pol(l+1,1)
     1		       + ((l+1)*mm(l) +  em(l))*pol(l-1,1)
	  endif
	enddo

	do l = 1, wave_L
	  ff2 = ff2 + ((l+1)*mp(l) + l*mm(l))*pol(l,1)
	  if(l.lt.2) then
	     ff3 = ff3 + (	ep(l) -   mp(l))*pol(l+1,2)
	  else
	     ff3 = ff3 + (	ep(l) -   mp(l))*pol(l+1,2) +
     1			 (	em(l) +   mm(l))*pol(l-1,2)
	  endif
	enddo

	do l = 2, wave_L
	  ff4 = ff4 + (mp(l) - ep(l) -mm(l) - em(l))* pol(l,2)
	enddo

	do l = 0, wave_L
	  if(l.lt.2) then
	     ff5 = ff5 + (l+1)*sp(l)*pol(l+1,1)
	  else
	     ff5 = ff5 + ((l+1)*sp(l)*pol(l+1,1) - l*sm(l)*pol(l-1,1))
	  endif
	enddo

	do l = 1, wave_L
	  ff6 = ff6 + (l*sm(l) - (l+1)*sp(l))*pol(l,1)
	enddo

	else
	do l = 0, max_L+1
	  ap(l)  = 0.0
	  am(l)  = 0.0
	  bp(l)  = 0.0
	  bm(l)  = 0.0
	  cp(l)  = 0.0
	  cm(l)  = 0.0
	enddo

	do l = 0, max_L
	  ap(l)  = 0.5*(l*mp(l) + (l+2)*ep(l))
	  am(l+1)= 0.5*((l+2)*mm(l+1) - l*em(l+1))
	  bp(l)  = ep(l)   - mp(l)
	  bm(l+1)= em(l+1) + mm(l+1)
	  cp(l)  =  (l+1)*sp(l)
	  cm(l+1)= -(l+1)*sm(l+1)
	enddo

	endif

      return
      end

      subroutine helicity_amps

      implicit none

      include 'mpintp.inc'
      include 'spp.inc'

      real   root2,s,s2,c,c2
      real   theta_cm
      integer l

      data root2/1.41421/

      call cgln_amps

      theta_cm = acos(csthcm)

      s   = sin(theta_cm)
      c   = cos(theta_cm)
      s2  = sin(theta_cm/2.0)
      c2  = cos(theta_cm/2.0)

      hh1 = 0.0
      hh2 = 0.0
      hh3 = 0.0
      hh4 = 0.0
      hh5 = 0.0
      hh6 = 0.0

      if(method_helicity.eq.1) then
c	 hh1 = -s*(ff3 + ff4*c)/root2
c	 hh2 = -(2*ff1 - 2*ff2*c + ff4*s2)/root2
c	 hh3 = -ff4*s2/root2
c	 hh4 = s*(2*ff2 + ff3 + ff4*c)/root2
c	 hh5 = ff5 + ff6*c
c	 hh6 = ff6*s
	hh1 = -s*c2*(ff3 + ff4)/root2
	hh2 =  c2*((ff2 - ff1) + 0.5*(1 - c)*(ff3 - ff4))*root2
	hh3 =  s*s2*(ff3 - ff4)/root2
	hh4 =  s2*((ff2 + ff1) + 0.5*(1 + c)*(ff3 + ff4))*root2
	hh5 =  c2*(ff5 + ff6)
	hh6 =  s2*(ff6 - ff5)

      else

	do l = 0, wave_L
	  hh1 = hh1 + s*c2*(bp(l)-bm(l+1))*(pol(l,2)-pol(l+1,2))/root2
	  hh2 = hh2 + c2*(ap(l)-am(l+1))*(pol(l,1)-pol(l+1,1))*root2
	  hh3 = hh3 + s*s2*(bp(l)+bm(l+1))*(pol(l,2)+pol(l+1,2))/root2
	  hh4 = hh4 + s2*(ap(l)+am(l+1))*(pol(l,1)+pol(l+1,1))*root2
	  hh5 = hh5 - c2*(cp(l)-cm(l+1))*(pol(l,1)-pol(l+1,1))
	  hh6 = hh6 - s2*(cp(l)+cm(l+1))*(pol(l,1)+pol(l+1,1))
	enddo

      endif

      end



C------------------------------------------------------------------------------
C INTERP
C
C Does interpolartion of database given two variables V1 and V2

      subroutine interp(v1,v2,vec_f)
C
      implicit none

      include 'mpintp.inc'

C Passed variables:

      real  vec_f(max_mp)
      real  v1,v2

C Local variables:
      real X(2)
      real fint
      integer mpp
C
      if (v1.gt.var1_max) then
c	write(6,'(a,2(1x,e15.8))')' FATAL ERROR: var1 .gt. max - ',v1,var1_max
c	STOP

      endif

      if (v1.lt.var1_min) then
c	write(6,'(a,2(1x,e15.8))')' FATAL ERROR: var1 .lt. min - ',v1,var1_min
c	STOP
	     goto 1
      endif

      if (v2.gt.var2_max) then
c	write(6,'(a,2(1x,e15.8))')' FATAL ERROR: var2 .gt. max - ',v2,var2_max
c	STOP
	 goto 1
      endif

      if (v2.lt.var2_min) then
c	write(6,'(a,2(1x,e15.8))')' FATAL ERROR: var2 .lt. min - ',v2,var2_min
c	STOP
	     goto 1
      endif

      if (method_spline.eq.1) then
C Get the 62 multipoles (F's) at the grid point (V1,V2)
C by 2-dimensional natural cubic spline interpolation.

	call splin2(var1,var2,sf1 ,d2sf1 ,nvar1,nvar2,v1,v2,vec_f(1) )
	call splin2(var1,var2,sf2 ,d2sf2 ,nvar1,nvar2,v1,v2,vec_f(2) )
	call splin2(var1,var2,sf3 ,d2sf3 ,nvar1,nvar2,v1,v2,vec_f(3) )
	call splin2(var1,var2,sf4 ,d2sf4 ,nvar1,nvar2,v1,v2,vec_f(4) )
	call splin2(var1,var2,sf5 ,d2sf5 ,nvar1,nvar2,v1,v2,vec_f(5) )
	call splin2(var1,var2,sf6 ,d2sf6 ,nvar1,nvar2,v1,v2,vec_f(6) )
	call splin2(var1,var2,sf7 ,d2sf7 ,nvar1,nvar2,v1,v2,vec_f(7) )
	call splin2(var1,var2,sf8 ,d2sf8 ,nvar1,nvar2,v1,v2,vec_f(8) )
	call splin2(var1,var2,sf9 ,d2sf9 ,nvar1,nvar2,v1,v2,vec_f(9) )
	call splin2(var1,var2,sf10,d2sf10,nvar1,nvar2,v1,v2,vec_f(10))
	call splin2(var1,var2,sf11,d2sf11,nvar1,nvar2,v1,v2,vec_f(11))
	call splin2(var1,var2,sf12,d2sf12,nvar1,nvar2,v1,v2,vec_f(12))
	call splin2(var1,var2,sf13,d2sf13,nvar1,nvar2,v1,v2,vec_f(13))
	call splin2(var1,var2,sf14,d2sf14,nvar1,nvar2,v1,v2,vec_f(14))
	call splin2(var1,var2,sf15,d2sf15,nvar1,nvar2,v1,v2,vec_f(15))
	call splin2(var1,var2,sf16,d2sf16,nvar1,nvar2,v1,v2,vec_f(16))
	call splin2(var1,var2,sf17,d2sf17,nvar1,nvar2,v1,v2,vec_f(17))
	call splin2(var1,var2,sf18,d2sf18,nvar1,nvar2,v1,v2,vec_f(18))
	call splin2(var1,var2,sf19,d2sf19,nvar1,nvar2,v1,v2,vec_f(19))
	call splin2(var1,var2,sf20,d2sf20,nvar1,nvar2,v1,v2,vec_f(20))
	call splin2(var1,var2,sf21,d2sf21,nvar1,nvar2,v1,v2,vec_f(21))
	call splin2(var1,var2,sf22,d2sf22,nvar1,nvar2,v1,v2,vec_f(22))
	call splin2(var1,var2,sf23,d2sf23,nvar1,nvar2,v1,v2,vec_f(23))
	call splin2(var1,var2,sf24,d2sf24,nvar1,nvar2,v1,v2,vec_f(24))
	call splin2(var1,var2,sf25,d2sf25,nvar1,nvar2,v1,v2,vec_f(25))
	call splin2(var1,var2,sf26,d2sf26,nvar1,nvar2,v1,v2,vec_f(26))
	call splin2(var1,var2,sf27,d2sf27,nvar1,nvar2,v1,v2,vec_f(27))
	call splin2(var1,var2,sf28,d2sf28,nvar1,nvar2,v1,v2,vec_f(28))
	call splin2(var1,var2,sf29,d2sf29,nvar1,nvar2,v1,v2,vec_f(29))
	call splin2(var1,var2,sf30,d2sf30,nvar1,nvar2,v1,v2,vec_f(30))
	call splin2(var1,var2,sf31,d2sf31,nvar1,nvar2,v1,v2,vec_f(31))
	call splin2(var1,var2,sf32,d2sf32,nvar1,nvar2,v1,v2,vec_f(32))
	call splin2(var1,var2,sf33,d2sf33,nvar1,nvar2,v1,v2,vec_f(33))
	call splin2(var1,var2,sf34,d2sf34,nvar1,nvar2,v1,v2,vec_f(34))
	call splin2(var1,var2,sf35,d2sf35,nvar1,nvar2,v1,v2,vec_f(35))
	call splin2(var1,var2,sf36,d2sf36,nvar1,nvar2,v1,v2,vec_f(36))
	call splin2(var1,var2,sf37,d2sf37,nvar1,nvar2,v1,v2,vec_f(37))
	call splin2(var1,var2,sf38,d2sf38,nvar1,nvar2,v1,v2,vec_f(38))
	call splin2(var1,var2,sf39,d2sf39,nvar1,nvar2,v1,v2,vec_f(39))
	call splin2(var1,var2,sf40,d2sf40,nvar1,nvar2,v1,v2,vec_f(40))
	call splin2(var1,var2,sf41,d2sf41,nvar1,nvar2,v1,v2,vec_f(41))
	call splin2(var1,var2,sf42,d2sf42,nvar1,nvar2,v1,v2,vec_f(42))
	call splin2(var1,var2,sf43,d2sf43,nvar1,nvar2,v1,v2,vec_f(43))
	call splin2(var1,var2,sf44,d2sf44,nvar1,nvar2,v1,v2,vec_f(44))
	call splin2(var1,var2,sf45,d2sf45,nvar1,nvar2,v1,v2,vec_f(45))
	call splin2(var1,var2,sf46,d2sf46,nvar1,nvar2,v1,v2,vec_f(46))
	call splin2(var1,var2,sf47,d2sf47,nvar1,nvar2,v1,v2,vec_f(47))
	call splin2(var1,var2,sf48,d2sf48,nvar1,nvar2,v1,v2,vec_f(48))
	call splin2(var1,var2,sf49,d2sf49,nvar1,nvar2,v1,v2,vec_f(49))
	call splin2(var1,var2,sf50,d2sf50,nvar1,nvar2,v1,v2,vec_f(50))
	call splin2(var1,var2,sf51,d2sf51,nvar1,nvar2,v1,v2,vec_f(51))
	call splin2(var1,var2,sf52,d2sf52,nvar1,nvar2,v1,v2,vec_f(52))
	call splin2(var1,var2,sf53,d2sf53,nvar1,nvar2,v1,v2,vec_f(53))
	call splin2(var1,var2,sf54,d2sf54,nvar1,nvar2,v1,v2,vec_f(54))
	call splin2(var1,var2,sf55,d2sf55,nvar1,nvar2,v1,v2,vec_f(55))
	call splin2(var1,var2,sf56,d2sf56,nvar1,nvar2,v1,v2,vec_f(56))
	call splin2(var1,var2,sf57,d2sf57,nvar1,nvar2,v1,v2,vec_f(57))
	call splin2(var1,var2,sf58,d2sf58,nvar1,nvar2,v1,v2,vec_f(58))
	call splin2(var1,var2,sf59,d2sf59,nvar1,nvar2,v1,v2,vec_f(59))
	call splin2(var1,var2,sf60,d2sf60,nvar1,nvar2,v1,v2,vec_f(60))
	call splin2(var1,var2,sf61,d2sf61,nvar1,nvar2,v1,v2,vec_f(61))
	call splin2(var1,var2,sf62,d2sf62,nvar1,nvar2,v1,v2,vec_f(62))
C
      elseif (method_spline.eq.2) then
C Get the 62 multipoles (F's) at the grid point (V1,V2)
C by 2-dimensional natural cubic spline interpolation.
C
	X(1) = v1
	X(2) = v2
C
	vec_f(1 ) = fint(2,X,NA,VAR,sf1 )
	vec_f(2 ) = fint(2,X,NA,VAR,sf2 )
	vec_f(3 ) = fint(2,X,NA,VAR,sf3 )
	vec_f(4 ) = fint(2,X,NA,VAR,sf4 )
	vec_f(5 ) = fint(2,X,NA,VAR,sf5 )
	vec_f(6 ) = fint(2,X,NA,VAR,sf6 )
	vec_f(7 ) = fint(2,X,NA,VAR,sf7 )
	vec_f(8 ) = fint(2,X,NA,VAR,sf8 )
	vec_f(9 ) = fint(2,X,NA,VAR,sf9 )
	vec_f(10) = fint(2,X,NA,VAR,sf10)
	vec_f(11) = fint(2,X,NA,VAR,sf11)
	vec_f(12) = fint(2,X,NA,VAR,sf12)
	vec_f(13) = fint(2,X,NA,VAR,sf13)
	vec_f(14) = fint(2,X,NA,VAR,sf14)
	vec_f(15) = fint(2,X,NA,VAR,sf15)
	vec_f(16) = fint(2,X,NA,VAR,sf16)
	vec_f(17) = fint(2,X,NA,VAR,sf17)
	vec_f(18) = fint(2,X,NA,VAR,sf18)
	vec_f(19) = fint(2,X,NA,VAR,sf19)
	vec_f(20) = fint(2,X,NA,VAR,sf20)
	vec_f(21) = fint(2,X,NA,VAR,sf21)
	vec_f(22) = fint(2,X,NA,VAR,sf22)
	vec_f(23) = fint(2,X,NA,VAR,sf23)
	vec_f(24) = fint(2,X,NA,VAR,sf24)
	vec_f(25) = fint(2,X,NA,VAR,sf25)
	vec_f(26) = fint(2,X,NA,VAR,sf26)
	vec_f(27) = fint(2,X,NA,VAR,sf27)
	vec_f(28) = fint(2,X,NA,VAR,sf28)
	vec_f(29) = fint(2,X,NA,VAR,sf29)
	vec_f(30) = fint(2,X,NA,VAR,sf30)
	vec_f(31) = fint(2,X,NA,VAR,sf31)
	vec_f(32) = fint(2,X,NA,VAR,sf32)
	vec_f(33) = fint(2,X,NA,VAR,sf33)
	vec_f(34) = fint(2,X,NA,VAR,sf34)
	vec_f(35) = fint(2,X,NA,VAR,sf35)
	vec_f(36) = fint(2,X,NA,VAR,sf36)
	vec_f(37) = fint(2,X,NA,VAR,sf37)
	vec_f(38) = fint(2,X,NA,VAR,sf38)
	vec_f(39) = fint(2,X,NA,VAR,sf39)
	vec_f(40) = fint(2,X,NA,VAR,sf40)
	vec_f(41) = fint(2,X,NA,VAR,sf41)
	vec_f(42) = fint(2,X,NA,VAR,sf42)
	vec_f(43) = fint(2,X,NA,VAR,sf43)
	vec_f(44) = fint(2,X,NA,VAR,sf44)
	vec_f(45) = fint(2,X,NA,VAR,sf45)
	vec_f(46) = fint(2,X,NA,VAR,sf46)
	vec_f(47) = fint(2,X,NA,VAR,sf47)
	vec_f(48) = fint(2,X,NA,VAR,sf48)
	vec_f(49) = fint(2,X,NA,VAR,sf49)
	vec_f(50) = fint(2,X,NA,VAR,sf50)
	vec_f(51) = fint(2,X,NA,VAR,sf51)
	vec_f(52) = fint(2,X,NA,VAR,sf52)
	vec_f(53) = fint(2,X,NA,VAR,sf53)
	vec_f(54) = fint(2,X,NA,VAR,sf54)
	vec_f(55) = fint(2,X,NA,VAR,sf55)
	vec_f(56) = fint(2,X,NA,VAR,sf56)
	vec_f(57) = fint(2,X,NA,VAR,sf57)
	vec_f(58) = fint(2,X,NA,VAR,sf58)
	vec_f(59) = fint(2,X,NA,VAR,sf59)
	vec_f(60) = fint(2,X,NA,VAR,sf60)
	vec_f(61) = fint(2,X,NA,VAR,sf61)
	vec_f(62) = fint(2,X,NA,VAR,sf62)
      endif
C
      return
C
 1    do mpp = 1,max_mp
	 vec_f(mpp) = 0.0
      enddo
      print *, 'v1= ',v1,' Q2 min= ',var1_min,' Q2 max= ',var1_max
      print *, 'v2= ',v2,' W  min= ',var2_min,' W  max= ',var2_max
      stop
C
      return
      end
c---------------------------------------------------------------------------
c	Calculatess the first and second derivatives of
c	Legendre Polynomials.
c---------------------------------------------------------------------------
      subroutine legendre

      implicit none

      include 'mpintp.inc'
      include 'spp.inc'

      real X
      integer l
c
      X = csthcm
c
      pol(0,1) = 0.
      pol(1,1) = 1.
      pol(2,1) = 3.*X
      pol(3,1) = (15.*X*X-3.)/2.
      pol(4,1) = (35.*X*X*X-15.*X)/2.
      pol(5,1) = (315.*X*X*X*X-210.*X*X+15.)/8.
      pol(6,1) = (693.*X*X*X*X*X-630.*X*X*X +105.*X)/8.
      pol(7,1) = (3003.*X*X*X*X*X*X-3465.*X*X*X*X+945*X*X-35.)/16.
      pol(0,2) = 0.
      pol(1,2) = 0.
      pol(2,2) = 3.
      pol(3,2) = 15.*X
      pol(4,2) = (105.*X*X-15.)/2.
      pol(5,2) = (315.*X*X*X-105.*X)/2.
      pol(6,2) = (3465.*X*X*X*X-1890.*X*X+105.)/8.
      pol(7,2) = (9009.*X*X*X*X*X-6930.*X*X*X +945.*X)/8.
c
c      do l = 0, 6
c	   print *, 'legedre', 'l=',l
c	   print *, pol(l,1), pol(l,2)
c      enddo
      return
      end

      subroutine maid_lee(iphy,a1,a2,a3,a4,a5,channel_opt,
     1		 sigma0,sig_t,sig_tt,sig_l,sig_lt,sig_ltp)
      implicit none
c
      real a1,a2,a3,a4,a5
c
      include 'mpintp.inc'
      include 'spp.inc'

      integer iphy, channel_opt, n_call
      real E_pi_cm
      real m_pi, alpi
      real sigma0, sig_t, sig_tt, sig_l, sig_lt,sig_ltp
c
      data n_call /0/
c
      m_pi  = m_pip
      if (channel_opt.eq.1) m_pi  = m_pi0
      e_hel = 0.0
c
      q2	  = a1
      w 	  = a2
      epsilon	  = a3
      csthcm	  = a4
      phicm	  = a5
c
      if (n_call.ne.0) goto 100
c
      method_helicity = 1
c
      mwave_L = 5
c
      if (iphy.eq.2) then
	if (channel_opt.eq.1) then
	  data_file = 'maid98-PPpi.tbl'
	elseif (channel_opt.eq.3) then
	  data_file = 'maid98-PNpi.tbl'
	endif
      endif
      if (iphy.eq.3) then
	if (channel_opt.eq.1) then
	  data_file = 'maid2000-PPpi.tbl'
	elseif (channel_opt.eq.3) then
	  data_file = 'maid2000-PNpi.tbl'
	endif
      endif
c
      print *, 'iphy=', iphy
      write(6,*) 'Reading multipoles from ',data_file

      call read_sf_file(data_file,2)

      n_call = n_call + 1

 100  E_pi_cm	 = 0.5*(W*W+m_pi**2-m_p**2)/W
      ppi_mag_cm = E_pi_cm**2 - m_pi**2
      ppi_mag_cm = sqrt(ppi_mag_cm)
      qv_mag_cm  = ((W*W+Q2+m_p**2)/2.0/W)**2-m_p**2
      qv_mag_cm  = sqrt(qv_mag_cm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     I found this is not good formula for nu_cm (K. Joo)
c     nu_cm  = qv_mag_cm**2-Q2
c      if (nu_cm.le.0) then
c	sigma0 = 0.0
c	sig_t  = 0.0
c	sig_l  = 0.0
c	sig_tt = 0.0
c	sig_lt = 0.0
c	return
c      endif
c      nu_cm	 = sqrt(nu_cm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nu_cm  = (W*W-m_p**2-Q2)/(2*W)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      fkt = 2.0*W*ppi_mag_cm/(W**2 - m_p**2)

      call xsection

      sigma0 = sigma_0
      sig_t  = sigma_t
      sig_l  = sigma_l
      sig_tt = sigma_tt
      sig_lt = sigma_lt
      sig_ltp = sigma_ltp

c      print *, 'MAID-LEE',sigma0,sigu,sigl,sigt,sigi

      return
      end
C------------------------------------------------------------------------------
C multipole_amps
C
C Fills multipole amplitudes using results from interpolation
C------------------------------------------------------------------------------
	subroutine multipole_amps

	implicit none

	include 'mpintp.inc'
	include 'spp.inc'

	real vec_f(max_mp)
	real factor
	integer i

	call interp(Q2,W,vec_f)
	factor = 0.141383e0	! Convert from 10**-3/m_pi+ to sqrt(microbarns)
c
c	nu_cm * S = qv_mag_cm * L
c
c We have to multiply nu_cm / qv_mag_cm because we are following Tiator's notation
c which uses L-multipole rather than S-multipole.
c Note that Sato-Lee uses S-multipole.
c
	    factor = 0.141383d0 * nu_cm / qv_mag_cm
c	 factor = nu_cm / qv_mag_cm
c	 factor = 0.141383d0
	sp(0) = cmplx(vec_F(1), vec_F(2) )*factor
	sp(1) = cmplx(vec_F(3), vec_F(4) )*factor
	sp(2) = cmplx(vec_F(5), vec_F(6) )*factor
	sp(3) = cmplx(vec_F(7), vec_F(8) )*factor
	sp(4) = cmplx(vec_F(9), vec_F(10))*factor
	sp(5) = cmplx(vec_F(11),vec_F(12))*factor

	sm(0) = 0.0
	sm(1) = cmplx(vec_F(13),vec_F(14))*factor
	sm(2) = cmplx(vec_F(15),vec_F(16))*factor
	sm(3) = cmplx(vec_F(17),vec_F(18))*factor
	sm(4) = cmplx(vec_F(19),vec_F(20))*factor
	sm(5) = cmplx(vec_F(21),vec_F(22))*factor

	factor = 0.141383d0
c
c	 factor = 1.0
	ep(0) = cmplx(vec_F(23),vec_F(24))*factor
	ep(1) = cmplx(vec_F(25),vec_F(26))*factor
	ep(2) = cmplx(vec_F(27),vec_F(28))*factor
	ep(3) = cmplx(vec_F(29),vec_F(30))*factor
	ep(4) = cmplx(vec_F(31),vec_F(32))*factor
	ep(5) = cmplx(vec_F(33),vec_F(34))*factor

	    em(0) = 0.0
		em(1) = 0.0
	em(2) = cmplx(vec_F(35),vec_F(36))*factor
	em(3) = cmplx(vec_F(37),vec_F(38))*factor
	em(4) = cmplx(vec_F(39),vec_F(40))*factor
	em(5) = cmplx(vec_F(41),vec_F(42))*factor

		mp(0) = 0.0
	mp(1) = cmplx(vec_F(43),vec_F(44))*factor
	mp(2) = cmplx(vec_F(45),vec_F(46))*factor
	mp(3) = cmplx(vec_F(47),vec_F(48))*factor
	mp(4) = cmplx(vec_F(49),vec_F(50))*factor
	mp(5) = cmplx(vec_F(51),vec_F(52))*factor

		mm(0) = 0.0
	mm(1) = cmplx(vec_F(53),vec_F(54))*factor
	mm(2) = cmplx(vec_F(55),vec_F(56))*factor
	mm(3) = cmplx(vec_F(57),vec_F(58))*factor
	mm(4) = cmplx(vec_F(59),vec_F(60))*factor
	mm(5) = cmplx(vec_F(61),vec_F(62))*factor
c
c	 print *, 'multipoles ------'
c	 print *, (sp(i),sm(i),ep(i),em(i),mp(i),mm(i),i=0,5)
	return
	end
C-------------------------------------------------------------------------------
C READ_DATA_FILE
C
C Reads data file for use in interpolation
C
CEV Both these routines operate in double precision arithmetic
C   to suppress round off errors in the interpolation procedure.
C   But INTERP returns single precision results.

      subroutine read_sf_file(rfile,iunit)

      implicit none

      include 'mpintp.inc'
      include 'spp.inc'
C Passed variables:

      character rfile*(*)
      integer iunit

C Local variables:

      integer jvar1,jvar2
      real dumvar1,dumvar2,dumvar3,dumvar4
      character*80 dummy
      character*2 dummyX
      character*5 dummyY
      real var1tmp1,var1tmp2
      real var2tmp1,var2tmp2

c      write(6,*) 'Input method of spline fitting (1=cubic, 2=linear)'
c      read(5,*) method_spline
      method_spline = 2

C First open the file.
      open(unit=iunit,file=rfile,form='formatted',status='old')

C Begin reading the file.

      do jvar1=1,nvar1	       ! Q2
	 do jvar2=1,nvar2      ! W
	   read(iunit,FMT=15,err=1000) dummyX,	var2(jvar2),
     1				       dummyY,	var1(jvar1)   ! w,q2
 15	   format(A8,f4.2,A7,f7.5)
	   if (jvar1.eq.1)     var1tmp1=var1(jvar1)
	   if (jvar1.eq.nvar1) var1tmp2=var1(jvar1)
	   if (jvar2.eq.1)     var2tmp1=var2(jvar2)
	   if (jvar2.eq.nvar2) var2tmp2=var2(jvar2)
C SL+
	   read(iunit,'(a)',err=1000) dummy
	   read(iunit,*,err=1000)  sf1(jvar1,jvar2), sf2(jvar1,jvar2),
     &				   sf3(jvar1,jvar2), sf4(jvar1,jvar2),
     &				   sf5(jvar1,jvar2), sf6(jvar1,jvar2)
	   read(iunit,*,err=1000)  sf7(jvar1,jvar2), sf8(jvar1,jvar2),
     &				   sf9(jvar1,jvar2),sf10(jvar1,jvar2),
     &				  sf11(jvar1,jvar2),sf12(jvar1,jvar2)
C SL-
	   read(iunit,'(a)',err=1000) dummy
	   read(iunit,*,err=1000) dumvar1	   ,dumvar2,
     &				  sf13(jvar1,jvar2),sf14(jvar1,jvar2),
     &				  sf15(jvar1,jvar2),sf16(jvar1,jvar2)
	   read(iunit,*,err=1000) sf17(jvar1,jvar2),sf18(jvar1,jvar2),
     &				  sf19(jvar1,jvar2),sf20(jvar1,jvar2),
     &				  sf21(jvar1,jvar2),sf22(jvar1,jvar2)
C EL+
	   read(iunit,'(a)',err=1000) dummy
	   read(iunit,*,err=1000) sf23(jvar1,jvar2),sf24(jvar1,jvar2),
     &				  sf25(jvar1,jvar2),sf26(jvar1,jvar2),
     &				  sf27(jvar1,jvar2),sf28(jvar1,jvar2)
	   read(iunit,*,err=1000) sf29(jvar1,jvar2),sf30(jvar1,jvar2),
     &				  sf31(jvar1,jvar2),sf32(jvar1,jvar2),
     &				  sf33(jvar1,jvar2),sf34(jvar1,jvar2)
C EL-
	   read(iunit,'(a)',err=1000) dummy
	   read(iunit,*,err=1000) dumvar1	  ,dumvar2,
     &				  dumvar3	  ,dumvar4,
     &				  sf35(jvar1,jvar2),sf36(jvar1,jvar2)
	   read(iunit,*,err=1000) sf37(jvar1,jvar2),sf38(jvar1,jvar2),
     &				  sf39(jvar1,jvar2),sf40(jvar1,jvar2),
     &				  sf41(jvar1,jvar2),sf42(jvar1,jvar2)
C ML+
	   read(iunit,'(a)',err=1000) dummy
	   read(iunit,*,err=1000) dumvar1	   ,dumvar2,
     &				  sf43(jvar1,jvar2),sf44(jvar1,jvar2),
     &				  sf45(jvar1,jvar2),sf46(jvar1,jvar2)
	   read(iunit,*,err=1000) sf47(jvar1,jvar2),sf48(jvar1,jvar2),
     &				  sf49(jvar1,jvar2),sf50(jvar1,jvar2),
     &				  sf51(jvar1,jvar2),sf52(jvar1,jvar2)
C ML-
	   read(iunit,'(a)',err=1000) dummy
	   read(iunit,*,err=1000) dumvar1	   ,dumvar2,
     &				  sf53(jvar1,jvar2),sf54(jvar1,jvar2),
     &				  sf55(jvar1,jvar2),sf56(jvar1,jvar2)
	   read(iunit,*,err=1000) sf57(jvar1,jvar2),sf58(jvar1,jvar2),
     &				  sf59(jvar1,jvar2),sf60(jvar1,jvar2),
     &				  sf61(jvar1,jvar2),sf62(jvar1,jvar2)
	 enddo
      enddo

C Close the input file.

      close(iunit)

C Determine minimum/maximum of variables

      if (var1tmp1.gt.var1tmp2) then
	     var1_max = var1tmp1
	     var1_min = var1tmp2
      else
	     var1_max = var1tmp2
	     var1_min = var1tmp1
      endif

      if (var2tmp1.gt.var2tmp2) then
	 var2_max = var2tmp1
	 var2_min = var2tmp2
      else
	 var2_max = var2tmp2
	 var2_min = var2tmp1
      endif

C Now construct natural cubic splines in the 2nd dimension (VAR2) and
C calculate the 2nd derivatives wrt. VAR2 (2nd entry of each D2SF) of
C each structure function (SF). This is done only once, so do it here.
      if(method_spline.eq.1) then
      call splie2(var2,sf1, nvar1,nvar2,d2sf1 )
      call splie2(var2,sf2, nvar1,nvar2,d2sf2 )
      call splie2(var2,sf3, nvar1,nvar2,d2sf3 )
      call splie2(var2,sf4, nvar1,nvar2,d2sf4 )
      call splie2(var2,sf5, nvar1,nvar2,d2sf5 )
      call splie2(var2,sf6, nvar1,nvar2,d2sf6 )
      call splie2(var2,sf7, nvar1,nvar2,d2sf7 )
      call splie2(var2,sf8, nvar1,nvar2,d2sf8 )
      call splie2(var2,sf9, nvar1,nvar2,d2sf9 )
      call splie2(var2,sf10,nvar1,nvar2,d2sf10)
      call splie2(var2,sf11,nvar1,nvar2,d2sf11)
      call splie2(var2,sf12,nvar1,nvar2,d2sf12)
      call splie2(var2,sf13,nvar1,nvar2,d2sf13)
      call splie2(var2,sf14,nvar1,nvar2,d2sf14)
      call splie2(var2,sf15,nvar1,nvar2,d2sf15)
      call splie2(var2,sf16,nvar1,nvar2,d2sf16)
      call splie2(var2,sf17,nvar1,nvar2,d2sf17)
      call splie2(var2,sf18,nvar1,nvar2,d2sf18)
      call splie2(var2,sf19,nvar1,nvar2,d2sf19)
      call splie2(var2,sf20,nvar1,nvar2,d2sf20)
      call splie2(var2,sf21,nvar1,nvar2,d2sf21)
      call splie2(var2,sf22,nvar1,nvar2,d2sf22)
      call splie2(var2,sf23,nvar1,nvar2,d2sf23)
      call splie2(var2,sf24,nvar1,nvar2,d2sf24)
      call splie2(var2,sf25,nvar1,nvar2,d2sf25)
      call splie2(var2,sf26,nvar1,nvar2,d2sf26)
      call splie2(var2,sf27,nvar1,nvar2,d2sf27)
      call splie2(var2,sf28,nvar1,nvar2,d2sf28)
      call splie2(var2,sf29,nvar1,nvar2,d2sf29)
      call splie2(var2,sf30,nvar1,nvar2,d2sf30)
      call splie2(var2,sf31,nvar1,nvar2,d2sf31)
      call splie2(var2,sf32,nvar1,nvar2,d2sf32)
      call splie2(var2,sf33,nvar1,nvar2,d2sf33)
      call splie2(var2,sf34,nvar1,nvar2,d2sf34)
      call splie2(var2,sf35,nvar1,nvar2,d2sf35)
      call splie2(var2,sf36,nvar1,nvar2,d2sf36)
      call splie2(var2,sf37,nvar1,nvar2,d2sf37)
      call splie2(var2,sf38,nvar1,nvar2,d2sf38)
      call splie2(var2,sf39,nvar1,nvar2,d2sf39)
      call splie2(var2,sf40,nvar1,nvar2,d2sf40)
      call splie2(var2,sf41,nvar1,nvar2,d2sf41)
      call splie2(var2,sf42,nvar1,nvar2,d2sf42)
      call splie2(var2,sf43,nvar1,nvar2,d2sf43)
      call splie2(var2,sf44,nvar1,nvar2,d2sf44)
      call splie2(var2,sf45,nvar1,nvar2,d2sf45)
      call splie2(var2,sf46,nvar1,nvar2,d2sf46)
      call splie2(var2,sf47,nvar1,nvar2,d2sf47)
      call splie2(var2,sf48,nvar1,nvar2,d2sf48)
      call splie2(var2,sf49,nvar1,nvar2,d2sf49)
      call splie2(var2,sf50,nvar1,nvar2,d2sf50)
      call splie2(var2,sf51,nvar1,nvar2,d2sf51)
      call splie2(var2,sf52,nvar1,nvar2,d2sf52)
      call splie2(var2,sf53,nvar1,nvar2,d2sf53)
      call splie2(var2,sf54,nvar1,nvar2,d2sf54)
      call splie2(var2,sf55,nvar1,nvar2,d2sf55)
      call splie2(var2,sf56,nvar1,nvar2,d2sf56)
      call splie2(var2,sf57,nvar1,nvar2,d2sf57)
      call splie2(var2,sf58,nvar1,nvar2,d2sf58)
      call splie2(var2,sf59,nvar1,nvar2,d2sf59)
      call splie2(var2,sf60,nvar1,nvar2,d2sf60)
      call splie2(var2,sf61,nvar1,nvar2,d2sf61)
      call splie2(var2,sf62,nvar1,nvar2,d2sf62)
      elseif(method_spline.eq.2) then
C****************************************************************
C the preparation for linear spline fitting
C
      do jvar1=1,nvar1	       ! Q2
	var(jvar1) = var1(jvar1)
      enddo
      do jvar2= nvar1+1,nvar1+nvar2	    ! W
	var(jvar2) = var2(jvar2-nvar1)
      enddo
C****************************************************************
      endif
      return

c 1000 STOP ' Error reading the structure function input file.'
 1000 print *, 'Error reading ', 'Q2= ',(0.2+0.05*(jvar1-1)),
     1	       'W= ', (1.1+0.01*(jvar2-1))
      stop
      end
C------------------------------------------------------------------------
CEV The following two routines concern natural cubic spline interpolation
C   of a function of two independent variables.
C------------------------------------------------------------------------

      SUBROUTINE splie2(x2a,ya,m,n,y2a)
C
C Given an m by n tabulated function ya(1:m,1:n), and tabulated independent
C variable x2a(1:n), this routine constructs 1-dimensional natural cubic
C splines of the 2nd dimension (1:n) of ya and returns the 2nd-derivatives in
C the array y2a(1:m,1:n).  From "Numerical Recipes" (FORTRAN, 1992 Ed.),
C Ch. 3.6, p. 121.  Call this routine only once, before a sequence of calls to
C SPLIN2.  The process is of order  m x n.
      IMPLICIT NONE
C Passed variables:
      INTEGER m,n
      REAL x2a(n),ya(m,n),y2a(m,n)
C Local variables:
      INTEGER i,j,NN
      REAL bound
      PARAMETER (NN=100)     ! Maximum expected value of n
      REAL ytmp(NN),y2tmp(NN)
C Boundary condition threshold for natural cubic spline:
      SAVE bound
      DATA bound/1.e30/
C
      do i=1,m
	 do j=1,n
	    ytmp(j)=ya(i,j)
	 enddo
C Natural spline:
	 call spline(x2a,ytmp,n,bound,bound,y2tmp)
	 do j=1,n
	    y2a(i,j)=y2tmp(j)
	 enddo
      enddo
C
      return
      END
C------------------------------------------------------------------------

      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
C
C Given an m by n tabulated function ya(1:m,1:n), and tabulated independent
C variables x1a(1:m), x2a(1:n), and tabulated 2nd-derivatives y2a(1:m,1:n) as
C produced by SPLIE2; and given a desired interpolating point x1, x2; this
C routine returns an interpolated function value y by 2-dimensional natural
C cubic spline interpolation.  From "Numerical Recipes" (FORTRAN, 1992 Ed.),
C Ch. 3.6, p. 121.  Call this routine sequentially after one call to SPLIE2.
C The process is of order  m x logn + m + logm + 1 ; i.e., roughly, of order
C m x logn.  Therefore, the optimum choice is n > m.
      IMPLICIT NONE
C Passed variables:
      INTEGER m,n
      REAL  x1a(m),x2a(n),ya(m,n),y2a(m,n),x1,x2,y
C Local variables:
      INTEGER i,j,NN
      REAL bound
      PARAMETER (NN=100)	   ! Maximum expected value of m, n
      REAL ytmp(NN),yytmp(NN),y2tmp(NN)
C Boundary condition threshold for natural cubic spline:
      SAVE bound
      DATA bound/1.e30/
C
      do i=1,m
	 do j=1,n
	    ytmp(j)=ya(i,j)
	    y2tmp(j)=y2a(i,j)
	 enddo
C Perform n evaluations of the splines of the (1:n)-dimension constructed by
C SPLIE2, using the 1-dimensional spline evaluator SPLINT.
	 call splint(x2a,ytmp,y2tmp,n,x2,yytmp(i))
      enddo
C Construct the 1-dimensional natural spline of the (1:m)-dimension and
C evaluate it.
      call spline(x1a,yytmp,m,bound,bound,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
C
      return
      END
C------------------------------------------------------------------------
CEV The subroutines SPLINE and SPLINT can be used for interpolating a
C   function of one independent variable.
C------------------------------------------------------------------------

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
C
C Given arrays y(1:n) and x(1:n) containing a tabulated function and its
C independent variable, respectively, with x(1) < x(2) < ... < x(n), and given
C values yp1 and ypn for the 1st derivative of the interpolating function at
C the points 1 and n, respectively, this routine returns an array y2(1:n) of
C length n which contains the 2nd derivative of the interpolating function at
C the tabulated points x(i).  If yp1 and/or ypn are equal to 10^30 or larger,
C the routine is signaled to set the corresponding boundary condition for a
C natural spline, with zero 2nd derivative on that boundary.  From "Numerical
C Recipes", Ch. 3.3, p. 109.
      IMPLICIT NONE
C Passed variables:
      INTEGER n
      REAL x(n),y(n),yp1,ypn,y2(n)
C Local variables:
      INTEGER i,j,NMAX
      PARAMETER (NMAX=500)	  ! Maximum expected value of n
      REAL p,qn,sig,un,u(NMAX)
C
      if (yp1.gt..99e30) then
C Natural lower boundary condition:
	 y2(1)=0.e0
	 u(1)=0.e0
      else
C Specified 1st derivative yp1:
	 y2(1)=-.5e0
	 u(1)=(3.e0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
C We need at least N+1 points to construct an interpolating polynomial
C of degree N:
      if (n.lt.4) stop ' Too few points for cubic spline'
C
C Decomposition loop of the tridiagonal algorithm.
C y2,u: temporary storage of the decomposed factors.
      do i=2,n-1
	 sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
	 p=sig*y2(i-1)+2.e0
	 y2(i)=(sig-1.e0)/p
	 u(i)=(6.e0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &	      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30) then
C Natural upper boundary condition:
	 qn=0.e0
	 un=0.e0
      else
C Specified 1st derivative ypn:
	 qn=.5e0
	 un=(3.e0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.e0)
C
C Backsubstitution loop of the tridiagonal algorithm.
      do j=n-1,1,-1
	 y2(j)=y2(j)*y2(j+1)+u(j)
      enddo
C
      return
      END
C------------------------------------------------------------------------

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
C
C Given the arrays ya(1:n) and xa(1:n) of length n, which tabulate a function
C and its independent variable (with the xa(i) in order, as in SPLINE),
C respectively, and given the array y2a(1:n), which is the output from SPLINE,
C and given a value of x, this routine returns a cubic-spline interpolated
C value y.  From "Numerical Recipes", Ch. 3.3, p. 110.
      IMPLICIT NONE
C Passed variables:
      INTEGER n
      REAL xa(n),ya(n),y2a(n),x,y
C Local variables:
      INTEGER k,klo,khi
      REAL a,b,h
C
      klo=1
      khi=n
C
C Locate the interpolation point by means of bisection.
C This is optimal if sequential calls to this routine are at random values
C of x. If sequential calls are in order, and closely spaced, it would be
C better to store previous values of klo and khi and test if they remain
C appropriate on the next call.
 1    if (khi-klo.gt.1) then
	 k=(khi+klo)/2
	 if (xa(k).gt.x) then
	    khi=k
	 else
	    klo=k
	 endif
	 go to 1
      endif
C klo,khi now bracket x.
      h=xa(khi)-xa(klo)
C The xa's must be distinct:
      if (h.eq.0.e0) stop ' Bad xa input in splint'
C Evaluate now the cubic spline polynomial:
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &	((a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi))*(h*h)/6.e0
C
      return
      END
      subroutine xsection

      implicit none

      include 'mpintp.inc'
      include 'spp.inc'

      real vl, vt, vlt, vtt, vltp
      real phicm_rad, ekin
      real sig_l,sig_t,sig_lt,sig_tt,sig_ltp

      call multipole_amps ! calc. multipole amplitudes
      call helicity_amps  ! calc. helicity amplitudes

      phicm_rad = phicm*pi/180.0

      vt   = 1.0
      vl   = epsilon
      vtt  = epsilon
      vlt  = sqrt(epsilon*(1+epsilon)/2)
      vltp = sqrt(epsilon*(1-epsilon)/2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this is relevant chagne due to nu_cm
c      ekin = Q2/nu_cm**2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ekin = sqrt(Q2) / nu_cm
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      fkt  = 2*W*ppi_mag_cm/(W**2 - m_p**2)

      sigma_t  = (cabs(hh1)**2+cabs(hh2)**2
     1		 +cabs(hh3)**2+cabs(hh4)**2)/2.
      sigma_l  =  cabs(hh5)**2+cabs(hh6)**2
      sigma_tt =  real(-conjg(hh1)*hh4+conjg(hh2)*hh3)
      sigma_lt =  real(conjg(hh5)*(hh1-hh4)  +
     1			conjg(hh6)*(hh2+hh3))*sqrt(2.0)
      sigma_ltp = -aimag(conjg(hh5)*(hh1-hh4) +
     1			conjg(hh6)*(hh2+hh3))*sqrt(2.0)

      sigma_l  = sigma_l*ekin**2
      sigma_t  = sigma_t
      sigma_tt = sigma_tt
      sigma_lt = sigma_lt*ekin
      sigma_ltp= sigma_ltp*ekin

      sigma_0  = fkt*( vt*sigma_t   + vl*sigma_l +
     1		 vtt*sigma_tt*cos(2*phicm_rad) +
     1		 vlt*sigma_lt*cos(  phicm_rad) ) +
     1	     e_hel*vltp*sigma_ltp*sin(phicm_rad)

c      print *, 'XSECTION',sigma_0,sigma_t,sigma_l,sigma_tt,sigma_lt

      sigma_0_5fold = gamma_v*sigma_0

      end


      subroutine d01fce(ndim, a, b, minpts, maxpts, functn, eps,
     * acc, lenwrk, wrkstr, finval, ifail)
      implicit real*8(a-h,o-z)
c     mark 8 release. nag copyright 1979.
c
c     adaptive multidimensional integration subroutine
c
c     *********  parameters for d01fce ****************************
c
c      input parameters
c
c     ndim    integer number of variables, must exceed 1 but
c	  not exceed 15.
c
c     a       real array of lower limits, with dimension ndim
c
c     b       real array of upper limits, with dimension ndim
c
c     minpts  integer minimum number of integrand values to be
c	  allowed, which must not exceed maxpts.
c
c     maxpts  integer maximum number of integrand values to be
c	  allowed, which must be at least
c	  2**ndim+2*ndim**2+2*ndim+1.
c
c     functn  externally declared user defined real function
c	  integrand. it must have parameters (ndim,z),
c	  where z is a real array of dimension ndim.
c
c     eps     real required relative accuracy, must be greater
c	  than zero
c
c     lenwrk  integer length of array wrkstr, must be at least
c	  2*ndim+4.
c
c     ifail   integer nag failure parameter
c	  ifail=0 for hard fail
c	  ifail=1 for soft fail
c
c      output parameters
c
c     minpts  integer number of integrand values used by the
c	  routine
c
c     wrkstr  real array of working storage of dimension (lenwrk).
c
c     acc     real estimated relative accuracy of finval
c
c     finval  real estimated value of integral
c
c     ifail   ifail=0 for normal exit, when estimated relative
c	    less integaccuracy rand values used.
c
c      ifail=1 if ndim.lt.2, ndim.gt.15, minpts.gt.maxpts,
c	    maxpts.lt.2**ndim+2*ndim*(ndim+1)+1, eps.le.0
c	    or lenwrk.lt.2*ndim+4.
c
c      ifail=2 if maxpts was too small for d01fce to obtain the
c	    required relative accuracy eps.  in this
c	    case d01fce returns a value of finval
c	    with estimated relative accuracy acc.
c
c      ifail=3 if lenwrk too small for maxpts integrand
c	    values.  in this case d01fce returns a
c	    value of finval with estimated accuracy
c	    acc using the working storage
c	    available, but acc will be greater
c	    than eps.
c
c     **************************************************************
c
c     .. scalar arguments ..
*     real eps, finval, acc
*     integer ifail, lenwrk, maxpts, minpts, ndim
c     .. array arguments ..
**    real a(ndim), b(ndim), wrkstr(lenwrk)
      dimension a(ndim), b(ndim), wrkstr(lenwrk)
c     .. function arguments ..
*     real functn
c     ..
c     .. local scalars ..
      character*8 srname
      double precision
     * abserr, df1, df2, difmax, f1, f2, f3, f4, half, lamda2,
     * lamda4, lamda5, one, ratio, rgncmp, rgnerr, rgnert, rgnval,
     * rgnvlt, rgnvol, rlndim, sum1, sum2, sum3, sum4, sum5, two,
     * twondm, weit1, weit2, weit3, weit4, weit5, weitp1, weitp2,
     * weitp3, weitp4, zero
      integer dvaxes, dvaxis, dvflag, funcls, ierror, j, k, maxaxs,
     * mxrgns, pointr, rgncls, rulcls, sbrgns, subrgn, subtmp,
     * tpontp, tpontr
c     .. local arrays ..
      dimension center(15), dif(15), oldcnt(15), width(15), z(15)
      integer dvcntl(15), dvcntr(15)
c     .. function references ..
*     real sqrt, x02aae
      integer p01aae, x02bbe
c     ..
      data srname /'  d01fce'/
      data zero, one, two, half /0.d0, 1.d0, 2.d0, 0.5d0/
c
c   subroutine initialisation and parameter checking
c
      if (ndim.lt.2 .or. ndim.gt.15) go to 560
      if (minpts.gt.maxpts) go to 560
      if (eps.le.zero) go to 560
      if (lenwrk.lt.2*ndim+4) go to 560
      funcls = 0
      finval = zero
      abserr = zero
      twondm = two**ndim
      rgnvol = twondm
      dvflag = 1
      fffff1 = float(x02bbe(one))
      fffff2 = 1.0/x02aae(0.0d0)
      maxaxs = int(dmin1(fffff1,fffff2))
c     maxaxs = int(amin1(float(x02bbe(one)),1.0/x02aae(0.0d0)))
      maxaxs = (maxaxs-ndim)/(ndim+1)
      mxrgns = lenwrk/(2*ndim+4)
      sbrgns = 0
      rgnvlt = zero
      rgnert = zero
      do 20 j=1,ndim
       center(j) = (a(j)+b(j))*half
       dif(j) = zero
       width(j) = (b(j)-a(j))*half
       dvcntl(j) = 1
       dvcntr(j) = 1
       oldcnt(j) = center(j)
       rgnvol = rgnvol*width(j)
   20 continue
c
c   end subroutine initialisation
c   basic rule initialisation
c
      rulcls = 2**ndim + 2*ndim*ndim + 2*ndim + 1
      funcls = rulcls
      if (maxpts.lt.rulcls) go to 560
      rlndim = ndim
      lamda2 = sqrt(9.0/70.0)
      lamda4 = sqrt(9.0/10.0)
      lamda5 = sqrt(9.0/19.0)
      weit1 = (12824.0-9120.0*rlndim+400.0*rlndim*rlndim)/19683.0
      weit2 = 980.0/6561.0
      weit3 = (1820.0-400.0*rlndim)/19683.0
      weit4 = 200.0/19683.0
      weit5 = 6859.0/19683.0/twondm
      weitp1 = (729.0-950.0*rlndim+50.0*rlndim**2)/729.0
      weitp2 = 245.0/486.0
      weitp3 = (265.0-100.0*rlndim)/1458.0
      weitp4 = 25.0/729.0
      ratio = (lamda2/lamda4)**2
c
c   end basic rule initialisation
      go to 100
c   divide subregion with largest error and prepare to use
c   basic rule on each portion
c
   40 subrgn = 1
      pointr = wrkstr(1)
      rgncls = rulcls
      rgnvol = twondm
      tpontr = pointr + 2
      do 60 j=1,ndim
       tpontr = tpontr + 2
       center(j) = wrkstr(tpontr-1)
       width(j) = wrkstr(tpontr)
       dvcntr(j) = 1
       dvcntl(j) = 1
       oldcnt(j) = center(j)
       rgnvol = rgnvol*width(j)
   60 continue
      dvaxes = wrkstr(pointr+2)
      if (dvaxes.lt.0) go to 600
   80 dvaxis = dvaxes
      dvaxes = dvaxis/(ndim+1)
      dvaxis = dvaxis - (ndim+1)*dvaxes
      dvcntl(dvaxis) = 2*dvcntl(dvaxis)
      rgncls = rgncls*2
      if (dvaxes.gt.0) go to 80
      if (funcls+rgncls.gt.maxpts) go to 580
      if (rgncls/rulcls+sbrgns-1.gt.mxrgns) dvflag = 2
      funcls = funcls + rgncls
c      print *,funcls
      abserr = abserr - wrkstr(pointr)
      finval = finval - wrkstr(pointr+1)
c
c   begin basic rule
  100 do 120 j=1,ndim
       z(j) = center(j)
  120 continue
      sum1 = functn(ndim,z)
      sum2 = zero
      sum3 = zero
      do 140 j=1,ndim
       z(j) = center(j) - lamda2*width(j)
       f1 = functn(ndim,z)
       z(j) = center(j) + lamda2*width(j)
       f2 = functn(ndim,z)
       z(j) = center(j) - lamda4*width(j)
       f3 = functn(ndim,z)
       z(j) = center(j) + lamda4*width(j)
       f4 = functn(ndim,z)
       sum2 = sum2 + f1 + f2
       sum3 = sum3 + f3 + f4
       df1 = f1 + f2 - two*sum1
       df2 = f3 + f4 - two*sum1
       dif(j) = dif(j) + abs(df1-ratio*df2)
       z(j) = center(j)
  140 continue
      sum4 = zero
      do 200 j=2,ndim
       z(j-1) = center(j-1) - lamda4*width(j-1)
       do 160 k=j,ndim
	  z(k) = center(k) - lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k) + lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k)
  160  continue
       z(j-1) = center(j-1) + lamda4*width(j-1)
       do 180 k=j,ndim
	  z(k) = center(k) - lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k) + lamda4*width(k)
	  sum4 = sum4 + functn(ndim,z)
	  z(k) = center(k)
  180  continue
       z(j-1) = center(j-1)
  200 continue
      sum5 = zero
      do 220 j=1,ndim
       z(j) = center(j) - lamda5*width(j)
  220 continue
  240 do 260 j=2,ndim
       if (z(j-1).lt.center(j-1)+width(j-1)) go to 280
       z(j-1) = center(j-1) - lamda5*width(j-1)
       z(j) = z(j) + two*lamda5*width(j)
  260 continue
      if (z(ndim).gt.center(ndim)+width(ndim)) go to 300
  280 sum5 = sum5 + functn(ndim,z)
      z(1) = z(1) + two*lamda5*width(1)
      go to 240
  300 rgnval = rgnvol*(weit1*sum1+weit2*sum2+weit3*sum3+weit4*
     * sum4+weit5*sum5)
      rgncmp = rgnvol*(weitp1*sum1+weitp2*sum2+weitp3*sum3+weitp4*
     * sum4)
      rgnerr = abs(rgnval-rgncmp)
c
c   end basic rule
c   store results of basic rule application
c
      rgnvlt = rgnvlt + rgnval
      rgnert = rgnert + rgnerr
      finval = finval + rgnval
      abserr = abserr + rgnerr
      if (dvflag.eq.0) go to 340
      if (dvflag.eq.2) go to 500
      pointr = mxrgns + sbrgns*(2*ndim+3) + 1
      sbrgns = sbrgns + 1
      wrkstr(sbrgns) = pointr
      subrgn = sbrgns
      tpontr = pointr + 2
      do 320 j=1,ndim
       tpontr = tpontr + 2
       wrkstr(tpontr-1) = center(j)
       wrkstr(tpontr) = width(j)
  320 continue
  340 wrkstr(pointr) = rgnert
      wrkstr(pointr+1) = rgnvlt
c   determine axis along which fourth difference is largest
      difmax = zero
      do 380 j=1,ndim
       if (difmax.gt.dif(j)) go to 360
       difmax = dif(j)
       dvaxis = j
  360	    dif(j) = zero
  380 continue
      tpontr = pointr + 2*(dvaxis+1)
      wrkstr(tpontr) = width(dvaxis)*half
      wrkstr(tpontr-1) = center(dvaxis) - wrkstr(tpontr)
      if (dvflag.ne.2) go to 400
      dvaxes = wrkstr(pointr+2)
      if (dvaxes.gt.maxaxs) dvaxes = -1
      dvaxis = dvaxis + (ndim+1)*dvaxes
  400 wrkstr(pointr+2) = dvaxis
      if (dvflag.eq.1) go to 460
c   determine the position in the parially ordered list of
c   the subregion which replaces most recently divided subregion
  420 subtmp = 2*subrgn
      if (subtmp.gt.sbrgns) go to 480
      tpontr = wrkstr(subtmp)
      if (subtmp.eq.sbrgns) go to 440
      tpontp = wrkstr(subtmp+1)
      if (wrkstr(tpontr).ge.wrkstr(tpontp)) go to 440
      subtmp = subtmp + 1
      tpontr = tpontp
  440 if (rgnert.ge.wrkstr(tpontr)) go to 480
      wrkstr(subtmp) = pointr
      wrkstr(subrgn) = tpontr
      subrgn = subtmp
      go to 420
c   when working storage is not used up, determine the
c   position in the partially ordered list for the description
c   of other portion(s) of most recently divided subregion
  460 subtmp = subrgn/2
      if (subtmp.lt.1) go to 480
      tpontr = wrkstr(subtmp)
      if (rgnert.le.wrkstr(tpontr)) go to 480
      wrkstr(subtmp) = pointr
      wrkstr(subrgn) = tpontr
      subrgn = subtmp
      go to 460
  480 rgnvlt = zero
      rgnert = zero
      if (dvflag.eq.2) go to 540
      dvflag = 1 - dvflag
c   count to determine the next part of the recently divided
c   subregion for application of the basic rule
  500 center(1) = center(1) + two*width(1)
      dvcntr(1) = dvcntr(1) + 1
      do 520 j=2,ndim
       if (dvcntr(j-1).le.dvcntl(j-1)) go to 100
       dvcntr(j-1) = 1
       center(j-1) = oldcnt(j-1)
       dvcntr(j) = dvcntr(j) + 1
       center(j) = center(j) + two*width(j)
  520 continue
      if (dvcntr(ndim).le.dvcntl(ndim)) go to 100
      center(ndim) = oldcnt(ndim)
      if (dvflag.eq.2) go to 340
c
c   end ordering of basic rule results
c   make checks for possible termination of routine
c
  540 acc = abserr/abs(finval)
      if (acc.gt.eps .or. funcls.lt.minpts) go to 40
c
c   loop back to apply basic rule
c
c   termination point, set ifail and return
c
      ierror = 0
      go to 620
  560 ierror = 1
      go to 620
  580 ierror = 2
      go to 620
  600 ierror = 3
  620 minpts = funcls
      ifail = p01aae(ifail,ierror,srname)
      return
      end

      double precision function x02aae(x)
      implicit real*8(a-h,o-z)
c     nag copyright 1975
c     mark 4.5 release
c+self,if=ibm.
cc     for ibm/360/370/3090
c      data z/z3380000000000000/
c      x02aae = z
c     for sun
      data z/1.1d-16/
      x02aae = z
c     * eps *
c     returns the value eps where eps is the smallest
c     positive
c     number such that 1.0 + eps > 1.0
c     the x parameter is not used
c     for icl 1900
c     x02aae = 2.0**(-37.0)
c+self,if=pc.
c     for pdp11
c      x02aae=2.d0**(-23.d0)
c+self.

      return
      end
c
      integer  function x02bbe(x)
      implicit real*8(a-h,o-z)
c     nag copyright 1975
c     mark 4.5 release
*     real x
c     * maxint *
c     returns the largest integer representable on the computer
c     the x parameter is not used
c     for icl 1900
c      x02bbe = 8388607
c     for ibm,sun,vax,ibm pc/386/486
       x02bbe = 2147483647
c   for pdp11
c     x02bbe=32767
      return
      end

      integer function p01aae(ifail, error, srname)
c     mark 1 release.  nag copyright 1971
c     mark 3 revised
c     mark 4a revised, ier-45
c     mark 4.5 revised
c     mark 7 revised (dec 1978)
c     returns the value of error or terminates the program.
      integer error, ifail, nout
      character*8 srname
c     test if no error detected
      if (error.eq.0) go to 20
c     determine output unit for message
      call x04aae (0,nout)
c     test for soft failure
      if (mod(ifail,10).eq.1) go to 10
c     hard failure
      write (nout,99999) srname, error
c     stopping mechanism may also differ
      stop
c     soft fail
c     test if error messages suppressed
   10 if (mod(ifail/10,10).eq.0) go to 20
      write (nout,99999) srname, error
   20 p01aae = error
      return
99999 format (1h0, 38herror detected by nag library routine , a8,
     * 11h - ifail = , i5//)
      end
      subroutine x04aae(i,nerr)
c     mark 7 release. nag copyright 1978
c     mark 7c revised ier-190 (may 1979)
c     if i = 0, sets nerr to current error message unit number
c     (stored in nerr1).
c     if i = 1, changes current error message unit number to
c     value specified by nerr.
c
c     *** note ***
c     this routine assumes that the value of nerr1 is saved
c     between calls.  in some implementations it may be
c     necessary to store nerr1 in a labelled common
c     block /ax04aa/ to achieve this.
c
c     .. scalar arguments ..
      integer i, nerr
c     ..
c     .. local scalars ..
      integer nerr1
c     ..
      data nerr1 /5/
      if (i.eq.0) nerr = nerr1
      if (i.eq.1) nerr1 = nerr
      return
      end
