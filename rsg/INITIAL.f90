!##############################################
!##############################################
!######                                ########
!######       PROGRAM RCPF.f90         ########
!######       Author: Luxin Zhang      ########
!######       Time: 04/25/2011         ########
!######       IGGCAS,CN                ########
!######       Ver 1.1-2                ########
!######                                ########
!##############################################
!##############################################	

program RCPF
implicit none
! RCPF : R - RSG, CP - CPML, F - Free surface
! 2D biot rotated-staggered-grid finite-difference code in velocity-stress formulation
! with convolution PML (CPML) absorbing conditions for an isotropic medium.

! For stress at grid point, and velocity at half grid

! x8t2
! for Free surface, it should be processed with specially attention besides setting the velocity to zero, density to small number close to zero.

! in isotropic medium, averaging the parameter is not needed. But in inhomogeneous medium, should averaging the parameter.

! parameter(nx=256,nz=256,nt=1000,nx=nx+40,nz=nz+40)
! input(rhos,rhof,kd,ks,kf,mu,ew,kappa,eta,tortuosity)
! convert the inputparameter to parameter which could be used directly in this code with another program.
! directly used inputparamer(rhob,rhof,rhom,lamda,mu,alpha,m,b)
!   rhob = ew*rhof+(1-ew)*rhos
!   rhom = relaxation_time*eta/kappa = tortuosity*rhof/ew
!   alpha = 1-kd/ks
!   m = 1 / ( (alpha-ew)/ks + ew/kf )
!   lamda = kd - 2/3*mu + alpha**2 * m
!   b = eta/kappa


!    parameter(nx=101,nz=101,nt=100,np=20,nx=nx-2*np,nz=nz-2*np,c1=9.0/8.0,c2=1.0/24.0)
!    integer,parameter :: nx=200,nz=200,nt=2000,np=20
!    real,parameter :: c1=1225.0/1024.0,c2=c1/15.0,c3=c1/125.0,c4=c1/1715.0
!     integer,parameter :: nx=nx-2*np,nz=nz-2*np
    
	character*77 parfile
	character*77 fname
!	�ⲿ���ǲ����ļ��и����Ĳ���	
    integer*4 nx,nz,nt,xs,zs,np,iw,dst,iwav !np���ձ߽�ĺ�ȣ��ı��ȣ��Ƿ���Ըı���������
	real*4 dx,dz,dt,vmax,r,chi0,a0,f0,ew,power_display,cutvect !dx ���񳤶�x dz ���񳤶�z Vmaxȷ���߽����ղ��� r���뷴��ϵ��
	integer*4 nr,fxr,fzr,dxr,dzr !�첨���������һ���첨��λ�ã��첨����
	integer*4 ifs,fslc !�����ڴ��ˣ�û��ֵ�����Լ��о�fslcӦ����2
	integer*4 ifa,sl
    real*4 pi
	real*4 c1,c2,c3,c4    
    real*4 dr,maxpar
	real*4 d0,alphamax
    real*4 ts,refdist,Courant_number
    integer*4 ix,iz,it,i
	
	real*4,allocatable,dimension(:,:) :: rho
!	real*4,allocatable,dimension(:,:) :: rhob,rhof,rhom,rho2 !�������飬�������������ôʹ�û���Ҫ�о����������b����f ��m ��2��! allocatable �Ƕ���ɱ��С�����飬���潫��ȷ������Ĵ�С��rhob,rhof,rhom,rho2����������	
!   real*4,allocatable,dimension(:,:) :: lamda,mu,alpha,m,b !�������飨p��ģ�� ����ģ�� Biotϵ�� ���ģ�� ��ɢϵ����
	real*4,allocatable,dimension(:,:) :: c11,c12,c13,c22,c23,c33
!	real*4,allocatable,dimension(:,:) :: lamda,mu
!	define the density  at half grid
!	real*4,allocatable,dimension(:,:) :: rhobhx,rhofhx,rhomhx,rho2hx����Щ�ܶȶ��ǰ�����λ�õģ����ǣ���ô����������hx�����x��hz�����z
!    real*4,allocatable,dimension(:,:) :: rhobhz,rhofhz,rhomhz,rho2hz����Щ�ܶȶ��ǰ�����λ�õģ����ǣ���ô���������ģ�����������Ϊ��������Ӧ���������˰�����λ��
!    real*4,allocatable,dimension(:,:) :: rhobxz,rhofxz,rhomxz,rho2xz !��Щ�ܶȶ��ǰ�����λ�õģ����ǣ���ô���������ģ������б�����ˣ�û��
    real*4,allocatable,dimension(:,:) :: rhoxz
!	define 'b' at half grid
!	real*4,allocatable,dimension(:,:) :: bhx,bhz
!	real*4,allocatable,dimension(:,:) :: bxz !ΪʲôBҲ�ǰ���������
    
! txz��xz�����Ӧ����vx��x�����ٶȣ�qx��Һ��x�����ٶ�
    real*4,allocatable,dimension(:,:) :: txx,txz,tzz,vx,vz
!	real*4,allocatable,dimension(:,:) :: txx,txz,tzz,p,vx,vz,qx,qz !�������飨txx txz p vx vz qx qz������������λ���ǲ�һ����
!	real*4,allocatable,dimension(:,:) :: F11,F21,F12,F22,F13,F23,F14,F24
!	real*4,allocatable,dimension(:,:) :: tmp
	
! dx* is derivative for x, and dz* is derivative for z
!	real*4,allocatable,dimension(:,:) :: dxtxx,dztxx,dxtxz,dztxz,dxtzz,dztzz,dxp,dzp !��ּ��㣬�õ��ǿռ�8�ף���Ҫ�õ�����c1 c2 c3 d4 ������λ��
!   real*4,allocatable,dimension(:,:) :: dxvx,dzvz,dxvz,dzvx,dxqx,dzqz,dxqz,dzqx !��ּ��㣬�õ��ǿռ�8�ף���Ҫ�õ�����c1 c2 c3 d4 ����λ��
	real*4,allocatable,dimension(:,:) :: dxtxx,dztxx,dxtxz,dztxz,dxtzz,dztzz !��ּ��㣬�õ��ǿռ�8�ף���Ҫ�õ�����c1 c2 c3 d4 ������λ��
    real*4,allocatable,dimension(:,:) :: dxvx,dzvz,dxvz,dzvx
!   parameter for CPML
	real*4,allocatable,dimension(:,:) :: ax,bx,az,bz,chix,chiz,cdx,cdz !CMPL��������Ӧ��Ӧ���Ǧ���d��kai������Ԥ��Ķ�һ��������������
	real*4,allocatable,dimension(:,:) :: ax_half,bx_half,az_half,bz_half,chix_half,chiz_half,cdx_half,cdz_half !CMPL��������Ӧ��Ӧ���Ǧ���d��kai������Ԥ��Ķ�һ������������
    
!    real*4,allocatable,dimension(:,:) :: psidxtxx,psidztzz,psidxtxz,psidztxz,psidxp,psidzp !pusi��CPML��΢�����ӣ�ͨ��pusai���Եõ���˥���ı߽�������ͨ��pusai����ȷ������Ĳ�������������Ӧ��CPML
!   real*4,allocatable,dimension(:,:) :: psidxvx,psidzvx,psidxvz,psidzvz,psidxqx,psidzqx,psidxqz,psidzqz !!pusi��CPML��΢�����ӣ�ͨ��pusai���Եõ���˥���ı߽�������ͨ��pusai����ȷ������Ĳ��������������ٶ�CPML
    real*4,allocatable,dimension(:,:) :: psidxtxx,psidztzz,psidxtxz,psidztxz !pusi��CPML��΢�����ӣ�ͨ��pusai���Եõ���˥���ı߽�������ͨ��pusai����ȷ������Ĳ�������������Ӧ��CPML
    real*4,allocatable,dimension(:,:) :: psidxvx,psidzvx,psidxvz,psidzvz
! output velocity and stress
!	real*4,allocatable,dimension(:,:) :: seisvx,seisvz,seisqx,seisqz !���Vx Vz qx qz���������в�����Һ����λ��ʱ��qx qzӦ��Ϊ��
!   real*4,allocatable,dimension(:,:) :: seistxx,seistzz,seistxz,seisp !���Ӧ��
    real*4,allocatable,dimension(:,:) :: seisvx,seisvz !���Vx Vz qx qz���������в�����Һ����λ��ʱ��qx qzӦ��Ϊ��
    real*4,allocatable,dimension(:,:) :: seistxx,seistzz,seistxz
    real*4,allocatable,dimension(:,:) :: wav !����ɶ��
    real*4,allocatable,dimension(:) :: wavelet,wlet !�Ӳ�
	integer :: status=0 !???��ɶ��

! cpu_time
	real time_begin,time_end !����ʱ��
    
!##########################################################################################################	
!##########################################################################################################
	write(*,*) 'Input the parfile name:'!�൱��load�� load��ı�����Ϊparfile
	read(*,*) parfile !�൱��load�� load��ı�����Ϊparfile
	!parfile='zeng2.wok'
	print *,parfile !���parfile
	
	open(99,file=parfile,form='formatted',access='sequential',iostat=status) !��parfile��һformatted��ʽ�򿪣���һ�������ݵĸ�����ʽ��ѧϰfortran�̳�
	write(*,*) status
	read(99,'(64x)')
	
	read(99,'(70x)')
	read(99,*) nx,nz,nt,dx,dz,dt !����x����ڵ�����z����ڵ�����ʱ��ڵ�����x������z������ʱ�䲽����CMPL�߽�������ô�ú���Ҫ����ib��ɢϵ������������һ��0/1�жϣ�
	write(*,*) nx,nz,nt,dx,dz,dt
!	dt=dt/1000.0   !��С1000��
	
	read(99,'(70x)')
	read(99,*) np,vmax,r,chi0,a0 !����chi0��a0����cpml�Ĳ��� һ��chi0=1��a0=��Ƶ
	write(*,*) np,vmax,r,chi0,a0
	
	read(99,'(70x)')
	read(99,*) f0,iw,dst,iwav,power_display,cutvect !������Ƶ���Ӳ����͡���������ʱ������power_display��cutvect�ǻ�ͼ�Ĳ��� ��windows����û��ʹ�û�ͼ���ӳ�����Ҫʹ������Ҫ��linux�������У����޸ĳ���ʹ֮���ø��ӳ���
	write(*,*) f0,iw,dst,iwav,power_display,cutvect
	
	read(99,'(70x)')
	read(99,*) xs,zs !������Դ����
	write(*,*) xs,zs
	
	read(99,'(70x)')
	read(99,*) nr,fxr,fzr,dxr,dzr !����첨���������һ���첨��λ�ã��첨����
	write(*,*) nr,fxr,fzr,dxr,dzr
	
	read(99,'(70x)')
	read(99,*) ifs,fslc,ifa
	write(*,*) ifs,fslc,ifa 	

        read(99,'(70x)')
	read(99,*) sl
	write(*,*) sl 
	close(99)		

!	initial the parameter 
    dr = sqrt(dx**2+dz**2)
    dx = 2.0*dx !���ں�������ת�������������������Ϊԭ����2������dxtxx dxtxz dxp dxvx dxvz dxqx dxqzʱ��ǰ��Ӧ�ü�1/2
    dz = 2.0*dz !ͬ�� ������Ϊԭ����2������dztxx dztxz dzp dzvx dzvz dzqx dzqzʱ��ǰ��Ӧ�ü�1/2
	!   modify the coefficient of spatial opertor    
    c1 = (1225.0/1024.0)
    c2 = c1/15.0 !Ӧ��Ϊ��
    c3 = c1/125.0 !Ϊ��
    c4 = c1/1715.0 !Ϊ��
    ts = 1.2/f0 ! time delay
	
	pi = 4.0*atan(1.0) !3.1415926
    d0 = -3*vmax*log(r)/(2.0*np*dx/2.0) !˥��ϵ���������ٿ�
    alphamax = pi*a0 !CPML��Ҫ�õ��Ķ���
	
!	allocate the arrays
!	allocate(rhob(nz,nx),rhof(nz,nx),rhom(nz,nx),rho2(nz,nx),lamda(nz,nx),mu(nz,nx),alpha(nz,nx),m(nz,nx),b(nz,nx)) !���������f ��m rho2�Ǹ�ɶ�� p��ģ�� ����ģ�� Biotϵ�� ���ģ�� ��ɢϵ����������㣩
    allocate(rho(nz,nx),c11(nz,nx),c12(nz,nx),c13(nz,nz),c22(nz,nx),c23(nz,nx),c33(nz,nx))
!	allocate(rho(nz,nx),lamda(nz,nx),mu(nz,nx))
	allocate(txx(nz,nx),txz(nz,nx),tzz(nz,nx),vx(nz,nx),vz(nz,nx))
	allocate(rhoxz(nz,nx))
!	allocate(txx(nz,nx),txz(nz,nx),tzz(nz,nx),p(nz,nx),vx(nz,nx),vz(nz,nx),qx(nz,nx),qz(nz,nx))	!��������txx txz tzz p vx vz
!	allocate(F11(nz,nx),F21(nz,nx),F12(nz,nx),F22(nz,nx),F13(nz,nx),F23(nz,nx),F14(nz,nx),F24(nz,nx))
!	allocate(rhobhx(nz,nx),rhofhx(nz,nx),rhomhx(nz,nx),rho2hx(nz,nx),rhobhz(nz,nx),rhofhz(nz,nx),rhomhz(nz,nx),rho2hz(nz,nx))
!    allocate(rhobxz(nz,nx),rhofxz(nz,nx),rhomxz(nz,nx),rho2xz(nz,nx))!���������f ��m rho2�Ǹ�ɶ�� p��ģ�� ����ģ�� Biotϵ�� ���ģ��
!	allocate(bhx(nz,nx),bhz(nz,nx))
!	allocate(bxz(nz,nx)) !��ɢϵ����������㣩
!	allocate(dxtxx(nz,nx),dztxx(nz,nx),dxtxz(nz,nx),dztxz(nz,nx),dxtzz(nz,nx),dztzz(nz,nx),dxp(nz,nx),dzp(nz,nx))!�����ֵdxtxx dztxx dxtxz dztxz dxtzz dztzz����Ϊ�õİ�����㣩
!	allocate(dxvx(nz,nx),dzvz(nz,nx),dxvz(nz,nx),dzvx(nz,nx),dxqx(nz,nx),dzqz(nz,nx),dxqz(nz,nx),dzqx(nz,nx))!��ֽ��dxvx dxvz dxqx dxqz
	allocate(dxtxx(nz,nx),dztxx(nz,nx),dxtxz(nz,nx),dztxz(nz,nx),dxtzz(nz,nx),dztzz(nz,nx))
	allocate(dxvx(nz,nx),dzvz(nz,nx),dxvz(nz,nx),dzvx(nz,nx))
	allocate(ax(nz,nx),bx(nz,nx),az(nz,nx),bz(nz,nx),chix(nz,nx),chiz(nz,nx),cdx(nz,nx),cdz(nz,nx))!CPML�Ĳ���ɶ�ģ���������Vx Vz qx qz �����������治֪��ѽ��
	allocate(ax_half(nz,nx),bx_half(nz,nx),az_half(nz,nx),bz_half(nz,nx),chix_half(nz,nx),chiz_half(nz,nx),cdx_half(nz,nx),cdz_half(nz,nx))!CPML�Ĳ���ɶ�ģ���������
!	allocate(psidxtxx(nz,nx),psidztzz(nz,nx),psidxtxz(nz,nx),psidztxz(nz,nx),psidxp(nz,nx),psidzp(nz,nx))!CPML�����µ�΢��Ӧ��
	allocate(psidxtxx(nz,nx),psidztzz(nz,nx),psidxtxz(nz,nx),psidztxz(nz,nx))
!	allocate(psidxvx(nz,nx),psidzvx(nz,nx),psidxvz(nz,nx),psidzvz(nz,nx),psidxqx(nz,nx),psidzqx(nz,nx),psidxqz(nz,nx),psidzqz(nz,nx))!CPML�����µ�΢���ٶ�
	allocate(psidxvx(nz,nx),psidzvx(nz,nx),psidxvz(nz,nx),psidzvz(nz,nx))
	allocate(wav(nz,nx))!��֪����ɶ
!	allocate(seisvx(nt,nr),seisvz(nt,nr),seisqx(nt,nr),seisqz(nt,nr),seistxx(nt,nr),seistzz(nt,nr),seistxz(nt,nr),seisp(nt,nr))!��������
	allocate(seisvx(nt,nr),seisvz(nt,nr),seistxx(nt,nr),seistzz(nt,nr),seistxz(nt,nr))
	allocate(wavelet(nt))!�Ӳ�����
	allocate(wlet(nt))
	
!	allocate(tmp(nz-8-2*np,nx-8-2*np))
	
! read parameter file
!    open(1,file='rhob.dat',access='stream')!�����b
 !   open(2,file='rhof.dat',access='stream')!�����f
  !  open(3,file='rhom.dat',access='stream')!�����m
!    open(4,file='lamda.dat',access='stream')!����P��ģ��
!    open(5,file='mu.dat',access='stream')!�������ģ��
   ! open(6,file='alpha.dat',access='stream')!����Biotϵ��
   ! open(7,file='m.dat',access='stream')!�������ģ��
   ! open(8,file='b.dat',access='stream')!�����ɢϵ��
    ! open(9,file='ew.dat',access='stream')
	open(1,file='rho.dat',access='stream')
    open(2,file='c11.dat',access='stream')
    open(3,file='c12.dat',access='stream')
    open(4,file='c13.dat',access='stream')
    open(5,file='c22.dat',access='stream')
    open(6,file='c23.dat',access='stream')
	open(7,file='c33.dat',access='stream')
    open(8,file='wlet.dat',access='stream')
	!�����Ǹ���ϵ������Ĺ��̣��������������ʽ
!	read(1)rhob
!	read(2)rhof
!	read(3)rhom
!	read(4)lamda
!	read(5)mu
!	read(6)alpha
!	read(7)m
!	read(8)b
	! read(9)ew
	read(1)rho
    read(2)c11
	read(3)c12
	read(4)c13
	read(5)c22
	read(6)c23
	read(7)c33
	read(8)wlet
	
!	rhob = rhob*1.0e3
!	rhof = rhof*1.0e3
!	rhom = rhom*1.0e3
!	lamda = lamda*1.0e9
!	mu = mu*1.0e9
!	m = m*1.0e9
!	b = b*1.0e9
	rho=rho*1.0e3
        c11=c11*1.0e9
	c12=c12*1.0e9
	c13=c13*1.0e9
	c22=c22*1.0e9
	c23=c23*1.0e9
	c33=c33*1.0e9
    wlet=wlet
	
	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	close(6)
	close(7)
	close(8)
	! close(9)

!   check if the free surface is needed (ifs .ne. 0) and set the density (very small, almost 0) and modules (0) upon the free surface
    if(ifs .eq. 1) then !1˵�������ɱ߽�����!�����ڴ��ˣ�û��ֵ�����Լ��о�fslcӦ����2
!        rhob(1:fslc-1,:) = rhob(1:fslc,:)*1.0e-3
 !       rhof(1:fslc-1,:) = rhof(1:fslc,:)*1.0e-3
  !      rhom(1:fslc-1,:) = rhom(1:fslc,:)*1.0e-3
 !       lamda(1:fslc,:) = 0.0
!		lamda(nz:nz,:) = 0.0
!		mu(1:fslc,:) = 0.0
!		mu(nz:nz,:) = 0.0
  !      alpha(1:fslc,:) = 0.0
!        m(1:fslc,:) = 0.0
 !       b(1:fslc,:) = 0.0
	    rho(1:fslc-1,:) = rho(1:fslc-1,:)*1.0e-3
        rho(nz:nz,:) = rho(nz:nz,:)*1.0e-3
       c11(1:fslc-1,:) = 0.0
        c12(1:fslc-1,:) = 0.0
        c13(1:fslc-1,:) = 0.0
		c22(1:fslc-1,:) = 0.0
		c23(1:fslc-1,:) = 0.0
		c33(1:fslc-1,:) = 0.0
		c11(nz:nz,:) = 0.0
        c12(nz:nz,:) = 0.0
		c13(nz:nz,:) = 0.0
		c22(nz:nz,:) = 0.0
		c23(nz:nz,:) = 0.0
		c33(nz:nz,:) = 0.0
    endif
	if(ifa .eq. 1) then
	rho(:,1:fslc-1) = rho(:,1:fslc-1)*1.0e-3
    rho(:,nx:nx) = rho(:,nx:nx)*1.0e-3
!	lamda(:,1:fslc-1)=0.0
!	lamda(:,nx:nx) = 0.0
!	mu(:,1:fslc-1)=0.0
!	mu(:,nx:nx) = 0.0
       c11(:,1:fslc-1)=0.0
       c12(:,1:fslc-1)=0.0
	   c13(:,1:fslc-1)=0.0
       c22(:,1:fslc-1)=0.0
	   c23(:,1:fslc-1)=0.0
	   c33(:,1:fslc-1)=0.0
		c11(:,nx:nx) = 0.0
		c12(:,nx:nx) = 0.0
		c13(:,nx:nx) = 0.0
		c22(:,nx:nx) = 0.0
		c23(:,nx:nx) = 0.0
		c33(:,nx:nx) = 0.0
    endif
 

!   ������ļ���С����nz*nx�ģ��Ѱ�����PML�߽�Ŀ�ȡ�����ɼ����������ʱ������ֻ����м��PML���֡�

!   compute rho2 = rhob*rhom - rhof**2
!    rho2 = rhob*rhom - rhof**2 !���Ǹ�ɶ�ܶȣ�û����
!	if (ib .eq. 0) b=0.0
	
! interpolate rho and b at half grid
    do iz=1,nz-1 !�ܶ�����Χ�ĸ����ܶȵĲ�ֵ��ƽ��ֵ��
        do ix=1,nx-1
!            rhobxz(iz,ix) = 0.25*(rhob(iz  ,ix  )+rhob(iz+1,ix  )+rhob(iz+1,ix+1)+rhob(iz  ,ix+1))
 !           rhomxz(iz,ix) = 0.25*(rhom(iz  ,ix  )+rhom(iz+1,ix  )+rhom(iz+1,ix+1)+rhom(iz  ,ix+1))
  !          rhofxz(iz,ix) = 0.25*(rhof(iz  ,ix  )+rhof(iz+1,ix  )+rhof(iz+1,ix+1)+rhof(iz  ,ix+1))
   !         rho2xz(iz,ix) = 0.25*(rho2(iz  ,ix  )+rho2(iz+1,ix  )+rho2(iz+1,ix+1)+rho2(iz  ,ix+1))
    !        bxz(iz,ix)    = 0.25*(b(iz  ,ix  )+b(iz+1,ix  )+b(iz+1,ix+1)+b(iz  ,ix+1))
        rhoxz(iz,ix) = 0.25*(rho(iz  ,ix  )+rho(iz+1,ix  )+rho(iz+1,ix+1)+rho(iz  ,ix+1))
		enddo
    enddo
	rhoxz(nz  ,1:nx-1) = 0.5*(rho(nz  ,1:nx-1)+rho(nz  ,2:nx  ))
    rhoxz(1:nz-1,nx) = 0.5*(rho(1:nz-1  ,nx  )+rho(2:nz  ,nx  ))
    rhoxz(nz,nx) = rho(nz,nx)
!    rhobxz(nz  ,1:nx-1) = 0.5*(rhob(nz  ,1:nx-1)+rhob(nz  ,2:nx  ))
 !   rhobxz(1:nz-1,nx) = 0.5*(rhob(1:nz-1  ,nx  )+rhob(2:nz  ,nx  ))
  !  rhobxz(nz,nx) = rhob(nz,nx)
    
!    rhomxz(nz  ,1:nx-1) = 0.5*(rhom(nz  ,1:nx-1)+rhom(nz  ,2:nx  ))
 !   rhomxz(1:nz-1,nx) = 0.5*(rhom(1:nz-1  ,nx  )+rhom(2:nz  ,nx  ))
  !  rhomxz(nz,nx) = rhom(nz,nx)
    
!    rhofxz(nz  ,1:nx-1) = 0.5*(rhof(nz  ,1:nx-1)+rhof(nz  ,2:nx  ))
 !   rhofxz(1:nz-1,nx) = 0.5*(rhof(1:nz-1  ,nx  )+rhof(2:nz  ,nx  ))
  !  rhofxz(nz,nx) = rhof(nz,nx)
    
!    rho2xz(nz  ,1:nx-1) = 0.5*(rho2(nz  ,1:nx-1)+rho2(nz  ,2:nx  ))
 !   rho2xz(1:nz-1,nx) = 0.5*(rho2(1:nz-1  ,nx  )+rho2(2:nz  ,nx  ))
  !  rho2xz(nz,nx) = rho2(nz,nx)
    
!    bxz(nz  ,1:nx-1) = 0.5*(b(nz  ,1:nx-1)+b(nz  ,2:nx  ))
 !   bxz(1:nz-1,nx) = 0.5*(b(1:nz-1  ,nx  )+b(2:nz  ,nx  ))
  !  bxz(nz,nx) = b(nz,nx)  	
	
!	if(ifs .eq. 1) then
!        rho2xz(1:fslc-1,:) = rho2xz(fslc+2,1)
 !   endif
    
	if(iw .eq. 1) then
		call dgauss(wavelet,nt,dt,f0,ts)
	elseif(iw .eq. 0) then
		call ricker(wavelet,nt,dt,f0,ts)
	elseif(iw.eq. 2) then
                wavelet=wlet
	else
		wavelet = 0.0
		wavelet(1) = 1.0
    endif

    open(55,file='wavelet.dat',access='stream')
    write(55)wavelet
    close(55)
    !write(*,*) wavelet(1:40)
    
!	/* let source contribute within limited distance */	
	if(iwav .eq. 1) then
		do ix = max(1,xs-3),min(nx,xs+3)
			do iz = max(1,zs-3),min(nz,zs+3)
				wav(iz,ix) = exp(-0.2*((real(iz)-zs)**2 + (real(ix)-xs)**2))
			end do
		end do
	else
		wav = 0.0
		wav(zs-2:zs+2,xs-2:xs+2) = 1.0
	endif
	open(56,file='wav.dat',access='stream')
    write(56)wav
    close(56)
!	if(iwav .eq. 1) then
!		do ix = max(1,xs-3),min(nx,xs+3)1,nx
!			do iz = 1,nz
!				wav(iz,ix) = exp(-0.2*((real(iz)-zs)**2 + (real(ix)-xs)**2))
!			end do
!		end do
!	else
!		wav = 0.0
!		wav(zs-2:zs+2,xs-2:xs+2) = 1.0
!	endif
! time begin
	call cpu_time(time_begin)

!---------------PML------------------------
! need vmax, so compute vmax_max first. Here assume vmax is known.
    
    ax = 0.0
    az = 0.0
    bx = 1.0
    bz = 1.0
    chix = 1.0
    chiz = 1.0
    cdx = 0.0
    cdz = 0.0

    ax_half = 0.0
    az_half = 0.0
    bx_half = 1.0
    bz_half = 1.0
    chix_half = 1.0
    chiz_half = 1.0
    cdx_half = 0.0
    cdz_half = 0.0

	do ix = 1,nx
	    if(ifs .eq. 0) then
	         do iz = 1,np+5
	              refdist = (np+5.0-(iz))/np
	              if(refdist > 1.0) refdist = 1.0
	              if(refdist < 0.0) refdist = 0.0
	              if(refdist >= 0.0) then
	                  chiz(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	                  cdz(iz,ix) = d0*(refdist)**2
	                  bz(iz,ix) = exp(-(cdz(iz,ix)/chiz(iz,ix)+ alphamax*(1.0-refdist)) * dt)
	                  az(iz,ix) = cdz(iz,ix) / (chiz(iz,ix)*(cdz(iz,ix) + chiz(iz,ix)*alphamax*(1.0-refdist))) * (bz(iz,ix) - 1.0)
	              endif
	              if(ix == 1) print *,'refdist at iz=',iz,'=',refdist

	              refdist = ((np+4.5-(iz))/np)
	              if(refdist > 1.0) refdist = 1.0
	              if(refdist < 0.0) refdist = 0.0
	              if(refdist >= 0.0) then
	                   chiz_half(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	                   cdz_half(iz,ix) = d0*(refdist)**2
	                   bz_half(iz,ix) = exp(-(cdz_half(iz,ix)/chiz_half(iz,ix)+ alphamax*(1.0-refdist)) * dt)
	                   az_half(iz,ix) = cdz_half(iz,ix) / (chiz_half(iz,ix)*(cdz_half(iz,ix) + chiz_half(iz,ix)*alphamax*(1.0-refdist))) * (bz_half(iz,ix) - 1.0)
	               endif
	               if(ix == 1) print *,'refdist half at iz=',iz,'=',refdist
	         enddo
        endif

     do iz = (nz-4-np),nz
          refdist = ((iz-0.5)-(nz-4.0-np))/np
         if(refdist > 1.0) refdist = 1.0
          if(refdist < 0.0) refdist = 0.0
          if(refdist >= 0.0) then
              chiz(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
              cdz(iz,ix) = d0*(refdist)**2
              bz(iz,ix) = exp(-(cdz(iz,ix)/chiz(iz,ix) + alphamax*(1.0-refdist)) * dt)
              az(iz,ix) = cdz(iz,ix) / (chiz(iz,ix)*(cdz(iz,ix) + chiz(iz,ix)*alphamax*(1.0-refdist))) * (bz(iz,ix) - 1.0)
         endif
         if(ix == 1) print *,'refdist at iz=',iz,'=',refdist

          refdist = (((iz-0.0)-(nz-4.0-np))/np)
         if(refdist > 1.0) refdist = 1.0
          if(refdist < 0.0) refdist = 0.0
          if(refdist >= 0.0) then
               chiz_half(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
               cdz_half(iz,ix) = d0*(refdist)**2
               bz_half(iz,ix) = exp(-(cdz_half(iz,ix)/chiz_half(iz,ix) + alphamax*(1.0-refdist)) * dt)
               az_half(iz,ix) = cdz_half(iz,ix) / (chiz_half(iz,ix)*(cdz_half(iz,ix) + chiz_half(iz,ix)*alphamax*(1.0-refdist))) * (bz_half(iz,ix) - 1.0)
          endif
          if(ix == 1) print *,'refdist half at iz=',iz,'=',refdist
     enddo
    enddo

	do iz = 1,nz
    if(ifa .eq. 0) then
	 do ix = 1,np+5
	  refdist = (np+5.0-(ix))/np
     if(refdist > 1.0) refdist = 1.0
	  if(refdist < 0.0) refdist = 0.0
	  if(refdist >= 0.0) then
	  chix(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	  cdx(iz,ix) = d0*(refdist)**2
	  bx(iz,ix) = exp(-(cdx(iz,ix)/chix(iz,ix)+ alphamax*(1.0-refdist)) * dt)
	  ax(iz,ix) = cdx(iz,ix) / (chix(iz,ix)*(cdx(iz,ix) + chix(iz,ix)*alphamax*(1.0-refdist))) * (bx(iz,ix) - 1.0)
	 endif
	 if(iz == 1) print *,'refdist at ix=',ix,'=',refdist

	  refdist = (np+4.5-(ix))/np
     if(refdist > 1.0) refdist = 1.0
	  if(refdist < 0.0) refdist = 0.0
	  if(refdist >= 0.0) then
	   chix_half(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	   cdx_half(iz,ix) = d0*(refdist)**2
	   bx_half(iz,ix) = exp(-(cdx_half(iz,ix)/chix_half(iz,ix)+ alphamax*(1.0-refdist)) * dt)
	   ax_half(iz,ix) = cdx_half(iz,ix) / (chix_half(iz,ix)*(cdx_half(iz,ix) + chix_half(iz,ix)*alphamax*(1.0-refdist))) * (bx_half(iz,ix) - 1.0)
	  endif
	  if(iz == 1) print *,'refdist half at ix=',ix,'=',refdist
	 enddo

	 do ix = (nx-4-np),nx
	  refdist = ((ix-0.5)-(nx-4-np))/np
     if(refdist > 1.0) refdist = 1.0
	  if(refdist < 0.0) refdist = 0.0
	  if(refdist >= 0.0) then
	  chix(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	  cdx(iz,ix) = d0*(refdist)**2
	  bx(iz,ix) = exp(-(cdx(iz,ix)/chix(iz,ix) + alphamax*(1.0-refdist)) * dt)
	  ax(iz,ix) = cdx(iz,ix) / (chix(iz,ix)*(cdx(iz,ix) + chix(iz,ix)*alphamax*(1.0-refdist))) * (bx(iz,ix) - 1.0)
	 endif
	 if(iz == 1) print *,'refdist at ix=',ix,'=',refdist

	  refdist = ((ix+0.0)-(nx-4-np))/np
     if(refdist > 1.0) refdist = 1.0
	  if(refdist < 0.0) refdist = 0.0
	  if(refdist >= 0.0) then
	   chix_half(iz,ix) = 1.0+(chi0-1.0)*(refdist)**2
	   cdx_half(iz,ix) = d0*(refdist)**2
	   bx_half(iz,ix) = exp(-(cdx_half(iz,ix)/chix_half(iz,ix) + alphamax*(1.0-refdist)) * dt)
	   ax_half(iz,ix) = cdx_half(iz,ix) / (chix_half(iz,ix)*(cdx_half(iz,ix) + chix_half(iz,ix)*alphamax*(1.0-refdist))) * (bx_half(iz,ix) - 1.0)
	  endif
	  if(iz == 1) print *,'refdist half at ix=',ix,'=',refdist

	 enddo
	 end if
	enddo

! print position of the source
	print *,'Position of the source:'
	print *
	print *,'x = ',xs
	print *,'z = ',zs
	print *

! define location of receivers
	print *,'There are ',nr,' receivers'
	print *

	print *,'the 1st receiver: (x_target,z_target) = ',fxr,fzr
	print *,'the last receiver: (x_target,z_target) = ',fxr+real(nr-1)*dxr,fzr+real(nr-1)*dzr
	print *,'distance between receivers: (dx,dz) = ',dxr,dzr

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
	Courant_number = vmax * dt / dr
	print *,'Courant number is ',Courant_number
	print *
	if(Courant_number > 1.0) then
		pause 'time step is too large, simulation will be unstable'
	endif

!-----------------------------------------------------------------------
!                   compute code
!-----------------------------------------------------------------------
! initial the initial values
	psidxtxx = 0.0
    psidztzz = 0.0
    psidxtxz = 0.0
    psidztxz = 0.0
!    psidxp = 0.0
!    psidzp = 0.0
	
    psidxvx = 0.0
    psidzvx = 0.0
    psidxvz = 0.0
    psidzvz = 0.0
	
!    psidxqx = 0.0
 !   psidzqx = 0.0
  !  psidxqz = 0.0
   ! psidzqz = 0.0
	
    vx = 0.0
    vz = 0.0
!    qx = 0.0
 !   qz = 0.0
	
	txx = 0.0
	txz = 0.0
	tzz = 0.0
!	p = 0.0

	write(*,*) 'Begin time step computing ...'
! time begin
!	call cpu_time(time_begin)	

	do it = 1,nt
        
        if(mod(it,100) .eq. 0) print *,'No. ',it,' of ',nt

! apply PML boundary
! compute the x ans z component of velocity
!!   compute the analytical solution of velocity in the stiff equations O(1)
!       vx = vx - rhofxz/rhobxz*qx * (exp(-rhobxz/rho2xz*bxz*dt) - 1.0)
!       vz = vz - rhofxz/rhobxz*qz * (exp(-rhobxz/rho2xz*bxz*dt) - 1.0)
!       qx = qx*exp(-rhobxz/rho2xz*bxz*dt)
!       qz = qz*exp(-rhobxz/rho2xz*bxz*dt)
       
!   compute the analytical solution of velocity in the stiff equations O(2)
!       vx = vx - rhofxz/rhobxz*qx * (exp(-rhobxz/rho2xz*bxz*dt/2.0) - 1.0)
 !      vz = vz - rhofxz/rhobxz*qz * (exp(-rhobxz/rho2xz*bxz*dt/2.0) - 1.0)
!       qx = qx*exp(-rhobxz/rho2xz*bxz*dt/2.0)
 !      qz = qz*exp(-rhobxz/rho2xz*bxz*dt/2.0)
       

!   compute the velocity of nonstiff equations
!   compute the spatial difference of stress and pressure
        dxtxx = 0.0
        dztxx = 0.0
        dxtxz = 0.0
        dztxz = 0.0
        dxtzz = 0.0
        dztzz = 0.0
!        dxp = 0.0
 !       dzp = 0.0
        
        do ix=4,nx-4
            do iz=4,nz-4
                dxtxx(iz,ix) = (c1*(txx(iz  ,ix+1)-txx(iz+1,ix  )) - c2*(txx(iz-1,ix+2)-txx(iz+2,ix-1)) + c3*(txx(iz-2,ix+3)-txx(iz+3,ix-2)) -c4*(txx(iz-3,ix+4)-txx(iz+4,ix-3)))
                dztxx(iz,ix) = (c1*(txx(iz+1,ix+1)-txx(iz  ,ix  )) - c2*(txx(iz+2,ix+2)-txx(iz-1,ix-1)) + c3*(txx(iz+3,ix+3)-txx(iz-2,ix-2)) -c4*(txx(iz+4,ix+4)-txx(iz-3,ix-3)))
                
                dxtzz(iz,ix) = (c1*(tzz(iz  ,ix+1)-tzz(iz+1,ix  )) - c2*(tzz(iz-1,ix+2)-tzz(iz+2,ix-1)) + c3*(tzz(iz-2,ix+3)-tzz(iz+3,ix-2)) -c4*(tzz(iz-3,ix+4)-tzz(iz+4,ix-3)))
                dztzz(iz,ix) = (c1*(tzz(iz+1,ix+1)-tzz(iz  ,ix  )) - c2*(tzz(iz+2,ix+2)-tzz(iz-1,ix-1)) + c3*(tzz(iz+3,ix+3)-tzz(iz-2,ix-2)) -c4*(tzz(iz+4,ix+4)-tzz(iz-3,ix-3)))
                
                dxtxz(iz,ix) = (c1*(txz(iz  ,ix+1)-txz(iz+1,ix  )) - c2*(txz(iz-1,ix+2)-txz(iz+2,ix-1)) + c3*(txz(iz-2,ix+3)-txz(iz+3,ix-2)) -c4*(txz(iz-3,ix+4)-txz(iz+4,ix-3)))
                dztxz(iz,ix) = (c1*(txz(iz+1,ix+1)-txz(iz  ,ix  )) - c2*(txz(iz+2,ix+2)-txz(iz-1,ix-1)) + c3*(txz(iz+3,ix+3)-txz(iz-2,ix-2)) -c4*(txz(iz+4,ix+4)-txz(iz-3,ix-3)))
                
!                dxp(iz,ix)   = (c1*(  p(iz  ,ix+1)-  p(iz+1,ix  )) - c2*(  p(iz-1,ix+2)-  p(iz+2,ix-1)) + c3*(  p(iz-2,ix+3)-  p(iz+3,ix-2)) -c4*(  p(iz-3,ix+4)-  p(iz+4,ix-3)))
!                dzp(iz,ix)   = (c1*(  p(iz+1,ix+1)-  p(iz  ,ix  )) - c2*(  p(iz+2,ix+2)-  p(iz-1,ix-1)) + c3*(  p(iz+3,ix+3)-  p(iz-2,ix-2)) -c4*(  p(iz+4,ix+4)-  p(iz-3,ix-3)))
            enddo
        enddo

		 dxtxx = (dztxx+dxtxx) / dx
		 dztzz = (dztzz-dxtzz) / dz
		 dxtxz = (dztxz+dxtxz) / dx
		 dztxz = (2.0*dztxz-dx*dxtxz) / dz
!		 dxp   = (dzp+dxp) / dx
!		 dzp   = (2.0*dzp-dx*dxp) / dz
!         dxtxx = (dxtxx) / dx
!		 dztzz = (dztzz) / dz
!		 dxtxz = (dxtxz) / dx
!		 dztxz = (dztxz) / dz
        
! update the psi_function of stress in the CPML
		psidxtxx = bx_half*psidxtxx + ax_half*dxtxx
		psidztzz = bz_half*psidztzz + az_half*dztzz
		psidxtxz = bx_half*psidxtxz + ax_half*dxtxz
		psidztxz = bz_half*psidztxz + az_half*dztxz
!		psidxp = bx_half*psidxp + ax_half*dxp
!		psidzp = bz_half*psidzp + az_half*dzp
		
! take dxtxx as dxtxx+psidxtxx
       dxtxx = dxtxx/chix_half + psidxtxx
        dztzz = dztzz/chiz_half + psidztzz
       dxtxz = dxtxz/chix_half + psidxtxz
        dztxz = dztxz/chiz_half + psidztxz
  !      dxp = dxp/chix_half + psidxp
   !     dzp = dzp/chiz_half + psidzp 
		   
! update the x and z component of solid and fluid velocity 
!	for velocity, the density should be averaged
!        qx = qx - dt/rho2xz*(rhofxz*(dxtxx+dztxz) + rhobxz*dxp)
 !   	qz = qz - dt/rho2xz*(rhofxz*(dxtxz+dztzz) + rhobxz*dzp)
 !       vx = vx + dt/rho2xz*(rhomxz*(dxtxx+dztxz) + rhofxz*dxp)
  !  	vz = vz + dt/rho2xz*(rhomxz*(dxtxz+dztzz) + rhofxz*dzp)
    	vx = vx + dt/rhoxz*((dxtxx+dztxz))
    	vz = vz + dt/rhoxz*((dxtxz+dztzz))
!   compute the analytical solution of velocity in the stiff equations O(2)
!       vx = vx - rhofxz/rhobxz*qx * (exp(-rhobxz/rho2xz*bxz*dt/2.0) - 1.0)
 !      vz = vz - rhofxz/rhobxz*qz * (exp(-rhobxz/rho2xz*bxz*dt/2.0) - 1.0)
  !     qx = qx*exp(-rhobxz/rho2xz*bxz*dt/2.0)
   !    qz = qz*exp(-rhobxz/rho2xz*bxz*dt/2.0)    	


!   compute the stress and pressure of nonstiff equations
!   compute the spatial difference of velocity
        dxvx = 0.0
        dzvx = 0.0
        dzvz = 0.0
        dxvz = 0.0
!        dxqx = 0.0
 !       dzqx = 0.0
  !      dzqz = 0.0
   !     dxqz = 0.0
		do ix = 5,nx-3
			do iz = 5,nz-3
			    dzvx(iz,ix) = (c1*(vx(iz  ,ix  )-vx(iz-1,ix-1)) - c2*(vx(iz+1,ix+1)-vx(iz-2,ix-2)) +c3*(vx(iz+2,ix+2)-vx(iz-3,ix-3)) -c4*(vx(iz+3,ix+3)-vx(iz-4,ix-4)))
				dxvx(iz,ix) = (c1*(vx(iz-1,ix  )-vx(iz  ,ix-1)) - c2*(vx(iz-2,ix+1)-vx(iz+1,ix-2)) +c3*(vx(iz-3,ix+2)-vx(iz+2,ix-3)) -c4*(vx(iz-4,ix+3)-vx(iz+3,ix-4)))
				
				dzvz(iz,ix) = (c1*(vz(iz  ,ix  )-vz(iz-1,ix-1)) - c2*(vz(iz+1,ix+1)-vz(iz-2,ix-2)) +c3*(vz(iz+2,ix+2)-vz(iz-3,ix-3)) -c4*(vz(iz+3,ix+3)-vz(iz-4,ix-4)))
				dxvz(iz,ix) = (c1*(vz(iz-1,ix  )-vz(iz  ,ix-1)) - c2*(vz(iz-2,ix+1)-vz(iz+1,ix-2)) +c3*(vz(iz-3,ix+2)-vz(iz+2,ix-3)) -c4*(vz(iz-4,ix+3)-vz(iz+3,ix-4)))
				
!				dzqx(iz,ix) = (c1*(qx(iz  ,ix  )-qx(iz-1,ix-1)) - c2*(qx(iz+1,ix+1)-qx(iz-2,ix-2)) +c3*(qx(iz+2,ix+2)-qx(iz-3,ix-3)) -c4*(qx(iz+3,ix+3)-qx(iz-4,ix-4)))
!				dxqx(iz,ix) = (c1*(qx(iz-1,ix  )-qx(iz  ,ix-1)) - c2*(qx(iz-2,ix+1)-qx(iz+1,ix-2)) +c3*(qx(iz-3,ix+2)-qx(iz+2,ix-3)) -c4*(qx(iz-4,ix+3)-qx(iz+3,ix-4)))
				
!				dzqz(iz,ix) = (c1*(qz(iz  ,ix  )-qz(iz-1,ix-1)) - c2*(qz(iz+1,ix+1)-qz(iz-2,ix-2)) +c3*(qz(iz+2,ix+2)-qz(iz-3,ix-3)) -c4*(qz(iz+3,ix+3)-qz(iz-4,ix-4)))
!				dxqz(iz,ix) = (c1*(qz(iz-1,ix  )-qz(iz  ,ix-1)) - c2*(qz(iz-2,ix+1)-qz(iz+1,ix-2)) +c3*(qz(iz-3,ix+2)-qz(iz+2,ix-3)) -c4*(qz(iz-4,ix+3)-qz(iz+3,ix-4)))
			enddo
		enddo	

   dxvx = (dzvx+dxvx) / dx
   dzvx = (2.0*dzvx-dx*dxvx) / dz
   dxvz = (dzvz+dxvz) / dx
   dzvz = (2.0*dzvz-dx*dxvz) / dz
!   dxqx = (dzqx+dxqx) / dx
 !  dzqx = (2.0*dzqx-dx*dxqx) / dz
  ! dxqz = (dzqz+dxqz) / dx
   !dzqz = (2.0*dzqz-dx*dxqz) / dz
 !   dxvx = (dxvx) / dx
  !  dzvx = (dzvx) / dz
   ! dxvz = (dxvz) / dx
   ! dzvz = (dzvz) / dz
		
! update the psi_function of the CPML
!	the CPML parameter here is at half grid
		psidxvx = bx*psidxvx + ax*dxvx
		psidzvx = bz*psidzvx + az*dzvx
		psidxvz = bx*psidxvz + ax*dxvz
		psidzvz = bz*psidzvz + az*dzvz
		
!		psidxqx = bx*psidxqx + ax*dxqx
!		psidzqx = bz*psidzqx + az*dzqx
!		psidxqz = bx*psidxqz + ax*dxqz
!		psidzqz = bz*psidzqz + az*dzqz
		
! take dxvx as dxvx+psidxvx	    
        dxvx = dxvx/chix + psidxvx
        dzvx = dzvx/chiz + psidzvx
        dxvz = dxvz/chix + psidxvz
        dzvz = dzvz/chiz + psidzvz
        
!        dxqx = dxqx/chix + psidxqx
 !       dzqx = dzqx/chiz + psidzqx
  !      dxqz = dxqz/chix + psidxqz
   !     dzqz = dzqz/chiz + psidzqz	        
			
		
! update the stress and pressure 
!        txx = txx + ((lamda+2.0*mu)*(dxvx) + lamda*(dzvz) + alpha*m*(dxqx+dzqz))*dt
 !       tzz = tzz + ((lamda+2.0*mu)*(dzvz) + lamda*(dxvx) + alpha*m*(dxqx+dzqz))*dt
  !      txz = txz + (mu*(dxvz+dzvx))*dt
   !     p = p + (-alpha*m*(dxvx+dzvz) - m*(dxqx+dzqz))*dt
       txx = txx + ((c11)*(dxvx) + c12*(dzvz)+c13*(dxvz+dzvx))*dt
        tzz = tzz + ((c22)*(dzvz) + c12*(dxvx)+c23*(dxvz+dzvx))*dt
        txz = txz + ((c13)*(dxvx) + c23*(dzvz)+c33*(dxvz+dzvx))*dt
              seisvx(it,1)=vx(404,410)
	!	seistxz(it,1)=txz(404,410)
		seisvz(it,1)=vz(404,450)

! store seismograms
		do i = 1,nr
			seisvx(it,i) = vx(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
			seisvz(it,i) = vz(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
!			seisqx(it,i) = qx(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
!			seisqz(it,i) = qz(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
			
			seistxx(it,i) = txx(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
			seistxz(it,i) = txz(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
			seistzz(it,i) = tzz(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
!			seisp(it,i) = p(fzr+(i-1)*dzr,fxr+(i-1)*dxr)
		end do
		if(sl.eq.1) then
		txx = txx + wavelet(it)*wav
!		tzz = tzz + wavelet(it)*wav
!		txz = txz + wavelet(it)*(1.0-ew)*wav/2.0
   
        endif
	    if(sl.eq.2) then
!		txx = txx + wavelet(it)*wav
		tzz = tzz+wavelet(it)*wav
!		txz = txz + wavelet(it)*(1.0-ew)*wav/2.0
  
        endif
        if(sl.eq.3) then
!		txx = txx + wavelet(it)*wav
!		tzz = tzz + wavelet(it)*wav
		txz = txz + wavelet(it)*(1.0-ew)*wav/2.0
!        p = p + wavelet(it)*wav
        endif
		
!		txx = txx + wavelet(it)*wav
!		tzz = tzz + wavelet(it)*wav
!		txz = txz + wavelet(it)*(1.0-ew)*wav/2.0
!        p = p + wavelet(it)*wav
		
! output the snapshot
		if(mod(it,dst) .eq. 0) then
			print *,'Time step # ',it
			print *,'Time: ',sngl((it-1)*dt),' seconds'
			print *,'Max norm velocity vector V (m/s) = ',max(maxval(sqrt(vx**2 + vz**2)),maxval(sqrt(vz**2)))
			print *	
			
			write(fname,"('snap_',i6.6,'_Vx.dat')") it
			maxpar = 1.0
			!maxpar = maxval(abs(vx))
			!maxpar = abs(vz(zs,xs))
!			tmp = vx(np+5:nz-4-np,np+5:nx-4-np)
			open(11,file=fname,access='stream')
			write(11) vx
			close(11)
!			call creat_2D_image(vx,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,0,power_display,cutvect)
            write(fname,"('snap_',i6.6,'_Vz.dat')") it
			!maxpar = maxval(abs(vz))
			!maxpar = 1.0
!			tmp = vz(np+5:nz-4-np,np+5:nx-4-np)
			open(12,file=fname,access='stream')
			write(12) vz
			close(12)
!			call creat_2D_image(vz,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,1,power_display,cutvect)
!            write(fname,"('snap_',i6.6,'_qx.dat')") it
			!maxpar = maxval(qx)
!			tmp = qx(np+5:nz-4-np,np+5:nx-4-np)
!			open(13,file=fname,access='stream')
!			write(13) qx
!			close(13)
!			call creat_2D_image(qx,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,2,power_display,cutvect)
!            write(fname,"('snap_',i6.6,'_qz.dat')") it
			!maxpar = maxval(qz)
!			tmp = qz(np+5:nz-4-np,np+5:nx-4-np)
!			open(14,file=fname,access='stream')
!			write(14) qz
!			close(14)
!			call creat_2D_image(qz,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,3,power_display,cutvect)

            write(fname,"('snap_',i6.6,'_txx.dat')") it
			!maxpar = maxval(txx)
!			tmp = txx(np+5:nz-4-np,np+5:nx-4-np)
			open(15,file=fname,access='stream')
			write(15) txx
			close(15)
!			call creat_2D_image(txx,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,4,power_display,cutvect)

            write(fname,"('snap_',i6.6,'_tzz.dat')") it
			!maxpar = maxval(tzz)
!			tmp = tzz(np+5:nz-4-np,np+5:nx-4-np)
			open(16,file=fname,access='stream')
			write(16) tzz
			close(16)
!			call creat_2D_image(tzz,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,5,power_display,cutvect)

            
            write(fname,"('snap_',i6.6,'_txz.dat')") it
!			!maxpar = maxval(txz)
!			tmp = txz(np+5:nz-4-np,np+5:nx-4-np)
			open(17,file=fname,access='stream')
			write(17) txz
			close(17)
!!			call creat_2D_image(txz,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np+4,6,power_display,cutvect)

			!maxpar = maxval(p)
!   		    write(fname,"('snap_',i6.6,'_p.dat')") it
 !  		    tmp = p(np+5:nz-4-np,np+5:nx-4-np)
!   		    open(18,file=fname,access='stream')
!   		    write(18) p
!   		    close(18)
!			call creat_2D_image(p,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np,7,power_display,cutvect)
		end if
		
	end do
	
	call cpu_time(time_end)
	write(*,*) 'Time of operation was ',time_end-time_begin,' seconds'


!   write the seismgram
    open(109,file='seisvz.dat',access='stream')
    write(109) seisvz
    close(109)

!	open(110,file='seisp.dat',access='stream')
 !   write(110) seisp
  !  close(110)

	open(111,file='seisvx.dat',access='stream')
    write(111) seisvx
    close(111)

!	open(112,file='seisqz.dat',access='stream')
 !   write(112) seisqz
  !  close(112)

!	open(113,file='seisqx.dat',access='stream')
 !   write(113) seisqx
  !  close(113)

	open(114,file='seistxx.dat',access='stream')
    write(114) seistxx
    close(114)

!	open(115,file='seistxz.dat',access='stream')
!    write(115) seistxz
!    close(115)

	open(116,file='seistzz.dat',access='stream')
    write(116) seistzz
    close(116)


	
!	deallocate(rhob,rhof,rhom,rho2,lamda,mu,alpha,m,b)
!    deallocate(rho,c11,c12,c13,c22,c23,c33)
!	deallocate(rho,lamda,mu)
!	deallocate(rhobhx,rhofhx,rhomhx,rho2hx,rhobhz,rhofhz,rhomhz,rho2hz)
!	deallocate(bhx,bhz)
!	deallocate(bxz)
    deallocate(rhoxz)	
	deallocate(txx,txz,tzz,vx,vz,seisvx,seisvz,seistxx,seistzz,seistxz)
	deallocate(dxtxx,dztxx,dxtxz,dztxz,dxtzz,dztzz,dxvx,dzvz,dxvz,dzvx)
	deallocate(psidxtxx,psidztzz,psidxtxz,psidztxz,psidzvx,psidxvz,psidzvz)
!	deallocate(rhobxz,rhofxz,rhomxz,rho2xz)	
!	deallocate(txx,txz,tzz,p,vx,vz,qx,qz,seisvx,seisvz,seisqx,seisqz,seistxx,seistzz,seistxz,seisp)
!	deallocate(dxtxx,dztxx,dxtxz,dztxz,dxtzz,dztzz,dxp,dzp,dxvx,dzvz,dxvz,dzvx,dxqx,dzqz,dxqz,dzqx)
!	deallocate(psidxtxx,psidztzz,psidxtxz,psidztxz,psidxp,psidzp,psidxvx,psidzvx,psidxvz,psidzvz,psidxqx,psidzqx,psidxqz,psidzqz)
	deallocate(ax,bx,az,bz,chix,chiz,cdx,cdz,ax_half,bx_half,az_half,bz_half,chix_half,chiz_half,cdx_half,cdz_half)
	deallocate(wavelet,wav)
	
!	deallocate(tmp)
    
	print *
	print *,'End of the simulation'
	pause
	print *
    
end program    
!*********************************
!*     END
!*********************************   



    
!*********************************
!*     Ricker Wavelet
!*********************************
      subroutine ricker(wv,nt,dt,f0,ts)
	  integer*4 nt
      real*4 wv(nt)
      real*4 aa,dt,f0,ts,pi
      integer ii

!cccccccccccccccccccccccccccccccccccccccc
    pi=4.*atan(1.0)
    do ii=1,nt
        aa=pi*f0*((ii-1.)*dt-ts)
        aa=aa*aa
        wv(ii)=(1.-2.*aa)*exp(-aa)
    enddo    
!c---------------------------------------
      return
      end subroutine

!*********************************
!*     dgauss wavelet
!*********************************
    subroutine dgauss(wv,nt,dt,f0,ts)
	integer*4 nt
    real*4 wv(nt)
    real*4 aa,dt,f0,ts,pi
    integer ii
!cccccccccccccccccccccccccccccccccccccccc
    pi = 4. * atan(1.0)
    do ii=1,nt
        aa = pi*f0*pi*f0
        
        wv(ii) = -2.0*aa*((ii-1)*dt-ts)*exp(-aa*((ii-1)*dt-ts)**2)
    end do
!cccccccccccccccccccccccccccccccccccccccc
    return
    end subroutine



!*********************************
!*     draw the snapshot
!*********************************
	subroutine creat_2D_image(sndat,nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,nr,np,ifield,power_display,cutvect)
	implicit none
!	real,parameter :: power_display=0.2
!	real,parameter :: cutvect=0.05
	logical,parameter :: white_background=.true.
	integer,parameter :: width_cross=5,thickness_cross=1,size_square=3

		integer nx,nz,it,xs,zs,fxr,fzr,dxr,dzr,np,nr
		real,dimension(nz,nx) :: sndat
		!integer,dimension(nr) :: xr,zr
		integer :: ix,iz,irec
		character(len=100) :: fname,system_command
		integer :: R,G,B
		real :: normalized,maxamp
		integer :: ifield
		real :: power_display,cutvect		

		if(ifield .eq. 0) then
			write(fname,"('image',i6.6,'_Vz.pnm')") it
			write(system_command,"('convert image',i6.6,'_Vz.pnm image',i6.6, &
				'_Vz.tif ; rm image',i6.6,'_Vz.pnm')") it,it,it
		else
			write(fname,"('image',i6.6,'_p.pnm')") it
			write(system_command,"('convert image',i6.6,'_p.pnm image',i6.6, &
				'_p.tif ; rm image',i6.6,'_p.pnm')") it,it,it
		endif

		open(unit=27, file=fname, status='unknown')

		write(27,"('P3')") ! write image in PNM P3 format

		write(27,*) nx,nz ! write image size
		write(27,*) '255' ! maximum value of each pixel color

! compute maximum amplitude
		maxamp = maxval(abs(sndat))

! image starts in upper-left corner in PNM format
		  do iz=nz,1,-1
			do ix=1,nx

! define data as vector component normalized to [-1:1] and rounded to nearest integer
! keeping in mind that amplitude can be negative
				normalized = sndat(iz,ix) / maxamp

! suppress values that are outside [-1:+1] to avoid small edge effects
				if(normalized < -1.0) normalized = -1.0
				if(normalized > 1.0) normalized = 1.0

! draw an orange cross to represent the source
				if((ix >= xs - width_cross .and. ix <= xs + width_cross .and. &
					iz >= zs - thickness_cross .and. iz <= zs + thickness_cross) .or. &
					(ix >= xs - thickness_cross .and. ix <= xs + thickness_cross .and. &
					iz >= zs - width_cross .and. iz <= zs + width_cross)) then
						 R = 255
						 G = 157
						 B = 0

! display two-pixel-thick black frame around the image
				else if(ix <= 2 .or. ix >= nx-1 .or. iz <= 2 .or. iz >= nz-1) then
						R = 0
						G = 0
						B = 0

! display edges of the PML layers
				else if((ix == np) .or. &
					(ix == nx - np) .or. &
					(iz == np) .or. &
					(iz == nz - np)) then
						R = 255
						G = 150
						B = 0

! suppress all the values that are below the threshold
				else if(abs(normalized) <= cutvect) then

! use a black or white background for points that are below the threshold
					if(white_background) then
						R = 255
						G = 255
						B = 255
					else
						R = 0
						G = 0
						B = 0
					endif

! represent regular image points using red if value is positive, blue if negative
				else if(normalized >= 0.0) then
					R = nint(255.0*normalized**power_display)
					G = 0
					B = 0
				else
					R = 0
					G = 0
					B = nint(255.0*abs(normalized)**power_display)
				endif

! draw a green square to represent the receivers
				do irec = 1,nr
					if((ix >= fxr+(irec-1)*dxr - size_square .and. ix <= fxr+(irec-1)*dxr  + size_square .and. &
						iz >= fzr+(irec-1)*dzr - size_square .and. iz <= fzr+(irec-1)*dzr + size_square) .or. &
						(ix >= fxr+(irec-1)*dxr  - size_square .and. ix <= fxr+(irec-1)*dxr  + size_square .and. &
						iz >= fzr+(irec-1)*dzr  - size_square .and. iz <= fzr+(irec-1)*dzr + size_square)) then
! use dark green color
							R = 30
							G = 180
							B = 60
					endif
				enddo

! write color pixel
				write(27,"(i3,' ',i3,' ',i3)") R,G,B

			enddo
		enddo

! close file
		close(27)

! call the system to convert image to GIF (can be commented out if "call system" is missing in your compiler)
		call system(system_command)
		print *,'convert over!'

		
		
	end subroutine 


