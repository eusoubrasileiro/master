program calor2D
!***********************************************************************
!************PROGRAMA QUE SOLUCIONA O PROBLEMA DE CALOR 2D**************
!*************NO REGIME TRANSIENTE ATRAVS DO MDF EXPLÖCITO*************
!**DISCIPLINA: "MODELAGEM NUMRICA DE ONDAS SÖSMICAS"(PPGG-ON)**********
!**PROFESSOR: LEANDRO DI BARTOLO****************************************
!**AULA 3 - DATA: 29/06/2015********************************************
!***********************************************************************
implicit none
!declara‡Æo de vari veis
integer :: nx,nz,nt !parametros numericos
real    :: dx,dz,dt,lbd !parametros numericos
real    :: tt1,tt2
real    :: x,z,t,alfa,Tesq,Tdir,Tsup,Tinf !parametros fisicos
real    :: start,finish !calcula tempo de execução
integer :: nsnap,ntsnap,count_snap=0 !numero de snapshots e contador
character(len=3) :: num_snap  !contador de snap para impressÆo
real, dimension(:,:), allocatable :: T1,T2 !vetores com as temperaturas
integer :: i,k,n !contadores
logical :: existe_arq
integer,parameter :: binpar=4

!call cpu_time(start)

call inicial() !Leitura de variaveis e calculos iniciais

!write(*,*) "Inicio da marcha no tempo..."
do n=1,nt!marcha no tempo
   !if(mod(n,50)==0) write(*,*) "passo",n !imprime na tela o avan‡o no tempo de 50 em 50 passos
   do i=2,nx-1 !loop do operador para c lculo da temperatura
       do k=2,nz-1
            T2(k,i)=T1(k,i)+lbd*(T1(k,i-1)-2*T1(k,i)+T1(k,i+1) + T1(k-1,i)-2*T1(k,i)+T1(k+1,i))
       enddo
   enddo
   !call imprime_snap() !imprime snapshot de temperatura
   T1(2:nz-1,2:nx-1)=T2(2:nz-1,2:nx-1) !atualiza as temperaturas para continuar a marcha
enddo
!call cpu_time(finish)
!write(*,*) "Fim da marcha no tempo! Tempo total de ",finish-start,"segundos."

!***********************************************************************
contains !subrotinas internas
!***********************************************************************

subroutine inicial()

implicit none

call CPU_time(tt1)

inquire(file="in.txt",exist=existe_arq) !testando se o arquivo est  na pasta
if(existe_arq) then
    !lendo parƒmetros de arquivo
    open(20,file="in.txt")
    read(20,*) x,z,t,nsnap
    read(20,*) dx,dz,dt
    read(20,*) alfa,Tesq,Tdir,Tsup,Tinf
    !c lculos iniciais
    lbd=alfa*dt/(dx*dx)
    nx=nint(x/dx)+1
    nz=nint(z/dz)+1
    nt=nint(t/dt)
    ntsnap = nt/nsnap
    !escreve na tela parametros lidos e calculados
    write(*,*) "Parametros lidos:"
    write(*,*) "x,z,t = ",x,",",z,",",t,", numero de snaps =",nsnap
    write(*,*) "dx,dz,dt = ",dx,",",dz,",",dt
    write(*,*) "alfa,Tesq,Tdir,Tsup,Tinf = ",alfa,",",Tesq,",",Tdir,",",Tsup,",",Tinf
    write(*,*) ""
    write(*,*) "Parametros numericos calculados:"
    write(*,*) "nx,nz,nt = ",nx,",",nz,",",nt,", ntsnap",ntsnap
    write(*,*) "lambda = ",lbd, "(ok! menor ou igual que 0.5)"
    write(*,*) ""
    if(lbd>0.5) then !teste de estabilidade
        write(*,*)"Lambda menor que 0.5: ESQUEMA INSTAVEL!"
        write(*,*)"Escolha parametros adequado e rode novamente."
        stop
    endif
    allocate(T1(nz,nx),T2(nz,nx)) !aloca‡Æo dinƒmica dos vetores T para a marcha no tempo
    T1=0.0 ; T2=0.0 !inicializa‡Æo dos vetores
    T1(:,1)=Tesq ; T1(:,nx)=Tdir ; T1(1,:)=Tsup ; T1(nz,:)=Tinf ; T2=T1 !aplica‡Æo das condi‡äes de contorno
else
    write(*,*)'arquivo de entrada "in.txt" nao encontrado!'
    stop
endif

call CPU_time(tt2)

Print*,"Tempo = ", tt2-tt1

end subroutine inicial

!***********************************************************************

subroutine imprime_snap()

implicit none

if (mod(n,ntsnap)==0) then !imprime snapshots de nsnap em nsnap passos
   count_snap=count_snap+1
   write(num_snap,"(I3.3)")count_snap
   write(*,*) "Imprimindo snap ", num_snap,"..."
   open(10,file="snap"//num_snap//".bin",status="replace",access="direct",form="unformatted",recl=binpar*nx*nz)
   write(10,rec=1) ((T2(k,i),k=1,nz),i=1,nx)
   close(10)
endif
end subroutine imprime_snap

end program calor2D
