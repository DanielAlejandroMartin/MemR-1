!***********************************************************************
!	Simple Code for Simulating The Memristor Metwork
! The program is organized as Follows:
! 1-Module Variables (includes All Variables needed for computation)
! 2- Main Program (Calls "Initial_Conditions and then LEARN_PROCESS
! 3-Subroutine Initial_Conditons
! 4- Subroutine LEARN_PROCESS: Most important one 
! 5- Subroutine MeasureProcedure
! 6- Subroutine DepressProcedure
! 7- Subroutine  Measure_Current
! 8- Subroutine Update_Resistence 
! 9- Subroutine Matrix_Make 
!10- Subroutine Update Variables



! The program solves a linear equation problem using  DGESV 
!From LAPACK (linear Algebra Package). LAPACK should be installed first. 
!Recommended operations From Terminal
! $ mkdir ResultsFolder
! $ gfortran ExampleCode.f90 -o Program -llapack
! $./Program

! RECOMEND TO OPEN IN GEANY



!***********************************************************************
!1- Module Vars  (In Module Vars I define all Variables)
!***********************************************************************



module vars
!********** NETWORK GEOMETRY **********
!Number of input (Ni) hiddden (Nj) and output (Nk) Neurons:
integer,parameter:: Ni=3
integer,parameter:: Nj=20
integer,parameter:: Nk=Ni
!********** Memristor & Nodes Voltabe/Resistente/Current **********
!Voltages at the input (Vi) Hidden and Ouptut (Nk) Neurons:
real(8),dimension(Ni)::Vi
real(8),dimension(Nj)::Vj
real(8),dimension(Nk)::Vk
!Resistences  between input and Hidden (Rij) and between Hidden and output (Rjk):
real(8),dimension(Ni,Nj)::Rij
real(8),dimension(Nj,Nk)::Rjk
!Currents  between input and Hidden (Iij) and between Hidden and output (Ijk):
real(8),dimension(Ni,Nj)::Iij
real(8),dimension(Nj,Nk)::Ijk

! Auxiliary integers that help writing the Linear Equiations Matrix:
integer,parameter:: Variables_Voltaje=Ni+Nj+Nk
integer,parameter:: Variables_CurrentA=Ni*Nj
integer,parameter::  Variables_CurrentB=Nj*Nk
integer,parameter:: Variables=Variables_Voltaje+Variables_CurrentA+Variables_CurrentB

!********** Learning Related Variables **********
real(8),dimension(Nk)::TotCorK ! Currents at each output Node
integer maxcorr !Output Node with Highest Current
INTEGER,DIMENSION(Ni):: truthtable !table to be learnt
real(8):: V0 !ApplyiedVoltages
real(8), parameter::  dt=1 !Time Discretization

!********** Learning-Algorithm Related Variables **********

integer step
integer Proces
integer To_Be_measured
INTEGER numberofDepressProcedures





!********** Solving-Linear-Equiatons Related Variables **********
!Variables For solving the linear algebra Problem
integer Equation_Number
real(8),dimension(Variables,variables):: MatrizEq
real(8),dimension(Variables,1):: Solutions_Vector
!Auxiliary variables for Linear Algebra subroutine 
integer:: pivot(Variables), ok


!**********  Memristor Parameters ********** 
! Memristor Parameters for input-hidden (ij) and for hidden-oputput (jk)
!memristors.
!beta: Resistence modification rate (See Equations)
!Ron: Lowest resistence	!Roff: Maximum resistence
!VT: Threshold Voltage	!This values are fixed at "Initial_Conditions" subroutine
real(8),dimension(Ni,Nj):: beta_ij 
real(8),dimension(Nj,Nk):: beta_jk 
real(8),dimension(Ni,Nj):: RON_ij 
real(8),dimension(Nj,Nk):: RON_jk
real(8),dimension(Ni,Nj):: ROFF_ij
real(8),dimension(Nj,Nk):: ROFF_jk 
real(8),dimension(Ni,Nj):: VT_ij 
real(8),dimension(Nj,Nk):: VT_jk

integer:: i,j,k !Auixilary Integers

end module

!***********************************************************************
!							MAIN PROGRAM		
!***********************************************************************
program MAIN !It Opens Evetything and calls LEARN_PROCESS
use vars
implicit none

call Initial_Conditions !subroutine for setting initial parameters for memrsitors

open (11,file="ResultsFolder/Large.dat",status="unknown") !Large Output
open (12,file="ResultsFolder/Short.dat",status="unknown") !Small Output

!***********************************************************************
!								DINAMICS		
!***********************************************************************

write(*,*) "Map to be learned"
do i=1,ni;
! I have chosen the Identity table for this example
truthtable(i)=i
write(*,*) i, truthtable(i)
enddo
write(*,*) "*********************"
write(*,*) "Step Error"
call LEARN_PROCESS

End Program


!***********************************************************************
!							Subroutines					
!***********************************************************************
! Initial Conditions

subroutine Initial_Conditions
use vars;! Variables of the whole system
implicit none
!Here we choose parameters for each memristor
!We also defien initial resistence

do i=1,ni; do j=1,nj
ROFF_ij(i,j)=10**6
beta_ij(i,j)=1.0-rand()*0.1 
VT_ij(i,j)= 10.**(-1.)*(1.+rand())/2.d0
Ron_ij(i,j)=100*(1+rand())/2.d0
enddo;enddo

do j=1,nj; do k=1,nk
ROFF_jk(j,k)=10**6
beta_jk(j,k)=1.0-rand()*0.1
VT_jk(j,k)= 10.**(-1.)*(1+rand())/2.d0
Ron_jk(j,k)=100*(1+rand())/2.d0
enddo;enddo
! We set each Memristor to its Highest Available Current
Rij=RON_ij
Rjk=RON_jk
end subroutine


!***********************************************************************
! Sobroutine LEARN_PROCESS
!It runds for up to s_max=50000 steps
!At each step it preforms up to C_max=80 Read prodedures
!***********************************************************************

subroutine LEARN_PROCESS
use vars!My sistem
implicit none
integer nerror

do step=1,50000 !s_max=50000
DO numberofDepressProcedures=1,80 !C_max=80

To_Be_measured=step+1-ni*int(step/ni)!Select which input Node we are "Writing"
V0=0.001
call MeasureProcedure
! Look For the ouptput with largest current.
maxcorr=1; 	do k=1,nk; 	if (totcork(k)>totcork(maxcorr)) maxcorr=k; enddo
!DepressProcedure if it is not the expected one.
if ((maxcorr-truthtable(To_Be_measured))**2>0) then
	V0=0.1
	call DepressProcedure 
	endif
ENDDO !numberofDepressProcedures 
1234 continue
Nerror=0


!CHECK
do To_Be_measured=1,ni
	V0=0.001; 	call MeasureProcedure
	maxcorr=1; 	do k=1,nk; 	if (totcork(k)>totcork(maxcorr)) maxcorr=k; enddo
if (maxcorr /= truthtable(To_Be_measured)) Nerror=Nerror+1

enddo
write(11,*) step, Nerror
write(*,*) step, Nerror
if (nerror==0) goto 5678
ENDDO !Experiment
call flush(11)
5678 continue
write(12,*) step

end subroutine

!***********************************************************************
!Subroutne DepressProcedure
!***********************************************************************


subroutine MeasureProcedure
use vars! Variables of the whole system
implicit none
!make matrix
	call Matrix_Make
!lapack solve linear equation
	call DGESV(variables, 1, MatrizEq, variables, pivot, Solutions_Vector, variables, ok)
!Update Variables, Update Memrisotors Memristence and Meassure Current on each output Node
	call Update_Variables
	call Update_Resistence
call Measure_Current
end subroutine


subroutine DepressProcedure
use vars! Variables of the whole system
implicit none
INTEGER ss !NUMBER OF UNIT TIMES WHERE i DepressProcedure

do ss=1,5
!dt=0.0025

	call Matrix_Make
	Solutions_Vector(ni+nj+maxcorr,1)=-V0
	call DGESV(variables, 1, MatrizEq, variables, pivot, Solutions_Vector, variables, ok)
	call Update_Variables
	call Update_Resistence
enddo
call Measure_Current
end subroutine


subroutine Measure_Current !measures current
use vars! Variables of the whole system
implicit none
	TotCorK=0
	do k=1,nk;	do j=1,nj;	TotCorK(k)=TotCorK(k)+Ijk(j,k); 	enddo; 		enddo; 

end subroutine	

subroutine Update_Resistence !Update Resistences
use vars! Variables of the whole system
implicit none
real(8) vm ! VM is the voltage Difference
real(8) F !F is the update function dR/dt= F Eq. 2 From text
do i=1,ni; do j=1,nj

Vm=Vi(i)-Vj(j)

F= beta_ij(i,j)*(Vm-0.5*(abs(Vm+VT_ij(i,j))-abs(Vm-VT_ij(i,j)) ))

Rij(i,j)=Rij(i,j)+f*dt

if (Rij(i,j)> ROFF_ij(i,j)) Rij(i,j)= ROFF_ij(i,j)
if (Rij(i,j)< RON_ij(i,j)) Rij(i,j)= RON_ij(i,j)
enddo; enddo

do j=1,nj; do k=1,nk
Vm=Vj(j)-Vk(k)

f= beta_jk(j,k)*(Vm-0.5*(abs(Vm+VT_jk(j,k))-abs(Vm-VT_jk(j,k)) ))


Rjk(j,k)=Rjk(j,k)+f*dt
if (Rjk(j,k)> ROFF_jk(j,k)) Rjk(j,k)= ROFF_jk(j,k)
if (Rjk(j,k)< RON_jk(j,k)) Rjk(j,k)= RON_jk(j,k)
enddo; enddo
end subroutine

subroutine Matrix_Make !Make the matrix
use vars! Variables of the whole system
implicit none
! variables are  ordered as follows
!Vi(1), Vi(2), .., Vi(Ni), Vj(1), .., Vj(Nj), Vk(1), .., Vk(Nk), 
!Iij(i=1,j=1),.., Iij(i=1,j=Nj), Iij(i=2,j=1), ..., Iij(i=Ni,j=nj),
!Ijk(j,k),..

! Set Matrix and Vector to 0
MatrizEq=0 !Matrix
Solutions_Vector=0 !SolutionVector


!********************** First Ni equations: Fix Vi=0 for all but one
do i=1,Ni
Matrizeq(i,i)=1
enddo
!Set To_Be_measured input neuron voltage as V0
i=To_Be_measured
Solutions_Vector(i,1)=V0


!*********************** Next Nj equations (Ni+1 to Ni+nj)
!Sum of currents over jth hiden node is 0
! Current Iij is variable number Variables_Voltaje+j+(i-1)*nj
do j=1,nj
do i=1,ni
Matrizeq(j+ni,Variables_Voltaje+j+(i-1)*nj)=1
enddo
! Current Ijk es variable number Variables_Voltaje+Variables_CurrentA+k+(j-1)*nk}
do k=1,nk
Matrizeq(j+ni,Variables_Voltaje+Variables_CurrentA+k+(j-1)*nk)=-1
enddo
enddo
!********************** !Next Nk equations (Ni+Nj+1 to Ni+nj+nk)
!Vk=0 !  Vk is variable Ni+Nj+k
do k=1,nk
Matrizeq(ni+nj+k,ni+nj+k)=1
enddo

!******************** !Next: Ni*Nj set of equations  (Eqs Ni+nj+nk to  Ni+nj+nk+Ni*Nj)
!Vi-Vj-I(i,j)*R(i,j)=0
! For input (i) to Hiden (j) Nodes
!Ecuation for i,j:
do i=1,ni
do j=1,nj
Equation_Number=Ni+nj+nk+j+(i-1)*nj
!Vj is elelment ni+j (of the vector of variables)
Matrizeq(Equation_Number,ni+j)=-1
!Vi is element i
Matrizeq(Equation_Number,i)=1
Matrizeq(Equation_Number,Variables_Voltaje+j+(i-1)*nj)=-1*Rij(i,j)
enddo
enddo

!*********************** !Next equiations (Ni+nj+nk+Ni*Nj+1 to the end)
!Vj-Vk-I(j,k)*R(j,k)
!Equiation  for j,k:
do j=1,nj
do k=1,nk
Equation_Number=Ni+nj+nk+Ni*Nj+k+(j-1)*nk
!Vk is elelment ni+nj+k
Matrizeq(Equation_Number,ni+nj+k)=-1
!Vj is elelment ni+j
Matrizeq(Equation_Number,ni+j)=1
! Current Ijk is variable number Variables_Voltaje+Variables_CurrentA+k+(j-1)*nk
Matrizeq(Equation_Number,Variables_Voltaje+Variables_CurrentA+k+(j-1)*nk)=-1*Rjk(j,k)
enddo
enddo


end subroutine

subroutine Update_Variables
use vars! Variables of the whole system
implicit none

!*** Voltage Variables
do i=1,ni
Vi(i)=Solutions_Vector(i,1)
enddo
do j=1,nj
Vj(j)=Solutions_Vector(ni+j,1)
enddo
do k=1,nk
Vk(k)=Solutions_Vector(ni+nj+k,1)
enddo

!*** Current Variables From Input to Hidden
do i=1,ni
do j=1,nj
Iij(i,j)=  Solutions_Vector(Variables_Voltaje+j+(i-1)*nj,1)
enddo
enddo
!*** Current Variables From Hidden  to Output
do j=1,nj
do k=1,nk
Ijk(j,k)= Solutions_Vector(Variables_Voltaje+Variables_CurrentA+k+(j-1)*nk,1)
enddo
enddo

end subroutine
