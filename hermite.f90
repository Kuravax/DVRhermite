program hermite
!Hermite polynomial basis centered around x=0, 1D stationary state simulation

    implicit none
    integer:: nodes, memstatus,i, io, iostatus,j
    double precision:: r_min, r_max, mass, length, corrfactor, omeganot,lambda
    double precision, allocatable:: xmatrix(:,:),xeigenvalues(:),veigenvalues(:),hmatrix(:,:),heigenvalues(:),vmatrix(:,:)

    open(newunit=io, file="input.txt", status="old", action="read", iostat = iostatus)
    if(iostatus /=0) then
        stop "***File open error***"
    end if
    read(io,*)
    read(io, *), r_min, r_max, nodes, mass, omeganot
    close(io)

    allocate(xmatrix(nodes,nodes), xeigenvalues(nodes),veigenvalues(nodes),hmatrix(nodes,nodes),&
    heigenvalues(nodes), vmatrix(nodes,nodes), STAT = memstatus)
    if(memstatus /= 0) then
    stop "*** Not enough memory for allocation***"
    end if
    hmatrix=0.0
    xmatrix=0.0
    xeigenvalues=0.0
    veigenvalues=0.0
    heigenvalues = 0.0
    vmatrix=0.0

    length = r_max-r_min
    !Compute tridiagonal X matrix
    call FBRXmatrix(nodes,xmatrix)
    !Diagonalize tridiagonal X matrix
    !Diagonalized X matrix serves as DVR-FBR matrix
    call DIAG2(nodes,nodes,xeigenvalues,xmatrix)
    !Scale up the X diag matrix
    corrfactor = length/(xeigenvalues(nodes)-xeigenvalues(1))
    !dvr potential
    call DVRPotential(nodes,xeigenvalues,veigenvalues,corrfactor)
    !Swap dvr to fbr
    call FBRpotential(nodes,veigenvalues,vmatrix,xmatrix)
    !FBR hamiltonian
    call FBRHamiltonian(nodes,hmatrix,vmatrix,omeganot)

    !Print 20 eigenvalues
    call DIAG2(nodes,nodes,heigenvalues,hmatrix)
    do i=1, min(20,nodes)
        print*,heigenvalues(i)
    end do
    !Save output
    open(newunit=io, file="output.txt", status="old", action="write", iostat = iostatus)
    if(iostatus /=0) then
        stop "***File open error***"
    end if
    write(io,*)"First 20 eigenvalues, correct eigenvalues, oscillator eigenvalues"
    do i=1, min(20,nodes)
        write(io,*) heigenvalues(i),(2*lambda*(i-0.5)-(i-0.5)**2)*4.5d3/(lambda**2),&
        (i-0.5)*omeganot
    end do
        write(io,*)"First 20 eigenvectors"
    do i=1, min(20,nodes)
        write(io,*) hmatrix(i,:)
    end do
    close(io)

end program
!generate the tridiagonal X matrix in FBR
subroutine FBRXmatrix(nodes,xmatrix)
    integer l, nodes
    double precision const1, const2,xmatrix(nodes,nodes)
    const1 = sqrt(acos(-1.0d0))
    const2= 2*sqrt(acos(-1.0d0))
    xmatrix(1,2) = 0
    !0,1 = 0. Then continue as would usually proceed.
    do l=2, nodes-1
        xmatrix(l,l-1) = const1
        xmatrix(l,l+1) = l*const2
        const1 = 2*const1
        const2 = 2*const2
    end do
    xmatrix(nodes,nodes-1) = 0.5*const1
return
end subroutine

!Generate Hamiltonian matrix in FBR, given diagonal V in DVR and Transform to DVR
subroutine FBRHamiltonian(nodes,hmatrix,FBRV,omeganot)
    integer i,j,nodes
    double precision hmatrix(nodes,nodes),FBRV(nodes,nodes)
    double precision omeganot
    do i=1,nodes
        hmatrix(i,i)=(i-0.5)*omeganot !numbering starts n=1 for matrix; n=0 for eigenenergy
        do j=1,i
            hmatrix(i,j) = hmatrix(i,j)+ FBRV(i,j)
            hmatrix(j,i) = hmatrix(j,i)+ FBRV(j,i)
        end do
    end do

return
end subroutine

!Diagonal DVR potential to square FBR potential.
subroutine FBRpotential (nodes, VDVR, FBRV, Tmatrix)
    integer i,j,k,nodes
    double precision VDVR(nodes), FBRV(nodes,nodes), Tmatrix(nodes,nodes), temp
    do i=1,nodes
        do j=1,i
            temp = 0.0
            do k=1,nodes
                temp = temp + Tmatrix(i,k)*VDVR(k)*Tmatrix(k,j)
            end do
            FBRV(i,j) = temp
            FBRV(j,i)= temp
        end do
    end do
return
end subroutine

subroutine DVRpotential(nodes, XDVR, VDVR,corrfactor)
    integer i, nodes
    double precision XDVR(nodes), VDVR(nodes), corrfactor
    do i=1, nodes
        VDVR(i) = 4.5d3*(1-exp((XDVR(i)-1.5d0)*-4.06d-9))**2 - 1/2*(XDVR(i))**2
    end do
return
end subroutine
