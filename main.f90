program sincdvr
    implicit none
    integer:: nodes, memstatus,i,io, iostatus,j
    double precision:: r_min, r_max, mass, length,Fconstant, r0
    double precision, allocatable:: xmatrix(:), vmatrix(:), tmatrix(:,:), kmatrix(:), hmatrix(:,:),eigenvalues(:),&
    sinpoints(:,:)

    open(newunit=io, file="input.txt", status="old", action="read", iostat = iostatus)
    if(iostatus /=0) then
        stop "***File open error***"
    end if
    read(io,*)
    read(io, *) r_min, r_max,nodes, mass
    read(io,*)
    read(io,*) Fconstant, r0
    close(io)

    allocate(xmatrix(nodes), vmatrix(nodes), tmatrix(nodes,nodes), kmatrix(nodes), &
    & hmatrix(nodes,nodes), eigenvalues(nodes), sinpoints(2,5*nodes+1),STAT = memstatus)
    if(memstatus /= 0) then
    stop "*** Not enough memory for allocation***"
    end if
    length = r_max - r_min
    call discretegrid(r_min,r_max,nodes,xmatrix)
    call DVRpotential(vmatrix,xmatrix,nodes,Fconstant,r0)
    call TransformMatrix(nodes,tmatrix)
    call FBRkinetic(mass,nodes,kmatrix,length)
    call FBRHamiltonian(nodes, tmatrix, vmatrix, kmatrix, hmatrix)
    !diagonalise hamiltonian matrix
    call DIAG2(nodes,nodes,eigenvalues,hmatrix)
    !save the first 20 eigenvalues and numeric error
    open(newunit=io, file="output.txt")
    if(iostatus /=0) then
        stop "***Output file error***"
    end if
    write(io,*) "Computed first 20 eigenvalues of Hamiltonian", "Nodes", nodes
    write(io,*) "n+1 ", "energy, hartree "
    do i=1,min(nodes,20)
        write(io,*) i, eigenvalues(i)
        print*, i, eigenvalues(i)
    end do
!    write(io,*) "Hamiltonian eigenvectors"
!    do i=1,min(nodes,20)
!        write(io,*) "Stationary State No. ", i
!        write(io,*) hmatrix(:,i)
!    end do

    !Print the scatter plot of first 4 stationary states.
    close(io)
    print*, "Output saved."
    open(newunit=io, file="plots.txt")
    if(iostatus /=0) then
        stop "***Output file error***"
    end if
    do j=1,4
        call FunctionPlot(nodes,sinpoints,hmatrix,r_min,r_max,j)
        write(io,*) "Wavefunction of state n=",j
        do i=1, 5*nodes+1
            write(io,*) sinpoints(:,i)
        end do
    end do
    print*,"Plot saved"
    close(io)

    deallocate(xmatrix, kmatrix, vmatrix, tmatrix,eigenvalues,hmatrix, STAT = memstatus)
    if(memstatus /= 0) then
    stop "*** Deallocation fail***"
    end if

end program


! generate the diagonal X matrix with (nodes) points
subroutine discretegrid(r_min,r_max,nodes,xmatrix)
    integer:: nodes, i
    double precision:: r_min, r_max, xmatrix(nodes), step
    step = (r_max-r_min)/(nodes+1)
    do i=1,nodes
    xmatrix(i) = r_min + step*i
    end do
return
end subroutine

! generate the diagonal kinetic energy matrix for FBR representation
! use basis of PIB functions
subroutine FBRkinetic(mass,nodes,kmatrix,step)
    integer:: nodes, i
    double precision:: mass, kmatrix(nodes), const, step
    const = (acos(-1.0)/step)**2 /(2*mass)
    do i=1,nodes
        kmatrix(i)= const * i*i
    end do
return
end subroutine

! generate the diagonal potential energy matrix for DVR representation
! 1. use harmonic oscillator potential
! 2. use coulomb field potential
subroutine DVRpotential(vmatrix,xmatrix,nodes, Fconstant, r0)
    integer:: nodes,i
    double precision:: vmatrix(nodes), xmatrix(nodes), Fconstant, r0, temp
    temp = Fconstant*0.5
    do i=1,nodes
        vmatrix(i)= temp * (xmatrix(i)-r0)**2
        !vmatrix(i) = temp* cosh(xmatrix(i)-r0)
    end do
return
end subroutine

! generate the transformation matrix from FBR to DVR
! T matrix is unitary and Hermitian
subroutine TransformMatrix (nodes, tmatrix)
    integer:: nodes,i ,j
    double precision:: tmatrix(nodes,nodes),tempconst, tempconst2,temp
    tempconst = sqrt(2.0/(nodes+1))
    tempconst2= acos(-1.0d0)/(nodes+1)
    do i=1,nodes
        do j=1,i
            temp = sin(i*j*tempconst2)*tempconst
            tmatrix(j,i) = temp
            tmatrix(i,j) = temp
        end do
    end do
end subroutine

! generate the Hamiltonian matrix in FBR
! needs to have diagonal V(DVR) and K(FBR) matrices to give
subroutine FBRHamiltonian(nodes, tmatrix, vmatrix, kmatrix, hmatrix)
    integer i,j,nodes,k
    double precision tmatrix(nodes,nodes), hmatrix(nodes,nodes), vmatrix(nodes),kmatrix(nodes), temp
    do i=1,nodes
        hmatrix(i,i) = kmatrix(i) !diagonal kinetic energy in FBR
        do j=1,i
            temp = 0.0d0
            do k=1,nodes
                temp = temp + (tmatrix(k,i)*vmatrix(k)*tmatrix(k,j)) !Swap basis of V from dvr to fbr
            end do
            !use the fact H is hermitian
            hmatrix(i,j) = temp + hmatrix(i,j)
            hmatrix(j,i) = temp + hmatrix(j,i)
        end do
    end do
return
end subroutine

!Plot the nth stationary state wavefunction on 5*nodes grid.
!Needs H eigenvectors matrix in FBR
subroutine FunctionPlot(nodes,sinpoints,hmatrix,r_min,r_max,eigen)
    double precision sinpoints(2,5*nodes), hmatrix(nodes,nodes)
    double precision r_min, r_max, length, const1, const2, step, temp
    integer nodes, eigen,i, j
    length = r_max - r_min
    const1 = sqrt(2/length)
    const2 = acos(-1.0)/length
    step = length/(5*nodes)
    do i=1, 5*nodes+1
        sinpoints(1,i) = r_min+(i-1)*step
    end do
    do i=1, 5*nodes+1
        temp = 0.0
        do j=1, nodes
            temp = temp + hmatrix(j,eigen)*const1*sin(j*const2*(sinpoints(1,i)-r_min))
        end do
        sinpoints(2,i) = temp
    end do
    !sqrt(2/L)* sin(npi(x-rmin)/L) - basis function of fbr
return
end subroutine
