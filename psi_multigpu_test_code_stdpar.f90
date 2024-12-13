!#######################################################################
!
! ****** PSI_MULTIGPU_TEST_CODE_STDPAR
!
!     This code mimics the basic MPI+DC tasks of PSI's
!     MAS Solar MHD code.
!
!     It sets up a Cartesian MPI topology, sets up a 3D grid
!     and tests a "seam" (point to point) MPI communication
!     using asyncronous Send and Recv calls.
!     This is used both on a allocatable sub-array, as well as
!     on a local buffer static array.
!
!     The code will automatically configure the topology based on
!     the number of MPI ranks it is called with.
!
!     If this code works on a multi-node GPU system,
!     than (most likely) so will MAS!
!
!     Author:  Ronald M. Caplan
!
!     Predictive Science Inc.
!     www.predsci.com
!     San Diego, California, USA 92121
!
!#######################################################################
! Copyright 2024 Predictive Science Inc.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!#######################################################################
!
module number_types
!
      use iso_fortran_env
!
      implicit none
!
      integer, parameter :: r_typ=REAL64
!
end module
!#######################################################################
module types
!
      use number_types
!
      type :: vvec
        real(r_typ), dimension(:,:,:), allocatable :: r !(nrm,nt,np)
        real(r_typ), dimension(:,:,:), allocatable :: t !(nr,ntm,np)
        real(r_typ), dimension(:,:,:), allocatable :: p !(nr,nt,npm)
      end type
!
end module
!#######################################################################
module mpidefs
!
      use mpi
!
      implicit none
!
! ****** Total number of processors.
      integer :: nproc
! ****** Total number of processors per node.
      integer :: nprocsh
! ****** Processor rank of this process in communicator
! ****** MPI_COMM_WORLD.
      integer :: iprocw
! ****** Processor rank of this process in communicator
! ****** comm_shared.
      integer :: iprocsh
! ****** Flag to designate that this is the processor with
! ****** rank 0 in communicator MPI_COMM_WORLD.
      logical :: iamp0
! ****** Communicator over all processors in the Cartesian topology.
      integer :: comm_all
! ****** Processor rank of this process in communicator
! ****** COMM_ALL.
      integer :: iproc
! ****** Processor rank in communicator COMM_ALL for the
! ****** processor that has rank 0 in MPI_COMM_WORLD.
      integer :: iproc0
! ****** Communicators over all processors in the phi dimension.
      integer :: comm_phi
! ****** Communicator over all shared processors on the node.
      integer :: comm_shared
! ****** Communicators over all processors in the theta and phi
! ****** dimensions.
      integer :: comm_tp
! ****** Communicators over all processors in the r dimension.
      integer :: comm_r
! ****** Processor rank in communicator COMM_R of the processor
! ****** that contains the lower radial boundary r=R0.
      integer :: iproc_rb0
! ****** Processor coordinate indices of this process
! ****** in the Cartesian topology.
      integer :: iproc_r,iproc_t,iproc_p
! ****** Processor coordinate indices of the neighboring
! ****** processors in the Cartesian topology.
      integer :: iproc_rm,iproc_rp
      integer :: iproc_tm,iproc_tp
      integer :: iproc_pm,iproc_pp
! ****** Number of processors along r, theta, and phi.
      integer :: nproc_r,nproc_t,nproc_p
! ****** Number of processors in 2D theta-phi plane.
      integer :: nproc_tp
! ****** Processor coordinate indices in 2D theta-phi plane.
      integer :: iproc2d_tp
! ****** Number type for REALs to be used in MPI calls.
      integer :: ntype_real
end module
!#######################################################################
module decomposition_params
!
!-----------------------------------------------------------------------
! ****** Input parameters that define the domain decomposition
! ****** among processors.
!-----------------------------------------------------------------------
!
      implicit none
!
! ****** Number of processors per dimension.
!
      integer, dimension(3) :: nprocs=(/-1,-1,-1/)
!
! ****** Number of mesh points per processor.
!
      integer, parameter :: mx_procs_per_dim=100
      integer, dimension(mx_procs_per_dim) :: mp_r=0
      integer, dimension(mx_procs_per_dim) :: mp_t=0
      integer, dimension(mx_procs_per_dim) :: mp_p=0
!
! ****** Mesh sizes for the "automatic" mesh decomposition.
!
      integer :: nr_auto
      integer :: nt_auto
      integer :: np_auto
!
end module
!#######################################################################
subroutine init_mpi
!
!-----------------------------------------------------------------------
!
! ****** Initialize MPI.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr,tcheck
!
!-----------------------------------------------------------------------
!
      call MPI_Init_thread (MPI_THREAD_FUNNELED,tcheck,ierr)
!
! ****** Get the total number of processors.
!
      call MPI_Comm_size (MPI_COMM_WORLD,nproc,ierr)
!
! ****** Get the index (rank) of the local processor in
! ****** communicator MPI_COMM_WORLD in variable IPROCW.
!
      call MPI_Comm_rank (MPI_COMM_WORLD,iprocw,ierr)
!
! ****** Create a shared communicator for all ranks in the node.
!
      call MPI_Comm_split_type (MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0, &
                                MPI_INFO_NULL,comm_shared,ierr)
!
! ****** Get the total number of processors in node.
!
      call MPI_Comm_size (comm_shared,nprocsh,ierr)
!
! ****** Get the index (rank) of the local processor in the local node.
!
      call MPI_Comm_rank (comm_shared,iprocsh,ierr)
!
! ****** Set the flag to designate whether this processor
! ****** has rank 0 in communicator MPI_COMM_WORLD.
!
      if (iprocw.eq.0) then
        iamp0=.true.
      else
        iamp0=.false.
      end if
!
      ntype_real=MPI_REAL8
!
end subroutine
!#######################################################################
subroutine check_proc_topology
!
!-----------------------------------------------------------------------
!
! ****** Check/set the requested processor topology.
!
!-----------------------------------------------------------------------
!
      use mpidefs
      use decomposition_params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,nreq
!
!-----------------------------------------------------------------------
!
! ****** Set the optimal values of the topology for unset dimensions.
!
      call set_proc_topology
!
! ****** Check that the number of processors available
! ****** matches the number requested.
!
      nreq=nprocs(1)*nprocs(2)*nprocs(3)
!
      if (nreq.ne.nproc) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_PROC_TOPOLOGY:'
          write (*,*) '### Processor topology specification error.'
          write (*,*) 'The number of processors requested does not'// &
                      ' equal the number available.'
          write (*,*) 'Number of processors requested = ',nreq
          write (*,*) 'Number of processors available = ',nproc
        end if
      end if
!
end subroutine
!#######################################################################
subroutine set_proc_topology
!
!-----------------------------------------------------------------------
!
! ****** Set the optimal values of the MPI rank topology
! ****** in dimensions not set by user.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use decomposition_params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1.0_r_typ
      real(r_typ), parameter :: zero=0.0_r_typ
      real(r_typ), parameter :: bigval=HUGE(1.0_r_typ)
!
!-----------------------------------------------------------------------
!
      integer, dimension(:), allocatable :: factors
      integer, dimension(:,:), allocatable :: rank_factors
      real(r_typ), dimension(:,:), allocatable :: nperrank
      real(r_typ), dimension(:), allocatable :: penalty
!
      integer :: i,j,k,fr,ft,fp,num_fac,num_rank_fac,best_idx
      real(r_typ) :: a12,a13,a23
!
!-----------------------------------------------------------------------
!
! ****** Extract nproc values.  A value of -1 indicates the dimension
! ****** should be autoset.
!
      nproc_r=nprocs(1)
      nproc_t=nprocs(2)
      nproc_p=nprocs(3)
!
! ****** If no dimensions are to be autoset, return.
!
      if(nproc_r.ne.-1.and.nproc_t.ne.-1.and.nproc_p.ne.-1) return
!
! ****** Get all factors of nproc and store them in factors array.
!
      i=1
      num_fac=0
      do while(i.le.nproc)
        if (MOD(nproc,i).eq.0) then
          num_fac=num_fac+1
        end if
        i=i+1
      enddo
      allocate (factors(num_fac))
      i=1
      num_fac=0
      do while(i.le.nproc)
        if (MOD(nproc,i).eq.0) then
          num_fac=num_fac+1
          factors(num_fac)=i
        end if
        i=i+1
      enddo
!
! ****** Set penalty function parameters and any fixed dimensions
! ****** based on which dimensions are to be autoset.
!
      a12=one
      a13=one
      a23=one
!
      if (nproc_r.ne.-1) then
        fr=nproc_r
        a12=zero
        a13=zero
      end if
      if (nproc_t.ne.-1) then
        ft=nproc_t
        a12=zero
        a23=zero
      end if
      if (nproc_p.ne.-1) then
        fp=nproc_p
        a13=zero
        a23=zero
      end if
!
! ****** Loop over all combinations of factors and save those that
! ****** yield the correct number of MPI ranks into rank_factors array.
!
      num_rank_fac=0
      do k=1,num_fac
        do j=1,num_fac
          do i=1,num_fac
            if(nproc_r.eq.-1) fr=factors(i)
            if(nproc_t.eq.-1) ft=factors(j)
            if(nproc_p.eq.-1) fp=factors(k)
            if (fr*ft*fp.eq.nproc) then
              num_rank_fac=num_rank_fac+1
            end if
          enddo
        enddo
      enddo
!
      if (num_rank_fac.eq.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in SET_PROC_TOPOLOGY:'
          write (*,*) '### Processor topology specification error.'
          write (*,*) 'No valid topologies found for selected options.'
          write (*,*) 'Number of MPI ranks = ',nproc
          write (*,*) 'NPROC_R = ',nproc_r
          write (*,*) 'NPROC_T = ',nproc_t
          write (*,*) 'NPROC_P = ',nproc_p
        end if
      end if
!
      allocate(rank_factors(num_rank_fac,3))
      allocate(nperrank(num_rank_fac,3))
      allocate(penalty(num_rank_fac))
!
      rank_factors(:,:)=-1
      penalty(:)=bigval
!
      num_rank_fac=0
      do k=1,num_fac
        do j=1,num_fac
          do i=1,num_fac
            if(nproc_r.eq.-1) fr=factors(i)
            if(nproc_t.eq.-1) ft=factors(j)
            if(nproc_p.eq.-1) fp=factors(k)
            if (fr*ft*fp.eq.nproc) then
              num_rank_fac=num_rank_fac+1
              rank_factors(num_rank_fac,1)=fr
              rank_factors(num_rank_fac,2)=ft
              rank_factors(num_rank_fac,3)=fp
            end if
          enddo
        enddo
      enddo
!
! ****** Get number of grid points per rank for each dimension.
!
!      nperrank(:,1)=real(nr_g)/rank_factors(:,1)
!      nperrank(:,2)=real(nt_g)/rank_factors(:,2)
!      nperrank(:,3)=real(np_g)/rank_factors(:,3)
!
! ****** Compute penalty function.
!
      penalty(:)=a12*(rank_factors(:,1)-rank_factors(:,2))**2 &
                +a23*(rank_factors(:,2)-rank_factors(:,3))**2 &
                +a13*(rank_factors(:,3)-rank_factors(:,1))**2
!
! ****** Eliminate any choices that yield less than a minimum number
! ****** of grid points per rank.
!
!      do i=1,num_rank_fac
!        if (nperrank(i,1).lt.4) penalty(i)=bigval
!        if (nperrank(i,2).lt.4) penalty(i)=bigval
!        if (nperrank(i,3).lt.3) penalty(i)=bigval
!      enddo
!
! ****** Find optimal topology.
!
      best_idx=MINLOC(penalty,1)
!
      if (penalty(best_idx).eq.bigval) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in SET_PROC_TOPOLOGY:'
          write (*,*) '### Processor topology specification error.'
          write (*,*) 'No valid topologies found for selected options'
          write (*,*) 'with selected grid.  '
          write (*,*) 'It is likely you are using too many MPI ranks.'
          write (*,*) 'Number of MPI ranks = ',nproc
          write (*,*) 'NPROC_R = ',nproc_r
          write (*,*) 'NPROC_T = ',nproc_t
          write (*,*) 'NPROC_P = ',nproc_p
!          write (*,*) 'NR = ',nr_g
!          write (*,*) 'NT = ',nt_g
!          write (*,*) 'NP = ',np_g
        end if
      end if
!
! ****** Set optimal topology.
!
      nprocs(1)=rank_factors(best_idx,1)
      nprocs(2)=rank_factors(best_idx,2)
      nprocs(3)=rank_factors(best_idx,3)
!
      deallocate(factors)
      deallocate(rank_factors)
      deallocate(nperrank)
      deallocate(penalty)
!
end subroutine
!#######################################################################
subroutine decompose_domain
!
!-----------------------------------------------------------------------
!
! ****** Decompose the domain into a Cartesian MPI topology.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use decomposition_params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      integer, parameter :: ndim=3
      integer, dimension(ndim) :: coords
      logical, dimension(ndim) :: periodic
      logical :: reorder
      logical, dimension(ndim) :: keep_dim
!
!-----------------------------------------------------------------------
!
! ****** Create a communicator over all processors, COMM_ALL,
! ****** that has a Cartesian topology.
!
! ****** Specify the periodicity of the coordinate system.
!
      periodic(1)=.false.
      periodic(2)=.false.
      periodic(3)=.true.
!
! ****** Allow re-ordering in the Cartesian topology.
!
      reorder=.true.
!
      call MPI_Cart_create (MPI_COMM_WORLD,ndim,nprocs, &
                            periodic,reorder,comm_all,ierr)
!
! ****** Get the index (rank) of the local processor in
! ****** communicator COMM_ALL in variable IPROC.
!
! ****** IMPORTANT NOTE:
! ****** If re-odering was allowed in the Cartesian topology
! ****** creation (above), then the rank of the local processor
! ****** in communicator COMM_ALL may be different from its rank
! ****** in communicator MPI_COMM_WORLD.
!
      call MPI_Comm_rank (comm_all,iproc,ierr)
!
! ****** Set the processor rank IPROC0 in communicator COMM_ALL
! ****** for the processor that has rank 0 in MPI_COMM_WORLD.
! ****** This value is broadcast to all the processors.
!
      if (iamp0) then
        iproc0=iproc
      end if
      call MPI_Bcast (iproc0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
! ****** Get the coordinate indices of this processor in the
! ****** Cartesian MPI topology.
!
      call MPI_Cart_coords (comm_all,iproc,ndim,coords,ierr)
!
      iproc_r=coords(1)
      iproc_t=coords(2)
      iproc_p=coords(3)
!
      nproc_r=nprocs(1)
      nproc_t=nprocs(2)
      nproc_p=nprocs(3)
!
      nproc_tp=nproc_t*nproc_p
!
! ****** Get the rank of the neighboring processors in the
! ****** Cartesian MPI topology.
!
      call MPI_Cart_shift (comm_all,0,1,iproc_rm,iproc_rp,ierr)
      call MPI_Cart_shift (comm_all,1,1,iproc_tm,iproc_tp,ierr)
      call MPI_Cart_shift (comm_all,2,1,iproc_pm,iproc_pp,ierr)
!
! ****** Create communicators for operations involving all
! ****** processors in the phi dimension.  These communicators
! ****** are stored in COMM_PHI (and generally represent different
! ****** communicators on different processors).
!
      keep_dim(1)=.false.
      keep_dim(2)=.false.
      keep_dim(3)=.true.
!
      call MPI_Cart_sub (comm_all,keep_dim,comm_phi,ierr)
!
! ****** Create communicators for operations involving
! ****** all processors in the theta and phi dimensions.
! ****** These communicators are stored in COMM_TP
! ****** (and generally represent different communicators on
! ****** different processors).
! ****** These communicators are used for operations that
! ****** involve radial planes.
!
      keep_dim(1)=.false.
      keep_dim(2)=.true.
      keep_dim(3)=.true.
!
      call MPI_Cart_sub (comm_all,keep_dim,comm_tp,ierr)
!
! ****** Get rank in the theta-phi communicator.
! ****** This is used for 2D IO.
!
      call MPI_Comm_rank (comm_tp,iproc2d_tp,ierr)
!
! ****** Create communicators for operations involving
! ****** all processors in the r dimension.
! ****** These communicators are stored in COMM_R
! ****** (and generally represent different communicators on
! ****** different processors).
!
      keep_dim(1)=.true.
      keep_dim(2)=.false.
      keep_dim(3)=.false.
!
      call MPI_Cart_sub (comm_all,keep_dim,comm_r,ierr)
!
end subroutine
!#######################################################################
subroutine seam_vvec (v)
!
!-----------------------------------------------------------------------
!
! ****** Seam the boundary points of a v vector between adjacent
! ****** processors.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types, ONLY : vvec
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(vvec) :: v
!
!-----------------------------------------------------------------------
!
      real(r_typ),dimension(:,:),allocatable :: sbuf1r,rbuf1r
      real(r_typ),dimension(:,:),allocatable :: sbuf2r,rbuf2r
      real(r_typ),dimension(:,:),allocatable :: sbuf1t,rbuf1t
      real(r_typ),dimension(:,:),allocatable :: sbuf2t,rbuf2t
      real(r_typ),dimension(:,:),allocatable :: sbuf1p,rbuf1p
      real(r_typ),dimension(:,:),allocatable :: sbuf2p,rbuf2p
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr
!
! ****** MPI tags for MPI_ISEND and MPI_IRECV.
!
      integer :: tagr=0
      integer :: tagt=1
      integer :: tagp=2
!
!-----------------------------------------------------------------------
!
      integer :: lbuf3r,lbuf3t,lbuf3p
      integer :: lbuf1r,lbuf1t,lbuf1p
      integer :: lbuf2r,lbuf2t,lbuf2p
      integer :: n1r,n2r,n3r,n1t,n2t,n3t,n1p,n2p,n3p
      integer :: req(12)
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
! ****** Get the dimensions of the arrays and buffer sizes:
!
      n1r=size(v%r,1);   n2r=size(v%r,2);   n3r=size(v%r,3)
      n1t=size(v%t,1);   n2t=size(v%t,2);   n3t=size(v%t,3)
      n1p=size(v%p,1);   n2p=size(v%p,2);   n3p=size(v%p,3)
!
      lbuf3r=n1r*n2r;    lbuf3t=n1t*n2t;    lbuf3p=n1p*n2p
      lbuf1r=n2r*n3r;    lbuf1t=n2t*n3t;    lbuf1p=n2p*n3p
      lbuf2r=n1r*n3r;    lbuf2t=n1t*n3t;    lbuf2p=n1p*n3p
!
! ****** Seam the third (periodic) dimension. Since seam data
!        is stride-1 in this case, no buffers are needed.
!
! ****** Launch async receives.
!
      call MPI_Irecv (v%r(:,:,  1),lbuf3r,ntype_real,iproc_pm,tagr, &
                      comm_all,req(1),ierr)
      call MPI_Irecv (v%r(:,:,n3r),lbuf3r,ntype_real,iproc_pp,tagr, &
                      comm_all,req(2),ierr)
      call MPI_Irecv (v%t(:,:,  1),lbuf3t,ntype_real,iproc_pm,tagt, &
                      comm_all,req(3),ierr)
      call MPI_Irecv (v%t(:,:,n3t),lbuf3t,ntype_real,iproc_pp,tagt, &
                      comm_all,req(4),ierr)
      call MPI_Irecv (v%p(:,:,  1),lbuf3p,ntype_real,iproc_pm,tagp, &
                      comm_all,req(5),ierr)
      call MPI_Irecv (v%p(:,:,n3p),lbuf3p,ntype_real,iproc_pp,tagp, &
                      comm_all,req(6),ierr)
!
! ****** Launch async sends.
!
      call MPI_Isend (v%r(:,:,n3r-1),lbuf3r,ntype_real,iproc_pp,tagr, &
                      comm_all,req(7),ierr)
      call MPI_Isend (v%r(:,:,    2),lbuf3r,ntype_real,iproc_pm,tagr, &
                      comm_all,req(8),ierr)
      call MPI_Isend (v%t(:,:,n3t-1),lbuf3t,ntype_real,iproc_pp,tagt, &
                      comm_all,req(9),ierr)
      call MPI_Isend (v%t(:,:,    2),lbuf3t,ntype_real,iproc_pm,tagt, &
                      comm_all,req(10),ierr)
      call MPI_Isend (v%p(:,:,n3p-1),lbuf3p,ntype_real,iproc_pp,tagp, &
                      comm_all,req(11),ierr)
      call MPI_Isend (v%p(:,:,    2),lbuf3p,ntype_real,iproc_pm,tagp, &
                      comm_all,req(12),ierr)
!
! ****** Wait for all seams to complete.
!
      call MPI_Waitall (12,req,MPI_STATUSES_IGNORE,ierr)
!
! ****** Seam the first dimension.
!
      if (nproc_r.gt.1) then
!
! ****** Load buffers.
!
        allocate (sbuf1r(n2r,n3r),rbuf1r(n2r,n3r), &
                  sbuf2r(n2r,n3r),rbuf2r(n2r,n3r), &
                  sbuf1t(n2t,n3t),rbuf1t(n2t,n3t), &
                  sbuf2t(n2t,n3t),rbuf2t(n2t,n3t), &
                  sbuf1p(n2p,n3p),rbuf1p(n2p,n3p), &
                  sbuf2p(n2p,n3p),rbuf2p(n2p,n3p))
!
        do concurrent (k=1:n3r,j=1:n2r)
          sbuf1r(j,k)=v%r(n1r-1,j,k)
          sbuf2r(j,k)=v%r(    2,j,k)
        enddo
!
        do concurrent (k=1:n3t,j=1:n2t)
          sbuf1t(j,k)=v%t(n1t-1,j,k)
          sbuf2t(j,k)=v%t(    2,j,k)
        enddo
!
        do concurrent (k=1:n3p,j=1:n2p)
          sbuf1p(j,k)=v%p(n1p-1,j,k)
          sbuf2p(j,k)=v%p(    2,j,k)
        enddo
!
        call MPI_Irecv (rbuf1r,lbuf1r,ntype_real,iproc_rm,tagr, &
                        comm_all,req(1),ierr)
        call MPI_Irecv (rbuf2r,lbuf1r,ntype_real,iproc_rp,tagr, &
                        comm_all,req(2),ierr)
        call MPI_Irecv (rbuf1t,lbuf1t,ntype_real,iproc_rm,tagt, &
                        comm_all,req(3),ierr)
        call MPI_Irecv (rbuf2t,lbuf1t,ntype_real,iproc_rp,tagt, &
                        comm_all,req(4),ierr)
        call MPI_Irecv (rbuf1p,lbuf1p,ntype_real,iproc_rm,tagp, &
                        comm_all,req(5),ierr)
        call MPI_Irecv (rbuf2p,lbuf1p,ntype_real,iproc_rp,tagp, &
                        comm_all,req(6),ierr)
!
! ****** Launch async sends.
!
        call MPI_Isend (sbuf1r,lbuf1r,ntype_real,iproc_rp,tagr, &
                        comm_all,req(7),ierr)
        call MPI_Isend (sbuf2r,lbuf1r,ntype_real,iproc_rm,tagr, &
                        comm_all,req(8),ierr)
        call MPI_Isend (sbuf1t,lbuf1t,ntype_real,iproc_rp,tagt, &
                        comm_all,req(9),ierr)
        call MPI_Isend (sbuf2t,lbuf1t,ntype_real,iproc_rm,tagt, &
                        comm_all,req(10),ierr)
        call MPI_Isend (sbuf1p,lbuf1p,ntype_real,iproc_rp,tagp, &
                        comm_all,req(11),ierr)
        call MPI_Isend (sbuf2p,lbuf1p,ntype_real,iproc_rm,tagp, &
                        comm_all,req(12),ierr)
!
! ****** Wait for all seams to complete.
!
        call MPI_Waitall (12,req,MPI_STATUSES_IGNORE,ierr)
!
! ****** Unload buffers.
!
        if (iproc_rm.ne.MPI_PROC_NULL) then
          do concurrent (k=1:n3r,j=1:n2r)
            v%r(1,j,k)=rbuf1r(j,k)
          enddo
          do concurrent (k=1:n3t,j=1:n2t)
            v%t(1,j,k)=rbuf1t(j,k)
          enddo
          do concurrent (k=1:n3p,j=1:n2p)
            v%p(1,j,k)=rbuf1p(j,k)
          enddo
        end if
!
        if (iproc_rp.ne.MPI_PROC_NULL) then
          do concurrent (k=1:n3r,j=1:n2r)
            v%r(n1r,j,k)=rbuf2r(j,k)
          enddo
          do concurrent (k=1:n3t,j=1:n2t)
            v%t(n1t,j,k)=rbuf2t(j,k)
          enddo
          do concurrent (k=1:n3p,j=1:n2p)
            v%p(n1p,j,k)=rbuf2p(j,k)
          enddo
        end if
!
        deallocate (sbuf1r,sbuf2r,sbuf1t,sbuf2t,sbuf1p,sbuf2p, &
                    rbuf1r,rbuf2r,rbuf1t,rbuf2t,rbuf1p,rbuf2p)
!
      end if
!
! ****** Seam the second dimension.
!
      if (nproc_t.gt.1) then
!
        allocate (sbuf1r(n1r,n3r),rbuf1r(n1r,n3r), &
                  sbuf2r(n1r,n3r),rbuf2r(n1r,n3r), &
                  sbuf1t(n1t,n3t),rbuf1t(n1t,n3t), &
                  sbuf2t(n1t,n3t),rbuf2t(n1t,n3t), &
                  sbuf1p(n1p,n3p),rbuf1p(n1p,n3p), &
                  sbuf2p(n1p,n3p),rbuf2p(n1p,n3p))
!
        do concurrent (k=1:n3r,j=1:n1r)
          sbuf1r(j,k)=v%r(j,n2r-1,k)
          sbuf2r(j,k)=v%r(j,    2,k)
        enddo
!
        do concurrent (k=1:n3t,j=1:n1t)
          sbuf1t(j,k)=v%t(j,n2t-1,k)
          sbuf2t(j,k)=v%t(j,    2,k)
        enddo
!
        do concurrent (k=1:n3p,j=1:n1p)
          sbuf1p(j,k)=v%p(j,n2p-1,k)
          sbuf2p(j,k)=v%p(j,    2,k)
        enddo
!
        call MPI_Irecv (rbuf1r,lbuf2r,ntype_real,iproc_tm,tagr, &
                        comm_all,req(1),ierr)
        call MPI_Irecv (rbuf2r,lbuf2r,ntype_real,iproc_tp,tagr, &
                        comm_all,req(2),ierr)
        call MPI_Irecv (rbuf1t,lbuf2t,ntype_real,iproc_tm,tagt, &
                        comm_all,req(3),ierr)
        call MPI_Irecv (rbuf2t,lbuf2t,ntype_real,iproc_tp,tagt, &
                        comm_all,req(4),ierr)
        call MPI_Irecv (rbuf1p,lbuf2p,ntype_real,iproc_tm,tagp, &
                        comm_all,req(5),ierr)
        call MPI_Irecv (rbuf2p,lbuf2p,ntype_real,iproc_tp,tagp, &
                        comm_all,req(6),ierr)
!
! ****** Launch async sends.
!
        call MPI_Isend (sbuf1r,lbuf2r,ntype_real,iproc_tp,tagr, &
                        comm_all,req(7),ierr)
        call MPI_Isend (sbuf2r,lbuf2r,ntype_real,iproc_tm,tagr, &
                        comm_all,req(8),ierr)
        call MPI_Isend (sbuf1t,lbuf2t,ntype_real,iproc_tp,tagt, &
                        comm_all,req(9),ierr)
        call MPI_Isend (sbuf2t,lbuf2t,ntype_real,iproc_tm,tagt, &
                        comm_all,req(10),ierr)
        call MPI_Isend (sbuf1p,lbuf2p,ntype_real,iproc_tp,tagp, &
                        comm_all,req(11),ierr)
        call MPI_Isend (sbuf2p,lbuf2p,ntype_real,iproc_tm,tagp, &
                        comm_all,req(12),ierr)
!
! ****** Wait for all seams to complete.
!
        call MPI_Waitall (12,req,MPI_STATUSES_IGNORE,ierr)
!
! ****** Unload buffers.
!
        if (iproc_tm.ne.MPI_PROC_NULL) then
          do concurrent (k=1:n3r,j=1:n1r)
            v%r(j,1,k)=rbuf1r(j,k)
          enddo
          do concurrent (k=1:n3t,j=1:n1t)
            v%t(j,1,k)=rbuf1t(j,k)
          enddo
          do concurrent (k=1:n3p,j=1:n1p)
            v%p(j,1,k)=rbuf1p(j,k)
          enddo
        end if
!
        if (iproc_tp.ne.MPI_PROC_NULL) then
          do concurrent (k=1:n3r,j=1:n1r)
            v%r(j,n2r,k)=rbuf2r(j,k)
          enddo
          do concurrent (k=1:n3t,j=1:n1t)
            v%t(j,n2t,k)=rbuf2t(j,k)
          enddo
          do concurrent (k=1:n3p,j=1:n1p)
            v%p(j,n2p,k)=rbuf2p(j,k)
          enddo
        end if
!
        deallocate (sbuf1r,sbuf2r,sbuf1t,sbuf2t,sbuf1p,sbuf2p, &
                    rbuf1r,rbuf2r,rbuf1t,rbuf2t,rbuf1p,rbuf2p)
!
      end if
!
end subroutine
!#######################################################################
program psi_multigpu_test_code_stdpar
!
      use number_types
      use types
      use mpidefs
      use decomposition_params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nr,nt,np,nr_2,nt_2,np_2
      integer :: i,j,k,ic,ierr
      integer, parameter :: n_per_dim=250
      integer, parameter :: n_cycles=100
      type(vvec), target :: v
!
      call init_mpi
!
! ****** Set the resolution for each MPI rank
!
      nr = n_per_dim
      nt = n_per_dim
      np = n_per_dim
!
      nprocs(1)=-1
      nprocs(2)=-1
      nprocs(3)=-1
!
      if (iamp0) then
        print*, 'Grid size per dimension per rank: ',n_per_dim
        print*, 'Grid size per rank: ',n_per_dim**3
        print*,' '
      endif
!
! ****** Check/set the processor topology.
!
      call check_proc_topology
!
      call decompose_domain
!
      if (iamp0) then
        print*, 'Number of ranks in DIM1: ',nproc_r
        print*, 'Number of ranks in DIM2: ',nproc_t
        print*, 'Number of ranks in DIM3: ',nproc_p
        print*, 'Total number of ranks: ',nproc
        print*,' '
      endif
!
      write(*,'(5(a,i3),a)') 'World rank ',iprocw, &
       ' has cart rank ',iproc, &
       ' and shared rank ',iprocsh,' (i.e. GPU device number ', &
       iprocsh+1,' out of ',nprocsh,')'
!
      nr_2=MAX(nr/2,1)
      nt_2=MAX(nt/2,1)
      np_2=MAX(np/2,1)
!
      allocate(v%r(nr-1,nt,np))
      allocate(v%t(nr,nt-1,np))
      allocate(v%p(nr,nt,np-1))
!
      do concurrent (k=1:np,j=1:nt,i=1:nr-1)
        v%r(i,j,k)=iproc
      enddo
      do concurrent (k=1:np,j=1:nt-1,i=1:nr)
        v%t(i,j,k)=iproc
      enddo
      do concurrent (k=1:np-1,j=1:nt,i=1:nr)
        v%p(i,j,k)=iproc
      enddo
!
!     Loop multiple times.
!
      if (iamp0) print*,"Starting seam cycles..."
!
      do ic=1,n_cycles
        call seam_vvec (v)
        if (iamp0) print*,"cycle ",ic," completed"
      enddo
!
      if (iamp0) then
        print*,"Run completed!"
        print*, "vr(:,:,1):",    v%r(nr_2,nt_2,1)
        print*, "vr(:,:,2):",    v%r(nr_2,nt_2,2)
        print*, "vr(:,:,np-1):", v%r(nr_2,nt_2,np-1)
        print*, "vr(:,:,np):",   v%r(nr_2,nt_2,np)
        print*, "vt(:,:,1):",    v%t(nr_2,nt_2,1)
        print*, "vt(:,:,2):",    v%t(nr_2,nt_2,2)
        print*, "vt(:,:,np-1):", v%t(nr_2,nt_2,np-1)
        print*, "vt(:,:,np):",   v%t(nr_2,nt_2,np)
        print*, "vp(:,:,1):",    v%p(nr_2,nt_2,1)
        print*, "vp(:,:,2):",    v%p(nr_2,nt_2,2)
        print*, "vp(:,:,np-2):", v%p(nr_2,nt_2,np-2)
        print*, "vp(:,:,np-1):", v%p(nr_2,nt_2,np-1)
!
        print*, "vr(:,1,:):",    v%r(nr_2,1,np_2)
        print*, "vr(:,2,:):",    v%r(nr_2,2,np_2)
        print*, "vr(:,nt-1,:):", v%r(nr_2,nt-1,np_2)
        print*, "vr(:,nt,:):",   v%r(nr_2,nt,np_2)
        print*, "vt(:,1,:):",    v%t(nr_2,1,np_2)
        print*, "vt(:,2,:):",    v%t(nr_2,2,np_2)
        print*, "vt(:,nt-2,:):", v%t(nr_2,nt-2,np_2)
        print*, "vt(:,nt-1,:):", v%t(nr_2,nt-1,np_2)
        print*, "vp(:,1,:):",    v%p(nr_2,1,np_2)
        print*, "vp(:,2,:):",    v%p(nr_2,2,np_2)
        print*, "vp(:,nt-1,:):", v%p(nr_2,nt-1,np_2)
        print*, "vp(:,nt,:):",   v%p(nr_2,nt,np_2)
      end if
!
      deallocate(v%r,v%t,v%p)
!
      call MPI_Finalize (ierr)
!
end program
