c#######################################################################
c
c ****** PSI_MULTIGPU_TEST_CODE
c
c     This code mimics the basic MPI+OpenACC tasks of PSI's
c     MAS Solar MHD code.
c
c     It sets up a Cartesian MPI topology, sets up a 3D grid
c     and tests a "seam" (point to point) MPI communication
c     using asyncronous Send and Recv calls, which use
c     OpenACC's 'host_data' to use CUDA-aware MPI.
c     This is used both on a allocatable sub-array, as well as
c     on a local buffer static array.
c
c     The code will automatically configure the topology based on
c     the number of MPI ranks it is called with.
c
c     The code uses an MPI shared communicator to set the GPU device
c     number using the 'set device' ACC pragma.
c
c     This code assumes you launch it with the number of MPI ranks per
c     node = number of GPUs per node
c     (e.g. for a dual socket node with 4 GPUs: mpiexec -npersocket 2)
c
c     If this code works on a multi-node GPU system,
c     than (most likely) so will MAS!
c
c     Author:  Ronald M. Caplan
c
c     Predictive Science Inc.
c     www.predsci.com
c     San Diego, California, USA 92121
c
c#######################################################################
c Copyright 2022 Predictive Science Inc.
c
c Licensed under the Apache License, Version 2.0 (the "License");
c you may not use this file except in compliance with the License.
c You may obtain a copy of the License at
c
c    http://www.apache.org/licenses/LICENSE-2.0
c
c Unless required by applicable law or agreed to in writing, software
c distributed under the License is distributed on an "AS IS" BASIS,
c WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
c implied.
c See the License for the specific language governing permissions and
c limitations under the License.
c#######################################################################
c
      module number_types
c
      use iso_fortran_env
c
      implicit none
c
      integer, parameter :: r_typ=REAL64
c
      end module
c#######################################################################
      module types
c
      use number_types
c
      type :: vvec
        real(r_typ), dimension(:,:,:), allocatable :: r !(nrm,nt,np)
        real(r_typ), dimension(:,:,:), allocatable :: t !(nr,ntm,np)
        real(r_typ), dimension(:,:,:), allocatable :: p !(nr,nt,npm)
      end type
c
      end module
c#######################################################################
      module mpidefs
c
      use mpi
c        
      implicit none
c
c ****** Total number of processors.
      integer :: nproc
c ****** Total number of processors per node.
      integer :: nprocsh
c ****** Processor rank of this process in communicator
c ****** MPI_COMM_WORLD.
      integer :: iprocw
c ****** Processor rank of this process in communicator
c ****** comm_shared.
      integer :: iprocsh
c ****** Flag to designate that this is the processor with
c ****** rank 0 in communicator MPI_COMM_WORLD.
      logical :: iamp0
c ****** Communicator over all processors in the Cartesian topology.
      integer :: comm_all
c ****** Processor rank of this process in communicator
c ****** COMM_ALL.
      integer :: iproc
c ****** Processor rank in communicator COMM_ALL for the
c ****** processor that has rank 0 in MPI_COMM_WORLD.
      integer :: iproc0
c ****** Communicators over all processors in the phi dimension.
      integer :: comm_phi
c ****** Communicator over all shared processors on the node.
      integer :: comm_shared
c ****** Communicators over all processors in the theta and phi
c ****** dimensions.
      integer :: comm_tp
c ****** Communicators over all processors in the r dimension.
      integer :: comm_r
c ****** Processor rank in communicator COMM_R of the processor
c ****** that contains the lower radial boundary r=R0.
      integer :: iproc_rb0
c ****** Processor coordinate indices of this process
c ****** in the Cartesian topology.
      integer :: iproc_r,iproc_t,iproc_p
c ****** Processor coordinate indices of the neighboring
c ****** processors in the Cartesian topology.
      integer :: iproc_rm,iproc_rp
      integer :: iproc_tm,iproc_tp
      integer :: iproc_pm,iproc_pp
c ****** Number of processors along r, theta, and phi.
      integer :: nproc_r,nproc_t,nproc_p
c ****** Number of processors in 2D theta-phi plane.
      integer :: nproc_tp
c ****** Processor coordinate indices in 2D theta-phi plane.
      integer :: iproc2d_tp
c ****** Number type for REALs to be used in MPI calls.
      integer :: ntype_real
      end module
c#######################################################################
      module decomposition_params
c
c-----------------------------------------------------------------------
c ****** Input parameters that define the domain decomposition
c ****** among processors.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Number of processors per dimension.
c
      integer, dimension(3) :: nprocs=(/-1,-1,-1/)
c
c ****** Number of mesh points per processor.
c
      integer, parameter :: mx_procs_per_dim=100
      integer, dimension(mx_procs_per_dim) :: mp_r=0
      integer, dimension(mx_procs_per_dim) :: mp_t=0
      integer, dimension(mx_procs_per_dim) :: mp_p=0
c
c ****** Mesh sizes for the "automatic" mesh decomposition.
c
      integer :: nr_auto
      integer :: nt_auto
      integer :: np_auto
c
      end module
c#######################################################################
      subroutine init_mpi
c
c-----------------------------------------------------------------------
c
c ****** Initialize MPI.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** MPI error return.
c
      integer :: ierr,tcheck
c
c-----------------------------------------------------------------------
c
      call MPI_Init_thread (MPI_THREAD_FUNNELED,tcheck,ierr)
c
c ****** Get the total number of processors.
c
      call MPI_Comm_size (MPI_COMM_WORLD,nproc,ierr)
c
c ****** Get the index (rank) of the local processor in
c ****** communicator MPI_COMM_WORLD in variable IPROCW.
c
      call MPI_Comm_rank (MPI_COMM_WORLD,iprocw,ierr)
c
c ****** Create a shared communicator for all ranks in the node.
c
      call MPI_Comm_split_type (MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,
     &                          MPI_INFO_NULL,comm_shared,ierr)
c
c ****** Get the total number of processors in node.
c
      call MPI_Comm_size (comm_shared,nprocsh,ierr)
c
c ****** Get the index (rank) of the local processor in the local node.
c
      call MPI_Comm_rank (comm_shared,iprocsh,ierr)
c
c ****** Set the flag to designate whether this processor
c ****** has rank 0 in communicator MPI_COMM_WORLD.
c
      if (iprocw.eq.0) then
        iamp0=.true.
      else
        iamp0=.false.
      end if
c
      ntype_real=MPI_REAL8
c
c ****** Set GPU device number for current rank.
c
!$acc set device_num(iprocsh)
c
      end subroutine
c#######################################################################
      subroutine check_proc_topology
c
c-----------------------------------------------------------------------
c
c ****** Check/set the requested processor topology.
c
c-----------------------------------------------------------------------
c
      use mpidefs
      use decomposition_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: i,nreq
c
c-----------------------------------------------------------------------
c
c ****** Set the optimal values of the topology for unset dimensions.
c
      call set_proc_topology
c
c ****** Check that the number of processors available
c ****** matches the number requested.
c
      nreq=nprocs(1)*nprocs(2)*nprocs(3)
c
      if (nreq.ne.nproc) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_PROC_TOPOLOGY:'
          write (*,*) '### Processor topology specification error.'
          write (*,*) 'The number of processors requested does not'//
     &                ' equal the number available.'
          write (*,*) 'Number of processors requested = ',nreq
          write (*,*) 'Number of processors available = ',nproc
        end if
      end if
c
      end subroutine
c#######################################################################
      subroutine set_proc_topology
c
c-----------------------------------------------------------------------
c
c ****** Set the optimal values of the MPI rank topology
c ****** in dimensions not set by user.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
      use decomposition_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1.0_r_typ
      real(r_typ), parameter :: zero=0.0_r_typ
      real(r_typ), parameter :: bigval=HUGE(1.0_r_typ)
c
c-----------------------------------------------------------------------
c
      integer, dimension(:), allocatable :: factors
      integer, dimension(:,:), allocatable :: rank_factors
      real(r_typ), dimension(:,:), allocatable :: nperrank
      real(r_typ), dimension(:), allocatable :: penalty
c
      integer :: i,j,k,fr,ft,fp,num_fac,num_rank_fac,best_idx
      real(r_typ) :: a12,a13,a23
c
c-----------------------------------------------------------------------
c
c ****** Extract nproc values.  A value of -1 indicates the dimension
c ****** should be autoset.
c
      nproc_r=nprocs(1)
      nproc_t=nprocs(2)
      nproc_p=nprocs(3)
c
c ****** If no dimensions are to be autoset, return.
c
      if(nproc_r.ne.-1.and.nproc_t.ne.-1.and.nproc_p.ne.-1) return
c
c ****** Get all factors of nproc and store them in factors array.
c
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
c
c ****** Set penalty function parameters and any fixed dimensions
c ****** based on which dimensions are to be autoset.
c
      a12=one
      a13=one
      a23=one
c
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
c
c ****** Loop over all combinations of factors and save those that
c ****** yield the correct number of MPI ranks into rank_factors array.
c
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
c
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
c
      allocate(rank_factors(num_rank_fac,3))
      allocate(nperrank(num_rank_fac,3))
      allocate(penalty(num_rank_fac))
c
      rank_factors(:,:)=-1
      penalty(:)=bigval
c
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
c
c ****** Get number of grid points per rank for each dimension.
c
c      nperrank(:,1)=real(nr_g)/rank_factors(:,1)
c      nperrank(:,2)=real(nt_g)/rank_factors(:,2)
c      nperrank(:,3)=real(np_g)/rank_factors(:,3)
c
c ****** Compute penalty function.
c
      penalty(:)=a12*(rank_factors(:,1)-rank_factors(:,2))**2
     &          +a23*(rank_factors(:,2)-rank_factors(:,3))**2
     &          +a13*(rank_factors(:,3)-rank_factors(:,1))**2
c
c ****** Eliminate any choices that yield less than a minimum number
c ****** of grid points per rank.
c
c      do i=1,num_rank_fac
c        if (nperrank(i,1).lt.4) penalty(i)=bigval
c        if (nperrank(i,2).lt.4) penalty(i)=bigval
c        if (nperrank(i,3).lt.3) penalty(i)=bigval
c      enddo
c
c ****** Find optimal topology.
c
      best_idx=MINLOC(penalty,1)
c
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
c          write (*,*) 'NR = ',nr_g
c          write (*,*) 'NT = ',nt_g
c          write (*,*) 'NP = ',np_g
        end if
      end if
c
c ****** Set optimal topology.
c
      nprocs(1)=rank_factors(best_idx,1)
      nprocs(2)=rank_factors(best_idx,2)
      nprocs(3)=rank_factors(best_idx,3)
c
      deallocate(factors)
      deallocate(rank_factors)
      deallocate(nperrank)
      deallocate(penalty)
c
      end subroutine
c#######################################################################
      subroutine decompose_domain
c
c-----------------------------------------------------------------------
c
c ****** Decompose the domain into a Cartesian MPI topology.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
      use decomposition_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      integer, parameter :: ndim=3
      integer, dimension(ndim) :: coords
      logical, dimension(ndim) :: periodic
      logical :: reorder
      logical, dimension(ndim) :: keep_dim
c
c-----------------------------------------------------------------------
c
c ****** Create a communicator over all processors, COMM_ALL,
c ****** that has a Cartesian topology.
c
c ****** Specify the periodicity of the coordinate system.
c
      periodic(1)=.false.
      periodic(2)=.false.
      periodic(3)=.true.
c
c ****** Allow re-ordering in the Cartesian topology.
c
      reorder=.true.
c
      call MPI_Cart_create (MPI_COMM_WORLD,ndim,nprocs,
     &                      periodic,reorder,comm_all,ierr)
c
c ****** Get the index (rank) of the local processor in
c ****** communicator COMM_ALL in variable IPROC.
c
c ****** IMPORTANT NOTE:
c ****** If re-odering was allowed in the Cartesian topology
c ****** creation (above), then the rank of the local processor
c ****** in communicator COMM_ALL may be different from its rank
c ****** in communicator MPI_COMM_WORLD.
c
      call MPI_Comm_rank (comm_all,iproc,ierr)
c
c ****** Set the processor rank IPROC0 in communicator COMM_ALL
c ****** for the processor that has rank 0 in MPI_COMM_WORLD.
c ****** This value is broadcast to all the processors.
c
      if (iamp0) then
        iproc0=iproc
      end if
      call MPI_Bcast (iproc0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
c ****** Get the coordinate indices of this processor in the
c ****** Cartesian MPI topology.
c
      call MPI_Cart_coords (comm_all,iproc,ndim,coords,ierr)
c
      iproc_r=coords(1)
      iproc_t=coords(2)
      iproc_p=coords(3)
c
      nproc_r=nprocs(1)
      nproc_t=nprocs(2)
      nproc_p=nprocs(3)
c
      nproc_tp=nproc_t*nproc_p
c
c ****** Get the rank of the neighboring processors in the
c ****** Cartesian MPI topology.
c
      call MPI_Cart_shift (comm_all,0,1,iproc_rm,iproc_rp,ierr)
      call MPI_Cart_shift (comm_all,1,1,iproc_tm,iproc_tp,ierr)
      call MPI_Cart_shift (comm_all,2,1,iproc_pm,iproc_pp,ierr)
c
c ****** Create communicators for operations involving all
c ****** processors in the phi dimension.  These communicators
c ****** are stored in COMM_PHI (and generally represent different
c ****** communicators on different processors).
c
      keep_dim(1)=.false.
      keep_dim(2)=.false.
      keep_dim(3)=.true.
c
      call MPI_Cart_sub (comm_all,keep_dim,comm_phi,ierr)
c
c ****** Create communicators for operations involving
c ****** all processors in the theta and phi dimensions.
c ****** These communicators are stored in COMM_TP
c ****** (and generally represent different communicators on
c ****** different processors).
c ****** These communicators are used for operations that
c ****** involve radial planes.
c
      keep_dim(1)=.false.
      keep_dim(2)=.true.
      keep_dim(3)=.true.
c
      call MPI_Cart_sub (comm_all,keep_dim,comm_tp,ierr)
c
c ****** Get rank in the theta-phi communicator.
c ****** This is used for 2D IO.
c
      call MPI_Comm_rank (comm_tp,iproc2d_tp,ierr)
c
c ****** Create communicators for operations involving
c ****** all processors in the r dimension.
c ****** These communicators are stored in COMM_R
c ****** (and generally represent different communicators on
c ****** different processors).
c
      keep_dim(1)=.true.
      keep_dim(2)=.false.
      keep_dim(3)=.false.
c
      call MPI_Cart_sub (comm_all,keep_dim,comm_r,ierr)
c
      return
      end
c#######################################################################
      subroutine seam_vvec (v)
c
c-----------------------------------------------------------------------
c
c ****** Seam the boundary points of a v vector between adjacent
c ****** processors.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types, ONLY : vvec
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(vvec) :: v
c
c-----------------------------------------------------------------------
c
      real(r_typ),dimension(:,:),allocatable :: sbuf1r,rbuf1r
      real(r_typ),dimension(:,:),allocatable :: sbuf2r,rbuf2r
      real(r_typ),dimension(:,:),allocatable :: sbuf1t,rbuf1t
      real(r_typ),dimension(:,:),allocatable :: sbuf2t,rbuf2t
      real(r_typ),dimension(:,:),allocatable :: sbuf1p,rbuf1p
      real(r_typ),dimension(:,:),allocatable :: sbuf2p,rbuf2p
c
c-----------------------------------------------------------------------
c
c ****** MPI error return.
c
      integer :: ierr
c
c ****** MPI tags for MPI_ISEND and MPI_IRECV.
c
      integer :: tagr=0
      integer :: tagt=1
      integer :: tagp=2
c
c-----------------------------------------------------------------------
c
      integer :: lbuf3r,lbuf3t,lbuf3p
      integer :: lbuf1r,lbuf1t,lbuf1p
      integer :: lbuf2r,lbuf2t,lbuf2p
      integer :: n1r,n2r,n3r,n1t,n2t,n3t,n1p,n2p,n3p
      integer :: req(12)
      integer :: i,j,k
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c
c ****** Get the dimensions of the arrays and buffer sizes:
c
      n1r=size(v%r,1);   n2r=size(v%r,2);   n3r=size(v%r,3)
      n1t=size(v%t,1);   n2t=size(v%t,2);   n3t=size(v%t,3)
      n1p=size(v%p,1);   n2p=size(v%p,2);   n3p=size(v%p,3)
c
      lbuf3r=n1r*n2r;    lbuf3t=n1t*n2t;    lbuf3p=n1p*n2p
      lbuf1r=n2r*n3r;    lbuf1t=n2t*n3t;    lbuf1p=n2p*n3p
      lbuf2r=n1r*n3r;    lbuf2t=n1t*n3t;    lbuf2p=n1p*n3p
c
c ****** Seam the third (periodic) dimension. Since seam data
c        is stride-1 in this case, no buffers are needed.
c
c ****** Launch async receives.
c
!$acc host_data use_device(v%r,v%t,v%p)
      call MPI_Irecv (v%r(:,:,  1),lbuf3r,ntype_real,iproc_pm,tagr,
     &                comm_all,req(1),ierr)
      call MPI_Irecv (v%r(:,:,n3r),lbuf3r,ntype_real,iproc_pp,tagr,
     &                comm_all,req(2),ierr)
      call MPI_Irecv (v%t(:,:,  1),lbuf3t,ntype_real,iproc_pm,tagt,
     &                comm_all,req(3),ierr)
      call MPI_Irecv (v%t(:,:,n3t),lbuf3t,ntype_real,iproc_pp,tagt,
     &                comm_all,req(4),ierr)
      call MPI_Irecv (v%p(:,:,  1),lbuf3p,ntype_real,iproc_pm,tagp,
     &                comm_all,req(5),ierr)
      call MPI_Irecv (v%p(:,:,n3p),lbuf3p,ntype_real,iproc_pp,tagp,
     &                comm_all,req(6),ierr)
c
c ****** Launch async sends.
c
      call MPI_Isend (v%r(:,:,n3r-1),lbuf3r,ntype_real,iproc_pp,tagr,
     &                comm_all,req(7),ierr)
      call MPI_Isend (v%r(:,:,    2),lbuf3r,ntype_real,iproc_pm,tagr,
     &                comm_all,req(8),ierr)
      call MPI_Isend (v%t(:,:,n3t-1),lbuf3t,ntype_real,iproc_pp,tagt,
     &                comm_all,req(9),ierr)
      call MPI_Isend (v%t(:,:,    2),lbuf3t,ntype_real,iproc_pm,tagt,
     &                comm_all,req(10),ierr)
      call MPI_Isend (v%p(:,:,n3p-1),lbuf3p,ntype_real,iproc_pp,tagp,
     &                comm_all,req(11),ierr)
      call MPI_Isend (v%p(:,:,    2),lbuf3p,ntype_real,iproc_pm,tagp,
     &                comm_all,req(12),ierr)
c
c ****** Wait for all seams to complete.
c
      call MPI_Waitall (12,req,MPI_STATUSES_IGNORE,ierr)
!$acc end host_data
c
c ****** Seam the first dimension.
c
      if (nproc_r.gt.1) then
c
c ****** Load buffers.
c
        allocate (sbuf1r(n2r,n3r),rbuf1r(n2r,n3r),
     &            sbuf2r(n2r,n3r),rbuf2r(n2r,n3r),
     &            sbuf1t(n2t,n3t),rbuf1t(n2t,n3t),
     &            sbuf2t(n2t,n3t),rbuf2t(n2t,n3t),
     &            sbuf1p(n2p,n3p),rbuf1p(n2p,n3p),
     &            sbuf2p(n2p,n3p),rbuf2p(n2p,n3p))
!$acc enter data create(sbuf1r,sbuf2r,sbuf1t,sbuf2t,sbuf1p,sbuf2p,
!$acc&                  rbuf1r,rbuf2r,rbuf1t,rbuf2t,rbuf1p,rbuf2p)
c
!$acc parallel default(present)
!$acc loop collapse(2)
        do k=1,n3r
          do j=1,n2r
            sbuf1r(j,k)=v%r(n1r-1,j,k)
            sbuf2r(j,k)=v%r(    2,j,k)
          enddo
        enddo
c
!$acc loop collapse(2)
        do k=1,n3t
          do j=1,n2t
            sbuf1t(j,k)=v%t(n1t-1,j,k)
            sbuf2t(j,k)=v%t(    2,j,k)
          enddo
        enddo
c
!$acc loop collapse(2)
        do k=1,n3p
          do j=1,n2p
            sbuf1p(j,k)=v%p(n1p-1,j,k)
            sbuf2p(j,k)=v%p(    2,j,k)
          enddo
        enddo
!$acc end parallel
c
!$acc host_data use_device(sbuf1r,sbuf2r,sbuf1t,
!$acc&                     sbuf2t,sbuf1p,sbuf2p,
!$acc&                     rbuf1r,rbuf2r,rbuf1t,
!$acc&                     rbuf2t,rbuf1p,rbuf2p)
        call MPI_Irecv (rbuf1r,lbuf1r,ntype_real,iproc_rm,tagr,
     &                  comm_all,req(1),ierr)
        call MPI_Irecv (rbuf2r,lbuf1r,ntype_real,iproc_rp,tagr,
     &                  comm_all,req(2),ierr)
        call MPI_Irecv (rbuf1t,lbuf1t,ntype_real,iproc_rm,tagt,
     &                  comm_all,req(3),ierr)
        call MPI_Irecv (rbuf2t,lbuf1t,ntype_real,iproc_rp,tagt,
     &                  comm_all,req(4),ierr)
        call MPI_Irecv (rbuf1p,lbuf1p,ntype_real,iproc_rm,tagp,
     &                  comm_all,req(5),ierr)
        call MPI_Irecv (rbuf2p,lbuf1p,ntype_real,iproc_rp,tagp,
     &                  comm_all,req(6),ierr)
c
c ****** Launch async sends.
c
        call MPI_Isend (sbuf1r,lbuf1r,ntype_real,iproc_rp,tagr,
     &                  comm_all,req(7),ierr)
        call MPI_Isend (sbuf2r,lbuf1r,ntype_real,iproc_rm,tagr,
     &                  comm_all,req(8),ierr)
        call MPI_Isend (sbuf1t,lbuf1t,ntype_real,iproc_rp,tagt,
     &                  comm_all,req(9),ierr)
        call MPI_Isend (sbuf2t,lbuf1t,ntype_real,iproc_rm,tagt,
     &                  comm_all,req(10),ierr)
        call MPI_Isend (sbuf1p,lbuf1p,ntype_real,iproc_rp,tagp,
     &                  comm_all,req(11),ierr)
        call MPI_Isend (sbuf2p,lbuf1p,ntype_real,iproc_rm,tagp,
     &                  comm_all,req(12),ierr)
c
c ****** Wait for all seams to complete.
c
        call MPI_Waitall (12,req,MPI_STATUSES_IGNORE,ierr)
!$acc end host_data
c
c ****** Unload buffers.
c
!$acc parallel default(present)
        if (iproc_rm.ne.MPI_PROC_NULL) then
!$acc loop collapse(2)
          do k=1,n3r
            do j=1,n2r
              v%r(1,j,k)=rbuf1r(j,k)
            enddo
          enddo
!$acc loop collapse(2)
          do k=1,n3t
            do j=1,n2t
              v%t(1,j,k)=rbuf1t(j,k)
            enddo
          enddo
!$acc loop collapse(2)
          do k=1,n3p
            do j=1,n2p
              v%p(1,j,k)=rbuf1p(j,k)
            enddo
          enddo
        end if
c
        if (iproc_rp.ne.MPI_PROC_NULL) then
!$acc loop collapse(2)
          do k=1,n3r
            do j=1,n2r
              v%r(n1r,j,k)=rbuf2r(j,k)
            enddo
          enddo
!$acc loop collapse(2)
          do k=1,n3t
            do j=1,n2t
              v%t(n1t,j,k)=rbuf2t(j,k)
            enddo
          enddo
!$acc loop collapse(2)
          do k=1,n3p
            do j=1,n2p
              v%p(n1p,j,k)=rbuf2p(j,k)
            enddo
          enddo
        end if
!$acc end parallel
c
!$acc exit data delete(sbuf1r,sbuf2r,sbuf1t,sbuf2t,sbuf1p,sbuf2p,
!$acc&                 rbuf1r,rbuf2r,rbuf1t,rbuf2t,rbuf1p,rbuf2p)
        deallocate (sbuf1r,sbuf2r,sbuf1t,sbuf2t,sbuf1p,sbuf2p,
     &              rbuf1r,rbuf2r,rbuf1t,rbuf2t,rbuf1p,rbuf2p)
c
      end if
c
c ****** Seam the second dimension.
c
      if (nproc_t.gt.1) then
c
        allocate (sbuf1r(n1r,n3r),rbuf1r(n1r,n3r),
     &            sbuf2r(n1r,n3r),rbuf2r(n1r,n3r),
     &            sbuf1t(n1t,n3t),rbuf1t(n1t,n3t),
     &            sbuf2t(n1t,n3t),rbuf2t(n1t,n3t),
     &            sbuf1p(n1p,n3p),rbuf1p(n1p,n3p),
     &            sbuf2p(n1p,n3p),rbuf2p(n1p,n3p))
!$acc enter data create(sbuf1r,sbuf2r,sbuf1t,sbuf2t,sbuf1p,sbuf2p,
!$acc&                  rbuf1r,rbuf2r,rbuf1t,rbuf2t,rbuf1p,rbuf2p)
c
!$acc parallel default(present)
!$acc loop collapse(2)
        do k=1,n3r
          do j=1,n1r
            sbuf1r(j,k)=v%r(j,n2r-1,k)
            sbuf2r(j,k)=v%r(j,    2,k)
          enddo
        enddo
c
!$acc loop collapse(2)
        do k=1,n3t
          do j=1,n1t
            sbuf1t(j,k)=v%t(j,n2t-1,k)
            sbuf2t(j,k)=v%t(j,    2,k)
          enddo
        enddo
c
!$acc loop collapse(2)
        do k=1,n3p
          do j=1,n1p
            sbuf1p(j,k)=v%p(j,n2p-1,k)
            sbuf2p(j,k)=v%p(j,    2,k)
          enddo
        enddo
!$acc end parallel
c
!$acc host_data use_device(sbuf1r,sbuf2r,sbuf1t,
!$acc&                     sbuf2t,sbuf1p,sbuf2p,
!$acc&                     rbuf1r,rbuf2r,rbuf1t,
!$acc&                     rbuf2t,rbuf1p,rbuf2p)
        call MPI_Irecv (rbuf1r,lbuf2r,ntype_real,iproc_tm,tagr,
     &                  comm_all,req(1),ierr)
        call MPI_Irecv (rbuf2r,lbuf2r,ntype_real,iproc_tp,tagr,
     &                  comm_all,req(2),ierr)
        call MPI_Irecv (rbuf1t,lbuf2t,ntype_real,iproc_tm,tagt,
     &                  comm_all,req(3),ierr)
        call MPI_Irecv (rbuf2t,lbuf2t,ntype_real,iproc_tp,tagt,
     &                  comm_all,req(4),ierr)
        call MPI_Irecv (rbuf1p,lbuf2p,ntype_real,iproc_tm,tagp,
     &                  comm_all,req(5),ierr)
        call MPI_Irecv (rbuf2p,lbuf2p,ntype_real,iproc_tp,tagp,
     &                  comm_all,req(6),ierr)
c
c ****** Launch async sends.
c
        call MPI_Isend (sbuf1r,lbuf2r,ntype_real,iproc_tp,tagr,
     &                  comm_all,req(7),ierr)
        call MPI_Isend (sbuf2r,lbuf2r,ntype_real,iproc_tm,tagr,
     &                  comm_all,req(8),ierr)
        call MPI_Isend (sbuf1t,lbuf2t,ntype_real,iproc_tp,tagt,
     &                  comm_all,req(9),ierr)
        call MPI_Isend (sbuf2t,lbuf2t,ntype_real,iproc_tm,tagt,
     &                  comm_all,req(10),ierr)
        call MPI_Isend (sbuf1p,lbuf2p,ntype_real,iproc_tp,tagp,
     &                  comm_all,req(11),ierr)
        call MPI_Isend (sbuf2p,lbuf2p,ntype_real,iproc_tm,tagp,
     &                  comm_all,req(12),ierr)
c
c ****** Wait for all seams to complete.
c
        call MPI_Waitall (12,req,MPI_STATUSES_IGNORE,ierr)
!$acc end host_data
c
c ****** Unload buffers.
c
!$acc parallel default(present)
        if (iproc_tm.ne.MPI_PROC_NULL) then
!$acc loop collapse(2)
          do k=1,n3r
            do j=1,n1r
              v%r(j,1,k)=rbuf1r(j,k)
            enddo
          enddo
!$acc loop collapse(2)
          do k=1,n3t
            do j=1,n1t
              v%t(j,1,k)=rbuf1t(j,k)
            enddo
          enddo
!$acc loop collapse(2)
          do k=1,n3p
            do j=1,n1p
              v%p(j,1,k)=rbuf1p(j,k)
            enddo
          enddo
        end if
c
        if (iproc_tp.ne.MPI_PROC_NULL) then
!$acc loop collapse(2)
          do k=1,n3r
            do j=1,n1r
              v%r(j,n2r,k)=rbuf2r(j,k)
            enddo
          enddo
!$acc loop collapse(2)
          do k=1,n3t
            do j=1,n1t
              v%t(j,n2t,k)=rbuf2t(j,k)
            enddo
          enddo
!$acc loop collapse(2)
          do k=1,n3p
            do j=1,n1p
              v%p(j,n2p,k)=rbuf2p(j,k)
            enddo
          enddo
        end if
!$acc end parallel
c
!$acc exit data delete(sbuf1r,sbuf2r,sbuf1t,sbuf2t,sbuf1p,sbuf2p,
!$acc&                 rbuf1r,rbuf2r,rbuf1t,rbuf2t,rbuf1p,rbuf2p)
        deallocate (sbuf1r,sbuf2r,sbuf1t,sbuf2t,sbuf1p,sbuf2p,
     &              rbuf1r,rbuf2r,rbuf1t,rbuf2t,rbuf1p,rbuf2p)
c
      end if
c
      end subroutine
c#######################################################################
      program psi_multigpu_test_code
c
      use number_types
      use types
      use mpidefs
      use decomposition_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nr,nt,np,nr_2,nt_2,np_2
      integer :: i,j,k,ic,ierr
      integer, parameter :: n_per_dim=250
      integer, parameter :: n_cycles=1000
      type(vvec), target :: v
c
      call init_mpi
c
c ****** Set the resolution for each MPI rank
c
      nr = n_per_dim
      nt = n_per_dim
      np = n_per_dim
c
      nprocs(1)=-1
      nprocs(2)=-1
      nprocs(3)=-1
c
      if (iamp0) then
        print*, 'Grid size per dimension per rank: ',n_per_dim
        print*, 'Grid size per rank: ',n_per_dim**3
        print*,' '
      endif
c
c ****** Check/set the processor topology.
c
      call check_proc_topology
c
      call decompose_domain
c
      if (iamp0) then
        print*, 'Number of ranks in DIM1: ',nproc_r
        print*, 'Number of ranks in DIM2: ',nproc_t
        print*, 'Number of ranks in DIM3: ',nproc_p
        print*, 'Total number of ranks: ',nproc
        print*,' '
      endif
c
      write(*,'(5(a,i3),a)'),'World rank ',iprocw,
     &' has cart rank ',iproc,
     &' and shared rank ',iprocsh,' (i.e. GPU device number ',
     & iprocsh+1,' out of ',nprocsh,')'
c
      nr_2=MAX(nr/2,1)
      nt_2=MAX(nt/2,1)
      np_2=MAX(np/2,1)
c
      allocate(v%r(nr-1,nt,np))
      allocate(v%t(nr,nt-1,np))
      allocate(v%p(nr,nt,np-1))
!$acc enter data create(v,v%r,v%t,v%p)
c
!$acc parallel default(present)
!$acc loop collapse(3)
      do k=1,np
        do j=1,nt
          do i=1,nr-1
            v%r(i,j,k)=iproc
          enddo
        enddo
      enddo
!$acc loop collapse(3)
      do k=1,np
        do j=1,nt-1
          do i=1,nr
            v%t(i,j,k)=iproc
          enddo
        enddo
      enddo
!$acc loop collapse(3)
      do k=1,np-1
        do j=1,nt
          do i=1,nr
            v%p(i,j,k)=iproc
          enddo
        enddo
      enddo
!$acc end parallel
c
c     Loop multiple times.
c
      if (iamp0) print*,"Starting seam cycles..."
c
      do ic=1,n_cycles
        call seam_vvec (v)
        if (iamp0) print*,"cycle ",ic," completed"
      enddo
c
      if (iamp0) then
        print*,"Run completed!"
!$acc update self(v%r,v%t,v%p)
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
c
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
c
!$acc exit data delete(v%r,v%t,v%p,v)
      deallocate(v%r,v%t,v%p)
c
      call MPI_Finalize (ierr)
c
      end program
