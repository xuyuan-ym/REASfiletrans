      module utils
      use reastype
      implicit none
      public :: mapper_file
      public :: mapper_emiss
      public :: getindices
      public :: linecnt
      public :: columncnt
      public :: reasmap
      public :: filemap
      public :: lolamaxmin
      type(reasspc), allocatable :: reasmap(:)  
      type(filespc), allocatable :: filemap(:)
      contains
!------>file spc count(public)
      subroutine mapper_file(nmax,spc_map,nspec)
      integer,intent(in) ::nmax
      character(len=*),intent(in) ::spc_map(nmax)
      integer,intent(out) ::nspec
!--------------
! local use varibales
!--------------
      integer :: i,istat 
      nspec=0
      do i =1,nmax
        if( spc_map(i)==' ') then
          nspec= i-1
          exit 
        else 
          write(*,*) 'mapper: spc_map =' ,trim(spc_map(i))
        end if 
      end do 
      if (nspec== 0) then
          write(*,*) 'the spc_map read is not success plz check '
          stop 'mapper'
      end if 
      write(*,*) 'max capable spc is : ',nmax
    
      endsubroutine mapper_file
!------------->mapper of emission ,considering of different mutipulor
      subroutine mapper_emiss(nmax,emiss_map,nespec)
      integer, intent(in) :: nmax
      character(len=*), intent(in) :: emiss_map(nmax)
      integer, intent(out):: nespec
      integer :: i, istat
      nespec=0
      do i=1, nmax
        if(emiss_map(i)==' ')then
          nespec=i-1
          exit
        else
          write(*,*) 'mapper emiss_map: ',trim(emiss_map(i))
        end if
      end do 
      print *, 'the number of emiss  map is ',nespec
     
     
      end subroutine mapper_emiss
!--------> assign /DIR/file to file%filename & file spc to file%file_spc
      subroutine assignfile(dir_reas,sub_dir,prefix,suffix,
     &spc_map,nspec) 
        integer,intent(in)  :: nspec
        character(len=*),intent(in) :: dir_reas
        character(len=*),intent(in) :: sub_dir
        character(len=*),intent(in) :: prefix
        character(len=*),intent(in) :: suffix
        character(len=*),intent(in) :: spc_map(nspec)
        
        !local v
        integer :: i,istat
      allocate(filemap(nspec),stat=istat)
      if (istat/= 0)then
      write(*,*) 'there is something wrong with allocate filemap'
      stop
      end if
      do i =1, nspec
       filemap(i)%filename=trim(dir_reas)//trim(sub_dir)//trim(prefix)
     &//trim(spc_map(i))//trim(suffix)
       filemap(i)%file_spc=trim(spc_map(i))
      end do 
      endsubroutine assignfile
!--------> assign reasspc mapped from file
      subroutine assignreas(nespec,nspec,emiss_map)
       integer,intent(in) :: nespec
       integer,intent(in) :: nspec
       character(len=*),intent(in) :: emiss_map(srcmax)
       !-------- >local 
       integer :: istat,i,j
       integer :: matchcnt    ! the number of '->'
       integer :: matcharr(1) ! array store pos of '->'
       integer :: fmatchidx   ! the idx tell witch file to open for read
       integer :: mlwidx      ! the pos for moweight, it is
       integer :: mlwidxarr(1)! array store pos of '(' for mol weight
       character(len=30) :: tempstr
       allocate(reasmap(nespec),stat=istat)
      if (istat/= 0)then
      write(*,*) 'there is something wrong with allocate reasmap'
      stop
      end if
      do i =1, nespec
      reasmap(i)%reasname=trim(emiss_map(i))
      fmatchidx=0
       write(*,*) 'emiss name for reas is :',trim(reasmap(i)%reasname) 
       call getindices('->',trim(reasmap(i)%reasname),matchcnt,matcharr)
       reasmap(i)%file_spc_r=trim(reasmap(i)%reasname(:matcharr(1)-1))
       do j =1, nspec
        if(reasmap(i)%file_spc_r == filemap(j)%file_spc) then
         fmatchidx= j
         write(*,*) 'find the first appear file to read',fmatchidx
         exit
        end if 
       end do 
       if (fmatchidx == 0) then
         write(*,*) 'error of find file for ',reasmap(i)%file_spc_r
         stop
       end if 
       reasmap(i)%open_file_name=filemap(fmatchidx)%filename
       call getindices('(',reasmap(i)%reasname,mlwidx,mlwidxarr)
       print *, mlwidxarr(1)
       reasmap(i)%vname=trim(reasmap(i)%reasname(matcharr(1)+2:
     &mlwidxarr(1)-1))
       write(*,*) 'vname in file is ',trim(reasmap(i)%vname)
       read(reasmap(i)%reasname(mlwidxarr(1)+1:),'(I4)') 
     &reasmap(i)%molecular_weight
       write(*,*) 'molecular_weight is :',reasmap(i)%molecular_weight
      end do 
      endsubroutine assignreas
!--------> count the indices of specific symbol(local)
      subroutine getindices(matchstr, str, match_cnt, match_ndx)
      integer , intent(out) :: match_cnt
      integer , optional, intent(out) :: match_ndx(:)
      character(len=*) , intent(in)   :: matchstr
      character(len=*) , intent(in)   :: str
!---------------
! local variables
!---------------
      integer :: pos,ndx
      logical :: set_ndx
      if(present(match_ndx)) then
        set_ndx=.true.
      else
        set_ndx=.false.
      end if        
      write(*,*) 'passed substr:',matchstr,'  passed str: ',trim(str)
     &,'  if use array to keep pos: ',set_ndx
      ndx=0
      match_cnt=0
      
      do 
        pos=index(str(ndx+1:),matchstr)
! test line for checking pos//write (*,*) 'the c at pos is :',str(pos:pos)
        if(pos==0) then
         exit
        end if 
        ndx=ndx+pos
        match_cnt=match_cnt+1
        if(set_ndx) then
          match_ndx(match_cnt)=ndx
        end if 
      end do 
      end subroutine getindices
!-----------------line cnt 
      subroutine linecnt(ifread,nlines)
      integer,intent(in) :: ifread
      integer,intent(out):: nlines
      !--------local
      integer            :: istat
      nlines=0
      do 
       read(ifread,*,iostat=istat)
       if (istat>0) then
        write(*,*) 'fail of read file on halfway',istat
        stop
       else if(istat<0) then
        write(*,*) 'line read success total is:', nlines
        exit
       else
        nlines=nlines+1
       end if 
      end do 
      rewind(ifread)
      end subroutine linecnt
!------------------col cnt
      subroutine columncnt(ifread,ncols)
      integer,intent(in) :: ifread
      integer,intent(out):: ncols
      ! local
      integer :: i
      character(len=1000)  :: sample
      ncols=0
      do i =1,9
      read(ifread,*) 
      end do 
      read(ifread,'(A)') sample
      do i =1,1000
      if(sample(i-1:i-1)==' '.and.
     &sample(i+1:i+1).ne.' '.and.
     &sample(i:i).ne.' ') then
      ncols=ncols+1
      end if 
      end do 
      print *, ncols
      rewind(ifread)
      end subroutine columncnt
      subroutine lolamaxmin(n,array,rmax,rmin)
      integer,intent(in) :: n
      real,intent(in)    :: array(n,1)
      real,intent(out)   :: rmax
      real,intent(out)   :: rmin
      integer            :: i
      real               :: temp
      write(*,*), 'test pass with the 3rd number plz check',array(3,1)
      rmax=array(1,1)
      do i =2,n
      if(array(i,1)>rmax) then
      rmax=array(i,1)
      end if 
      end do 
      rmin=array(1,1)
      do i =2,n
      if(array(i,1)<rmin) then
      rmin=array(i,1)
      end if 
      end do  
      end subroutine lolamaxmin
      end module utils
