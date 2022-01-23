      program reasv1
      
      use utils, only: mapper_file
      use utils, only: mapper_emiss
      use utils, only: getindices
      use utils, only:reasmap
      use utils, only:filemap
      use utils, only:assignfile
      use utils, only:assignreas
      use utils, only:linecnt,columncnt
      use utils, only:lolamaxmin
      use reastype
      include 'netcdf.inc'
!-------------------
!define the dummp variables in namelist.input 
!-------------------
      integer  :: stats  ! status for file iostat
      character(len=100) :: dir_reas     
      character(len=100) :: sub_dir
      character(len=100) :: prefix
      character(len=100) :: suffix
      character(len=namesize)  :: spc_map(srcmax)
      character(len=namesize)  :: emiss_map(srcmax)
      integer   :: year
      character(len=4) :: yearc
      integer :: nspec,nespec   
      integer :: reaspos
      integer :: reasposa(2)
!-------------------
!define the varibles for file write loop
!-------------------
      integer :: status,astatus  ! status for netcdf,astatus->gridmalloc
      integer :: ncid
      integer :: i,j
      integer :: ifread,fios,nlines,ncols,nlon,nlat
      integer :: latidx,lonidx
      real :: latmax,latmin,lonmax,lonmin
      character(len=100) ::filemeta(9)
      real, allocatable  ::reastxt(:,:)
      real, allocatable  ::reasnc(:,:,:)
      real, allocatable  ::lat(:),lon(:)
      integer :: time2010(12)
      integer :: time(12)
      integer :: dimlonid,dimlatid,dimtimeid
      integer :: dims(3),length(3),start(3)
      integer :: varid,vartimeid,varlatid,varlonid
!-------------------
      namelist /control/dir_reas,sub_dir,prefix,suffix,spc_map,
     &emiss_map,year
      sub_dir=' '
      spc_map(:)=' '
      emiss_map(:)=' '
      read(*,nml=control,iostat=stats)
      if(stats/=0) then
        write(*,*) 'can not read the namelist',stats
        stop
      end if      
      call mapper_file(srcmax,spc_map,nspec)  
      call mapper_emiss(srcmax,emiss_map,nespec)
      print *, 'the number of reas file input is :', nspec
      print *, 'the number of reas emiss is      :', nespec
      call assignfile(dir_reas,sub_dir,prefix,suffix,
     &spc_map,nspec)
      call assignreas(nespec,nspec,emiss_map)
      do ifread =1, nespec
       open(ifread,file=trim(reasmap(ifread)%open_file_name),
     &iostat=fios)
       if(fios/=0) then
        write(*,*) 'fail to read reasfile
     &:',reasmap(ifread)%open_file_name,'ios=',fios
        stop
       end if 
       do i =1,9
        read(ifread,'(A)'),filemeta(i)
       end do 
       call linecnt(ifread,nlines)
       call columncnt(ifread,ncols)
       allocate(reastxt(nlines,ncols),STAT=astatus)
       write(*,*) 'successful allocate array if =0 ->',astatus
       do i =1,9
        read(ifread,*)
       end do 
       do i =1,nlines
        read(ifread,'((2F8.2,12E14.7))')(reastxt(i,j),j=1,ncols)
       end do 
       rewind(ifread)
       call lolamaxmin(nlines,reastxt(:,1),lonmax,lonmin)
       call lolamaxmin(nlines,reastxt(:,2),latmax,latmin)
       write(*,*) '======> latmax:',latmax,'=======> latmin',latmin
       write(*,*) '======> lonmax:',lonmax,'=======> lonmin',lonmin
       nlon=(lonmax-lonmin)/0.25+1
       nlat=(latmax-latmin)/0.25+1
       write(*,*) 'londim is ',nlon,'latdim is ',nlat
       allocate(reasnc(nlon,nlat,12))
       allocate(lat(nlat))
       allocate(lon(nlon))
       do i =1,nlines
        do j =1,12
          latidx=(reastxt(i,2)-latmin)/0.25+1
          lonidx=(reastxt(i,1)-lonmin)/0.25+1
          reasnc(lonidx,latidx,j)=reastxt(i,j+2)
        end do 
       end do  
        
       do i =1,nlat
         lat(i)=latmin+0.25*(i-1)
       end do 
       do i =1,nlon
         lon(i)=lonmin+0.25*(i-1)
       end do 
       data time2010 / 40177, 40208, 40236, 40267, 40297, 40328,
     & 40358, 40389, 40420, 40450, 40481, 40511/
       do i =1,12
        if(year==2010) then
        time(i)=time2010(i)       
        else 
        time(i)=time2010(i)+(year-2010)*365
        end if 
       end do 
!!!!!!!!!!!ncfile write
       write(yearc,'(I4)') year     
       print *, yearc
       status=nf_create('REAS'//'_'//yearc//'_'//
     &trim(reasmap(ifread)%vname)
     &//'.nc',nf_netcdf4,ncid)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_redef(ncid)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_def_dim(ncid,'lat',nlat,dimlatid)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_def_dim(ncid,'lon',nlon,dimlonid)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_def_dim(ncid,'time',12,dimtimeid)
       if(status/=nf_noerr) call errhandle(status)

       dims(:)=(/dimlonid,dimlatid,dimtimeid/)
       status=nf_def_var(ncid,reasmap(ifread)%vname,nf_float,
     &3,dims,varid)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_def_var(ncid,'lat',nf_float,1,dimlatid,varlatid)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_def_var(ncid,'lon',nf_float,1,dimlonid,varlonid)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_def_var(ncid,'time',nf_int,1,dimtimeid,vartimeid)                        
       if(status/=nf_noerr) call errhandle(status)

       status=nf_put_att_text(ncid,varid,'Emissions',9,'this is gas ')
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_att_text(ncid,varid,'Unit',19,'mol km^-2 hr^-1')
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_att_int(ncid,varid,'molecular_weight',nf_int
     &,1,reasmap(ifread)%molecular_weight)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_att_text(ncid,varlonid,'long_name',18,'longitued')
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_att_text(ncid,varlonid,'comment',18,'centerofcell')
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_att_text(ncid,varlatid,'long_name',18,'latitued')
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_att_text(ncid,varlatid,'comment',18,'centerofcell')
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_att_text(ncid,vartimeid,'long_name',18,'2010:1-12')
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_att_text(ncid,vartimeid,'comment',18,'monthly data'
     &)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_att_text(ncid,vartimeid,'units',27,'days since
     &1900-01-01 00:00')
       start(:)=1      
       length(:)=(/nlon,nlat,12/)
       status=nf_put_vara_real(ncid,varid,start,length,reasnc)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_vara_real(ncid,varlatid,(/1/),nlat,lat)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_vara_real(ncid,varlonid,(/1/),nlon,lon)
       if(status/=nf_noerr) call errhandle(status)
       status=nf_put_vara_int(ncid,vartimeid,(/1/),12,time)
       if(status/=nf_noerr) call errhandle(status)

       deallocate(lat)
       deallocate(lon)
       deallocate(reasnc)
       deallocate(reastxt)
       
      end do 
      endprogram reasv1
      subroutine errhandle(status)
      integer :: status
      if(status.ne.nc_noerr) then
      print *, "error :",nf_strerror(status)
      stop
      end if 
      end subroutine errhandle
