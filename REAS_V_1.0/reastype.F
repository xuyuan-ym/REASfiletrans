      module reastype
      integer,parameter :: linemax=500 !namelist realated lenth
      integer,parameter :: namesize=32 !reasname realated lenth
      integer,parameter :: srcmax=50   !reasspc  volumn
      integer,parameter :: attlen=30   !netcdf attribute lenth
      type filespc
       character(len=namesize) ::file_spc
       character(len=linemax)  ::filename
      end type filespc
      type reasspc
       character(len=namesize) :: file_spc_r
       character(len=linemax)  :: reasname
       character(len=namesize) :: vname
       character(len=linemax)  :: open_file_name
       integer                 :: molecular_weight
      end type reasspc
      end module reastype
