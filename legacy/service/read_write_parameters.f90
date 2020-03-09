module rwpar_routes
implicit none

contains

!======================================================================
subroutine readpar(filename,comments,nofcom)

integer, parameter:: maxnumpar=256
character*(*) filename, comments(maxnumpar) ! filename - name of file with parameters of task
integer n, nofcom                           ! nofcom   - number of used lines in filename
                                            ! calculated in this subroutine
open (90,file=filename,status='old',err=190)
  do n=1,maxnumpar
     comments(n)=' '
     read (90,'(a)',end=19,err=191) comments(n)
  end do
  write(*,*)' Warning!'
  write(*,*)' Too many lines in: ',filename
  write(*,*)' Their number is greater than',maxnumpar
19 continue
 close(90)

 nofcom=n-1

 write(*,'(a,a)')' Input ocean task parameters from ', filename(1:len_trim (filename))
 do n=1,nofcom
  write(*,'(1x,a)') comments(n)(1:len_trim (comments(n)))
 end do
 write(*,*)

 return

190   write(*,*) ' error in open file: ',  & 
                filename(1:len_trim (filename))
      stop 1

191   write(*,*) ' error in reading file: ', &
                filename(1:len_trim (filename))
      stop 2

endsubroutine readpar

!======================================================================
subroutine writepar(filename,comments,nofcom)

 integer       nofcom         !number of used lines in filename (external)
 character*(*) filename       !name of file with parameters of task
 character*(*) comments(nofcom)
 integer n, l, lonc

 open (90,file=filename,err=190)
  lonc=len(comments(1))
   do n=1,nofcom
    l=lonc
    do while(comments(n)(l:l)==' '.and.l>1)
     l=l-1
    end do
    write (90,'(a)',err=191) comments(n)(1:l)
   end do
 close(90)

 write(*,*)' output of task parameters in ', filename(1:len_trim (filename))
 do n=1,nofcom
  write(*,'(1x,a)') comments(n)(1:78)
 end do
return

190   write(*,*) ' error in open file: ',  &
                 filename(1:len_trim (filename))
      stop 1

191   write(*,*) ' error in reading file: ',   &
                 filename(1:len_trim (filename))
      stop 2

endsubroutine writepar
!=============================================================================================
!      This function finds the first lexeme in the string. 
!      We call lexeme the character set separated by space or tab.
!      Memory for OUT_STRING must be allocated by caller.
!      Rusakov Noida. 

subroutine get_first_lexeme(in_string, out_string)
implicit none
character(*) in_string
character(*) out_string

!REMOVE LEADING AND TRIM BLANKS
 out_string = adjustl(in_string)
 out_string = trim   (out_string)
 out_string = out_string(1 : index(out_string, ' ')) 
endsubroutine get_first_lexeme

endmodule rwpar_routes