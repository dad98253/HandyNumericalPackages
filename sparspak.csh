#!/bin/csh
#  sparspak.csh  17 June 2000
#
setenv DEBUG TRUE
set echo
#
cd ~/src/sparspak
mkdir temp
cd temp
rm *
f90split ../sparspak.f90
#
foreach FILE (`ls *.f90`)
#
  if ( $DEBUG == TRUE ) then
    f90 -c -C -fullwarn -g -static -u $FILE
    if ( $status != 0 ) exit
  else
    f90 -c $FILE
    if ( $status != 0 ) exit
  endif
  rm $FILE
end
#
ar qc libsparspak.a *.o
rm *.o
#
mv libsparspak.a ~/lib/$ARCH
cd ~/src/sparspak
rmdir temp
#
echo "A new version of sparspak has been created."
