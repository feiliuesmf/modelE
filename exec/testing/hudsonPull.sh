cd /discover/nobackup/ccruz/devel/modelE.clones/master
original_head=$(git rev-parse HEAD) || exit 
git pull || exit 
updated_head=$(git rev-parse HEAD) || exit 
if test "$updated_head" = "$original_head" ; then 
   echo "simplex repository has not changed" 
   exit 1 
else 
   echo "*** If you see this message then it means that simplex.giss.nasa.gov:/giss/gitrepo/modelE.git changed ***"
   exit 0 
fi
