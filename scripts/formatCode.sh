#!/bin/sh
command="astyle -A3 -t -p  -n" 
for file in $@
do
    case $file in
	*.c | *.cc) $command $file;;
	*.h | *.hh) $command $file;;
	*.cpp)      $command $file;;
    esac
#  /home/robgon/astyle/build/gcc/bin/astyle --indent=force-tab --pad-oper --pad-paren --delete-empty-lines --suffix=none --indent-namespaces --indent-col1-comments -n ${FICHEROS}_copy
done

