#!/bin/bash
DOCCONF=fulldoc
mainDir=$PWD
docDir=doc/doxygen/${DOCCONF}/
if [ ! -d "${docDir}" ];then
    mkdir -p ${docDir}
fi

if [ ! -d "${docDir}/html" ];then
    cd ${docDir}
    git clone -b gh-pages  git@github.com:CMS-HGCal/TestBeam.git html
    cd ${mainDir}
else
    cd ${docDir}/html
    if [ "`git branch | grep -c gh-pages`" == "0" ];then
	cd ${mainDir}
	rm ${docDir}/html/ -Rf
	cd ${docDir}/
	git clone -b gh-pages  git@github.com:CMS-HGCal/TestBeam.git html
	cd ${mainDir}
    fi
fi


cd ${mainDir}
commit=`git log | head -1 | awk '{print $2}'`

cat ${DOCCONF} | sed "s/PROJECT_NUMBER.*/PROJECT_NUMBER = ${commit}/" | doxygen -


cd ${docDir}/html/
ls
git remote -v 
git branch 
git pull
git add *.html
git add *.css *.js
git add -f *.gif *.png
git add *.map
git add *.md5
git add search
git add search/*.html search/*.js
git commit -m "updated documentation" -a
git commit -m "updated documentation" -a
git push origin gh-pages:gh-pages
cd -

