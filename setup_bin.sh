#!/bin/bash

## creates sym links for all executable files in ~/src with -type f and -perm +111
## in the $HOME/bin directory

for exec in $(find $HOME/src/* -type f -perm +111)
do 
    if [ -x $exec ]
	then cp -rs $exec ~/bin/
    fi
done


## since creation is for all executable files, some are not necesary, and thus removed manually to have a relatively clean bin file
## perl, shell and python scripts were also removed.
cd $HOME/bin

rm *.fa *.sample *.bed Makefile *.rst *genome *.sh *.html *.css *.fna *.ebwt AUTH* *README* LICENSE* MANUAL* NEWS* TUTORIAL VERSION *.fq *.raw *.pl md5* test.t *.lua *.py *.so *.conf *.bat rtd.css_t *~ *.bt2 *.jar *.cfg *.txt *.svg *.h *.src *.bfc *.png *.pxd *.properties *.out *.c

cd ../

## creates sym links for all shell and perl scripts in $HOME/src
## in the $HOME/scripts directory
if [ -f $HOME/scripts ]
    then echo found scripts
    else mkdir scripts
fi

for exec in $(find $HOME/src/* -name "*.sh")
do 
    if [ -x $exec ]
	then cp -rs $exec $HOME/scripts/
    fi
done

for exec in $(find $HOME/src/* -name "*.pl")
do 
    if [ -x $exec ]
	then cp -rs $exec $HOME/scripts/
    fi
done

