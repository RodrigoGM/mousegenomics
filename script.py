#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author : Leonor Palmeira

import sys

inFile=sys.argv[1]
outFile=sys.argv[2]

f=open(inFile,"r")
all=f.readlines()
f.close()

s=""
rand=1
for i in all:
    tmp=i.split()
    if tmp[2]!=".":
        s+=i
    else:
        tmp[2]="uk"+str(rand)
        s+="\t".join(tmp)+"\n"
        rand+=1

f=open(outFile,"w")
f.write(s)
f.close()
