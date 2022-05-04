#!/usr/bin/python
from socket import gethostname
from sys import argv
print("\n\n########\nHello world from ",gethostname())
print("\n")

try:
    print("testing Arguments...\t",argv[1],'...\tsuccess!')
except:
    print("Either no argument supplied or task failed", argv)


