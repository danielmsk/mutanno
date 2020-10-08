#!/usr/bin/env python
# -*- coding: utf-8 -*-
#### search_test.py
#### made by Min-Seok Kwon
#### 2019-06-04 14:08:53
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path="/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path="/ms1/bin/python_lib"
else:
    sys_path="/home/mk446/bin/python_lib"
sys.path.append(sys_path)
import file_util
import proc_util
import time

def search_test(output):
    f = open(output, 'w')
    i = 0
    for line in open('search_term.txt'):
        line = line.strip()
        if len(line)>0:
            i += 1
            arr = line.split('\t')
            cmd = 'curl -X GET "localhost:9200/test/'+arr[0]+'" -H \'Content-Type: application/json\' -d\''+arr[1]+'\''
            t1 = time.time()
            proc_util.run_cmd(cmd)
            t2 = time.time()
            cont = str(i) + '\t' + str(t2-t1) + '\t' + cmd
            # print (cont)
            f.write(cont+'\n')
    f.close()

if __name__ == "__main__":
    search_test(sys.argv[1])
