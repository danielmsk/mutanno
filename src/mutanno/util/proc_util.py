#!/usr/bin/env python
import os


def cpu_num():
    import psutil
    return psutil.NUM_CPUS


def mem_usage():
    import psutil
    return psutil.virtual_memory()


def cpu_usage():
    import psutil
    return psutil.cpu_percent(interval=1, percpu=False)


def disk_usage(path="/"):
    import psutil
    return psutil.disk_usage(path)


def back_run(scmd):
    import file_util
    svrname = os.uname()[1]
    fname = "./tmp_" + svrname + ".sh"
    # file_util.fileSave (fname, scmd.strip() + " & \n", "w")
    file_util.fileSave(fname, scmd.strip() + " > /dev/null 2>&1", "w")
    run_cmd("chmod 755 " + fname)
    run_cmd(fname + " & ")


def run_cmd(scmd, flag=False):
    if flag:
        print(scmd)
    rst = os.popen(scmd)
    rst_cont = rst.read()
    # rst = subprocess.Popen([scmd], stdout=subprocess.PIPE, shell=True)
    # (rst_cont, error) = rst.communicate()
    # sys.stdout = old_stdout
    return rst_cont


def proc_list(v=""):
    import str_util
    cmd = "ps -ef"
    if v != "":
        cmd += " | grep " + v
    cont = run_cmd(cmd)
    arr = cont.strip().split("\n")
    plist = []
    for line in arr:
        p = {}
        arr2 = str_util.strip_space(line).split(" ")
        p['user'] = arr2[0]
        p['pid'] = arr2[1]
        p['pid2'] = arr2[2]
        p['id'] = arr2[3]
        p['start_time'] = arr2[4]
        p['job_id'] = arr2[5]
        p['run_time'] = arr2[6]
        p['proc'] = " ".join(arr2[7:])
        plist.append(p)

    return plist


def in_proc_list(proc_name):
    plist = proc_list()
    flag = False
    for p in plist:
        if proc_name in p['proc']:
            flag = True
            break
    return flag


def cnt_proc(proc_name):
    plist = proc_list()
    cnt = 0
    for p in plist:
        if proc_name in p['proc']:
            cnt += 1
            break
    return cnt


def get_bjob_list():
    cmd = "bjobs -aw"
    rst = run_cmd(cmd)
    bsubjoblist = []
    i = 0
    for line in rst.split("\n"):
        # print line
        if i > 0 and len(line) > 20:
            for kk in range(10):
                line = line.replace("  ", " ")
            arr = line.split(" ")
            bjob = {}
            bjob['id'] = arr[0]
            bjob['status'] = arr[2]
            bjob['sh'] = " ".join(arr[6:])
            bsubjoblist.append(bjob)
        i += 1
    return bsubjoblist


def get_bjob(tag, bjob_list=[]):
    if len(bjob_list) == 0:
        bjob_list = get_bjob_list()
    tjob = {'status': 'NA', 'id': ''}
    for bjob in bjob_list:
        if tag in bjob['sh']:
            tjob = bjob
            break
    return tjob


'''
# alternative approach
def run_cmd (u_cmd):
    p = Popen([u_cmd], shell=True, stdout=PIPE)
    p.wait()
    return p.stdout.read()
'''


class RemoteSvr():

    def __init__(self, server_name=""):
        if server_name == "":
            import os
            self.server_name = os.uname()[1]
        else:
            self.server_name = server_name

    def update_status(self, path=""):
        import web_util
        import db_util
        import file_util
        self.logpath = path
        d = {}
        d['cpu_usage'] = cpu_usage()
        d['server_name'] = self.server_name
        d['check_time'] = "NOW()"
        sql = db_util.insertSQL(d, "server_status")
        if path == "":
            web_util.dbgate(sql, "pcaso")
        else:
            import time
            file_util.fileSave(path + "/" + self.server_name + ".log", str(time.time()) + "\t" + sql+"\n", "w")

    def log_updater(self):
        import file_util
        import web_util

        backupkey = []
        for line in open(self.logpath + "/prev.log"):
            backupkey.append(line.strip().split("\t")[0])

        cont = ""
        for fname in file_util.walk(self.logpath):
            for line in open(fname):
                arr = line.strip().split("\t")
                if not arr[0] in backupkey:
                    web_util.dbgate(arr[1], "pcaso")
                    cont += line

        file_util.fileSave(self.logpath + "/prev.log", cont, "w")

    def job_run(self, jobfile):
        import file_util
        if file_util.is_exist(jobfile) and cpu_usage() < 90:
            i = 0
            cont = ""
            scmd = ""
            for line in open(jobfile):
                if i == 0:
                    scmd = line.strip()
                else:
                    cont += line
                i += 1
            file_util.fileSave(jobfile, cont, "w")
            if len(scmd) > 0:
                back_run(scmd)
                if i > 1:
                    self.job_run(jobfile)


if __name__ == "__main__":
    run_cmd("ls")
