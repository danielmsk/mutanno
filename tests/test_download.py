import sys
import shlex
sys.path.append('..')
from src import mutanno

# import mutanno
prog = "mutanno"

cmdlist = []
cmdlist.append("""
    download \
    -source_path /Users/pcaso/db/MUTANNO/TESTDATASOURCE \
    -source CLINVAR \
    -version latest \
    -refversion hg38
""")

def test_run_download():
    for cmd in cmdlist:
        cmd = prog + " " + cmd.strip()
        sys.argv = shlex.split(cmd)
        print(' '.join(sys.argv))
        # print(cmd)
        # print(shlex.quote(sys.argv))
        mutanno.cli()

if __name__ == "__main__":
    test_run_download()
