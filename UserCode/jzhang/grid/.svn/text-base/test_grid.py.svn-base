import sys
import re
from cStringIO import StringIO
import jobsub_submit as js

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

def condorout_to_jobid(condorout):
    jobid = []
    for line in condorout:
        m = re.search(r"JobsubJobId.*: (.*)", line)
        if m:
            jobid.append(m.group(1))
    return  jobid

print(sys.argv)
sys.argv=['jobsub_submit','--group=coupp', '-N', '1', '--resource-provides=usage_model=DEDICATED',
    'file:///bluearc/storage/recon/devel/code/SBCcode/UserCode/jzhang/test.sh']
with Capturing() as output:
    print(sys.argv)
    # a = js.main(sys.argv)
    # print(a)
print(output)

# condor output is
condorout = ['/fife/local/scratch/uploads/coupp/jzhang/2017-07-05_132523.716114_3821', '', '/fife/local/scratch/uploads/coupp/jzhang/2017-07-05_132523.716114_3821/test.sh_20170705_132527_155349_0_1_.cmd', '', 'submitting....', '', 'Submitting job(s).', '', '1 job(s) submitted to cluster 18685985.', '', 'JobsubJobId of first job: 18685985.0@fifebatch1.fnal.gov', '', 'Use job id 18685985.0@fifebatch1.fnal.gov to retrieve output', '']
jobid = condorout_to_jobid(condorout)
print(jobid)

#  jobsub_q --jobid 18685985.0@fifebatch1.fnal.gov
