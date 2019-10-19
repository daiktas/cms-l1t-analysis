import unittest
from cmsl1t.batch.condor import _parse_condor_submit_output

NJOBS = 558
CLUSTER_ID = 695608
test_out = "Submitting job(s)" + "." * 557 + "{0} job(s) submitted to cluster {1}."
test_out = test_out.format(NJOBS, CLUSTER_ID)


class TestLSFBatch(unittest.TestCase):

    def test_parse_condor_submit_output(self):
        job_ids = list(_parse_condor_submit_output(test_out))
        self.assertEqual(len(job_ids), NJOBS)
        for i in range(NJOBS):
            job_id = str(CLUSTER_ID) + '.' + str(i)
            self.assertEqual(job_ids[i], job_id)
