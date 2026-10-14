#!/usr/bin/env python3
"""benchmark_min_dim.check_backends (review 20261013-180049 G3): a backend=cpu run may profile cpu and cpu-batch (A30);
a cuda run only cuda; anything else stops the measurement. No binary, no GPU: python3 test_backend_check.py"""
import os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_min_dim import check_backends


class BackendCheck(unittest.TestCase):
    def test_accepted(self):
        check_backends('d', 'cpu', {'cpu'})
        check_backends('d', 'cpu', {'cpu', 'cpu-batch'})
        check_backends('d', 'cpu', {'cpu-batch'})
        check_backends('d', 'cuda', {'cuda'})

    def test_refused(self):
        for backend, seen in (('cpu', {'cpu', 'cuda'}), ('cuda', {'cuda', 'cpu'}), ('cuda', {'cpu-batch'}), ('cpu', {'sparse'})):
            with self.assertRaises(SystemExit):
                check_backends('d', backend, seen)


if __name__ == '__main__':
    unittest.main()
