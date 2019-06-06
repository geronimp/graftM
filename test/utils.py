import unittest


class T(unittest.TestCase):
    def assertEqualFastaUnsorted(self, expected_str, observed_file_path):
        expected_splits = []
        expected_lines = list([s.strip() for
                               s in expected_str.split('\n') if s != ''])
        for i, l in enumerate(expected_lines):
            if i % 2 == 0 and l.strip() != '':
                expected_splits.append(
                    [expected_lines[i].strip(), expected_lines[i+1].strip()])
        observed_splits = []
        with open(observed_file_path) as f:
            lines = list([l.strip() for l in f.readlines()])
            for i, line in enumerate(lines):
                if i % 2 == 0:
                    observed_splits.append([lines[i], lines[i+1]])
        self.assertEqual(sorted(expected_splits), sorted(observed_splits))
