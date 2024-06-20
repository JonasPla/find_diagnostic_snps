import sys
import io
import unittest
from contextlib import redirect_stdout
from contextlib import contextmanager
import diagnostic_snps # type: ignore

@contextmanager
def redirect_stdin(new_stdin):
    old_stdin = sys.stdin
    sys.stdin = new_stdin
    try:
        yield
    finally:
        sys.stdin = old_stdin

def test_script(args, stdin_input=None):
    sys.argv = args
    f_stdout = io.StringIO()

    if stdin_input is not None:
        f_stdin = io.StringIO(stdin_input)
        stdin_context = redirect_stdin(f_stdin)
    else:
        stdin_context = redirect_stdin(sys.stdin)

    try:
        with redirect_stdout(f_stdout), stdin_context:
            diagnostic_snps.main()
    except SystemExit:
        # Catch SystemExit to prevent the test script from exiting
        pass

    return f_stdout.getvalue()

class TestDiagnosticSNPs(unittest.TestCase):

    def test_simple(self):
        stdout = test_script([
            "diagnostic_snps.py",
            "-i",
            "test.aln",
            "-ref",
            "A",
        ])
        expected_output = "pos\tA\tB\tC\tD\t\n5\tA\t-\tN\tN\t\n6\t-\tC\tC\tN\t\n"
        self.assertEqual(stdout, expected_output)

    def test_strict_neg(self):
        stdout = test_script([
            "diagnostic_snps.py",
            "-i",
            "test2.aln",
            "-ref",
            "A",
            "--strict"
        ])
        expected_output = "pos\tA\tB\tC\tD\tA2\tB2\tC2\tD2\t\n"
        self.assertEqual(stdout, expected_output)

    def test_strict_pos(self):
        stdout = test_script([
            "diagnostic_snps.py",
            "-i",
            "test.aln",
            "-ref",
            "A",
            "--strict"
        ])
        expected_output = "pos\tA\tB\tC\tD\t\n5\tA\t-\tN\tN\t\n6\t-\tC\tC\tN\t\n"
        self.assertEqual(stdout, expected_output)

    def test_ignore(self):
        stdout = test_script([
            "diagnostic_snps.py",
            "-i",
            "test.aln",
            "-ref",
            "A",
            "--ignore",
            "B"
        ])
        expected_output = "pos\tA\tC\tD\t\n6\t-\tC\tN\t\n"
        self.assertEqual(stdout, expected_output)

    def test_stdin(self):
        stdin_input = ">A\nACGTAA-GT\n>B\nACGT-ACGT\n>C\nACGTNACGT\n>D\nACGTNANGT"
        stdout = test_script(["diagnostic_snps.py", "-ref", "A"], stdin_input)
        expected_output = "pos\tA\tB\tC\tD\t\n5\tA\t-\tN\tN\t\n6\t-\tC\tC\tN\t\n"
        self.assertEqual(stdout, expected_output)

if __name__ == "__main__":
    unittest.main()
