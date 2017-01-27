
import os

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(ESPRESO_TESTS)

if __name__ == '__main__':
    print "\nBENCHMARKS TESTS:"
    os.system("python " + os.path.join(ROOT, "tests", "benchmarks.py"))
    print "\nINPUT TESTS:"
    os.system("python " + os.path.join(ROOT, "tests", "input.py"))
    print "\nSOLVER TESTS:"
    os.system("python " + os.path.join(ROOT, "tests", "solver.py"))