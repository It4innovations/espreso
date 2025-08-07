#!/usr/bin/env python3

import os, argparse
import numpy as np

def load_matrices(root_dir="results/last/DEBUG"):
    matrices = {}

    def read_sparse(file_path, sym=True):
        with open(file_path, 'r') as f:
            [nrows, ncols, nonzeros] = map(int, f.readline().split())
            A = np.zeros((nrows, ncols), dtype=float)
            for i, line in enumerate(f):
                [r, c, v] = line.split()
                A[int(r), int(c)] = float(v)
                if sym:
                    A[int(c), int(r)] = float(v)
        return A

    def read_dense(file_path):
        with open(file_path, 'r') as f:
            [nrows, ncols] = map(int, f.readline().split())
            A = np.zeros((nrows, ncols), dtype=float)
            for r, line in enumerate(f):
                for c, v in enumerate(line.split()):
                    A[int(r), int(c)] = float(v)
        return A

    def name(file):
        if path[5].removesuffix(".txt")[-1].isdigit():
            return path[5].removesuffix(".txt")
        else:
            return path[5].removesuffix(".txt") + "0"

    for dirpath, dirnames, filenames in sorted(os.walk(root_dir)):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            path = os.path.relpath(file_path, root_dir).split(os.sep)
            key = ":".join(path[:4])
            if key not in matrices:
                matrices[key] = {}

            if path[4] == "system" and path[5].startswith("K"):
                matrices[key][name(path[5])] = read_sparse(file_path)
            if path[4] == "system" and path[5].startswith("RegMat"):
                matrices[key][name(path[5])] = read_sparse(file_path)
            elif path[4] == "system" and path[5].startswith("R"):
                matrices[key][name(path[5])] = read_dense(file_path)

    return matrices

def analyze(args, matrices):
    invalid_kernels = 0
    invalid_regularization = 0
    for path, feti in matrices.items():
        if path.endswith("0"):
            print("\n===============================\n")
        print(path)
        for name, matrix in sorted(feti.items()):
            if name.startswith("K"):
                domain = name[1:]
                K = matrix
                reg_mat = feti["RegMat" + domain]
                N = feti["R" + domain]

                eps = np.finfo(float).eps
                tol = max(K.shape) * eps
                KN = np.linalg.norm(np.matmul(K, N.T), ord=2, axis=0)
                if args.svd:
                    s = np.linalg.svd(K, compute_uv=False)
                    s_reg = np.linalg.svd(K + reg_mat, compute_uv=False)

                print("  DOMAIN " + domain)
                if args.verbose:
                    e_k = sorted(np.linalg.eigvalsh(K))
                    e_reg = sorted(np.linalg.eigvalsh(K + reg_mat))
                    # U, s, V = np.linalg.svd(K, full_matrices=False)
                    print("  eig(K)       ", ", ".join(map("{:+.2e}".format, e_k[:8])), " ... ", ", ".join(map("{:+.2e}".format, e_k[-4:])))
                    print("  eig(K+RegMat)", ", ".join(map("{:+.2e}".format, e_reg[:8])), " ... ", ", ".join(map("{:+.2e}".format, e_reg[-4:])))
                    print("  max(K)       ", "{:+.2e}".format(K.max()))
                    print("  norm(K*R)    ", ", ".join(map("{:+.2e}".format, KN)))
                    print("  norm(K*R)/max", ", ".join(map("{:+.2e}".format, KN / K.max())))
                    # print("  SVD(K)       ", "U:", U.shape, "V:", V.shape, "s:", ", ".join(map("{:+.2e}".format, s[:4])), " ... ", ", ".join(map("{:+.2e}".format, s[-8:])))
                    if args.svd:
                        print("  SVD(K)       ", "s:", ", ".join(map("{:+.2e}".format, s[:4])), " ... ", ", ".join(map("{:+.2e}".format, s[-8:])))

                print("  ----------------------------")
                print("  R.rows                     = {}".format(N.shape[0]))
                print("  max(|K*R|)/max(K)          = {:+.2e}".format(max(KN)/K.max()))
                print("  COND                       = {:+.2e}".format(s_reg[0] / s_reg[-1]))
                if args.svd:
                    print("  zero sing. vals (K)        = {}".format(np.sum(s <= tol * np.max(s))))
                    print("  zero sing. vals (K+RegMat) = {}".format(np.sum(s_reg <= tol * np.max(s_reg))))
                    if np.sum(s <= tol * np.max(s)) != N.shape[0]:
                        print("  INVALID NUMBER OF KERNELS")
                        invalid_kernels += 1
                    if np.sum(s_reg <= tol * np.max(s_reg)):
                        print("  INVALID REGULARIZATION")
                        invalid_regularization += 1
                print("  ----------------------------")
    
    print("INVALID KERNELS {}".format(invalid_kernels))
    print("INVALID REGULARIZATION {}".format(invalid_regularization))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script for analysation of matrices printed by espreso.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose logging." )
    parser.add_argument('-s', '--svd', action='store_true', help="Compute SVD decomposition." )
    args = parser.parse_args()

    analyze(args, load_matrices())








