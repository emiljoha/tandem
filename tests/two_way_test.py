import numpy as np
from sympy.utilities.iterables import multiset_permutations
from sympy.combinatorics.permutations import Permutation
import pdb
num_examples = 1
num_particles = 4
num_orbitals = 10


def orthonormal_matrix(num_particles, num_orbitals):
    X, R = np.linalg.qr(np.random.uniform(size=(num_orbitals, num_particles)))
    return X


def save_HF_states(file_name_HF, file_name_one_mat, num_matrices, num_particles, num_orbitals):
    file = open(file_name_HF, "w")
    file_mat = open(file_name_one_mat, "w")
    for m in range(num_matrices):
        X = orthonormal_matrix(num_particles, num_orbitals)
        # test X
        for i in range(num_particles):
            for j in range(num_particles):
                if i == j:
                    assert(np.dot(X[:, i], X[:, j]) - 1 < 1e-4)
                else:
                    assert(np.dot(X[:, i], X[:, j]) < 1e-4)
        assert(np.abs(np.trace(np.matmul(np.transpose(X), X))) - num_particles < 1e-5)
        init_state = np.concatenate((np.zeros(num_orbitals - num_particles),
                                     np.ones(num_particles)))
        assert(len(init_state) == num_orbitals)
        assert(np.shape(np.nonzero(init_state))[1] == num_particles)
        for state in multiset_permutations(init_state):
            value = 0
            index = np.sort(num_orbitals - 1 - np.nonzero(state)[0])
            for p in multiset_permutations(np.arange(num_particles)):
                D = 1
                for i in range(num_particles):
                    D *= X[index[p[i]], i]
                sign = (-1)**Permutation(p).parity()
                value += D * sign
            file.write(str(value) + " ")
        file.write("\n")
        one_mat = np.matmul(X, np.transpose(X))
        for i in range(num_orbitals):
            for j in range(num_orbitals):
                file_mat.write(str(one_mat[i, j]) + " ")
            file_mat.write("\n")
        file_mat.write("\n")
    file.close()
    file_mat.close()

    
if __name__ == "__main__":
    save_HF_states("HF.wf", "one_mat.m", num_examples,
                   num_particles, num_orbitals)
