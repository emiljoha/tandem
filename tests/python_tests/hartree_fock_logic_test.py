import numpy as np

def wf_from_D(D):
    num_particles, num_orbitals = np.shape(D)
    assert(num_particles == 4)
    wf_dict = {}
    for j1 in range(num_orbitals):
        for j2 in range(num_orbitals):
            for j3 in range(num_orbitals):
                for j4 in range(num_orbitals):
                    if len(set((j1, j2, j3, j4))) == 4: # No doubles
                        inversions, result = merge_sort([j1, j2, j3, j4])
                        result = tuple(result)
                        if result in wf_dict:
                            wf_dict[result] += D[0][j1]*D[1][j2]*D[2][j3]*D[3][j4] * (-1)**inversions
                        else:
                            wf_dict[result] = D[0][j1]*D[1][j2]*D[2][j3]*D[3][j4] * (-1)**inversions
    state = [False]*num_orbitals
    state[-num_particles:] = [True]*num_particles
    def indices(state):
        index = []
        for orb, occupation in enumerate(reversed(state)):
            if occupation is True:
                index.append(orb)
        return tuple(index)
    wf = []
    while True:
        wf.append(wf_dict[indices(state)])
        if not next_permutation(state):
            break
    return wf

def merge_sort(A):
    if len(A) <= 1:
        return 0, A
    middle = len(A)/2
    left_inversions, left = merge_sort(A[:middle])
    right_inversions, right = merge_sort(A[middle:])
    merge_inversions, merged = merge(left, right)
    inversions = left_inversions + right_inversions + merge_inversions
    return inversions, merged

def merge(left, right):
    result = []
    i, j, inversions = 0, 0, 0
    while i < len(left) and j < len(right):
        if left[i] <= right[j]:
            inversions += j
            result.append(left[i])
            i += 1
        else:
            result.append(right[j])
            j += 1
    inversions += j*(len(left)-i)
    result += left[i:]
    result += right[j:]
    return inversions, result
                    
def next_permutation(a):
    """Generate the lexicographically next permutation inplace.

    https://en.wikipedia.org/wiki/Permutation#Generation_in_lexicographic_order
    Return false if there is no next permutation.
    """
    # Find the largest index i such that a[i] < a[i + 1]. If no such
    # index exists, the permutation is the last permutation
    for i in reversed(range(len(a) - 1)):
        if a[i] < a[i + 1]:
            break  # found
    else:  # no break: not found
        return False  # no next permutation

    # Find the largest index j greater than i such that a[i] < a[j]
    j = next(j for j in reversed(range(i + 1, len(a))) if a[i] < a[j])

    # Swap the value of a[i] with that of a[j]
    a[i], a[j] = a[j], a[i]

    # Reverse sequence from a[i + 1] up to and including the final element a[n]
    a[i + 1:] = reversed(a[i + 1:])
    return True

def D_file_to_wf_file(D_file_name):
    D = np.loadtxt(D_file_name)
    import pdb; pdb.set_trace()
    D = np.transpose(D)
    wf = wf_from_D(D)
    wf_file_name = '%s_wf_py_%s' % (D_file_name[:-6], D_file_name[-5:])
    np.savetxt(wf_file_name, wf, newline=' ')
    return wf_file_name

class HamiltonianBucket:
    def __init__(self, hamiltonian_names, num_particles, interaction_scaling, Interaction):
        self._H = []
        for file_name in hamiltonian_names:
            H_touple = np.loadtxt('%s.dat' % file_name)
            H1 = H_touple[0]
            H2 = H_touple[1]
            # Factor two for mysterius reasons and the minus
            # sign is a difference in convention. nuclear
            # scaling vrenorm equally god given.
            if interaction_scaling == 'harmonic':
                particle_scale_factor = (4.0 - 1.0) / (float(num_particles) - 1.0)
                vrenorm = 2
            elif interaction_scaling == 'CA':
                particle_scale_factor = 1 / (float(num_particles) - 1.0)
                vrenorm = -(42.0 / (40.0 + float(num_particles)))**0.3 / 4.0
            else:
                raise ValueError('No %s interaction scaling found')
            for h in [-(H1 * particle_scale_factor + vrenorm * V_in * H2) for V_in in Interaction]:
                self._H.append(h)
        self._H = np.array(self._H)

    def get_all(self):
        return self._H

def Energy(D, H):
    rho = np.matmul(D, np.transpose(D))
    num_orbitals = np.shape(rho)[0]
    two_matrix = np.zeros((num_orbitals, num_orbitals, num_orbitals, num_orbitals))
    for i in range(num_orbitals):
        for j in range(num_orbitals):
            for k in range(num_orbitals):
                for l in range(num_orbitals):
                    two_matrix[i][j][k][l] = rho[i][l]*rho[j][k] - rho[i][k]*rho[j][l]
    flat_2rdm = np.reshape(two_matrix, [-1])
    # print max(flat_2rdm), max(two_rdm)
    # print min(flat_2rdm), min(two_rdm)
    return [-np.sum(h*flat_2rdm) for h in H]

def D_to_2rdm(D_file_name, two_matrix_file_name, I):
    D = np.loadtxt(D_file_name)
    two_rdm = np.loadtxt(two_matrix_file_name)
    H = HamiltonianBucket(['H_ca_20_edges'], 4, 'CA', [I]).get_all()
    E = Energy(D, H)
    E2 = [np.sum(h*two_rdm) for h in H]
    assert(abs(E[0] - E2[0]) < 0.1)


def test_hartree_fock_conversion():
    for D, twordm, I in [('ca_0_D.txt', 'two_matrix_0.txt', 0), ('ca_0.5_D.txt', 'two_matrix_1.txt', 0.5),
                      ('ca_1_D.txt', 'two_matrix_2.txt', 1), ('ca_1.5_D.txt', 'two_matrix_3.txt', 1.5),
                      ('ca_2_D.txt', 'two_matrix_4.txt', 2)]:
        D_to_2rdm(D, twordm, I)
    # for D_file_name in ['ca_0_D.txt', 'ca_0.5_D.txt', 'ca_1_D.txt', 'ca_1.5_D.txt', 'ca_2_D.txt']:
    #     print D_file_to_wf_file(D_file_name)
