import PyTandem as pt

def lists_are_equal(list1, list2):
    if not len(list1) == len(list2):
        return False
    for i in range(len(list1)):
        if not list1[i] == list2[i]:
            return False
    return True

def test_is_spin_zero():
    num_orbitals = 4
    num_particles = 2
    spins = [-3, -1, 1, 3]
    expected_result = [False, False, True, True, False, False]
    actual_result = pt.is_spin_zero(num_orbitals, num_particles, spins)
    assert(lists_are_equal(expected_result, actual_result))

def test_basis_energies():
    num_orbitals = 4
    num_particles = 2
    energies = [1, 2, 3, 4]
    expected_result = [3, 4, 5, 5, 6, 7]
    actual_result = pt.basis_energies(num_orbitals, num_particles, energies)
    assert(lists_are_equal(expected_result, actual_result))

