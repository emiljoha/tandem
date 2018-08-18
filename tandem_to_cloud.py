from multiprocessing import Pool
from npy_to_cloud import upload, already_exists
import numpy as np
import PyTandem as pt
import gc
import tqdm


def generate_and_save_to_cloud(args):
    distribution = args['distribution']
    num_particles = args['num_particles']
    num_orbitals = args['num_orbitals']
    num_examples = args['num_examples']
    i = args['i']
    name = '%s_%s_%s_%s_%s.cho.npy' % (distribution, num_particles,
                                       num_orbitals, num_examples, i)
    cloud_path = 'data/%s/%s' % (num_orbitals, name)
    if already_exists(cloud_path, bucket_name='tandem_10'):
        print '%s already exists' % name
        return None
    data = pt.tandem(num_particles, num_examples, num_orbitals, distribution)
    gc.collect()
    upload(data, cloud_path, bucket_name='tandem_10')
    # print(name)
    gc.collect()


def wf_variance_from_energy_and_spin(num_orbitals, num_particles,
                                     one_basis_spins,
                                     one_basis_energies, temperature):
    zero_spin_bool_list = pt.is_spin_zero(num_orbitals, num_particles,
                                          one_basis_spins)
    # We want to do numpy broadcasting magic
    zero_spin_bool_list = np.array(zero_spin_bool_list)
    no_interaction_energies = pt.basis_energies(num_orbitals,
                                                num_particles,
                                                one_basis_energies)
    # Note the numpy broadcasting.
    energy_diff_from_ground = no_interaction_energies[0] - \
        np.array(no_interaction_energies)
    return np.exp(energy_diff_from_ground / temperature) * zero_spin_bool_list


def sample_around_HF_wf(temperature, num_examples, wave_function,
                        num_orbitals, num_particles, one_basis_spins,
                        one_basis_energies):
    wf_variance = wf_variance_from_energy_and_spin(num_orbitals,
                                                   num_particles,
                                                   one_basis_spins,
                                                   one_basis_energies,
                                                   temperature)
    sampled_wf = np.random.normal(wave_function, wf_variance,
                                  (num_examples, len(wave_function)))
    normalized_wf = sampled_wf / np.reshape(np.linalg.norm(sampled_wf,
                                                           axis=1),
                                            (len(sampled_wf),
                                             1))
    return normalized_wf


def HF_generate_and_save_to_cloud_from_wf(args):
    wave_function = args['wave_function']
    temperature = args['temperature']
    num_particles = args['num_particles']
    num_orbitals = args['num_orbitals']
    num_examples = args['num_examples']
    i = args['i']
    one_basis_spins = args['one_basis_spins']
    one_basis_energies = args['one_basis_energies']
    wf_hash = str(hash(tuple(wave_function)))[:5]
    name = 'HF-ES-%s_%s_%s_%s_%s_%s.cho.npy' % (temperature, wf_hash,
                                                num_particles,
                                                num_orbitals,
                                                num_examples, i)
    cloud_path = 'data/%s/%s' % (num_orbitals, name)
    if already_exists(cloud_path, bucket_name='tandem_10'):
        # print '%s already exists' % name
        return None
    wf_around_hf = sample_around_HF_wf(temperature, num_examples,
                                       wave_function, num_orbitals,
                                       num_particles, one_basis_spins,
                                       one_basis_energies)
    data = pt.tandem_on_wf(wf_around_hf, num_particles, num_orbitals)
    upload(data, cloud_path, bucket_name='tandem_10')
    del data
    del wave_function
    gc.collect()


def HF_generation(batch_size=10, num_files=10, num_orbitals=20,
                  num_cpus=2, D_files=[('../data/CA/HF/ca_1_D.txt', 4)],
                  temperatures=[0.1],
                  one_basis_spins=[-7, -5, -3, -1, 1, 3, 5, 7,
                                   -3, -1, 1, 3,
                                   -1, 1,
                                   -5, -3, -1, 1, 3, 5],
                  one_basis_energies=[-8.6240, -8.6240, -8.6240, -8.6240, -8.6240, -8.6240, -8.6240, -8.6240,
                                      -5.6793, -5.6793, -5.6793, -5.6793,
                                      -4.1370, -4.1370,
                                      -1.3829, -1.3829, -1.3829, -1.3829, -1.3829, -1.3829]):
    arguments = []
    for D_file_name, num_particles in D_files:
        assert(len(temperatures))
        D = np.transpose(np.loadtxt(D_file_name))
        wave_function = pt.hf_wf_from_D(D, num_orbitals, num_particles)
        for T in temperatures:
            for i in range(num_files):
                arguments.append(
                    {'wave_function': wave_function,
                     'temperature': T,
                     'num_particles': num_particles,
                     'num_orbitals': num_orbitals,
                     'num_examples': batch_size,
                     'i': i,
                     'one_basis_spins': one_basis_spins,
                     'one_basis_energies': one_basis_energies})
    pool = Pool(num_cpus)
    # HF_generate_and_save_to_cloud_from_wf(arguments[0])
    for res in tqdm.tqdm(pool.imap_unordered(
            HF_generate_and_save_to_cloud_from_wf,
            arguments), total=len(arguments)):
        if res is not None:
            print res.get()


def tandem_generation(num_orbitals, num_examples):
    arguments = []
    pool = Pool(16)
    for distribution in ['zero_spin']:
        for num_particles in [4]:
            for i in range(10):
                arguments.append({'distribution': distribution,
                                  'num_particles': num_particles,
                                  'num_orbitals': num_orbitals,
                                  'num_examples': num_examples, 'i': i})
        for _ in tqdm.tqdm(pool.imap_unordered(
                generate_and_save_to_cloud, arguments), total=len(arguments)):
            pass
        arguments = []
        print(distribution)


if __name__ == '__main__':
    HF_generation(batch_size=2, num_files=2, num_orbitals=20,
                  num_cpus=2, D_files=[('../data/CA/HF/ca_1_D.txt', 4),
                                       ('../data/CA/HF/ca_0.5_D.txt', 4),
                                       ('../data/CA/HF/ca_1.5_D.txt', 4),
                                       ('../data/CA/HF/ca_2_D.txt', 4)],
                  temperatures=[0.1],
                  one_basis_spins=[-7, -5, -3, -1, 1, 3, 5, 7,
                                   -3, -1, 1, 3, -1, 1, -5, -3, -1, 1, 3, 5],
                  one_basis_energies=[-8.6240, -8.6240, -8.6240, -8.6240, -8.6240,
                                      -8.6240, -8.6240, -8.6240, -5.6793, -5.6793, -5.6793, -5.6793,
                                      -4.1370, -4.1370, -1.3829, -1.3829, -1.3829, -1.3829, -1.3829,
                                      -1.3829])

