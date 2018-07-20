from multiprocessing import Pool
import numpy as np
import time
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

def sample_around_HF_wf(wave_function, temperature, num_examples):
    wave_function = np.array(wave_function)
    sampled_wf = \
        np.random.normal(wave_function,
                        np.array([(np.abs(wave_function) * temperature).tolist()] * num_examples))
    normalized_wf = sampled_wf / np.reshape(np.linalg.norm(sampled_wf, axis=1), (len(sampled_wf), 1))
    return normalized_wf


def HF_generate_and_save_to_cloud_from_wf(args):
    wave_function = args['wave_function']
    temperature = args['temperature']
    num_particles = args['num_particles']
    num_orbitals = args['num_orbitals']
    num_examples = args['num_examples']
    i = args['i']
    wf_hash = str(hash(tuple(wave_function)))[:5]
    name = 'HF-%s_%s_%s_%s_%s_%s.cho.npy' % (temperature, wf_hash, num_particles,
                                       num_orbitals, num_examples, i)
    cloud_path = 'data/%s/%s' % (num_orbitals, name)
    if already_exists(cloud_path, bucket_name='tandem_10'):
        # print '%s already exists' % name
        return None
    wf_around_hf = sample_around_HF_wf(wave_function, temperature, num_examples)
    data = pt.tandem_on_wf(wf_around_hf, num_particles, num_orbitals)
    upload(data, cloud_path, bucket_name='tandem_10')
    del data
    del normalized_wf
    del wave_function
    gc.collect()


def HF_generate_and_save_to_cloud_from_D(args):
    D_file_name = args['D_file_name']
    num_orbitals = args['num_orbitals']
    num_particles = args['num_particles']
    temperatures = args['temperatures']
    num_files = args['num_files']
    batch_size = args['batch_size']
    num_cpus = args['num_cpus']
    assert(len(temperatures))
    D = np.transpose(np.loadtxt(D_file_name))
    wave_function = pt.hf_wf_from_D(D, num_orbitals, num_particles)
    arguments = []
    for T in temperatures:
        for i in range(num_files):
            arguments.append(
                {'wave_function': wave_function,
                 'temperature': T,
                 'num_particles': num_particles,
                 'num_orbitals': num_orbitals,
                 'num_examples': batch_size,
                 'i': i})
    # HF_generate_and_save_to_cloud_from_wf(arguments[0])
    pool = Pool(num_cpus)
    for res in tqdm.tqdm(pool.imap_unordered(HF_generate_and_save_to_cloud_from_wf,
                                             arguments), total=len(arguments)):
        if res is not None:
            print res.get()

def HF_generation(batch_size=10, num_files=10, num_orbitals=20,
                  num_cpus=2, D_files=[('../data/CA/HF/ca_1_D.txt', 4)], temperatures=[0.1]):
    arguments = []
    for D_file_name, num_particles in D_files:
        arguments.append({'D_file_name': D_file_name,
                          'temperatures': temperatures,
                          'num_files': num_files,
                          'batch_size': batch_size,
                          'num_particles': num_particles,
                          'num_orbitals': num_orbitals,
                          'num_cpus': num_cpus})
    
    for args in arguments:
        HF_generate_and_save_to_cloud_from_D(args)


def tandem_generation():
    arguments = []
    pool = Pool(16)
    for distribution in ['zero_spin']:
        for num_particles in [4]:
            for i in range(10):
                arguments.append({'distribution': distribution,
                                  'num_particles': num_particles,
                                  'num_orbitals': num_orbitals,
                                  'num_examples': num_examples, 'i': i})
        for _ in tqdm.tqdm(pool.imap_unordered(generate_and_save_to_cloud, arguments), total=len(arguments)):
            pass
        arguments = []
        print(distribution)


if __name__ == '__main__':
    HF_generation(batch_size=10, num_files=1000, num_orbitals=20,
                  num_cpus=64, D_files=[('../data/CA/HF/ca_0_D.txt', 4),
                                        ('../data/CA/HF/ca_0.5_D.txt', 4),
                                        ('../data/CA/HF/ca_1.5_D.txt', 4),
                                        ('../data/CA/HF/ca_2_D.txt', 4)],
                  temperatures=[0.1])
