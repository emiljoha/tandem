from tempfile import TemporaryFile
from multiprocessing import Pool
import numpy as np
import time
from google.cloud import storage
from PyTandem import tandem
import gc
import tqdm

def download(cloud_path, bucket):
    blob2 = bucket.get_blob(cloud_path)
    fh = TemporaryFile()
    blob2.download_to_file(fh)
    fh.seek(0)
    D = np.load(fh)
    fh.close()
    return D

def upload(data, cloud_path, bucket_name):
    client = storage.Client()
    bucket = client.get_bucket(bucket_name)
    fh = TemporaryFile()
    np.save(fh, data)
    del data
    gc.collect()
    fh.seek(0)
    blob = storage.Blob(cloud_path, bucket)
    blob.upload_from_file(fh)
    fh.close()
    del fh
    del bucket
    del client
    gc.collect()

def already_exists(cloud_path, bucket_name):
    client = storage.Client()
    bucket = client.get_bucket(bucket_name)
    blob = storage.Blob(cloud_path, bucket)
    return blob.exists(client)
    

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
    data = tandem(num_particles, num_examples, num_orbitals, distribution)
    gc.collect()
    upload(data, cloud_path, bucket_name='tandem_10')
    # print(name)
    gc.collect()

num_examples = 10
num_orbitals = 10
arguments = []
pool = Pool(2)
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

