from google.cloud import storage
from tempfile import TemporaryFile
import numpy as np
import gc

def download(cloud_path, bucket_name):
    client = storage.Client()
    bucket = client.get_bucket(bucket_name)
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
