import numpy as np
import npy_to_cloud

def local_test_npy_to_cloud():
    bucket_name = 'tandem_10'
    cloud_path = 'test_data.npy'
    np.random.seed(42)
    test_data = np.random.uniform(size=(10, 20))
    assert not npy_to_cloud.already_exists(cloud_path, bucket_name)
    npy_to_cloud.upload(test_data, cloud_path, bucket_name)
    test_data_from_cloud = npy_to_cloud.download(cloud_path, bucket_name)
    np.testing.assert_array_equal(test_data, test_data_from_cloud)
    assert npy_to_cloud.already_exists(cloud_path, bucket_name)
    print 'npy_to_cloud tests PASSED'

if __name__ == '__main__':
    local_test_npy_to_cloud()
