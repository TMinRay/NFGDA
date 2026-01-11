import os

def get_nf_input_name(l2_file, path_config):
    nf_input_file = l2_file.split('.')[0]+'.npz'
    return os.path.join(path_config.nf_dir, nf_input_file)

def get_nf_detection_name(l2_file, path_config):
    npzout = 'nf_pred'+l2_file+'.npz'
    return os.path.join(path_config.nf_preds_dir, npzout)

def get_nf_forecast_name(l2_file, path_config):
    npzout = 'nf_forecast'+l2_file+'.npz'
    return os.path.join(path_config.nf_forecast_dir, npzout)