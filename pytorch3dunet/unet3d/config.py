import torch
import yaml
from pathlib import Path
from argparse import ArgumentParser
import os
from pytorch3dunet.unet3d import utils
import logging
from pytorch3dunet.unet3d.utils import set_default_log_level, set_filename

checkpointname = "checkpoint"

def load_config(runconfigPath, nworkers, pdb_workers, device_str):
    runconfig = yaml.safe_load(open(runconfigPath, 'r'))
    runFolder = Path(runconfig.get('runFolder', Path(runconfigPath).parent))
    train_config = runFolder / 'train_config.yml'

    config = yaml.safe_load(open(train_config, 'r'))

    from pytorch3dunet.datasets.config import RunConfig
    class_config = RunConfig(runFolder=runFolder, runconfig=runconfig, nworkers=nworkers, pdb_workers=pdb_workers,
                          loaders_config=config['loaders'])

    logger = utils.get_logger('ConfigLoader')
    logger.info(f'Read config:\n{class_config.pretty_format()}')

    config['dry_run'] = runconfig.get('dry_run', False)
    config['dump_inputs'] = runconfig.get('dump_inputs', False)

    os.makedirs(class_config.loaders_config.tmp_folder, exist_ok=True)

    config['trainer']['checkpoint_dir'] = str(runFolder / checkpointname)


    if device_str is not None:
        logger.info(f"Device specified in config: '{device_str}'")
        if device_str.startswith('cuda') and not torch.cuda.is_available():
            logger.warn('CUDA not available, using CPU')
            device_str = 'cpu'
    else:
        device_str = "cuda:0" if torch.cuda.is_available() else 'cpu'
        logger.info(f"Using '{device_str}' device")

    device = torch.device(device_str)
    config['device'] = device
    return config, class_config

def parse_args():

    parser = ArgumentParser()
    parser.add_argument("-r", "--runconfig", dest='runconfig', type=str, required=True,
                        help=f"The run config yaml file")
    parser.add_argument("-p", "--pdbworkers", dest='pdbworkers', type=int, required=True,
                        help=f"Number of workers for the pdb data generation. Typically this can (and should) be "
                             f"higher than numworkers")
    parser.add_argument("-n", "--numworkers", dest='numworkers', type=int, required=True,
                        help=f"Number of workers")
    parser.add_argument("-d", "--device", dest='device', type=str, required=False,
                        help=f"Device")
    parser.add_argument("--debug", dest='debug', default=False, action='store_true')
    parser.add_argument("--profile", dest='profile', default=False, action='store_true')
    parser.add_argument("--logfile", dest='logfile', type=str, required=False)

    args = parser.parse_args()

    if args.debug:
        set_default_log_level(logging.DEBUG)
    if args.logfile:
        set_filename(args.logfile)

    runconfig = args.runconfig
    nworkers = int(args.numworkers)
    pdbworkers = int(args.pdbworkers)

    config, class_config = load_config(runconfig, nworkers, pdbworkers, args.device)
    return args, config, class_config