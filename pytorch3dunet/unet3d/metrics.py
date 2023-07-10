import importlib
import torch
from pytorch3dunet.unet3d.losses import compute_per_channel_dice
from pytorch3dunet.unet3d.utils import get_logger, get_attr, expand_as_one_hot
import numpy as np
from openbabel import pybel

logger = get_logger('EvalMetric')


class DiceCoefficient:
    """Computes Dice Coefficient.
    Generalized to multiple channels by computing per-channel Dice Score
    (as described in https://arxiv.org/pdf/1707.03237.pdf) and theTn simply taking the average.
    Input is expected to be probabilities instead of logits.
    This metric is mostly useful when channels contain the same semantic class (e.g. affinities computed with different offsets).
    DO NOT USE this metric when training with DiceLoss, otherwise the results will be biased towards the loss.
    """

    def __init__(self, epsilon=1e-6, **kwargs):
        self.epsilon = epsilon

    def __call__(self, input, target, pdbData=None):
        # Average across channels in order to get the final score
        return torch.mean(compute_per_channel_dice(input, target, epsilon=self.epsilon))


class MeanIoU:
    """
    Computes IoU for each class separately and then averages over all classes.
    """

    def __init__(self, skip_channels=(), ignore_index=None, thres=0.5, **kwargs):
        """
        :param skip_channels: list/tuple of channels to be ignored from the IoU computation
        :param ignore_index: id of the label to be ignored from IoU computation
        """
        self.ignore_index = ignore_index
        self.skip_channels = skip_channels
        self.thres = thres

    def __call__(self, input, target, pdbObj=None):
        """
        :param input: 5D probability maps torch float tensor (NxCxDxHxW)
        :param target: 4D or 5D ground truth torch tensor. 4D (NxDxHxW) tensor will be expanded to 5D as one-hot
        :return: intersection over union averaged over all channels
        """
        assert input.dim() == 5

        n_classes = input.size()[1]

        if target.dim() == 4:
            target = expand_as_one_hot(target, C=n_classes, ignore_index=self.ignore_index)

        assert input.size() == target.size()

        per_batch_iou = []
        for _input, _target in zip(input, target):
            binary_prediction = self._binarize_predictions(_input, n_classes, self.thres)

            if self.ignore_index is not None:
                # zero out ignore_index
                mask = _target == self.ignore_index
                binary_prediction[mask] = 0
                _target[mask] = 0

            # convert to uint8 just in case
            binary_prediction = binary_prediction.byte()
            _target = _target.byte()

            per_channel_iou = []
            for c in range(n_classes):
                if c in self.skip_channels:
                    continue

                per_channel_iou.append(self._jaccard_index(binary_prediction[c], _target[c]))

            assert per_channel_iou, "All channels were ignored from the computation"
            mean_iou = torch.mean(torch.tensor(per_channel_iou))
            per_batch_iou.append(mean_iou)

        return torch.mean(torch.tensor(per_batch_iou))

    def _binarize_predictions(self, input, n_classes, thres):
        """
        Puts 1 for the class/channel with the highest probability and 0 in other channels. Returns byte tensor of the
        same size as the input tensor.
        """
        if n_classes == 1:
            # for single channel input just threshold the probability map
            result = input >= thres
            return result.long()

        _, max_index = torch.max(input, dim=0, keepdim=True)
        return torch.zeros_like(input, dtype=torch.uint8).scatter_(0, max_index, 1)

    def _jaccard_index(self, prediction, target, pdbObj=None):
        """
        Computes IoU for a given target and prediction tensors
        """
        return torch.sum(prediction & target).float() / torch.clamp(torch.sum(prediction | target).float(), min=1e-8)


def get_evaluation_metric(config):
    """
    Returns the evaluation metric function based on provided configuration
    :param config: (dict) a top level configuration object containing the 'eval_metric' key
    :return: an instance of the evaluation metric
    """

    modules = ['pytorch3dunet.unet3d.metrics', 'pytorch3dunet.unet3d.pdb_metrics']

    assert 'eval_metric' in config, 'Could not find evaluation metric configuration'
    metric_config = config['eval_metric']
    metric_class = get_attr(metric_config['name'], modules)
    return metric_class(**metric_config)


def get_log_metrics(config):

    modules = ['pytorch3dunet.unet3d.metrics', 'pytorch3dunet.unet3d.pdb_metrics']

    ret = []
    if 'log_metrics' in config:
        metric_configs = config['log_metrics']
        for metric_config in metric_configs:
            metric_class = get_attr(metric_config['name'], modules)
            ret.append(metric_class(**metric_config))
    return ret


# Metrics from PURESNET
class BindingPocket:
    def __init__(self, occupied_cells):
        self.occupied_cells = occupied_cells


def create_3d_grid(pocket, resolution):
    min_coords = np.min(pocket.occupied_cells, axis=0)
    shifted_coords = pocket.occupied_cells - min_coords

    max_coords = np.max(shifted_coords, axis=0)
    grid_shape = np.ceil((max_coords) / resolution).astype(int) + 1
    grid = np.zeros(grid_shape, dtype=bool)

    for cell in shifted_coords:
        cell_idx = np.floor(cell / resolution).astype(int)
        grid[tuple(cell_idx)] = True
    return grid


def intersection_over_union(pocket1, pocket2, resolution):
    grid1 = create_3d_grid(pocket1, resolution)
    grid2 = create_3d_grid(pocket2, resolution)

    common_grid_shape = np.maximum(grid1.shape, grid2.shape)

    grid1_padded = np.pad(
        grid1,
        [(0, int(common_grid_shape[i] - grid1.shape[i])) for i in range(3)],
        mode='constant',
        constant_values=False
    )

    grid2_padded = np.pad(
        grid2,
        [(0, int(common_grid_shape[i] - grid2.shape[i])) for i in range(3)],
        mode='constant',
        constant_values=False
    )

    intersection_grid = np.logical_and(grid1_padded, grid2_padded)
    union_grid = np.logical_or(grid1_padded, grid2_padded)

    intersection_volume = np.sum(intersection_grid) * resolution ** 3
    union_volume = np.sum(union_grid) * resolution ** 3

    return intersection_volume / union_volume


def intersection_over_lig(lig, pocket, resolution):
    grid1 = create_3d_grid(lig, resolution)
    grid2 = create_3d_grid(pocket, resolution)

    common_grid_shape = np.maximum(grid1.shape, grid2.shape)

    grid1_padded = np.pad(
        grid1,
        [(0, int(common_grid_shape[i] - grid1.shape[i])) for i in range(3)],
        mode='constant',
        constant_values=False
    )

    grid2_padded = np.pad(
        grid2,
        [(0, int(common_grid_shape[i] - grid2.shape[i])) for i in range(3)],
        mode='constant',
        constant_values=False
    )

    intersection_grid = np.logical_and(grid1_padded, grid2_padded)

    intersection_volume = np.sum(intersection_grid) * resolution ** 3
    lig_volume = np.sum(grid1_padded) * resolution ** 3

    return intersection_volume / lig_volume


def coordinates(molecule):
    ligand_coords = [atom.coords for atom in molecule.atoms]
    return np.array(ligand_coords)


def get_DVO(pocket1_molecule, pocket2_molecule, resolution=1):
    pocket1_coords = coordinates(pocket1_molecule)
    pocket2_coords = coordinates(pocket2_molecule)
    pocket1 = BindingPocket(pocket1_coords)
    pocket2 = BindingPocket(pocket2_coords)
    return intersection_over_union(pocket1, pocket2, resolution)


def get_PLI(ligand_molecule, pocket_molecule, resolution=1):
    lig_coords = coordinates(ligand_molecule)
    pkt_coords = coordinates(pocket_molecule)
    ligand = BindingPocket(lig_coords)
    pocket = BindingPocket(pkt_coords)
    return intersection_over_lig(ligand, pocket, resolution)