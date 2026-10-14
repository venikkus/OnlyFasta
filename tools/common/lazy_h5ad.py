"""Read obs, var and one raw counts matrix from an .h5ad file without loading other data.

anndata's backed mode (read_h5ad(backed="r")) keeps only X on disk and reads every entry of
adata.layers fully into memory, so a large atlas with a raw counts layer does not fit in RAM.
open_counts() leaves the counts matrix on disk; slicing rows reads only those rows
(fast for CSR-stored matrices, which is the AnnData default for cells x genes).

    handle, adata = open_counts("atlas.h5ad", "counts")
    try:
        rows = adata.layers["counts"][np.array([0, 5, 9]), :]   # scipy sparse matrix
    finally:
        handle.close()

The returned object exposes obs, var, var_names, n_obs, n_vars and either X (counts_layer "X")
or layers[counts_layer]; row indices must be sorted increasingly.
"""
from types import SimpleNamespace

import h5py


def _anndata_io():
    try:
        from anndata.io import read_elem, sparse_dataset
    except ImportError:  # anndata < 0.11
        from anndata.experimental import read_elem, sparse_dataset
    return read_elem, sparse_dataset


def open_counts(path, counts_layer):
    """Return (open h5py file, AnnData-like view); close the file when done."""
    read_elem, sparse_dataset = _anndata_io()
    handle = h5py.File(path, "r")
    try:
        key = "X" if counts_layer == "X" else f"layers/{counts_layer}"
        if key not in handle:
            layers = list(handle["layers"].keys()) if "layers" in handle else []
            raise KeyError(f"{path}: no {key!r}; available layers: {layers}")
        element = handle[key]
        matrix = sparse_dataset(element) if isinstance(element, h5py.Group) else element
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
    except BaseException:
        handle.close()
        raise
    view = SimpleNamespace(
        obs=obs, var=var, var_names=var.index, n_obs=len(obs), n_vars=len(var),
        X=matrix if counts_layer == "X" else None,
        layers={} if counts_layer == "X" else {counts_layer: matrix},
    )
    return handle, view
