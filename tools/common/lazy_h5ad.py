"""Read obs, var and one raw counts matrix from an .h5ad file without loading other data.

anndata's backed mode (read_h5ad(backed="r")) keeps only X on disk and reads every entry of
adata.layers fully into memory, so a large atlas with a raw counts layer does not fit in RAM.
open_counts() leaves the counts matrix on disk; selecting rows reads them in large contiguous
blocks of bounded size. Reading row by row is avoided on purpose: files written by R/rhdf5 often
store each array as a single HDF5 chunk, and then every small read re-reads the whole array.

    handle, adata = open_counts("atlas.h5ad", "counts")
    try:
        rows = adata.layers["counts"][np.array([0, 5, 9]), :]   # scipy CSR matrix
    finally:
        handle.close()

The returned object exposes obs, var, var_names, n_obs, n_vars and either X (counts_layer "X")
or layers[counts_layer]; row indices must be sorted increasingly.
"""
from types import SimpleNamespace

import h5py
import numpy as np
import scipy.sparse as sp

# Max stored values read per block: ~12 bytes each (float64 data + int32 index) -> ~240 MB.
BLOCK_NNZ = 20_000_000
BLOCK_DENSE = 20_000_000


def _read_elem():
    try:
        from anndata.io import read_elem
    except ImportError:  # anndata < 0.11
        from anndata.experimental import read_elem
    return read_elem


def _check_rows(key, shape):
    rows, cols = key if isinstance(key, tuple) else (key, slice(None))
    if cols != slice(None):
        raise IndexError("Only row selection with all columns is supported")
    rows = np.asarray(rows, dtype=np.int64)
    if rows.ndim != 1 or (len(rows) and (rows[0] < 0 or rows[-1] >= shape[0])):
        raise IndexError("Row indices must be a 1-D array within the matrix")
    if np.any(np.diff(rows) <= 0):
        raise IndexError("Row indices must be sorted increasingly without duplicates")
    return rows


class CSRRows:
    """Row access to a CSR matrix stored in an HDF5 group (data, indices, indptr)."""

    def __init__(self, group, shape):
        self.data, self.indices = group["data"], group["indices"]
        self.indptr = group["indptr"][:].astype(np.int64)
        self.shape = tuple(int(x) for x in shape)
        self.dtype = self.data.dtype

    def __getitem__(self, key):
        rows = _check_rows(key, self.shape)
        starts, ends = self.indptr[rows], self.indptr[rows + 1]
        lengths = ends - starts
        data, indices, i = [], [], 0
        while i < len(rows):
            # One contiguous read covering as many selected rows as fit in the block budget.
            j = max(i + 1, int(np.searchsorted(ends, starts[i] + BLOCK_NNZ, side="right")))
            lo, hi = starts[i], ends[j - 1]
            block_data, block_indices = self.data[lo:hi], self.indices[lo:hi]
            n = lengths[i:j]
            offsets = np.repeat(starts[i:j] - lo, n)
            within = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
            data.append(block_data[offsets + within])
            indices.append(block_indices[offsets + within])
            i = j
        indptr = np.concatenate([[0], np.cumsum(lengths)])
        empty = np.array([], dtype=self.dtype)
        return sp.csr_matrix(
            (np.concatenate(data) if data else empty,
             np.concatenate(indices) if indices else np.array([], dtype=np.int32), indptr),
            shape=(len(rows), self.shape[1]))


class DenseRows:
    """Row access to a dense matrix stored as an HDF5 dataset."""

    def __init__(self, dataset):
        self.dataset, self.shape, self.dtype = dataset, dataset.shape, dataset.dtype

    def __getitem__(self, key):
        rows = _check_rows(key, self.shape)
        per_block = max(1, BLOCK_DENSE // max(1, self.shape[1]))
        parts, i = [], 0
        while i < len(rows):
            j = max(i + 1, int(np.searchsorted(rows, rows[i] + per_block, side="left")))
            block = self.dataset[rows[i]:rows[j - 1] + 1]
            parts.append(block[rows[i:j] - rows[i]])
            i = j
        return np.concatenate(parts) if parts else np.empty((0, self.shape[1]), self.dtype)


def open_counts(path, counts_layer):
    """Return (open h5py file, AnnData-like view); close the file when done."""
    read_elem = _read_elem()
    handle = h5py.File(path, "r")
    try:
        key = "X" if counts_layer == "X" else f"layers/{counts_layer}"
        if key not in handle:
            layers = list(handle["layers"].keys()) if "layers" in handle else []
            raise KeyError(f"{path}: no {key!r}; available layers: {layers}")
        element = handle[key]
        if isinstance(element, h5py.Group):
            encoding = element.attrs.get("encoding-type", "")
            if encoding != "csr_matrix":
                raise ValueError(f"{path}: {key} is stored as {encoding!r}; row access needs csr_matrix "
                                 "(cells x genes). Convert it with adata.X = adata.X.tocsr() and save again.")
            matrix = CSRRows(element, element.attrs["shape"])
        else:
            matrix = DenseRows(element)
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
