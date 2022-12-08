Example Data
============

Some examples, especially from the `misc` directory, required their input data
to be in MATLAB or GNU Octave compatible HDF5 files. These files are generated in
MATLAB using
```matlab
  save('-v7.3','filename.mat',...);
```
or in GNU Octave by calling
```matlab
  save('-hdf5', 'filename.h5',...);
```

The required datasets inside the files are mentioned in the description of the
examples.

The following data files are available:

 - `glyap-40x40-demo.h5`: Example data of dimension 40x40 for the generalized
                          Lyapunov solvers in `misc`.

