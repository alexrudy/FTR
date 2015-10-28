#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def main():
    """Try some array unpacking. Accept shape from arguments."""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("shape", type=int, nargs=2, help="Shape to pack/unpack.")
    parser.add_argument('-r','--random', action='store_true', help="Use random numbers.")
    opt = parser.parse_args()
    
    shape = tuple(opt.shape)
    nft = int(shape[-1] / 2) + 1
    shape_ft = tuple(list(shape)[:-1] + [nft])
    print("Packed Shape = ({}x{}) -> ({}x{}) = Unpacked Shape".format(*(shape_ft + shape)))
    
    if opt.random:
        in_array = np.random.randn(*shape_ft) + 1j * np.random.randn(*shape_ft)
    else:
        in_array = np.arange(shape_ft[0])[:,None] + 1j * np.arange(shape_ft[1])[None, :]
    in_array.shape = shape_ft
    

    output_array_c = unpack_array_c(in_array, shape)
    output_array_py = unpack_array_py(in_array, shape)
    
    success = np.allclose(output_array_c, output_array_py)
    if success:
        print("Unpacking Success!")
    else:
        print("Array data does not match!")
        print("Input Array ({}x{}):".format(*shape_ft))
        print(in_array)
        print("(c) Output Array ({}x{}):".format(*shape))
        print(output_array_c)
        print("(py) Output Array ({}x{}):".format(*shape))
        print(output_array_py)
    
def _prepare_arrays(input_data, output_shape):
    """Prepare input and output arrays."""
    input_data  = np.asarray(input_data)
    input_shape = input_data.shape
    output_data = np.empty(output_shape, dtype=input_data.dtype)
    
    # Total number of elements in each type of array.
    nft = np.prod(input_shape)
    nn  = np.prod(output_shape)
    
    # Extract single dimension parameters, where relevant.
    ny = output_shape[0]
    nf = input_shape[-1]
    nx = output_shape[-1]
    if ny != input_shape[0]:
        raise ValueError("Shape mismatch: ({}x ) -> ({}x )".format(input_shape[0], ny))
    if (nx // 2) + 1 != nf:
        raise ValueError("Shape mismatch: ( x{}) -> ( x{}), expected ( x{})".format(nf, nx, (nx // 2) + 1))
    return output_data
    
def unpack_array_c(input_data, output_shape):
    """Unpack an array."""
    input_data  = np.asarray(input_data)
    output_data = _prepare_arrays(input_data, output_shape)
    input_shape = input_data.shape
    
    # Total number of elements in each type of array.
    nft = np.prod(input_shape)
    nn  = np.prod(output_shape)
    
    # Extract single dimension parameters, where relevant.
    ny = output_shape[0]
    nf = input_shape[-1]
    nx = output_shape[-1]
    
    # Determine the central element.
    nfo = nf - 1 if nx % 2 == 0 else nf
    for i in range(nft):
        
        # Extract the cartesian coordinates in the original grid.
        y = i // nf
        x = i %  nf
        
        # re-map onto the wider grid.
        io = (y * nx) + x
        ii = (y * nf) + x
        output_data.flat[io] = input_data.flat[ii]
        
        # Check our indexing math with numpy's version.
        assert np.allclose([y,x], np.unravel_index(ii, input_shape))
        assert np.allclose([io], np.ravel_multi_index([y, x], output_shape))
        
        # If we are at y=0, flip around center.
        if y == 0 and 0 < x < nfo:
            yo = 0
            xo = nx - x
            io = (yo * nx) + xo
            
            # Check our math!
            assert np.allclose([io], np.ravel_multi_index([yo, xo], output_shape))
            assert np.allclose([yo, xo], np.unravel_index(io, output_shape))
            
            # Insert data.
            output_data.flat[io] = np.conj(input_data.flat[ii])
            
        # If we are beyond y=0, flip in both axes.
        if y > 0 and 0 < x < nfo:
            yo = ny - y
            xo = nx - x
            io = (yo * nx) + xo
            
            # Check our math!
            assert np.allclose([io], np.ravel_multi_index([yo, xo], output_shape))
            assert np.allclose([yo, xo], np.unravel_index(io, output_shape))
            
            # Insert data.
            output_data.flat[io] = np.conj(input_data.flat[ii])
    
    return output_data
    
def unpack_array_py(input_data, output_shape):
    """Try unpacking the IDL way, and compare the results."""
    input_data  = np.asarray(input_data)
    output_data = _prepare_arrays(input_data, output_shape)
    input_shape = input_data.shape
    
    # Extract single dimension parameters, where relevant.
    ny = output_shape[0]
    nf = input_shape[-1]
    nx = output_shape[-1]
    o = 1 if nx % 2 == 1 else 2
    output_data[ :  , 0:nf] = input_data
    output_data[0   ,nf:nx] = np.conj(input_data[0,nf-o:0:-1])
    output_data[1:ny,nf:nx] = np.conj(input_data[ny-1:0:-1,nf-o:0:-1])
    return output_data

if __name__ == '__main__':
    main()