import mylib
import ctypes
import numpy as np

print("Try test_empty:")
mylib.test_empty()
 
print("\nTry test_add:")
print(mylib.test_add(34.55, 23))

print("\nTry test_add_double:")
print(mylib.test_add_double(34.55, 43.0))

# Create a 25 elements array
numel = 25
data = (ctypes.c_int * numel)(*[x for x in range(numel)])

# Pass the above array and the array length to C:
print("\nTry passing an array of 25 integers to C:")
mylib.test_passing_array(data, numel)
 
print("data from Python after returning from C:")
for indx in range(numel):
    print(data[indx], end=" ")
print("")

x = [0.0, 0.0, 0.0, 0.0]
y = [1.0, 2.0, 3.0, 4.0]
z = [3.2, 4.2, 5.2, 6.2]
xarr = (ctypes.c_double * len(x))(*x)
yarr = (ctypes.c_double * len(y))(*y)
zarr = (ctypes.c_double * len(z))(*z)

mylib.test_add_array(xarr, yarr, zarr, 4)

print("data from Python after returning from C:")
for indx in range(4):
    print(xarr[indx], end=" ")
print("")

np_arr = np.ctypeslib.as_array(xarr)
print(np_arr)

npx = np.array(x).astype(np.float64)
npy = np.array(y).astype(np.float64)
npz = np.array(z).astype(np.float64)
xa = npx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
ya = npy.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
za = npz.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

mylib.test_add_array(xa, ya, za, 4)
npx = np.ctypeslib.as_array(xa, shape=(4,))
print(npx)

mylib.test_axpby_double(xa, ya, za, 2.0, 2.0, 4)
npa = np.ctypeslib.as_array(xa, shape=(4,))
print(npa)

npx = np.array(x).astype(np.float32)
npy = np.array(y).astype(np.float32)
npz = np.array(z).astype(np.float32)
xb = npx.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
yb = npy.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
zb = npz.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
print(npx)
mylib.test_axpby_float(xb, yb, zb, 2.0, 2.0, 4)
npb = np.ctypeslib.as_array(xb, shape=(4,))
print(x)
print(npb)
print(npx)
































