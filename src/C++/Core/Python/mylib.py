""" Python wrapper for the C shared library mylib"""
import sys, platform
import ctypes, ctypes.util
 
# Find the library and load it
mylib_path = ctypes.util.find_library("./mylib")
if not mylib_path:
    print("Unable to find the specified library.")
    sys.exit()

try:
    mylib = ctypes.CDLL(mylib_path)
except OSError:
    print("Unable to load the system C library")
    sys.exit()

# Make the function names visible at the module level and add types
test_empty = mylib.test_empty
 
test_add = mylib.test_add
test_add.argtypes = [ctypes.c_float, ctypes.c_float]
test_add.restype = ctypes.c_float

test_add_double = mylib.test_add_double
test_add_double.argtypes = [ctypes.c_double, ctypes.c_double]
test_add_double.restype = ctypes.c_double

test_passing_array = mylib.test_passing_array
test_passing_array.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.c_int]
test_passing_array.restype = None

test_add_array = mylib.test_add_array
test_add_array.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int]
test_add_array.restype = None

test_axpby_double = mylib.test_axpby_double
test_axpby_double.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double, ctypes.c_int]
test_axpby_double.restype = None

test_axpby_float = mylib.test_axpby_float
test_axpby_float.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.c_float, ctypes.c_float, ctypes.c_int]
test_axpby_float.restype = None

