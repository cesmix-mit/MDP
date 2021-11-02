inline void Unmap(dstype *y, dstype *x, dstype *h, int *image, int triclinic, int dim, int n, int backend)
{
	if (backend == 1)
		cpuUnmap(y, x, h, image, triclinic, dim, n);
#ifdef USE_OMP
	if (backend == 4)
		ompUnmap(y, x, h, image, triclinic, dim, n);
#endif
#ifdef USE_HIP
	if (backend == 3)
		hipUnmap(y, x, h, image, triclinic, dim, n);
#endif
#ifdef USE_CUDA
	if (backend == 2)
		gpuUnmap(y, x, h, image, triclinic, dim, n);
#endif
}
