#define N_REDUCETHREADS 1024

unsigned int nextPow2( unsigned int x )
{
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}


__device__ real rng_uni(unsigned int *state)
{
        unsigned int x = *state;

	    x = x ^ (x >> 13);
		x = x ^ (x << 17);
		x = x ^ (x >> 5);

        *state = x;

        return ((real)x + (real)1.0) / 4294967296.0f;
}

template <class Basis>
__global__ void kernelWC(
						unsigned int *rng_state,
						const Basis * basis,
						const SimuParamsWC * pars,
						real *x, 
						real *y,
						real *z,
						real *phase,
						real *Gx,
						real *Gy,
						real *Gz
					)
					
{
  // thread index
  int idx = threadIdx.x + blockDim.x*blockIdx.x;
  // In the original WC code this was held in CONSTANT memory. But for any 
  real d_sigma = sqrt(2.0f * basis->getD() * pars->timestep);
	// RNG stuff
  unsigned int lrng_state = rng_state[idx];
  
  unsigned int t;
  volatile real dx, dy, dz, x1, y1, z1, x2, y2, z2, r, phi, phase1, kahanC, kahanT, kahanY;
  
  x1 = x[idx]; // temporary coordinates in register memory
  y1 = y[idx];
  z1 = z[idx];
  phase1 = 0.f; // accumulated particle phase, in register memory
  kahanC = 0.f; // accumulators for Kahan summation of phase
  kahanT = 0.f;
  kahanY = 0.f;
  
  for (t=0; t < pars->steps; t++)
    {

      // sample random displacements
      r = sqrt(-2.0f * log(rng_uni(&lrng_state)));
      phi = 2.0f * PI * (-0.5f + rng_uni(&lrng_state));
      dx = r * cos(phi);
      dy = r * sin(phi);
      r = sqrt(-2.0f * log(rng_uni(&lrng_state)));
      phi = 2.0f * PI * (-0.5f + rng_uni(&lrng_state));
      dz = r * cos(phi);
      
      // increment temporary coordinates
      x2 = x1 + dx * d_sigma;
      y2 = y1 + dy * d_sigma;
      z2 = z1 + dz * d_sigma;

      // Rejection sampling for boundary conditions
      if(!basis->inside(x2,y2,z2))
	{
	  // update coordinates
	  x2 = x1;
	  y2 = y1;
	  z2 = z1;
	}
      
	  kahanY = Gx[t]*x2 + Gy[t]*y2 + Gz[t]*z2 - kahanC;
      kahanT = phase1 + kahanY;
      kahanC = (kahanT - phase1) - kahanY;
      phase1 = kahanT;
	  
	  x1 = x2;
	  y1 = y2;
	  z1 = z2;
    }
  
	  // rng_state[idx] = lrng_state;
	  phase[idx] = phase1*GAMMA*pars->timestep;   // store accumulated phase
	  // x[idx] = x2;           // store particle positions
	  // y[idx] = y2;
	  // z[idx] = z2;
}



__global__ void reducePhaseKernel(real *d_idata,
				  real *d_odata,
				  unsigned int n,
				  unsigned int blockSize,
				  unsigned int creal)
{
  extern __shared__ real sdata[];
  
  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  // COMPUTE EFFECTIVE PHASE AT THIS POINT: phi = cos(G * phase)
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
  
  //real mySum = (i < n) ? cosf(G * d_idata[i]) : 0.f;
  real mySum = (creal==1 ?
		 ((i < n) ? cosf(d_idata[i]) : (real)0) :
		 ((i < n) ? sinf(d_idata[i]) : (real)0));

  if (i + blockSize < n) 
    mySum += (creal==1 ?
	      cosf(d_idata[i+blockSize]) :
	      sinf(d_idata[i+blockSize]));
  
  sdata[tid] = mySum;
  __syncthreads();
  
  // do reduction in shared mem
  for(unsigned int s=blockDim.x/2; s>32; s>>=1) 
    {
      if (tid < s)
	{
	  sdata[tid] = mySum = mySum + sdata[tid + s];
	}
      __syncthreads();
    }
  
  if (tid < 32)
    {
      // now that we are using warp-synchronous programming (below)
      // we need to declare our shared memory volatile so that the compiler
      // doesn't reorder stores to it and induce incorrect behavior.
      volatile real *smem = sdata;
      if (blockSize >=  64) { smem[tid] = mySum = mySum + smem[tid + 32]; }
      if (blockSize >=  32) { smem[tid] = mySum = mySum + smem[tid + 16]; }
      if (blockSize >=  16) { smem[tid] = mySum = mySum + smem[tid +  8]; }
      if (blockSize >=   8) { smem[tid] = mySum = mySum + smem[tid +  4]; }
      if (blockSize >=   4) { smem[tid] = mySum = mySum + smem[tid +  2]; }
      if (blockSize >=   2) { smem[tid] = mySum = mySum + smem[tid +  1]; }
    }
  
  // write result for this block to global mem 
  if (tid == 0) d_odata[blockIdx.x] = sdata[0];
}

void reducePhase(int size, int threads, int blocks, real *d_idata, real *d_odata,  char creal)
{
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  
  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = (threads <= 32) ? 2 * threads * sizeof(real) : threads * sizeof(real);
  
  reducePhaseKernel<<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size, threads, (unsigned int)creal);
}

__global__ void reduceSumKernel(real *d_idata,
			real *d_odata,
			unsigned int n,
			unsigned int blockSize)
{
  extern __shared__ real sdata[];

//  - takes log(n) steps for n input elements
//  - uses n threads
//  - only works for power-of-2 arrays
//			
//  This version unrolls the last warp to avoid synchronization where it 
//  isn't needed.
//  
//  Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory. 
//  In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.  
//  If blockSize > 32, allocate blockSize*sizeof(T) bytes.
  
  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
  
  real mySum = (i < n) ? d_idata[i] : 0.f;
  if (i + blockSize < n) 
    mySum += d_idata[i+blockSize];  
  
  sdata[tid] = mySum;
  __syncthreads();
  
  // do reduction in shared mem
  for(unsigned int s=blockDim.x/2; s>32; s>>=1) 
    {
      if (tid < s)
	{
	  sdata[tid] = mySum = mySum + sdata[tid + s];
	}
      __syncthreads();
    }
  
  if (tid < 32)
    {
      // now that we are using warp-synchronous programming (below)
      // we need to declare our shared memory volatile so that the compiler
      // doesn't reorder stores to it and induce incorrect behavior.
      volatile real *smem = sdata;
      if (blockSize >=  64) { smem[tid] = mySum = mySum + smem[tid + 32]; }
      if (blockSize >=  32) { smem[tid] = mySum = mySum + smem[tid + 16]; }
      if (blockSize >=  16) { smem[tid] = mySum = mySum + smem[tid +  8]; }
      if (blockSize >=   8) { smem[tid] = mySum = mySum + smem[tid +  4]; }
      if (blockSize >=   4) { smem[tid] = mySum = mySum + smem[tid +  2]; }
      if (blockSize >=   2) { smem[tid] = mySum = mySum + smem[tid +  1]; }
    }
  
  // write result for this block to global mem 
  if (tid == 0) d_odata[blockIdx.x] = sdata[0];
}


void reduceSum(int size, int threads, int blocks, real *d_idata, real *d_odata)
{
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  
  // when there is only one warp per block, we need to allocate two warps 
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = (threads <= 32) ? 2 * threads * sizeof(real) : threads * sizeof(real);
  
  reduceSumKernel<<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size, threads);
}


real computePhaseSx(int n, real* d_idata, char creal)
{
  real Sx = (real)0;
  char needReadBack = 1;
  
  int threads, blocks;
  int maxThreads = N_REDUCETHREADS;
  
  // calculate number of threads and blocks for initial reduction
  threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
  blocks = (n + (threads * 2 - 1)) / (threads * 2);
  
  // allocate mem for the result on host side
  real* h_odata = (real*) malloc(blocks*sizeof(real));
  for(int i=0; i<blocks; i++)
    {
      h_odata[i] = (real)0;
    }

  // allocate device memory and data
  real* d_odata = NULL;
  safe_cuda( cudaMalloc((void**) &d_odata, blocks*sizeof(real)) );
  
  // copy data directly to device memory
  safe_cuda( cudaMemcpy(d_odata, h_odata, blocks*sizeof(real), cudaMemcpyHostToDevice) );
  
  // execute the kernel
  //printf("Reduce threads: %i\t Blocks: %i\t n: %i\n", threads, blocks, n);
  reducePhase(n, threads, blocks, d_idata, d_odata,  creal);
  
  // sum partial block sums on GPU
  int s=blocks;
  while(s > 1) 
    {
      int threads = 0, blocks = 0;
      threads = (s < maxThreads*2) ? nextPow2((s + 1)/ 2) : maxThreads;
      blocks = (s + (threads * 2 - 1)) / (threads * 2);
      
      reduceSum(s, threads, blocks, d_odata, d_odata);
      
      s = (s + (threads*2-1)) / (threads*2);
    }
  
  if (s > 1)
    {
      // copy result from device to host
      safe_cuda( cudaMemcpy( h_odata, d_odata, s * sizeof(real), cudaMemcpyDeviceToHost) );
      
      for(int i=0; i < s; i++) 
	{
	  Sx += h_odata[i];
	}
      
      needReadBack = 0;
    }
  
  cudaThreadSynchronize();
  
  if (needReadBack == 1)
    {
      // copy final sum from device to host
      safe_cuda( cudaMemcpy( &Sx, d_odata, sizeof(real), cudaMemcpyDeviceToHost) );
    }
  
  // cleanup
  free(h_odata);
  safe_cuda(cudaFree(d_odata));
  
  return Sx;
}

