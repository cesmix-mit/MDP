/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

/* ----------------------------------------------------------------------
   uniform RN
------------------------------------------------------------------------- */

double cpuRandomUniform(int *seed)
{
  int k = seed[0]/IQ;
  seed[0] = IA*(seed[0]-k*IQ) - IR*k;
  if (seed[0] < 0) seed[0] += IM;
  double ans = AM*seed[0];
  return ans;
}

/* ----------------------------------------------------------------------
   gaussian RN
------------------------------------------------------------------------- */

double cpuRandomGaussian(int *seed, int *save, double *second)
{
  double first,v1,v2,rsq,fac;

  if (!save[0]) {
    do {
      v1 = 2.0*cpuRandomUniform(seed)-1.0;
      v2 = 2.0*cpuRandomUniform(seed)-1.0;
      rsq = v1*v1 + v2*v2;
    } while ((rsq >= 1.0) || (rsq == 0.0));
    fac = sqrt(-2.0*log(rsq)/rsq);
    second[0] = v1*fac;
    first = v2*fac;
    save[0] = 1;
  } else {
    first = second[0];
    save[0] = 0;
  }
  return first;
}

void cpuRandomResetSeed(int *seed, int *save, int ibase, double *coord)
{
  int i;

  char *str = (char *) &ibase;
  int n = sizeof(int);

  unsigned int hash = 0;
  for (i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }

  str = (char *) coord;
  n = 3 * sizeof(double);
  for (i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }

  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);

  seed[0] = hash & 0x7ffffff;
  if (!seed[0]) seed[0] = 1;
  
  for (i = 0; i < 5; i++) cpuRandomUniform(seed);
  
  save[0] = 0;    
}

double cpuGamdev(const int ia, int *seed)
{
  int j;
  double am,e,s,v1,v2,x,y;

  if (ia < 1) return 0.0;
  if (ia < 6) {
    x=1.0;
    for (j=1; j<=ia; j++)
      x *= cpuRandomUniform(seed);

    // make certain, that -log() doesn't overflow.
    if (x < 2.2250759805e-308)
      x = 708.4;
    else
      x = -log(x);
  } else {
  restart:
    do {
      do {
        do {
          v1 = cpuRandomUniform(seed);
          v2 = 2.0*cpuRandomUniform(seed) - 1.0;
        } while (v1*v1 + v2*v2 > 1.0);

        y=v2/v1;
        am=ia-1;
        s=sqrt(2.0*am+1.0);
        x=s*y+am;
      } while (x <= 0.0);

      if (am*log(x/am)-s*y < -700 || v1<0.00001) {
        goto restart;
      }

      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    } while (cpuRandomUniform(seed) > e);
  }
  return x;
}

