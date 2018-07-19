/*

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

double Error(double *F,int N) {
  double Fss = 0.0;
  int i;
  for (i=1;i<N-1;i++) {
      Fss = Fss + F[i]*F[i];
  }
  Fss = pow(Fss,0.5)/((double) N - 2);
  return Fss;
}

double Force(double *x, double *F, int N) {
  int i,j,a,b;
  double dx_ij;
  double Fss  = 0.0;
  double Fr,Fl,Fm,dF;
  int priv_nloops,thread_id;

  #pragma omp parallel private(priv_nloops, thread_id)
  {

  #pragma omp for
  for (i=0;i<N;i++) {

      F[i] = 0.0;
      for (j=0;j<N;j++) {
        if (i != j) {
          dx_ij = x[i]-x[j];
          if (dx_ij == 0.0) {
            printf("Error.\n");
          }
          dF = 1.0/(dx_ij*dx_ij);
          //printf("i = %d; j = %d; Fm = %.2g; dF = %.2g; Fm/dF = %.1g\n",i,j,Fm,dF,Fm/dF);
          if (dx_ij > 0.0) {
            F[i] = F[i] + dF;
          } else {
            F[i] = F[i] - dF;
          }
        }
      }
      if (i > 0 && i < (N-1)) {
        //Fss = Fss + F[i]*F[i];
      }
    }

  }
  //Fss = pow(Fss,0.5)/((double) N - 2);
  return Fss;
}

write(double *x, int N, int t, double Fss, double alpha) {
  int i;
  FILE *fp;
  fp = fopen("pushcode.txt", "w");
  fprintf(fp, "%d\n",t);
  fprintf(fp, "%.16f\n",Fss);
  fprintf(fp, "%.16f\n",alpha);
  for (i=0;i<N;i++) {
    fprintf(fp, "%.16f\n",x[i]);
  }
  fclose(fp);
}

int read(double *x) {
  FILE *stream;
  char *line = NULL;
  size_t len = 0;
  int read;
  int i = 0;
  stream = fopen("x.txt", "r");
  if (stream == NULL) {
    return 1;
  }
 
  while ((read = getline(&line, &len, stream)) != -1) {
    //printf("Retrieved line of length %d :\n", read);
    //printf("%s", line);
    x[i] = atof(line);
    //printf("i = %d; %f\n",i,x[i]);
    i = i+1;
  }
  free(line);
  fclose(stream);
  return i;
}

int main(int argc, char * argv[]) {

  int i,j;
  int N = atoi(argv[1]);
  int cont = 1;
  int t = 0;
  int t_stop = atoi(argv[2]);

  double dx_ij = 0.0;
  double Fss = 0.0;
  double Fsso = 0.0;
  double alpha = 1.0;
  double alphao = 1.0;
  double F_stop = 0.0001;

  double *F,*Fo;
  double *x,*xo;
  F = (double *) malloc(N*sizeof(double));
  if (F == NULL) {
    printf("pushcode.c: Could not allocate F.");
    return 1;
  }
  Fo = (double *) malloc(N*sizeof(double));
  if (Fo == NULL) {
    printf("pushcode.c: Could not allocate Fo.");
    return 1;
  }
  x = (double *) malloc(N*sizeof(double));
  if (x == NULL) {
    printf("pushcode.c: Could not allocate x.");
    return 1;
  }
  xo = (double *) malloc(N*sizeof(double));
  if (xo == NULL) {
    printf("pushcode.c: Could not allocate xo.");
    return 1;
  }

  read(x);

  //for (i=0;i<N;i++) {
    //x[i] = -1.0 + 2.0*((double) i)/(N-1);
    //printf("pushcode.c: x[%d] = %f\n",i,x[i]);
  //}
 
  Fss = Force(x,F,N); // Update F

  while (t < t_stop && Fss > F_stop) {

    //printf("pushcode.c: t = %03d\n",t);
    for (i=0;i<N;i++) {
      //printf("pushcode.c: F[%d] = %f; x[%d] = %f\n",i,F[i],i,x[i]);
    }

    Fss = Force(x,F,N); // Update F

    for (i=0;i<N;i++) {
      xo[i] = x[i];
      Fo[i] = F[i];
    }
    Fsso = Fss;

    alpha  = 2*alpha;
    alphao = 2*alpha;

    //alpha  = 1.0;
    //alphao = 1.0;

    while (1) {

      for (i=0;i<N;i++) {
        //printf("pushcode.c: x[%d] = %f\n",i,x[i]);
      }
      //printf("pushcode.c: t = %03d; Fss = %.6g; alpha = %.6g\n",t,Fss,alpha);

      for (i=0;i<N-1;i++) {
        x[i] = xo[i] + alpha*Fo[i];

        if (i == 0   && x[i] < -1.0) {x[i] = -1.0;}
        if (i == N-1 && x[i] >  1.0) {x[i] =  1.0;}

        if (i > 0 && x[i] < x[i-1]) {
          alpha = alpha/2.0;
          break;
        }
      }

      if (i == N-1) { // No change in alpha yet.
        Fss = Force(x,F,N);
        if (Fss >= Fsso) {
          alpha = alpha/2.0;
        }
      }

      if (alpha == alphao || alpha < 2.3e-16) {
        break;
      }
      alphao = alpha;
    }
    //printf("pushcode.c: n = %d; t = %03d; Fss = %.2f; alpha = %.6g\n",N,t,Fss,alpha);
    //printf("n = %03d; t = %03d; alpha = %.2e; Fss = %.2e;\n",N,t,alpha,Fss);

    t = t+1;
  }

  for (i=0;i<N;i++) {
    //printf("pushcode.c: x[%d] = %f\n",i,x[i]);
  }
  //printf("pushcode.c: N = %d; t = %03d; Fss = %.6g; alpha = %.6g\n",N,t,Fss,alpha);

  write(x,N,t,Fss,alpha);

  free(x);
  free(xo);
  free(F);
  free(Fo);

  return 0;

}
