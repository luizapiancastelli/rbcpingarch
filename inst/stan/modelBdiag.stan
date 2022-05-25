data {
  int<lower=1> n;
  matrix[n,2] Y;
  real med_y1;
  real med_y2;
}
parameters {
  matrix<lower = 0, upper = 1>[2,2] A;
  vector<lower = 0, upper = 1-max( [A[1,1]+A[2,1], A[1,2]+A[2,2]] )>[2] B;
  matrix<lower = 0>[2,1] omega;
  real phi;
}
transformed parameters{
 matrix<lower = 0>[n,2] lambda;
 matrix[2,2] I;
 matrix[2,2] B_mat;
 I = diag_matrix(rep_vector(1.0, 2));
 B_mat = diag_matrix(B);

 lambda[1,1] = omega[1,1] + A[1,1]*mean(Y[,1]) + B_mat[1,1]*med_y1;
 lambda[1,2] = omega[2,1] + A[2,2]*mean(Y[,2]) + B_mat[2,2]*med_y2;

  for(t in 2:n){
   lambda[t,] = to_row_vector( to_vector(omega') + A*(lambda[t-1,]') + B_mat*(Y[t-1, ]'));
  }
}
model {
  vector[n-1] y1; vector[n-1] y2;
  vector[n-1] l1; vector[n-1] l2;
  y1 = Y[2:n,1]; y2 = Y[2:n,2];
  l1 = lambda[2:n, 1]; l2 = lambda[2:n, 2];

  target += dot_product(y1,log(l1)) +  dot_product(y2,log(l2));
  target += -dot_product(l1, (1 + y2*(exp(phi) -1)) );
  target += -dot_product(l2, exp( -l1*(exp(phi)-1) + phi*y1) );
  target += phi*dot_product(y1, y2);
}
