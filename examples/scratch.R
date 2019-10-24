bcf_out = readRDS('examples/my_data.rds', refhook = NULL)

z = c(1,0,0,1,0,0,1,1,0,1)
z = z[bcf_out$perm]

allfit = (bcf_out$yhat[,bcf_out$perm] - bcf_out$muy)/bcf_out$sdy

allfit_con = (bcf_out$mu[,bcf_out$perm] - bcf_out$muy)/bcf_out$sdy

Tm = bcf_out$tau[,bcf_out$perm] * (1.0/ (bcf_out$sdy*(bcf_out$b1 - bcf_out$b0)) ) 


b_diff_matrix = t(matrix(1, 10, 1)%*%(bcf_out$b1 - bcf_out$b0))

b0_matrix = t(matrix(1, 10, 1)%*%bcf_out$b0)

b_scale = t(t(b_diff_matrix)*z) + b0_matrix

all_fit_mod = Tm*b_scale