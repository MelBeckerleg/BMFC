#Generate a rank k approximation of a database with missing data

import sys
sys.path.append('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_PYTHON_mwe/models/')
sys.path.append('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_PYTHON_mwe/solvers/')


import numpy as np
import data_init as di
import rank_k
import tiling_error
from GenerateOpts import GenerateOpts 
from SolverOpts import SolverOpts 



data_opts=GenerateOpts()
A=data_opts.generate_function(data_opts)
M=data_opts.subsample_function(A,data_opts)

solver_opts=SolverOpts()


[X,Y]=rank_k_solve(solver_opts,M)

err=calculate_error(X,Y,A,X_true,Y_true)

print('RMSE error from method: '+ solver_opts.tag+str(err.rmse))
