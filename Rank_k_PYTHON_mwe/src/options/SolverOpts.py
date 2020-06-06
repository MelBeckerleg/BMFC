import rank_k

class SolverOpts:
    def __init__(self):
        self.tag='TBMC_partition'
        self.solver = rank_k.TBMC
        self.rank_one_solver=rank_k.init_y_partition
        self.p2 = 0
        self.kmax = 0
        #rank 1 options
        self.mu = 1
        self.seed = 'NA'
        #self.rank_1_opts={'mu':1,'seed':'NA'}
