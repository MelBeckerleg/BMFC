import data_init as di

class GenerateOpts:
    def __init__(self):
        self.generate_function=di.generate_row_clusters
        self.m = 500
        self.n = 500
        self.ratio = 0.3
        self.noise = 0.03
        self.kmax = 20
        self.rho=0.7 #subsampling factor
        self.subsample_function=di.generate_sub


