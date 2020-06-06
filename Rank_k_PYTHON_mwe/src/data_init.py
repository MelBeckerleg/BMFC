import numpy as np
import copy
import csv
import pandas as pd
import scipy.sparse

class Data:
    def __init__(self, args1): #args1=generate_opts; args2=sub_sample_opts
        for key, vals in args1.iteritems():
            setattr(self, key, vals)

    def generate(self):
        self.data = self.generator_function(self)
        if self.transpose_generated:
            self.data= np.transpose(self.data)

    def generate_root(self,fn_root):
        data = self.generator_root_function(self)
        if self.transpose_generate_root:
                data= np.transpose(data)
        self.root = temp.values

    def sub_sample(self, data_values):
        self.sub, self.test,self.train = self.subsampling_function(self,data_values)
        #DATA.subsampling_function(M,DATA.ssargs) works!!!!


#####Generator functions

def generate_k_planted_tiles(self):
    n=self.n;m=self.m;noise=self.noise;kmax=self.kmax;ratio=self.ratio
    try:
        randseed=self.randseed
        np.random.seed(randseed)
    except(AttributeError):
        print('no rand seed set')
    M = 1. * np.zeros([m, n])
    unit_size = n * m / (kmax)*ratio#######want a negative skew!
    k=0
    while k<kmax:
        x = np.zeros([1, m])[0]
        y = np.zeros([1, n])[0]
        y[np.random.permutation(n)[0:int(np.ceil(np.sqrt(unit_size)))]] = 1
        x[np.random.permutation(m)[0:int(np.ceil(np.sqrt(unit_size) ))]] = 1
        M = M + np.outer(x, y)
        k=k+1
    M = 2. * (M > 0) - 1
    noises=np.random.rand(m,n)
    M[np.where(noises>1-noise)]=-M[np.where(noises>1-noise)]
    del noises
    return M

def generate_block_diagonal(self):
    n=self.n;m=self.m;noise=self.noise;kmax=self.kmax;ratio=self.ratio
    row_clusters=np.random.permutation(m)-1;U=np.zeros([m,kmax]);
    col_clusters=np.random.permutation(n)-1;V=np.zeros([n,kmax]);
    for i in range(kmax):
        U[row_clusters[int(i*m/kmax):int((i+1)*m/kmax)],i] =1
        V[row_clusters[int(i*n/kmax):int((i+1)*n/kmax)],i] =1
    M=2.*np.dot(U,V.T)-1.
    noises=np.random.rand(m,n)
    M[np.where(noises>1-noise)]=-M[np.where(noises>1-noise)]
    del noises
    return M

def generate_row_clusters(self,keep_factors=1):
    n=self.n;m=self.m;noise=self.noise;kmax=self.kmax;ratio=self.ratio
    row_clusters=np.random.permutation(m)-1;U=np.zeros([m,kmax]);
    for i in range(kmax):
        U[row_clusters[int(i*m/kmax):int((i+1)*m/kmax)],i] =1
    Vt=np.reshape([(np.random.rand(n)>(1-ratio))*1. for _ in range(kmax)],[kmax,n])
    M=2.*np.dot(U,Vt)-1.
    noises=np.random.rand(m,n)
    M[np.where(noises>1-noise)]=-M[np.where(noises>1-noise)]
    if keep_factors:
        self.X_orig=U;self.Y_orig=Vt.T
    del noises
    return M

def generate_checkerboard(self):
    n=self.n;m=self.m;noise=self.noise;kmax=self.kmax;ratio=self.ratio
    try:
        randseed=self.randseed
        np.random.seed(randseed)
    except(AttributeError):
        print('no rand seed set')
    M = 1. * np.zeros([m, n])
    row_clusters=np.random.permutation(m)
    col_clusters=np.random.permutation(n)
    for i in range(kmax):
        for j in range(kmax):
            M[int(i*m/kmax):int((i+1.)*m/kmax),int(j*n/kmax):int((j+1.)*n/kmax)]=1.*(np.random.rand(1)>(1-ratio))
    M=M[row_clusters][col_clusters]#retrieve related clusters with [np.ravel(np.where(row_clusters==i)) for i in range(len(row_clusters))]
    M = 2. * (M > 0) - 1
    noises=np.random.rand(m,n)
    M[np.where(noises>1-noise)]=-M[np.where(noises>1-noise)]
    return M


def generate_random(self):
    n=self.n;m=self.m;ratio=self.ratio;
    try:
        randseed=self.randseed
        np.random.seed(randseed)
    except(NameError):
        print('no rand seed set')
    data = 2 * (np.random.rand(m, n) > ratio) * 1. - 1.
    return data
def load_chembl(self):
    #fn = '/home/user/Documents/Ethera/MP2/Data/Chembl/denseSubset.50k.chemblonly.csv'
    #vals=[]
    fout =open('/home/beckerleg/Ethera/MP2/Data/read_ch/mats/ch_mat.txt', 'r')
    #fout =open('/home/user/Documents/Ethera/FirstYear/DataTrail/ch_mat_vals.txt', 'r')
    vals=np.transpose([np.ravel(rows.strip('\n').split('  ')).astype(int) for rows in fout])
    fout.close()
    idx=np.where(pd.DataFrame(vals[0]).isin(np.where(np.histogram(vals,bins=vals[0].max()+1)[0]>20)[0]))[0]
    data=scipy.sparse.coo_matrix((vals[2][idx]*(-1.),(pd.factorize(vals[0][idx])[0],pd.factorize(vals[1][idx])[0])),shape=(np.max(pd.factorize(vals[0][idx])[0])+1,np.max(pd.factorize(vals[1][idx])[0])+1))
    return np.asarray(data.todense())*1.

def load_RC(self):
    fn='/home/beckerleg/Ethera/Datasets/RCdata/rating_final.csv'
    tmp=pd.read_csv(fn)
    x=tmp[u'userID'].factorize()[0];y=tmp[u'placeID'].factorize()[0];z=(tmp[u'rating'].factorize()[0]>1)*2.-1
    data=scipy.sparse.coo_matrix((z,(x,y)))
    return np.asarray(data.todense())

def load_gene_expression(self):
    #target_url=self.fn
    target_url = ('https://www.ebi.ac.uk/biostudies/files/S-EPMC384712/pnas_101_12_4164__08531DataSet1.txt')
    df = pd.read_csv(target_url,delimiter='\t')
    df=(df-df.mean())/df.std()
    df=(df>0)*2.-1.
    return np.asarray(df)


def load_fn_preds(self):
    fn=self.fn
    ##fn_orig='/home/user/Documents/Ethera/FirstYear/Clusters/Results/SQallp956c'
    ##fn_root= '/home/user/Documents/Ethera/FirstYear/DataTrail/SQ956.csv'
    with open(fn, 'rb') as csvfile:
        temp = pd.read_csv(csvfile, sep=',', header=None)
    return temp

def load_fn_bin(obj):
    fn=obj.fn
    ##fn_orig='/home/user/Documents/Ethera/FirstYear/Clusters/Results/SQallp956c'
    ##fn_root= '/home/user/Documents/Ethera/FirstYear/DataTrail/SQ956.csv'
    with open(fn, 'rb') as csvfile:
        temp = pd.read_csv(csvfile, sep=',', header=None)
        data = 2. * (temp.values >= 0.5) - 1
    return data
def load_movie_lens(self):
    ratings=pd.read_csv('/home/beckerleg/Ethera/Datasets/ml-100k/u.data',sep='\t',header=None)
    ratings=ratings.rename(columns= {0:'userId',1:'movieId',2:'rating',3:'timestamp'})
    ratings['rating']=(ratings['rating']>4)*2.-1 #could
    data=scipy.sparse.coo_matrix((ratings['rating'],(pd.factorize(ratings['userId'])[0],pd.factorize(ratings['movieId'])[0])),shape=(len(np.unique(ratings['userId'])),len(np.unique(ratings['movieId']))))
    try:
        np.random.seed(self.randseed)
    except:
        print('no rand seed set')
    try:
        rho=self.frac
        data=data.tolil()[np.ravel(np.random.random_integers(len(np.unique(ratings['userId']))-1,size=int(np.floor(rho*len(np.unique(ratings['userId'])))))),:]
        data=data.tolil()[:,np.ravel(np.random.random_integers(len(np.unique(ratings['movieId']))-1,size=int(np.floor(rho*len(np.unique(ratings['movieId']))))))]
    except:
        print('using full ML data')
    data=np.asarray(data.todense())
    data=data[[i for i,x in enumerate(data) if np.any(x)],:]
    data=data[:,[i for i,x in enumerate(np.transpose(data)) if np.any(x)]]
    #data=np.transpose(data)
    return data

def download_target(self):
    target_url=self.fn
    target_url = ('https://archive.ics.uci.edu/ml/machine-learning-databases/00233/CNAE-9.data')
    data = pd.read_csv(target_url)

    #https://towardsdatascience.com/machine-learning-nlp-text-classification-using-scikit-learn-python-and-nltk-c52b92a7c73a

##########Subsampling functions
def generate_sub(self,data_values):
    rho=self.rho;
    try:
        randseed=self.randseed
        np.random.seed(randseed)
    except(AttributeError):
        print('no rand seed set')
    I, J = np.where(data_values)
    dp = len(I)
    rnd = np.random.permutation(dp)
    test = rnd[int(np.ceil(rho * dp)):]
    train=rnd[:int(np.ceil(rho*dp))]
    sub = copy.deepcopy(data_values)
    sub[I[test], J[test]] = 0
    return sub, test,train


def load_sub(self,data_values):
    fn_sub=self.fn_sub;rho=self.rho;
    I, J = np.where(data_values)
    rnd = np.asarray(np.load(fn_sub))
    test = rnd[int(np.ceil(rho * (len(rnd)))):]
    train=rnd[:int(np.ceil(rho*len(rnd)))]
    sub = copy.deepcopy(data_values)
    sub[I[test], J[test]] = 0
    return sub, test,train


