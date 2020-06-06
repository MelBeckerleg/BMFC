import copy
import numpy as np

#Rank_k_solve class for generating a rank k binary factorisation of mxk database
#A=XY' where X is mxk binary matrux and Y is nxk binary matrix,
# Author: Mel Beckerleg
# Date 10/2019

class Rank_k_solve:
    def __init__(self, solve_opts):
        for key,vals in solve_opts.iteritems():
            setattr(self,key,vals)
    def solve(self,*M):
        [X, Y] = self.solver(self,*M)
        m,k=np.shape(X)
        try:
            tile_size=[sum(X)[i]*sum(Y)[i] for i in xrange(k)][::-1] #flipped order
            X=X[:,np.argsort(tile_size)];Y=Y[:,np.argsort(tile_size)]
        except(IndexError):
            print('not sorted')
        return X, Y

###### options for rank_k_solve
    '''solve_dict = {'best': rank_k.init_noisyblock, 'lp': rank_k.init_y_lp,
                  'avg': rank_k.init_y_avg, 'partition': rank_k.init_y_partition, 'CGIHT': rank_k.init_y_CGIHT,
                  'svd': rank_k.init_y_svd,'lp_weighted':rank_k.init_y_lpweighted}
    '''

########
def stop_at_1(opts,M):
    #perform a rank-one binary approximation
    fun1 = opts.rank_one_solver
    x,y = fun1(M,opts.rank_one_opts)
    Xfinal=np.reshape([x,np.zeros(len(x))],[2,len(x)])
    Yfinal=np.reshape([y,np.zeros(len(y))],[2,len(y)])
    return Xfinal, Yfinal
########
def load_solved(opts):
    transpose = opts.rank_k_opts['transpose_solve']
    fn = opts.load_fn
    temp = np.load(fn)
    if transpose:
        Xfinal=temp['arr_1']
        Yfinal=temp['arr_0']
    else:
        Xfinal=temp['arr_0']
        Yfinal=temp['arr_1']
    return Xfinal,Yfinal
########
def iterative_rank1(opts,M):
    fun1=opts.rank_one_solver
    x,y=fun1(M,opts.rank_one_opts)
    x,y=iterative_update(M,x, y, max_iters=20, p2=opts.rank_k_opts['p2'])
    Xfinal=np.reshape([x,np.zeros(len(x))],[2,len(x)]).T
    Yfinal=np.reshape([y,np.zeros(len(y))],[2,len(y)]).T
    return Xfinal, Yfinal
    #return x,y
############### Rank_k fcts
def iterative_update(M,x, y, max_iters=20, p2=1):
    convergence = 0
    epsilon = 0.005
    iter = 0
    while iter < max_iters and convergence != 1 and p2:
        x = np.ravel((np.dot(M, np.transpose(y)) >= 0) * 1.)
        y = np.ravel(((np.dot(np.transpose(x), M)) >= 0) * 1.)
        iter += 1
        if np.count_nonzero(1. * M - (2. * (np.outer(x, y)) - 1.))/ (1. * len(np.where(M))) < epsilon:
            convergence = 1
    return x, y

def partition_check(Msub, xsub, y):
    #check whether to continue to partition
    flag = 0
    tol = 0.05
    i_vals= np.where(sum(np.transpose(Msub)))
    #check agreement
    val=np.max([np.linalg.norm(((Msub[i,:]+1)/2.-y)[np.where(Msub[i,:])])/np.linalg.norm(Msub[i,:]) for i in i_vals])
    if sum(xsub) and sum(1 - xsub) and val > tol:
        flag = 1
    return flag


def save_function(x1, Msub, xsub, ysub, Xfinal, Yfinal, X, Y):
    #update X,Y,Xfinal,Yfinal
    if partition_check(Msub, xsub, ysub):
        if not sum(np.all(X == np.zeros([1, len(x1)]), axis=1)): #check there's room left to store:
            X = np.concatenate((X, np.zeros([5, len(x1)])), axis=0)
            Y = np.concatenate((Y, np.zeros([5, len(ysub)])), axis=0)
        zero_sub = np.where(np.all(X == np.zeros([1, len(x1)]), axis=1))
        first_zero_sub = zero_sub[0][0]
        X[first_zero_sub, :] = copy.deepcopy(x1)
        Y[first_zero_sub, :] = copy.deepcopy(ysub)
    elif sum(xsub): #no split required
        if not sum(np.all(Xfinal == np.zeros([1, len(x1)]), axis=1)):
            Xfinal = np.concatenate((Xfinal, np.zeros([5, len(x1)])), axis=0)
            Yfinal = np.concatenate((Yfinal, np.zeros([5, len(ysub)])), axis=0)
        zero_sub = np.where(np.all(Xfinal == np.zeros([1, len(x1)]), axis=1))
        first_zero_sub = zero_sub[0][0]
        if sum(x1) and sum(ysub):
            Xfinal[first_zero_sub, :] = copy.deepcopy(x1)
            Yfinal[first_zero_sub, :] = copy.deepcopy(ysub)
    return Xfinal, Yfinal, X, Y


def TBMC(opts,A):
#implementation of the algorithm outlined in 
#Beckerleg, M. and Thompson, A., 2020. A divide-and-conquer algorithm for binary matrix completion. Linear Algebra and its Applications.
    try:
        np.random.seed(opts.rand_seed_solve) #
    except(AttributeError):
        print('norandseedsetforTBMC')
    data = A
    fun1 = opts.rank_one_solver
    p2=opts.rank_k_opts['p2']
    max_iters = opts.rank_k_opts['kmax']
    # set_up:
    max_number_of_squares = 2*max_iters
    m,n = np.shape(data)
    X = np.zeros([max_number_of_squares, m]) #working with Xt and Yt
    Y = np.zeros([max_number_of_squares, n])
    Xfinal = np.zeros([max_number_of_squares, m])
    Yfinal = np.zeros([max_number_of_squares, n])
    # put everything in to start with
    X[0, :] = 1.
    Y[0, :] = 1.
    iters=0
    not_converged=1
    while not_converged and iters < max_iters:
        print(iters)
        x = X[0, :]
        data_chunk = data[np.ravel(np.where(x)), :]
        xsub, ysub = fun1(data_chunk,opts.rank_one_opts)
        xsub, ysub = iterative_update(data_chunk, xsub, ysub,p2)
        x1 = copy.deepcopy(x)
        x1[np.ravel(np.where(x))] = xsub
        [Xfinal, Yfinal, X, Y] = save_function(x1, data_chunk, xsub, ysub, Xfinal, Yfinal, X, Y)
        x1[np.ravel(np.where(x))] = 1 - xsub
        [Xfinal, Yfinal, X, Y] = save_function(x1, data_chunk, xsub, ysub, Xfinal, Yfinal, X, Y)
        X = np.delete(X, (0), axis=0)
        Y = np.delete(Y, (0), axis=0)
        tmp1, tmp2 = np.nonzero(X)
        not_converged = (tmp1 + 1).any()
        iters += 1
    if not_converged:
        print('max number of iterations reached')
        Xfinal = np.concatenate((Xfinal, X), axis=0)
        Yfinal = np.concatenate((Yfinal, Y), axis=0)
    return Xfinal.T, Yfinal.T

###########################################################################
def binary_rescale(W,H):
    m,k=np.shape(W);n,k=np.shape(H)
    Dw=np.diag([np.max(W[:,i]) for i in xrange(k)])
    Dwinv=np.diag([1/np.max(W[:,i]) for i in xrange(k)])
    Dh=np.diag([np.max(H[:,i]) for i in xrange(k)])
    Dhinv=np.diag([1/np.max(H[:,i]) for i in xrange(k)])
    W=np.dot(W,np.dot(np.sqrt(Dwinv),np.sqrt(Dh)))
    H=np.dot(H,np.dot(np.sqrt(Dw),np.sqrt(Dhinv)))
    return W,H

def binarise(W,H):
    W=(W>0.5)*1.
    H=(H>0.5)*1.
    return W,H

def NMF_rank_k(gopts,N):
#BMF using NMF, as in Zhang, Z., Li, T., Ding, C. and Zhang, X., 2007, October. Binary matrix factorization with applications. In Seventh IEEE International Conference on Data Mining (ICDM 2007) (pp. 391-400). IEEE.
    import sklearn.decomposition.nmf_missing as nmf
    M=copy.deepcopy(N)
    M[M==0]=np.nan
    M=(M+1.)/2.
    k=gopts.rank_k_opts['kmax']
    cut_off=0.5
    H=H.T
    W,H=binary_rescale(W,H)
    W,H=binarise(W,H)
    rank=len(H)
    if rank<=1:
        W=np.reshape([np.ravel(W),np.zeros(len(W))],[2,len(W)])
        H=np.reshape([np.ravel(H),np.zeros(len(np.transpose(H)))],[2,len(np.transpose(H))])
    return W,H

def NMF_rank_one(M,*empty):
    import sklearn.decomposition.nmf as nmf
    M[M==0]=np.nan
    M=(M+1.)/2.
    cut_off=0.5
    w,h,n_iter=nmf.non_negative_factorization(M,n_components=1,solver='mu')
    #check dimension
    w=np.ravel((w>cut_off)*1.)
    h=np.ravel((h>cut_off)*1.)
    del M
    return w,h
##########################################################################
def pn_project(X):
    X=2.*(X>0)-1
    return X

def sv_thresh(X,gamma):
    m,n=np.shape(X)
    print('type is' +str(type(X)))
    if not str(type(X))=='float':
        X=X.astype(float)
    U,S,V=np.linalg.svd(X,full_matrices=False)
    Smat=np.zeros((m,n))
    Smat[:len(S),:len(S)]=((np.diag(S)-gamma)>0)*1.
    D=np.dot(U*S*((S-gamma)>0),V)
    print(['D is', D])
    return D

def rank_r_approx(A,r):
    #rank r approximation using SVD
    u,s,v = np.linalg.svd(A, full_matrices=False)
    Ar = np.zeros(np.shape(A))
    for i in xrange(r):
        Ar += s[i] * np.outer(u.T[i], v[i])
    return Ar

def approx_cluster(A,r,thresh):
    #cluster based on proximity to randomly selected rows
    m,n=np.shape(A)
    to_cluster=np.ones(m)
    clusters=np.zeros(m)
    vals=np.zeros((r,n))
    iter=0
    while iter<r and np.sum(to_cluster):
        print(iter)
        #rows
        chosen=A[np.random.choice(np.where(to_cluster)[0]),:]
        indices=[i for i in xrange(m) if np.linalg.norm(A[i,:]-chosen,ord=2)<thresh and to_cluster[i]]
        to_cluster[indices]=0
        #columns
        clusters[indices]=iter
        iter=iter+1
    return clusters

def spectral_method(gopts,M):
    #implementation of method outlined in
    #Xu, J., Wu, R., Zhu, K., Hajek, B., Srikant, R. and Ying, L., 2014, June. Jointly clustering rows and columns of binary matrices: Algorithms and trade-offs. In The 2014 ACM international conference on Measurement and modeling of computer systems (pp. 29-41).
    m,n=np.shape(M)
    k=gopts.rank_k_opts['kmax']
    #erasure probability
    epsilon=1.-np.count_nonzero(M)/(1.*m*n)
    #threshold for distance from user.
    thresh=np.sqrt((1-epsilon))*k*np.log(m)*12
    #assign subsets
    #probability of assignment
    delta=(1-epsilon)/4.
    omega=np.where(M)
    rho_approx=1.*len(omega[0])/m/n
    assign=np.random.random(len(omega[0]))
    omega_one=np.where(assign<(1/2.))
    omega_two=np.where((assign>(1/2.-delta))*(assign<(1.-delta)))
    R_one=np.zeros((m,n));R_one[omega[0][omega_one],omega[1][omega_one]]=M[omega[0][omega_one],omega[1][omega_one]]
    R_two=np.zeros((m,n));R_two[omega[0][omega_two],omega[1][omega_two]]=M[omega[0][omega_two],omega[1][omega_two]]
    #cluster using rank one approximations
    PR_one=1.*(rank_r_approx(R_one,k)>0.0)
    PR_two=1.*(rank_r_approx(R_two,k)>0.0)
    #cluster
    row_clusters=approx_cluster(PR_one,k,thresh)
    col_clusters=approx_cluster(np.transpose(PR_one),k,thresh)
    #block assign values
    block_vals=[[(1/(rho_approx*len(np.where(col_clusters==cval)[0])*len(np.where(row_clusters==rval)[0]))*np.sum(PR_two[np.where(row_clusters==rval)[0],:][:,np.where(col_clusters==cval)[0]])>0.5)*1. for cval in xrange(k) ]for rval in xrange(k)]
    # recluster rows.
    Y=np.asarray([[block_vals[rval][int(col_clusters[j])] for j in xrange(n)] for rval in xrange(k)])
    #distance of each row from each of the centres
    vals=np.asarray([np.linalg.norm((PR_two-Y[cval,:]),axis=1) for cval in xrange(k)])
    #reassigned clusters
    row_clusters=np.asarray([np.where(vals.T[i]==vals.min(axis=0)[i])[0][0] for i in xrange(m)])
    #return X,Y give row and oclumn membership of clusters
    X=np.asarray([1.*(row_clusters==i+1) for i in xrange(k)])
    Y=1.*(Y>0)
    return X.T,Y.T




def convex_method(gopts,M):
    #implementation of method outlined in
    #Xu, J., Wu, R., Zhu, K., Hajek, B., Srikant, R. and Ying, L., 2014, June. Jointly clustering rows and columns of binary matrices: Algorithms and trade-offs. In The 2014 ACM international conference on Measurement and modeling of computer systems (pp. 29-41).
    k=gopts.rank_k_opts['kmax']
    m,n=np.shape(M)

    epsilon=1.-np.count_nonzero(M)/(1.*m*n)
    print(epsilon)
    penalty=3.*np.sqrt((1-epsilon)*np.sqrt(n*m))
    print(penalty)
    penalty_min=np.sqrt((1-epsilon)*np.sqrt(n*m))
    print(penalty_min)
    mu=1;Y_c=0;Y_p=0;alpha_c=1;alpha_p=1;num_iters=5
    gamma=0.1
    for iter in xrange(num_iters):
        Z_c=Y_c+((alpha_p-1)/alpha_c)*(Y_c-Y_p)
        print(Z_c)
        Y_p=Y_c
        Y_c=pn_project(sv_thresh(Z_c+M/mu,penalty/mu))
        print(['Y_c is',Y_c])
        alpha_p=alpha_c
        alpha_c=(1+np.sqrt(1+4*pow(alpha_p,2)))/2.
        penalty=max(gamma*penalty,penalty_min)
    #binarise?
    Y_c=(Y_c>0)*1.
    print(Y_c)
    C=np.unique(Y_c,axis=0)
    vals=[np.all(np.equal(row,C),axis=1) for row in Y_c]
    vals=np.transpose(np.reshape(vals,np.shape(vals))*1.)
    return vals,C
#########################################################################
# intialisation functions# to use used as fun1, must return x,y

def glpk_tile_ip(M,*empty):
    import cvxopt
    import cvxopt.glpk
    import numpy as np
    import scipy.sparse as sps
    n, m = np.shape(M)
    numzspos = len(np.nonzero(M == 1.)[0])
    numzsneg = len(np.nonzero(M == -1.)[0])
    ##  Create constraint matrix
    colpos = np.asarray(np.concatenate((np.linspace(0, numzspos - 1, numzspos),
                                        np.where(M == 1)[0] + numzspos + numzsneg,
                                        np.where(M == 1)[1] + numzspos + numzsneg + n),
                                       axis=0))  # revisit this to make it more streamlined!
    colneg = np.asarray(np.concatenate((np.linspace(numzspos, numzspos + numzsneg - 1, numzsneg),
                                        np.where(M == -1)[0] + numzspos + numzsneg,
                                        np.where(M == -1)[1] + numzspos + numzsneg + n),
                                       axis=0))  # revisit this to make it more streamlined!
    # colupper=np.asarray([i for i in xrange(n+m+numzspos+numzsneg)])
    # collower=np.asarray([i for i in xrange(n+m+numzspos+numzsneg)])
    # col=np.concatenate((colpos,colneg,colupper,collower),axis=0)
    col = np.concatenate((colpos, colneg), axis=0)
    rowpos = np.asarray(np.concatenate((np.linspace(0, numzspos - 1, numzspos), np.linspace(0, numzspos - 1, numzspos),
                                        np.linspace(0, numzspos - 1, numzspos)), axis=0))
    rowneg = np.asarray(np.concatenate((np.linspace(numzspos, numzspos + numzsneg - 1, numzsneg),
                                        np.linspace(numzspos, numzspos + numzsneg - 1, numzsneg),
                                        np.linspace(numzspos, numzspos + numzsneg - 1, numzsneg)), axis=0))
    # rowupper=np.asarray([i+numzspos+numzsneg for i in xrange(n+m+numzspos+numzsneg)])
    # rowlower=np.asarray([i+n+m+2*numzspos+2*numzsneg for i in xrange(n+m+numzspos+numzsneg)])
    # row=np.asarray(np.concatenate((rowpos,rowneg,rowupper,rowlower),axis=0))
    row = np.asarray(np.concatenate((rowpos, rowneg), axis=0))
    # values=np.concatenate((2*np.ones([numzspos,]),-np.ones([2*numzspos,]),-np.ones([numzsneg,]),np.ones([2*numzsneg,]),np.ones([n+m+numzspos+numzsneg,]),-np.ones([n+m+numzspos+numzsneg,])),axis=0)
    values = np.concatenate(
        (2 * np.ones([numzspos, ]), -np.ones([2 * numzspos, ]), -np.ones([numzsneg, ]), np.ones([2 * numzsneg, ])),
        axis=0)
    # Q=sps.coo_matrix((values,(row,col)),shape=(numzspos+numzsneg+2*(numzsneg+numzspos+n+m),(numzsneg+numzspos+n+m)))
    Q = sps.coo_matrix((values, (row, col)), shape=(numzspos + numzsneg, (numzsneg + numzspos + n + m)))
    Q2 = cvxopt.spmatrix(Q.data.tolist(), Q.row.tolist(), Q.col.tolist(), size=Q.shape)
    # c=np.concatenate((np.zeros([numzspos,]),np.ones([numzsneg,]),np.ones([n+m+numzspos+numzsneg,]),np.zeros([n+m+numzspos+numzsneg,])))
    c = np.concatenate((np.zeros([numzspos, ]), np.ones([numzsneg, ])))
    f = np.concatenate((np.ones([numzspos, ]), -np.ones([numzsneg, ]), np.zeros([n + m, ])), axis=0)
    # format
    f = cvxopt.matrix(-f, tc='d')
    G = cvxopt.matrix(Q2, tc='d')
    h = cvxopt.matrix(c, tc='d')
    # solve
    (status, sln) = cvxopt.glpk.ilp(f, G, h, B=set([i for i in xrange(n + m + numzspos + numzsneg)]))
    x = np.ravel(np.asarray(sln))[numzspos + numzsneg:numzspos + numzsneg + n]
    y = np.ravel(np.asarray(sln))[numzspos + numzsneg + n:]
    z = np.ravel(np.asarray(sln))[:numzspos + numzsneg]
    # dont want x,y empty
    if not sum(x > 0):
        x = np.zeros([n, ])
    if not sum(y > 0):
        y = np.zeros([m, ])
    return z, x, y



def cylp_tile_lp(M,mu=1,*empty):
    from cylp.cy.CyClpSimplex import CyClpSimplex
    from cylp.py.modeling.CyLPModel import CyLPArray
    import numpy as np
    import scipy.sparse as sps
    ###########
    # np.
    # delete when as function
    # A=1.*(np.random.rand(1500,5000)>0.5)
    # A=np.zeros((5,5))
    # A[0:100,0:300]=1
    # A[0:2,0:2]=1
    # M=2*A-1
    ###
    n, m = np.shape(M)
    numzs = len(np.nonzero(M == -1.)[0])
    ########
    ##  Create constraint matrix
    col = np.asarray(np.concatenate(
        (np.linspace(0, numzs - 1, numzs), np.where(M == -1)[0] + numzs, np.where(M == -1)[1] + numzs + n), axis=0))
    row = np.asarray(np.concatenate(
        (np.linspace(0, numzs - 1, numzs), np.linspace(0, numzs - 1, numzs), np.linspace(0, numzs - 1, numzs)), axis=0))
    values = np.concatenate((-np.ones([numzs, ]), np.ones([2 * numzs, ])), axis=0)
    Q = sps.coo_matrix((values, (row, col)), shape=(numzs, numzs + n + m))
    ## constraint vector
    c = np.ones([numzs, ])
    ## objective vector
    f = np.concatenate(
        (-mu*np.ones([numzs, ]), .5 * sum(np.transpose(M) == 1) * np.ones([n, ]), .5 * sum(M == 1) * np.ones([m, ])),
        axis=0)
    #print("set up")
    #########
    ## Create simplex model
    s = CyClpSimplex()  ###don't double define on the clusters, size allocation error
    #print("simplex")
    xyz = s.addVariable('xyz', n + m + numzs)
    s += Q * xyz <= c
    s += 0 <= xyz <= 1
    s.objective = CyLPArray(-f) * xyz
    #print("s built")
    #########j
    ## Solve
    # s.primal()
    s.initialSolve()
    #print("s solve")
    ###########
    ## Recover tile
    output = s.primalVariableSolution["xyz"]
    z = output[0:numzs]
    x = output[numzs:numzs + n]
    y = output[numzs + n:]
    #print('shape is')
    #print(np.shape(y))
    del s,col,row,Q,values,output
    return z, x, y



def init_y_partition(a, *empty):
    (n, m) = np.shape(a)
    select = int(np.floor(np.random.rand(1) * m))  # separator column
    # take the average of all rows that have a positive entry in 'select'
    y = 1. * np.ravel(1. * (np.mean(np.transpose(a[np.where((a[:, select]) > 0), :]), axis=1)) > 0)
    x = 1. * np.ravel(a[:, select] > 0)
    return x, y


def init_y_avg(a, *empty):
    m,n = np.shape(a)
    y = 1. * (sum(a) / (1. * m) > 0.0)
    x = np.zeros([1, m])[0]
    x[np.where((sum(np.transpose(a[:, np.where(y)[0]])) > 0) * 1.)] = 1
    return x, y



def init_y_ip(M, *empty):
    z, x, y = glpk_tile_ip(M)
    return x, y

def init_y_lpweighted(M,opts, *empty):
    mu=opts['mu']
    z, x, y = cylp_tile_lp(M,mu)
    return x, y




