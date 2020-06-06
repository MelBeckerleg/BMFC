import numpy as np


class Tiling_error:
    def __init__(self, DATA, X, Y, error_opts):  # opts:
        for key,vals in DATA.__dict__.iteritems():
            setattr(self, key, vals)
        self.X = X
        self.Y = Y
        m,k=np.shape(X)
        self.tile_size=[sum(X)[i]*sum(Y)[i] for i in xrange(k)]
        self.numtiles=np.count_nonzero(self.tile_size)
        for key, vals in error_opts.iteritems():
            setattr(self, key, vals)

    def tiled_error(self):
        #err = [[metric(self.data[self.idxfun(self, k)],
         #      np.dot(np.transpose(self.X[self.row0(k):k + 1]), self.Y[self.row0(k):k + 1])[self.idxfun(self, k)]) for k
          #     in self.numtilesrange(self.numtiles)] for metric,idxfun in self.metrics]
        if self.numtiles:
            err = [self.metric(self.data[self.idxfun(self, k)],
                   2.*np.dot(self.X[:,self.row0(k):k + 1], self.Y[:,self.row0(k):k + 1].T)[self.idxfun(self, k)]-1.,self,k) for k
                   in xrange(len(self.tile_size))]
        else:
            err = [self.metric(self.data[self.idxfun(self, 0)],
                   2.*np.dot(self.X[:,self.row0(0): 1], self.Y[:,self.row0(0): 1].T)[self.idxfun(self, 0)]-1.,self,0)]

        return err
###########

###########options for row0
def bytile(k):  # row0 for single tile
    return k


def alltiles(k):  # row0 for all tiles/ numtilesrange for justfirsttile
    return 0

##########options for numtilesrange
def xrangenumtiles(numtiles):
    return xrange(numtiles)


def numtiles(numtiles):
    return [numtiles]


############options for idxfun
def idxfun_sample(obj,*k):
    I, J = np.where(obj.data)
    idxI = I[obj.test]
    idxJ = J[obj.test]
    return idxI, idxJ

def idxfun_all(obj,*k):
    I, J = np.where(obj.data)
    return I,J

def idxfun_train(obj,*k):
    I, J = np.where(obj.data)
    idxI = I[obj.train]
    idxJ = J[obj.train]
    return idxI,idxJ

def idxfun_tile_only(obj,k):
    Iidx,Jidx=obj.tile_only_sub(obj)
    #print('Xwhere is', np.where(np.ravel(obj.X[k:k + 1, :])))
    idxI = np.isin(np.ravel(Iidx),np.where(np.ravel(obj.X[k:k + 1, :]))[0])
    idxJ = np.isin(np.ravel(Jidx),np.where(np.ravel(obj.Y[k:k + 1, :]))[0])
    sel = np.where(np.multiply(idxI, idxJ))[0]
    Iidx = Iidx[sel]
    Jidx = Jidx[sel]
    return Iidx, Jidx

######################## options for error types

def tfpn_alt(actual, predicted,*obj):
    #TN = np.count_nonzero([(np.ravel(actual) > 0)[i] == (np.ravel(predicted) > 0)[i] for i in xrange(len(np.ravel(predicted)))])
    #FN = np.count_nonzero([(np.ravel(actual) < 0)[i] == (np.ravel(predicted) > 0)[i] for i in xrange(len(np.ravel(predicted)))])
    #FP = np.count_nonzero([(np.ravel(actual) > 0)[i] == (np.ravel(predicted) < 0)[i] for i in xrange(len(np.ravel(predicted)))])
    #TP = np.count_nonzero([(np.ravel(actual) < 0)[i] == (np.ravel(predicted) < 0)[i] for i in xrange(len(np.ravel(predicted)))])
    actual=-1.*actual;predicted=-1.*predicted
    TP = np.count_nonzero((np.ravel(actual) > 0) == (np.ravel(predicted) > 0))
    FP = np.count_nonzero((np.ravel(actual) < 0) == (np.ravel(predicted) > 0))
    FN = np.count_nonzero((np.ravel(actual) > 0) == (np.ravel(predicted) < 0))
    TN = np.count_nonzero((np.ravel(actual) < 0) == (np.ravel(predicted) < 0))
    if (2. * TP + FN + FP):
        fscore = 2. * TP / (2. * TP + FN + FP)
    elif TN:
        fscore=np.nan
        print('Fscore undefined for single class')
    else:
        fscore=np.nan
        print('No values to assess Fscore')
    return fscore

def tfpn(actual, predicted,*obj):
    TP = np.count_nonzero((np.ravel(actual) > 0) == (np.ravel(predicted) > 0))
    FP = np.count_nonzero((np.ravel(actual) < 0) == (np.ravel(predicted) > 0))
    FN = np.count_nonzero((np.ravel(actual) > 0) == (np.ravel(predicted) < 0))
    TN = np.count_nonzero((np.ravel(actual) < 0) == (np.ravel(predicted) < 0))
    if (2. * TP + FN + FP):
        fscore = 2. * TP / (2. * TP + FN + FP)
    elif TN:
        fscore=np.nan
        print('Fscore undefined for single class')
    else:
        fscore=np.nan
        print('No values to assess Fscore')
    return fscore#, TP, FP, FN, TN


def error_l2(actual, predicted,*obj):
    err = np.linalg.norm(actual - predicted)

    return err


def error_l1(actual, predicted,*vars):
    #err = np.count_nonzero(actual - predicted)
    err=sum(abs(actual - predicted))
    return err


def error_l0_percent(actual, predicted,*vars):
    err = np.divide(np.count_nonzero(actual - predicted), len(actual)*1.)

    return err


def error_l0(actual, predicted,*vars):
    err = np.count_nonzero(actual - predicted)
    return err


def density(actual, predicted,obj,k):# requires tile_only for idxfun to make anysense
    tile_size=obj.tile_size[k]*1.
    density_pos = np.divide(np.count_nonzero(actual > 0), tile_size)
    density_neg = np.divide(np.count_nonzero(actual < 0), tile_size)
    density_scale = 1.*len(actual)
    return density_pos#,density_neg

def proportion(actual, predicted,*vars):# requires tile_only for idxfun to make anysense
    proportion_pos = np.divide(np.count_nonzero(actual > 0), 1.*len(actual))
    proportion_neg = np.divide(np.count_nonzero(actual < 0), 1.*len(actual))
    density_scale = 1.*len(actual)
    return proportion_pos#,proportion_neg

def capture_error(metrics_opts):
    for metric in metrics_opts.metrics:
        val = []
        fn = metric_opts.metrics()
        val.append(fn(opts))

def tile_size(actual,predicted,obj,k):
    tile_size=obj.tile_size[k]
    return tile_size
# fairly_standard=tiling_error({IDX},{METRICS})

def tile_support(actual,predicted,obj,k):
    tile_support=np.count_nonzero(actual > 0)
    return tile_support

def tile_pndensity(actual,predicted,obj,k):
    tile_size=obj.tile_size[k]*1.
    tile_support=np.divide(np.count_nonzero(actual),tile_size)
    return tile_support
