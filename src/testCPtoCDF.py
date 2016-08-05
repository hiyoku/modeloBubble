import cPickle as cp
import den

denList = cp.load(open('output/v0', 'r'))
den.create_netcdf(denList)
