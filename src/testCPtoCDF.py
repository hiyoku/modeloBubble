import cPickle as cp
import BubbleModelStorage as BMS

denList = cp.load(open('v2', 'r'))
BMS.create_netcdf(denList, name="dataset_semt")
