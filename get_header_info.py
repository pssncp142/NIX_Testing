###############################################################################
# Yigit Dallilar 14/10/2019
# MPE
#
# List selected header values of given fits files
###############################################################################

from astropy.io import fits
import re
import os
from NIX_Testing import NIX_Image_List

keywords = ['MJD-OBS', 'HIERARCH ESO DET READ CURNAME']

data_dir = '/home/ydallilar/Documents/NIX/nixDetBackup'
test_ids = ['PER-085-05-2', 'PER-085-05-3']
tbl_fmt = '%20s,%50s,%20s,%40s'

table_header =  ['TEST_ID', 'FILENAME']

config = {'test_ids' : test_ids, 'data_dir' : data_dir}

if __name__=='__main__':

    NIX_Data = NIX_Image_List(config)
    NIX_Data.printTable(keywords, tbl_fmt)


