# converts tif => mjpeg or mjpeg => tif
import cv2
import numpy as np
from skimage import io
import matplotlib.pyplot as plt
def main():
    inputtif = "/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/compression/compression/testin.tif"
    mj2file = '/groups/mousebrainmicro/home/base/CODE/MOUSELIGHT/compression/compression/output5.avi'
    tif2mj2(inputtif,mj2file)

def tif2mj2(inputtif,mj2file):
    im = io.imread(inputtif)
    # fourcc = cv2.VideoWriter_fourcc(*"MJPG")
    fourcc = cv2.VideoWriter_fourcc(*'WMV1')
    print(fourcc)
    out = cv2.VideoWriter(mj2file, fourcc, 20.0, (im.shape[2], im.shape[1]))
    tr=np.zeros((im.shape[1], im.shape[2],3), dtype="uint8")
    print im.shape
    # convert to mj2
    for idx in range(im.shape[0]):
        imidx = im[idx,:,:]/256
        # print 'working on: '+str(idx)
        tr[:,:,0] = imidx
        tr[:, :, 1] = imidx
        tr[:, :, 2] = imidx
        out.write(tr)
    out.release()



if __name__ == '__main__':
   main()