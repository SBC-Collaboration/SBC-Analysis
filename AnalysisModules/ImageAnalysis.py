import os
import re
import cv2
import numpy as np
import numpy.matlib as mat
# from matplotlib import pyplot as plt
import time
import collections
# import ipdb

class Bubble(object):

    def __init__(self):
        self.bubbles = collections.OrderedDict(
            runid=[],
            ev=[],
            cam=[],
            trigFrame=[],
            frame=[],
            iPostTrig=[],
            ibubimage=[],
            nbubimage=[],
            ipix=[],
            jpix=[],
            xc=[],
            yc=[],
            dx=[],
            dy=[],
            theta=[],
            mass=[],
            ADCThresh=[]
        )

    def emptyBubble(self):
        for key in self.bubbles.keys():
            self.bubbles[key].append(0)

    def addRunEV(self, runid, ev):
        self.bubbles['runid'] = [runid] * len(self.bubbles['cam'])
        self.bubbles['ev'] = [ev] * len(self.bubbles['cam'])

    def convertDType(self):
        self.bubbles['runid'] = np.int32(self.bubbles['runid'])
        self.bubbles['ev'] = np.int32(self.bubbles['ev'])
        self.bubbles['cam'] = np.uint8(self.bubbles['cam'])
        self.bubbles['trigFrame'] = np.uint8(self.bubbles['trigFrame'])
        self.bubbles['iPostTrig'] = np.uint8(self.bubbles['iPostTrig'])
        self.bubbles['frame'] = np.uint8(self.bubbles['frame'])
        self.bubbles['ibubimage'] = np.uint8(self.bubbles['ibubimage'])
        self.bubbles['nbubimage'] = np.uint8(self.bubbles['nbubimage'])
        self.bubbles['ipix'] = np.float32(self.bubbles['ipix'])
        self.bubbles['jpix'] = np.float32(self.bubbles['jpix'])
        self.bubbles['xc'] = np.float32(self.bubbles['xc'])
        self.bubbles['yc'] = np.float32(self.bubbles['yc'])
        self.bubbles['dx'] = np.float32(self.bubbles['dx'])
        self.bubbles['dy'] = np.float32(self.bubbles['dy'])
        self.bubbles['theta'] = np.float32(self.bubbles['theta'])
        self.bubbles['mass'] = np.float32(self.bubbles['mass'])
        self.bubbles['ADCThresh'] = np.float32(self.bubbles['ADCThresh'])


def BuildImageList(runDir, ev):
    out = []
    # all bmp files
    evDir = os.path.join(runDir, str(ev))
    if os.path.exists(evDir):
        imageList = list(filter(lambda x: re.search('cam.*\d+.*image.*\d+.*\.bmp', x),
                               os.listdir(evDir)))
    else:
        return out

    # extract camera index and frame index
    camIx = np.zeros(len(imageList), dtype=np.int)
    frameIx = np.zeros(len(imageList), dtype=np.int)

    for i in range(0, len(imageList)):
        m = re.match(r"cam.*?(\d+).*?image.*?(\d+).*?\.bmp", imageList[i])
        camIx[i] = np.int(m.group(1))
        frameIx[i] = np.int(m.group(2))

    # sort cameras
    imageIx = np.argsort(camIx)
    frameIx = frameIx[imageIx]
    imageList = [imageList[x] for x in imageIx]

    # sort frames for each camera
    camIds = np.unique(camIx)
    k = 0
    for i in range(0, len(camIds)):
        slice = np.arange(k, k + np.sum(camIx == camIds[i]))
        out.append([os.path.join(evDir, imageList[x]) for x in slice[np.argsort(frameIx[slice])]])
        k += slice.shape[0]

    return out


# de-noise the images and compute the mean and std for pretrigger frames
def ComputeRefFrame(imageList, nPreTrigger):
    # read images and filter with a gaussian filter, and convert to 3D array for each camera
    img = [np.dstack(
        tuple([cv2.GaussianBlur(cv2.imread(imageList[iCam][iFrame], cv2.IMREAD_GRAYSCALE), (7, 7), 3, 3) for
               iFrame in
               range(0, np.amin([len(imageList[iCam]), nPreTrigger]))]))
        for iCam in range(0, len(imageList))]

    # plt.imshow(img[0][0])
    # plt.show()

    # mean and std of pretrigger frames for each camera
    meanFrames = [np.mean(cam, axis=2) for cam in img]
    stdFrames = [np.std(cam, axis=2) for cam in img]

    return meanFrames, stdFrames


def FindBubbles(imageList, nPreTrigger, nBubbleFrame, meanFrames, stdFrames, ADCThresh, nBubbleSize):
    # Load non-pretrigger frames
    img = [np.dstack(
        tuple([cv2.GaussianBlur(cv2.imread(imageList[iCam][iFrame], cv2.IMREAD_GRAYSCALE), (7, 7), 3, 3) for
               iFrame in
               range(nPreTrigger, len(imageList[iCam]))]))
        for iCam in range(0, len(imageList))]

    # take the difference between a frame and the reference frame
    img = [(cam - refFrame[:, :, None]) for cam, refFrame in zip(img, meanFrames)]

    # minimum value in frame diff
    #- minDiff = [np.amin(cam, axis=(0, 1)) for cam in img]
    # shift frame diff
    #- img = [(cam-x0) for cam, x0 in zip(img, minDiff)]

    # thresholding the frame diff, convert to bw image
    bwImg = [np.uint8(cam <= -ADCThresh) for cam in img]
    nOverThresh = [np.sum(cam, axis=(0, 1)) for cam in bwImg]
    isCamTrig = [np.any(x >= nBubbleSize) for x in nOverThresh]

    # detect bubbles on trigger frame
    # bubbles = collections.OrderedDict(
    #     cam=[],
    #     trigFrame=[],
    #     frame=[],
    #     iPostTrig=[],
    #     ibubimage=[],
    #     nbubimage=[],
    #     ipix=[],
    #     jpix=[],
    #     xc=[],
    #     yc=[],
    #     dx=[],
    #     dy=[],
    #     theta=[],
    #     mass=[],
    # )
    Bub = Bubble()
    for iCam in range(0, len(img)):
        if isCamTrig[iCam]:
            # trigger frame index
            iTrigger = np.nonzero(nOverThresh[iCam] >= nBubbleSize)[0][0]
            for iFrame in range(iTrigger, np.amin([iTrigger + nBubbleFrame, bwImg[iCam].shape[2]])):
                # -- based on cv2.connectedComponentsWithStats
                n_label, labels, stats, centroids = cv2.connectedComponentsWithStats(bwImg[iCam][:, :, iFrame].copy())
                is_Bub = stats[:,4] >= nBubbleSize
                n_Bub = np.sum(is_Bub[1:]) # first element is background
                ix_Bub = np.nonzero(is_Bub)[0]

                for i_Bub in range(1, n_Bub + 1):
                    Bub.bubbles['cam'].append(iCam)
                    Bub.bubbles['trigFrame'].append(nPreTrigger + iTrigger)
                    Bub.bubbles['iPostTrig'].append(iFrame - iTrigger)
                    Bub.bubbles['frame'].append(nPreTrigger + iFrame)
                    Bub.bubbles['ibubimage'].append(i_Bub)
                    Bub.bubbles['nbubimage'].append(n_Bub)
                    Bub.bubbles['ipix'].append(centroids[ix_Bub[i_Bub]][0])
                    Bub.bubbles['jpix'].append(centroids[ix_Bub[i_Bub]][1])
                    Bub.bubbles['xc'].append(centroids[ix_Bub[i_Bub]][0])
                    Bub.bubbles['yc'].append(centroids[ix_Bub[i_Bub]][1])
                    Bub.bubbles['dx'].append(stats[ix_Bub[i_Bub], 2])
                    Bub.bubbles['dy'].append(stats[ix_Bub[i_Bub], 3])
                    Bub.bubbles['theta'].append(0)
                    Bub.bubbles['mass'].append(0)
                    Bub.bubbles['ADCThresh'].append(ADCThresh)
                # --
                # -- based on cv2.findContours
                # modImg, contours, hierarchy = cv2.findContours(bwImg[iCam][:, :, iFrame].copy(),
                #                                                cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
                # masses = np.array([cv2.contourArea(x) for x in contours])
                # goodBlob = np.nonzero(masses > nBubbleSize)[0]
                #
                # # plt.imshow(bwImg[iCam][:, :, iFrame])
                # # plt.show()
                #
                # for i_Bub in range(0, len(goodBlob)):
                #     center, (dx, dy), theta = cv2.minAreaRect(contours[goodBlob[i_Bub]])
                #     # fitEllipse need at least 5 points
                #     center, (dx, dy), theta = cv2.minAreaRect(contours[i_Bub])
                #     Bub.bubbles['cam'].append(iCam)
                #     Bub.bubbles['trigFrame'].append(nPreTrigger + iTrigger)
                #     Bub.bubbles['iPostTrig'].append(iFrame - iTrigger)
                #     Bub.bubbles['frame'].append(nPreTrigger + iFrame)
                #     Bub.bubbles['ibubimage'].append(i_Bub + 1)
                #     Bub.bubbles['nbubimage'].append(len(goodBlob))
                #     Bub.bubbles['ipix'].append(center[0])
                #     Bub.bubbles['jpix'].append(center[1])
                #     Bub.bubbles['xc'].append(center[0])
                #     Bub.bubbles['yc'].append(center[1])
                #     Bub.bubbles['dx'].append(dx)
                #     Bub.bubbles['dy'].append(dy)
                #     Bub.bubbles['theta'].append(theta)
                #     Bub.bubbles['mass'].append(masses[goodBlob[i_Bub]])
                # --
        if (not isCamTrig[iCam]) or (len(Bub.bubbles['cam']) < 1):
            Bub.bubbles['cam'].append(iCam)
            Bub.bubbles['trigFrame'].append(0)
            Bub.bubbles['iPostTrig'].append(0)
            Bub.bubbles['frame'].append(0)
            Bub.bubbles['ibubimage'].append(0)
            Bub.bubbles['nbubimage'].append(0)
            Bub.bubbles['ipix'].append(0)
            Bub.bubbles['jpix'].append(0)
            Bub.bubbles['xc'].append(0)
            Bub.bubbles['yc'].append(0)
            Bub.bubbles['dx'].append(0)
            Bub.bubbles['dy'].append(0)
            Bub.bubbles['theta'].append(0)
            Bub.bubbles['mass'].append(0)
            Bub.bubbles['ADCThresh'].append(0)

    # bubbles['cam']=np.uint8(bubbles['cam'])
    # bubbles['trigFrame']=np.uint8(bubbles['trigFrame'])
    # bubbles['iPostTrig']=np.uint8(bubbles['iPostTrig'])
    # bubbles['frame']=np.uint8(bubbles['frame'])
    # bubbles['ibubimage']=np.uint8(bubbles['ibubimage'])
    # bubbles['nbubimage']=np.uint8(bubbles['nbubimage'])
    # bubbles['ipix']=np.float32(bubbles['ipix'])
    # bubbles['jpix']=np.float32(bubbles['jpix'])
    # bubbles['xc']=np.float32(bubbles['xc'])
    # bubbles['yc']=np.float32(bubbles['yc'])
    # bubbles['dx']=np.float32(bubbles['dx'])
    # bubbles['dy']=np.float32(bubbles['dy'])
    # bubbles['theta']=np.float32(bubbles['theta'])
    # bubbles['mass'] = np.float32(bubbles['mass'])

    return Bub


def BubbleFinder(runDir, ev, nPreTrigger, nFrames, ADCThresh, nBubbleSize):
    imageList = BuildImageList(runDir, ev)
    if len(imageList) > 0:
        meanFrames, stdFrames = ComputeRefFrame(imageList, nPreTrigger)
        Bub = FindBubbles(imageList, nPreTrigger, nFrames, meanFrames, stdFrames, ADCThresh, nBubbleSize)
    else:
        Bub = Bubble()
        Bub.emptyBubble()

    runname = os.path.basename(runDir)
    runid_str = runname.split('_')
    runid = np.int32(runid_str)

    Bub.addRunEV(runid, ev)
    Bub.convertDType()

    return Bub


def test():
    t0 = time.time()
    imageList = BuildImageList('/mnt/XENON_DAQ/SBC-17-data/20170628_0', 5)
    meanFrames, stdFrames = ComputeRefFrame(imageList, 10)
    a = FindBubbles(imageList, 10, 3, meanFrames, stdFrames, 15, 9)
    print(a)
    print(time.time() - t0)

def test2():
    a= BubbleFinder('/mnt/XENON_DAQ/SBC-17-data/20170628_0', 5,
                 10, 3, 15, 4)
    print(a)

# test2()
