# event list to process

# import jobs2grid
#
# jobs2grid.SubmitJobs(evlist, run/event)
#
# python_dir = '/coupp/app/home/coupp/anaconda2/bin'
# datadir
# outdir
# execdir = '/coupp/app/condor-exec/coupp'
# runseries

def append_to(element, to=[]):
    to.append(element)
    return to

my_list = append_to(12)
print my_list

my_other_list = append_to(42)
print my_other_list

my_other_list = append_to(42, to=[])
print my_other_list


hist = cv2.calcHist([bwImg[iCam][:, :, iFrame]],[0],None,[256],[0,256])


gim = cv2.GaussianBlur(img[0][:,:,iFrame], (7, 7), 3, 3)

# from matplotlib import pyplot as plt
print 'l'

plt.subplot(1,2,1); plt.imshow(img[0][:,:,iFrame], cmap='gray'); plt.subplot(1,2,2); plt.imshow(np.uint8(img[0][:,:,iFrame] <= -14), cmap='gray'); plt.show()

plt.imshow(np.uint8(img[0][:,:,iFrame] <= -14), cmap='gray'); plt.show()

plt.imshow(stdFrames[0], cmap='gray'); plt.show()

plt.imshow(np.uint8(img[0][:,:,iFrame] <= -4*stdFrames[0]), cmap='gray'); plt.show()

# adaptive threshold
ii=iFrame-1; aa=np.uint8((img[0][:,:,ii]-np.amin(img[0][:,:,ii]))*255.0/(np.amax(img[0][:,:,ii])-np.amin(img[0][:,:,ii]))); plt.imshow(aa, cmap='gray'); plt.show()

th3 = cv2.adaptiveThreshold(aa,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,11,0); plt.imshow(th3, cmap='gray'); plt.show()

bb=np.abs(aa);np.amax(np.tril(np.tile(bb,(bb.size,1))), axis=1)

plt.imshow(img[0][:,:,iFrame], cmap='gray'); plt.show()




th3 = cv2.adaptiveThreshold(zz,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,11,0); plt.imshow(th3, cmap='gray'); plt.show()

xx=img[0][:,:,iFrame-1]; plt.imshow(xx); plt.show()
xx=img[0][:,:,iFrame-1]; yy=cv2.GaussianBlur(xx, (9, 9), 3, 3); zz=np.uint8((yy-np.amin(yy))*255.0/(np.amax(yy)-np.amin(yy))); plt.imshow(zz, cmap='gray'); plt.show()
yy=dd; zz=np.uint8((yy-np.amin(yy))*255.0/(np.amax(yy)-np.amin(yy))); plt.imshow(zz, cmap='gray'); plt.show()
ret2,th2 = cv2.threshold(zz,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU); plt.imshow(th2); plt.show()

aa=img[0][:,:,iFrame-1]
(mm,nn)=aa.shape; ii = np.transpose(np.tile(np.arange(1,mm+1),(nn,1))); jj = np.tile(np.arange(1,nn+1),(mm,1)); kk=np.float32(jj)/ii;
uu=aa.copy()
uu[kk>np.float_(nn)/mm]=np.mean(aa)
dd=aa.copy(); dd[kk<np.float_(nn)/mm]=np.mean(aa)

plt.imshow(xx); plt.show()

se1 = cv2.getStructuringElement(cv2.MORPH_RECT, (5,5))
se2 = cv2.getStructuringElement(cv2.MORPH_RECT, (2,2))
mask = cv2.morphologyEx(img_bw, cv2.MORPH_CLOSE, se1)
mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, se2)

plt.subplot(1,2,1); plt.imshow(binImg, cmap='gray'); plt.subplot(1,2,2); plt.imshow(mask, cmap='gray'); plt.show()
modImg, contours, hierarchy = cv2.findContours(binImg.copy(), cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

ii=iFrame;plt.subplot(1,2,1); plt.imshow(img[0][:,:,ii], cmap='gray'); plt.subplot(1,2,2); plt.imshow(binImg, cmap='gray'); plt.show()
ii=iFrame;plt.subplot(1,2,1); plt.imshow(img[0][:,:,ii], cmap='gray'); plt.subplot(1,2,2); plt.imshow(np.uint8(img[0][:,:,ii] <= -14), cmap='gray'); plt.show()

bb=img[0][:,:,ii].copy()
yy=bb; zz=np.uint8((yy-np.amin(yy))*255.0/(np.amax(yy)-np.amin(yy))); plt.imshow(zz, cmap='gray'); plt.show()



