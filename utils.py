from __future__ import print_function
import os, glob
import numpy as np
from numpy.lib.stride_tricks import as_strided as ast
import scipy.misc
from scipy.optimize import curve_fit
from math import exp
from datetime import datetime
import xml.etree.ElementTree as ET
from sklearn.metrics import mean_squared_error
from skimage.morphology import opening, disk,closing,ball,erosion, dilation
from skimage.filters import gaussian, threshold_otsu, sobel, rank
from sklearn.metrics import mean_squared_error
import nibabel as nib

def paramap(img, xmask, ymask, zmask, res, time, tf, compressfactor, windSize_x, windSize_y, windSize_z):
    print('*************************** Starting Parameteric Map *****************************')
    # print('Prep For Loop:');print(str(datetime.now()));
    start_time = datetime.now()
    #1a. Windowing and image info
    global windSize, voxelscale, compression, imgshape, timeconst, times, xlist, ylist, zlist, windows, typefit;
    windSize = (windSize_x, windSize_y, windSize_z);
    voxelscale = res[0]*res[1]*res[2];
    compression = compressfactor; 
    imgshape = img.shape;
    typefit = tf;
    #img = img - np.mean(img[:,0:4,:,:,:,:],axis=1);img[img < 1]=0;

    # Make expected calculation time

    #1b. Creat time point and position lists
    timeconst = time;#time/(img.shape[1]+1);
    times = [];#times = np.arange(1,img.shape[3]+1);
    times = [i*time for i in range(1, img.shape[3]+1)];

    xlist = np.arange(min(xmask), max(xmask)+windSize[0], windSize[0])
    ylist = np.arange(min(ymask), max(ymask)+windSize[1], windSize[1])
    zlist = np.arange(min(zmask), max(zmask)+windSize[2], windSize[2])
    final_map = np.empty([img.shape[0], img.shape[1], img.shape[2]], dtype=list)
    first_loop = True
    for x_base in range(len(xlist)):
        for y_base in range(len(ylist)):
            for z_base in range(len(zlist)):
                cur_mask = np.zeros([img.shape[0], img.shape[1], img.shape[2]])
                indices = []
                for x in range(windSize[0]):
                    cur_index = []
                    cur_index.append(xlist[x_base]+x)
                    for y in range(windSize[1]):
                        cur_index.append(ylist[y_base]+y)
                        for z in range(windSize[2]):
                            cur_index.append(zlist[z_base]+z)
                            indices.append(cur_index.copy())
                            cur_index.pop()
                        cur_index.pop()
                    cur_index.pop()
                sig_indices = False
                for i in indices:
                    if max(img[i[0],i[1],i[2]]) != 0:
                        cur_mask[i[0],i[1],i[2]] = 1
                        sig_indices = True
                if not sig_indices:
                    continue
                cur_TIC = generate_TIC(img, cur_mask, times, 24.9,  voxelscale)
                normalizer = np.max(cur_TIC[:,1]);
                cur_TIC[:,1] = cur_TIC[:,1]/normalizer;

                # Bunch of checks
                if np.isnan(np.sum(cur_TIC[:,1])):
                    print('STOPPED:NaNs in the VOI')
                    return;
                if np.isinf(np.sum(cur_TIC[:,1])):
                    print('STOPPED:InFs in the VOI')
                    return;

                # Do the fitting
                try:
                    params, popt, RMSE, wholecurve = data_fit(cur_TIC,'BolusLognormal',normalizer, time);
                except RuntimeError:
                    print('RunTimeError')
                    #params = np.array([np.max(TIC[:,1])*normalizer, np.max(TIC[:,0])*np.max(TIC[:,1])*normalizer, np.max(TIC[:,0]), np.max(TIC[:,0])*4, 0]);
                    params = np.array([np.max(cur_TIC[:,1])*normalizer, np.trapz(cur_TIC[:,1]*normalizer, x=cur_TIC[:,0]), cur_TIC[np.argmax(cur_TIC[:,1]),0], np.max(cur_TIC[:,0])*2, 0]);

                for i in indices:
                    final_map[i[0], i[1],i[2]] = [popt[0], params[0], params[2], params[3]]

                # if first_loop:
                #     # first_loop_end_time = datetime.now()
                #     # print("Estimated time till completion:")
                #     # estimate = (first_loop_end_time.second-start_time.second)
                #     # estimate = estimate*len(xlist)*len(ylist)*len(zlist)
                #     # minutes = estimate/60
                #     # print(str(str(int(minutes))+" minutes, " + str(estimate-(int(minutes)*60))+" seconds"))
                #     first_loop = False

    print('Paraloop ended:')#;print(str(datetime.now()));
    return final_map;

def generate_TIC(window, mask, times, compression, voxelscale):
    TICtime=times;TIC=[]; 
    bool_mask = np.array(mask, dtype=bool)
    for t in range(0,window.shape[3]):
        tmpwin = window[:,:,:,t];      
        TIC.append(np.exp(tmpwin[bool_mask]/compression).mean()*voxelscale);
        # TIC.append(np.around((tmpwin[bool_mask]/compression).mean()*voxelscale, decimals=1)); 
    TICz = np.array([TICtime,TIC]).astype('float64'); TICz = TICz.transpose();
    TICz[:,1]=TICz[:,1]-np.mean(TICz[0:2,1]);#Substract noise in TIC before contrast.
    if TICz[np.nan_to_num(TICz)<0].any():#make the smallest number in the TIC 0.
        TICz[:,1]=TICz[:,1]+np.abs(np.min(TICz[:,1]));
    else:
        TICz[:,1]=TICz[:,1]-np.min(TICz[:,1]);
    return TICz;

def data_fit(TIC,model,normalizer, timeconst):
    #Fitting function
    #Returns the parameters scaled by normalizer
    #Beware - all fitting - minimization is done with data normalized 0 to 1. 
    if model == 'BolusLognormal':
        #kwargs = {"max_nfev":5000}
        popt, pcov = curve_fit(bolus_lognormal, TIC[0], TIC[1], p0=(1.0,3.0,0.5,0.1),bounds=([0., 0., 0., -1.], [np.inf, np.inf, np.inf, 10.]),method='trf')#p0=(1.0,3.0,0.5,0.1) ,**kwargs
        popt = np.around(popt, decimals=1);
        auc = popt[0]; rauc=normalizer*popt[0]; mu=popt[1]; sigma=popt[2]; t0=popt[3]; mtt=timeconst*np.exp(mu+sigma*sigma/2);
        tp = timeconst*exp(mu-sigma*sigma); wholecurve = bolus_lognormal(TIC[0], popt[0], popt[1], popt[2], popt[3]); pe = normalizer*np.max(wholecurve);
        rt0 = timeconst*t0;# + tp;
        
        # Filters to block any absurb numbers based on really bad fits. 
        if tp > 220: tp = 220; #pe = 0.1; rauc = 0.1; rt0 = 0.1; mtt = 0.1;
        if rt0 > 160: rt0 = 160; #pe = 0.1; rauc = 0.1; tp = 0.1; mtt = 0.1;
        if mtt > 2000: mtt = 2000; #pe = 0.1; rauc = 0.1; tp = 0.1; rt0 = 0.1;
        if pe > 1e+07: pe = 1e+07;
        if rauc > 1e+08: rauc = 1e+08;
        params = np.array([pe, rauc, tp, mtt, rt0]);
        # Get error parameters=
        residuals = TIC[1] - bolus_lognormal(TIC[0], popt[0], mu, sigma, t0);
        ss_res = np.sum(residuals[~np.isnan(residuals)]**2);# Residual sum of squares
        ss_tot = np.sum((TIC[1]-np.mean(TIC[1]))**2);# Total sum of squares
        r_squared = 1 - (ss_res / ss_tot);# R squared
        RMSE = (scipy.sum(residuals[~np.isnan(residuals)]**2)/(residuals[~np.isnan(residuals)].size-2))**0.5;#print('RMSE 1');print(RMSE);# RMSE
        rMSE = mean_squared_error(TIC[1], bolus_lognormal(TIC[0], popt[0], mu, sigma, t0))**0.5; wholecurve *= normalizer;#print('RMSE 2');print(rMSE);
        return params, popt, RMSE, wholecurve;

def bolus_lognormal(x, auc, mu, sigma, t0):      
    curve_fit=(auc/(2.5066*sigma*(x-t0)))*np.exp(-1*(((np.log(x-t0)-mu)**2)/(2*sigma*sigma))) 
    return np.nan_to_num(curve_fit)

def read_xmlraw_image_func(filename):    
    # get .raw filename
    filename_raw=filename[0:len(filename)-3]+('0.raw');
    fff = open(filename_raw,'rb')

    # parsing xml file
    tree = ET.parse(filename);
    root = tree.getroot();

    # ADD MAX FRAMES TO LOAD
    numfiles = len(root); # Comment this out if using max num frames on next lines

    for i in range(0, numfiles):
        if  root[i].tag=='Columns':
          M=int(root[i].text);
        if  root[i].tag=='Rows':
          N=int(root[i].text);  
        if  (root[i].find('Geometry') == None) == False:
          geometry = root[i].find('Geometry')
          if not geometry:
              continue
          layer = geometry.find('Layers')
          if not layer:
              continue
          layer = layer.find('Layer')
          if not layer:
              continue
          regionLocationMaxz1 = layer.find('RegionLocationMaxz1')
          physicalDeltaX = layer.find('PhysicalDeltaX')
          physicalDeltaY = layer.find('PhysicalDeltaY')
          physicalDeltaZ = layer.find('PhysicalDeltaZ')
          if regionLocationMaxz1:
              P=int(regionLocationMaxz1.text)+1
          if physicalDeltaX:
              voxelX = float(physicalDeltaX.text)
          if physicalDeltaY:
              voxelY = float(physicalDeltaY.text)
          if physicalDeltaZ:
              voxelZ = float(physicalDeltaZ.text)
          try:
              voxel=[float(physicalDeltaX.text)*10, float(physicalDeltaY.text)*10, float(physicalDeltaZ.text)*10] # cm to mm
              P=int(regionLocationMaxz1.text)+1
          except:
              continue
          
        if root[i].tag=='AcquisitionDateTime':
          tval=root[i].text;
          dateStr=tval[0:4]+'-'+tval[4:6]+'-'+tval[6:8]+' '+tval[8:10]+':'+tval[10:12]+':'+tval[12:len(tval)];
          time=float(tval[12:len(tval)])+float(tval[10:12])*60+float(tval[8:10])*3600;
    
    x = np.fromfile(fff,dtype=np.uint8)

    try:
        P
    except NameError:
        P = int(x.size / M / N)
        voxel=[(P-1)*10, (P-1)*10, (P-1)*10]

    shapes = (M,N,P);
    img = np.reshape(x, (P,N,M))

    return img, voxel, time, shapes, dateStr

def read3D(data, newres, cut):#=[[-1,-1,5,5],[-1,-1,25,5],[-1,-1,5,5]]):
    # cut=[[10,15,5,5],[50,5,20,4],[15,15,5,5]] #Size reduce with user selected caps
    # cut=[[-1,-1,5,5],[-1,-1,5,5],[-1,-1,5,5]] #Automatic size reduce
    # cut=[[0,0,5,5],[0,0,20,4],[0,0,5,5]] # Keep original size      
 
    N_lines_z_axis_cut=cut[0][0:2] #10,10   
    N_lines_z_axis_cut_limit=cut[0][2:4]#5,5
    N_lines_y_axis_cut=cut[1][0:2]#50,5  
    N_lines_y_axis_cut_limit=cut[1][2:4]#20,5
    N_lines_x_axis_cut=cut[2][0:2]##[15, 15]    
    N_lines_x_axis_cut_limit=cut[2][2:4]#5,5
    xmldir = data+('/*.xml');
    xmlnamedir = sorted(glob.glob(xmldir));

    # Ignore .mevis.xml files
    for file in xmlnamedir:
        if (len(file) >= 10 and file[-10:] == ".mevis.xml"):
            xmlnamedir.remove(file)

    if not len(xmlnamedir):
        print("No XML files found in folder")
        return

        
    imarray = []#np.zeros((len(xmlnamedir),shapes[2],shapes[1],shapes[0]),dtype='uint8')
    timeinitial = -1000; ix=-1;
    imi_mid, res, timelast, shapes, dateStr = read_xmlraw_image_func(xmlnamedir[np.uint16(len(xmlnamedir)/2)]); 
    imi10, res, timelast, shapes, dateStr = read_xmlraw_image_func(xmlnamedir[10]); 
    imi_mid, res, timelast, shapes, dateStr = read_xmlraw_image_func(xmlnamedir[1]);

    finalRes = res
    
    finalImi = []
    prevX = -1
    for xmlname in xmlnamedir:
        #imarray[0,xmlnamedir.index(xmlname),0,:,:,:], res, timelast, shapes, dateStr = read_xmlraw_image_func(xmlname)
        imi, res, timelast, shapes, dateStr = read_xmlraw_image_func(xmlname); 
        sz=imi.shape

        if timeinitial==-1000 and np.sum(cut) != 0:             
            timeinitial = timelast
            mp=np.mean(imi10,0)+np.mean(imi_mid,0);
            mpb=mp>0.3*threshold_otsu(mp); 
            lc=rank.otsu(mp.astype('uint16'), disk(10)); mpb2=mp>0.7*np.mean(lc); 
            mpb=(mpb+mpb2)>0;
            mp1=np.sum(mpb[:,40:mp.shape[1]-20],1);
            ix = np.where(mp1>0)
            mp0=np.sum(mpb,0)
            ix0 = np.where(mp0>0)          
            mpz=np.mean(imi10,1)+np.mean(imi_mid,1)+np.mean(imi,1);
            mpzb=mpz>1.2*threshold_otsu(mpz)           
            mpz1=np.sum(mpzb,1)
            iz = np.where(mpz1>0)     
            #ix[0][-1]=sz[1];ix[0][0]=0; #for test
            #ix0[0][-1]=sz[2];ix0[0][0]=0; #for test          
            #iz[0][-1]=sz[0];iz[0][0]=0; #for test   

            if N_lines_y_axis_cut[0]!=0: 
                N_lines_y_axis_cut[0]=max(0,max(N_lines_y_axis_cut[0]*(ix[0][0]>0),ix[0][0]-N_lines_y_axis_cut_limit[0])); 
            if N_lines_y_axis_cut[1]!=0:          
                N_lines_y_axis_cut[1]=sz[1]-min(sz[1],max((sz[1]-N_lines_y_axis_cut[1])*(ix[0][-1]<sz[1]),ix[0][-1]+N_lines_y_axis_cut_limit[1]))

            if N_lines_x_axis_cut[0]!=0:          
                N_lines_x_axis_cut[0]=max(0,max(N_lines_x_axis_cut[0]*(ix0[0][0]>0),ix0[0][0]-N_lines_x_axis_cut_limit[0])) 
            if N_lines_x_axis_cut[1]!=0:          
                N_lines_x_axis_cut[1]=sz[2]-min(sz[2],max((sz[2]-N_lines_x_axis_cut[1])*(ix0[0][-1]<sz[2]),ix0[0][-1]+N_lines_x_axis_cut_limit[1])) 

            if N_lines_z_axis_cut[0]!=0: 
                N_lines_z_axis_cut[0]=max(0,max(N_lines_z_axis_cut[0]*(iz[0][0]>0),iz[0][0]-N_lines_z_axis_cut_limit[0])) 
            if N_lines_z_axis_cut[1]!=0:          
                N_lines_z_axis_cut[1]=sz[0]-min(sz[0],max((sz[0]-N_lines_z_axis_cut[1])*(iz[0][-1]<sz[1]),iz[0][-1]+N_lines_z_axis_cut_limit[1]))

            #print(N_lines_z_axis_cut,N_lines_y_axis_cut,N_lines_x_axis_cut)

            print('Orginal image volume size',imi.shape)


        if imi.shape[0] == 1:
            continue

        if np.sum(cut) == 0:
            print('Orginal image volume size - no change',imi.shape);
        
        else:
            # reduce matrix size      
            imi=imi[N_lines_z_axis_cut[0]:sz[0]-N_lines_z_axis_cut[1],N_lines_y_axis_cut[0]:sz[1]- N_lines_y_axis_cut[1],N_lines_x_axis_cut[0]:sz[2]-N_lines_x_axis_cut[1]]

        if prevX != -1 and prevX != imi.shape[0]:
            print(imi.shape)
            continue
        prevX = imi.shape[0]
        imarray.append(imi)
        finalImi = imi
        
    if not len(finalImi) or not finalImi.shape[0] or not finalImi.shape[1] or not finalImi.shape[2]:
        print("Inputted image uses 2d data. Cannot parse into 3d data")
        exit(1)

    time = timelast - timeinitial; #total time of cine in seconds.
    if newres!=0:
        print('voxel size is changed from ', res, 'to voxel size of ', newres)
    
    print('Reduced image volume size',np.array(finalImi).shape)

    sh1_og=np.shape(imi_mid);sh1_og=(sh1_og[0]*sh1_og[1]*sh1_og[2])/(1024**2)
    sh1=np.shape(finalImi);sh1=(sh1[0]*sh1[1]*sh1[2])/(1024**2)
    imarray1 = np.zeros((1,len(imarray),1,finalImi.shape[0],finalImi.shape[1],finalImi.shape[2]),dtype='uint8')
    imarray = np.array(imarray)
    imarray1[0,:,0,:,:,:] = np.asarray(imarray)

    return imarray1, finalRes, time;

def xml2nifti(folderPath, fileDestination):
    print('Started xml2nii:')
    print(folderPath)
    print(str(datetime.now()))

    path = os.path.normpath(folderPath)
    splitpath = path.split(os.sep)
    name = splitpath[-1]

    # Read data
    # For clarification about magic array, see read3D function
    imarray_org, orgres, time = read3D(folderPath, 0, [[-1,-1,5,5],[-1,-1,25,5],[-1,-1,5,5]])
    timeconst = time/(imarray_org.shape[1]+1)
    print('Done 3D to 4D:')
    print(str(datetime.now()))

    imarray_tosub = np.min(imarray_org[:,:,:,:,:,:],axis=1)
    imarray_org = (imarray_org - imarray_tosub).astype('int16')
    del imarray_tosub   
    imarray_org[imarray_org < 0] = 0
    imarray_org = imarray_org.astype('uint8')
    imarray_org2 = np.squeeze(imarray_org)
    imarray_org2 = imarray_org2.swapaxes(0,3)
    imarray_org2 = imarray_org2.swapaxes(1,2)
    print('Saving cleaned up 4D:')
    print(str(datetime.now()))

    affine = np.eye(4)
    niiarray = nib.Nifti1Image(imarray_org2.astype('uint8'), affine)
    niiarray.header['pixdim'] = [4., orgres[0], orgres[1], orgres[2], timeconst, 0., 0., 0.]
    outputPath = os.path.join(fileDestination, str(name+'.nii.gz'))
    nib.save(niiarray, outputPath)
    return outputPath