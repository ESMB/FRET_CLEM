#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:57:38 2021

@author: Mathew
"""


from skimage.io import imread
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from skimage import filters
import czifile

# Path to the images:

pathlist=[]

pathlist.append("/Users/Mathew/Documents/Current analysis/FRET Heatmap/")


# Filename contains

filename=".czi"


# Folders to analyse:
    
def load_image(toload):
    
    img = czifile.imread(toload)
    
    return img

def z_project(image):
    
    mean_int=np.mean(image,axis=0)
  
    return mean_int

# Subtract background:
def subtract_bg(image):
    background = threshold_local(image, 11, offset=np.percentile(image, 1), method='median')
    bg_corrected =image - background
    return bg_corrected

def threshold_image_otsu(input_image):
    threshold_value=filters.threshold_otsu(input_image)  
    
    # threshold_value=input_image.mean()+3*input_image.std()
    print(threshold_value)
    binary_image=input_image>threshold_value

    return threshold_value,binary_image

def threshold_image_standard(input_image,thresh):
     
    binary_image=input_image>thresh

    return binary_image

# Threshold image using otsu method and output the filtered image along with the threshold value applied:
    
def threshold_image_fixed(input_image,threshold_number):
    threshold_value=threshold_number   
    binary_image=input_image>threshold_value

    return threshold_value,binary_image

# Label and count the features in the thresholded image:
def label_image(input_image):
    labelled_image=measure.label(input_image)
    number_of_features=labelled_image.max()
 
    return number_of_features,labelled_image
    
# Function to show the particular image:
def show(input_image,color=''):
    if(color=='Red'):
        plt.imshow(input_image,cmap="Reds")
        plt.show()
    elif(color=='Blue'):
        plt.imshow(input_image,cmap="Blues")
        plt.show()
    elif(color=='Green'):
        plt.imshow(input_image,cmap="Greens")
        plt.show()
    else:
        plt.imshow(input_image)
        plt.show() 
    
        
# Take a labelled image and the original image and measure intensities, sizes etc.
def analyse_labelled_image(labelled_image,original_image):
    measure_image=measure.regionprops_table(labelled_image,intensity_image=original_image,properties=('area','perimeter','centroid','orientation','major_axis_length','minor_axis_length','mean_intensity','max_intensity'))
    measure_dataframe=pd.DataFrame.from_dict(measure_image)
    return measure_dataframe

# This is to look at coincidence purely in terms of pixels

def coincidence_analysis_pixels(binary_image1,binary_image2):
    pixel_overlap_image=binary_image1&binary_image2         
    pixel_overlap_count=pixel_overlap_image.sum()
    pixel_fraction=pixel_overlap_image.sum()/binary_image1.sum()
    
    return pixel_overlap_image,pixel_overlap_count,pixel_fraction

# Look at coincidence in terms of features. Needs binary image input 

def feature_coincidence(binary_image1,binary_image2):
    number_of_features,labelled_image1=label_image(binary_image1)          # Labelled image is required for this analysis
    coincident_image=binary_image1 & binary_image2        # Find pixel overlap between the two images
    coincident_labels=labelled_image1*coincident_image   # This gives a coincident image with the pixels being equal to label
    coinc_list, coinc_pixels = np.unique(coincident_labels, return_counts=True)     # This counts number of unique occureences in the image
    # Now for some statistics
    total_labels=labelled_image1.max()
    total_labels_coinc=len(coinc_list)
    fraction_coinc=total_labels_coinc/total_labels
    
    # Now look at the fraction of overlap in each feature
    # First of all, count the number of unique occurances in original image
    label_list, label_pixels = np.unique(labelled_image1, return_counts=True)
    fract_pixels_overlap=[]
    for i in range(len(coinc_list)):
        overlap_pixels=coinc_pixels[i]
        label=coinc_list[i]
        total_pixels=label_pixels[label]
        fract=1.0*overlap_pixels/total_pixels
        fract_pixels_overlap.append(fract)
    
    
    # Generate the images
    coinc_list[0]=1000000   # First value is zero- don't want to count these. 
    coincident_features_image=np.isin(labelled_image1,coinc_list)   # Generates binary image only from labels in coinc list
    coinc_list[0]=0
    non_coincident_features_image=~np.isin(labelled_image1,coinc_list)  # Generates image only from numbers not in coinc list.
    
    return coinc_list,coinc_pixels,fraction_coinc,coincident_features_image,non_coincident_features_image,fract_pixels_overlap


for path in pathlist:       
    # :pad the images:
    for root, dirs, files in os.walk(path):
      for name in files:
              # print(name)
              if filename in name:
                  if 'pdf' not in name:
                      resultsname = name
                      print(resultsname)
    
    image=load_image(path+resultsname)
    
    lipid=image[:,:,0,:,:,:, :,:]
    donor=image[:,:,2,:,:,:, :,:]
    FRET=image[:,:,3,:,:,:, :,:]
    acceptor=image[:,:,4,:,:,:, :,:]
    
  
    
    donor_ave=np.mean(donor,axis=3)
    acceptor_ave=np.mean(acceptor,axis=3)
    FRET_ave=np.mean(FRET,axis=3)
    lipid_ave=np.mean(lipid,axis=3)
    
    donor_flat=donor_ave[0, 0, 0, :, :, 0]
    acceptor_flat=acceptor_ave[0, 0, 0, :, :, 0]
    FRET_flat=FRET_ave[0, 0, 0, :, :, 0]
    lipid_flat=lipid_ave[0, 0, 0, :, :, 0]
    
    
    value,donor_binary=threshold_image_otsu(donor_flat)
    
    value,acceptor_binary=threshold_image_otsu(acceptor_flat)
    
    
    both_binary=donor_binary*acceptor_binary
    
    
    fret_im=acceptor_flat/(acceptor_flat+donor_flat)
    
    
    fret_im_thresh=both_binary*fret_im
    
    imsr2 = Image.fromarray(lipid_flat)
    imsr2.save(path+'Lipid.tif')
    
    imsr2 = Image.fromarray(fret_im_thresh)
    imsr2.save(path+'FRET.tif')
    
    
   


