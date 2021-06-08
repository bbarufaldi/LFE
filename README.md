# LFE
Laplacian fractional entropy (LFE) measures evaluate higher-order statistical properties of x-ray images objectively [1-3]. The LFE measure is calculated from the relative entropy of the response histogram of an image compared to that of a Gaussian histogram matched for mean and variance.[1] These measures have shown to be an effective method to evaluate phantom realism and to select simulation parameters [2-3].

```diff
+ If you use this software to develop your work, cite manuscript [3] and our most recent publications for future references

[1]
[2]
[3]

```

This open-source project provides a Jave implementation of the LFE methods.

## Requirement

- JRE 8+ 

Make sure that you have install the Jave runtime environment on your local machine.

## Run

You can run the software via command line using XML files. As example, we provide a template of a phantom image for testing and configurable XML file ("config.xml").

Use the the following command to execute the LFE software:

```
java -jar LFE_Measurements.jar config.xml
```

## Config XML

The XML file is divided in three sections: Gabor, LFE, and Images:

### Images
- Mask 

Image path of binary mask (Gray 8-bits, format .RAW).

- Original

X-ray image path of binary mask (Gray 16-bits Unsigned, format .RAW).

- Width and height

Image and mask dimensions (pixels).


### Gabor

- Width and height 

Window dimensions of your Gabor kernels used for 2D convolution.

- DX and DY 

Pixel pitch of the x-ray image.

- FREQ_CEN


- A_RATIO


- PSHIFT


- angles


- BW_OCT



### LFE

- NH


- EXC_FRAC


- Optmizer_Gauss_ID and Optmizer_Lap_ID


  
  - 1: SIMPLEX
  - 2: BOBYQA 
  - 3: POWELL
  
- Save_Plot

Flag that writes the images of the LFE plots (.PNG). Set this flag to "false" if your system does not support a graphical interface.

The datapoints are saved in a separate .CSV file.

