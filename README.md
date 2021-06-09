# Laplacian fractional entropy
Laplacian fractional entropy (LFE) measures evaluate higher-order statistical properties of x-ray images objectively [1-3]. The LFE measure is calculated from the relative entropy of the response histogram of an image compared to that of a Gaussian histogram matched for mean and variance.[1] These measures have shown to be an effective method to evaluate phantom realism and to select simulation parameters [2-3].

```diff
+ If you use this software to develop your work, cite our most recent publications for future references

[1] Abbey, Craig K et al. “Non-Gaussian statistical properties of breast images.” Medical physics vol. 39,11 (2012): 7121-30. doi:10.1118/1.4761869
[2] Abbey, Craig K et al. “Evaluation of non-Gaussian statistical properties in virtual breast phantoms.” J Med Imaging vol. 6,2 (2019):2329-4302. 10.1117/1.JMI.6.2.025502
[3] Barufaldi, B et al. “Computational Breast Anatomy Simulation Using Multi-scale Perlin Noise.” IEEE Trans Med Imaging vol. (2021):

```

This open-source project provides a Java implementation of the LFE methods.

## Requirement

- JRE 8+ 

Make sure that you have install the Java runtime environment on your local machine.

## Instructions

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

Window dimensions of your Gabor kernels used for 2D convolution (spatial domain).

- DX and DY 

Pixel pitch of the x-ray image.

- FREQ_CEN

Frequency of of Gabor kernel.

- A_RATIO

Aspect ratio of Gabor kernel.

- PSHIFT

Phase shift of Gabor kernel.

- angles

Number of rotation angles (PI/n) of Gabor kernel. Ex. n=4 (0°, 45°, 90°, and 180°).

- BW_OCT

Octave band width of Gabor kernel.


### LFE

- NH

The number of histogram bins.

- EXC_FRAC

Histogram exclusion fraction. This parameter controls how many points fall in the last histogram bin which is outside the range of consideration. This is done to avoid problems with small counts in extreme values, which LFE will regard as non-gaussian behavior.

- Optmizer_Gauss_ID and Optmizer_Lap_ID

IDs of minimization algorithms used to minimize non-linear funtions to the optimal solution. These "solvers" do not require compute derivatives. For more information visit [Commons Math](https://commons.apache.org/proper/commons-math/) website.
  
  1: SIMPLEX: Simplex-based direct search optimization.
  
  2: BOBYQA: Powell's BOBYQA algorithm.
  
  3: POWELL: Powell's algorithm.
  
- Save_Plot

Flag that writes the images of the LFE plots (.PNG). The plots include the histogram and fitted probability functions. Set this flag to "false" if your system does not handle graphics.

The datapoints are saved in a separate .CSV file.

