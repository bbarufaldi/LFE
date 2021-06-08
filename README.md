# LFE
Laplacian fractional entropy (LFE) measures evaluate higher-order statistical properties of x-ray images objectively [1-3]. The LFE measure is calculated from the relative entropy of the response histogram of an image compared to that of a Gaussian histogram matched for mean and variance.[1] These measures have shown to be an effective method to evaluate phantom realism and to select simulation parameters [2-3].

This open-source project provides a Jave implementation of the LFE methods. The sogtware is compatible to any operational system*

*Requirement*:

JRE 8+ - Make sure that you have install the Jave runtime environment on your local machine

*Run*:

You can run the software via command line using XML files. As example, we provide a template of a configurable XML file ("config.xml").

Use the the following command to execute the LFE software:

java -jar LFE_Measurements.jar config.xml

*Config XML*:

The XML file is divided in three sections: Gabor, LFE, and Images.

Gabor parameters:

1) Width and Height: Window size of your Gabor kernels used for 2D convolution.
2) DX and DY: pixel pitch of the x-ray image.
3) FREQ_CEN: 
4) A_RATIO:
5) PSHIFT:
6) angles:
7) BW_OCT:

