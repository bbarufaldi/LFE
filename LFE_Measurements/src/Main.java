
/*

* Laplacian fractional entropy (LFE) is a measure of higher-order statistical properties of images [1-3]. 
* The LFE measure is calculated from the relative entropy of the response histogram of an image compared to that of a Gaussian histogram matched for mean and variance [1]. 
* The approach accounts for entropy effects from binning into a histogram. The Laplacian distribution is used as a scaling constant, since Laplacian response histograms are a model for natural scenes. 
* The entropy relative to a Gaussian is normalized to the relative entropy of a Laplacian distribution relative to the Gaussian (and also accounting for binning effects). 
* Thus 100% LFE can be interpreted as having non-Gausian "structure" that is equivalent to a Laplacian distribution. 
* These measures have shown to be an effective method to evaluate the impact of anatomy and noise in mammography [1], and virtual breast phantom realism [2-3].

* If you use this software to develop your work, cite our most recent publications for future references

[1] Abbey, Craig K et al. “Non-Gaussian statistical properties of breast images.” Med Phys vol. 39,11 (2012): 7121-30. doi:10.1118/1.4761869
[2] Abbey, Craig K et al. “Evaluation of non-Gaussian statistical properties in virtual breast phantoms.” J Med Imaging vol. 6,2 (2019):2329-4302. doi:10.1117/1.JMI.6.2.025502
[3] Barufaldi, B et al. “Computational Breast Anatomy Simulation Using Multi-scale Perlin Noise.” IEEE Trans Med Imaging vol. 0 (2021):1-10 doi:10.1109/TMI.2021.3087958

 */

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.io.FileOpener;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.plugin.FFT;
import ij.plugin.FFTMath;
import ij.plugin.ImageCalculator;
import ij.plugin.filter.Convolver;
import ij.plugin.filter.FFTCustomFilter;
import ij.plugin.filter.FFTFilter;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Container;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.AbstractSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.jtransforms.fft.FloatFFT_2D;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

/**
 *
 * @author Bruno Barufaldi
 */
public class Main {
    
    //Options
    static final private int SIMPLEX = 1;
    static final private int BOBYQA = 2;
    static final private int POWELL = 3;
    
    static final private int SPATIAL = 0;
    static final private int FOURIER = 1;
    
    // Default parameters
    static private int NH = 99; 
    static private double EXC_FRAC = 0.01;
    static private int opt_gauss = BOBYQA;       
    static private int opt_lap = POWELL;         
    static private double[] H_LO, H_HI;
    static private double[] P_H;
    static private boolean flag_plot = false;
    
    // Gabor parameters
    static private int opt_gabor = FOURIER;
    static private double A_RATIO = 1.0;
    static private double PSHIFT = Math.PI/2;
    static private double FREQ_CEN = 1;
    static private int nAngles = 6;
    static private double[] angle;
    static private double BW_OCT = 1.4;
        
    // Kernell parameters
    static private int NX = 31; 
    static private int NY = 31;
    static private double DX = 0.07; //pixel pitch
    static private double DY = 0.07;
    
    // Image Parameters
    static private String imageName, maskName;
    static private ImagePlus original, mask;
    
     /**
     * @param args the command line arguments
     */
    
    public static void main(String[] args) {
        
        //read config file
        File xml = new File(args[0]);
        //File xml = new File("config.xml");
        readConfig(xml);
        
        ImageStack is = calculate_Gabors();
        
        //LFE method [Abbey et al., Med. Phys vol. 39,11(2012):7121-30. doi:10.1118/1.4761869]
        for(int k=1 ; k <= is.getSize(); k++){
            //make_resp_hist.pro
            double[][] H = make_resp_hist(is.getProcessor(k)); 
            //lfe_fit_histo.pro
            lfe_fit_histo(H, k-1);
        }

    }
    
    public static void readConfig(File xml) {
        
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();  
        
        try {
            DocumentBuilder db = dbf.newDocumentBuilder();
            Document doc = db.parse(xml);
            doc.getDocumentElement().normalize();
            
            //Lists
            NodeList gaborList = doc.getElementsByTagName("Gabor");
            NodeList imageList = doc.getElementsByTagName("Images");
            NodeList lfeList = doc.getElementsByTagName("LFE");
            
            Node node;
            for (int i = 0; i < gaborList.getLength(); i++){  
                node = gaborList.item(i);  
                if (node.getNodeType() == Node.ELEMENT_NODE){  
                    Element eElement = (Element) node;
                    
                    opt_gabor = Integer.parseInt(eElement.getElementsByTagName("domain").item(0).getTextContent());
                    NX = Integer.parseInt(eElement.getElementsByTagName("width").item(0).getTextContent());
                    NY = Integer.parseInt(eElement.getElementsByTagName("height").item(0).getTextContent());
                    DX = Double.parseDouble(eElement.getElementsByTagName("DX").item(0).getTextContent());
                    DY = Double.parseDouble(eElement.getElementsByTagName("DY").item(0).getTextContent());
                    FREQ_CEN = Double.parseDouble(eElement.getElementsByTagName("FREQ_CEN").item(0).getTextContent());
                    A_RATIO = Double.parseDouble(eElement.getElementsByTagName("A_RATIO").item(0).getTextContent());
                    PSHIFT = Math.PI/Double.parseDouble((eElement.getElementsByTagName("PSHIFT").item(0).getTextContent().replace("PI/", "")));
                    nAngles = Integer.parseInt(eElement.getElementsByTagName("angles").item(0).getTextContent());
                    BW_OCT = Double.parseDouble(eElement.getElementsByTagName("BW_OCT").item(0).getTextContent());
                    angle = new double[nAngles];
                    
                }
            }
            
            //Window checking
            if(NX%2 == 0){
                NX = NX-1;
                if(NY%2 == 0)
                    NY = NY-1;
                System.out.println("WARNING: Window dimensions should be created using 'odd' values.\nThe window dimensions were updated. Check the folder 'kernel'.");
            }
            
            for (int i = 0; i < lfeList.getLength(); i++){  
                node = lfeList.item(i);  
                if (node.getNodeType() == Node.ELEMENT_NODE){  
                    Element eElement = (Element) node;
                    
                    NH = Integer.parseInt(eElement.getElementsByTagName("NH").item(0).getTextContent());
                    EXC_FRAC = Double.parseDouble(eElement.getElementsByTagName("EXC_FRAC").item(0).getTextContent());
                    opt_gauss = Integer.parseInt(eElement.getElementsByTagName("Optmizer_Gauss_ID").item(0).getTextContent());
                    opt_lap = Integer.parseInt(eElement.getElementsByTagName("Optmizer_Lap_ID").item(0).getTextContent());
                    flag_plot = Boolean.parseBoolean(eElement.getElementsByTagName("Save_Plot").item(0).getTextContent());
                }
            }
            
            FileInfo info = new FileInfo();
            FileOpener opener;
            for (int i = 0; i < imageList.getLength(); i++){  
                node = imageList.item(i);  
                 
                if (node.getNodeType() == Node.ELEMENT_NODE){  
                    Element eElement = (Element) node;
                    
                    info.fileName = eElement.getElementsByTagName("original").item(0).getTextContent();
                    info.fileType = FileInfo.GRAY16_UNSIGNED; 
                    info.width = Integer.parseInt(eElement.getElementsByTagName("width").item(0).getTextContent());
                    info.height = Integer.parseInt(eElement.getElementsByTagName("height").item(0).getTextContent());
                    info.intelByteOrder = false;
                    opener = new FileOpener(info);
                    original = opener.open(false);
                    
                    info.fileName = eElement.getElementsByTagName("mask").item(0).getTextContent();
                    info.fileType = FileInfo.GRAY8; 
                    opener = new FileOpener(info);
                    mask = opener.open(false);
                    
                    imageName = original.getTitle();
                    maskName = mask.getTitle();
                    //original.show();
                    //mask.show();
                }
            }
        } catch (ParserConfigurationException | SAXException | IOException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public static ImageStack calculate_Gabors(){
        
        int width = original.getWidth();
        int height = original.getHeight();
        ImageStack is = new ImageStack(width, height);
        
        //Parameters used for kernel creation
        double DEL_F = (Math.pow(2.0, BW_OCT) - 1.0)*FREQ_CEN/(Math.pow(2.0, BW_OCT) + 1.0);
        double SD_X_FREQ = DEL_F/Math.sqrt(2.0*Math.log(2.0));
        double SD_X = 1.0/(2.0*Math.PI*SD_X_FREQ);
        double SD_Y = A_RATIO*SD_X;
        
        // Rotate kernel from 0 to 180 degrees in nAngles
        double rotationAngle = Math.PI/(double)nAngles;
        double theta;
        int X_CEN, Y_CEN;
        FloatProcessor filter;
        ImageStack kernels;
        
        kernels = new ImageStack(NX, NY);
        X_CEN = (int) Math.round(NX / 2);
        Y_CEN = (int) Math.round(NY / 2);
        
        for (int i=0; i<nAngles; i++){   
            
            theta = rotationAngle * i;
            angle[i] = theta;
            filter = new FloatProcessor(NX, NY);

            for (int x=-X_CEN; x<=X_CEN; x++){
                for (int y=-Y_CEN; y<=Y_CEN; y++){

                    double XR = (double)x * Math.cos(theta) * DX + (double)y * Math.sin(theta) * DY;
                    double YR = (double)y * Math.cos(theta) * DY - (double)x * Math.sin(theta) * DX;

                    double ETERM = -0.5*(Math.pow(XR/SD_X,2) + Math.pow(YR/SD_Y,2));
                    double gb = Math.exp(ETERM)*Math.cos(2.0*Math.PI*FREQ_CEN*XR + PSHIFT);
                    filter.setf(x+X_CEN, y+Y_CEN, (float)gb );
                }
            }
            kernels.addSlice("kernel angle = " + theta, filter);
        }
        
        // Save kernels
        FileSaver fs;
        File k_folder = new File("kernel");
        if(!k_folder.exists())
            k_folder.mkdir();
        
        for(int i=1 ; i <= kernels.getSize(); i++){
            fs = new FileSaver(new ImagePlus("_"+(i-1),kernels.getProcessor(i)));
            fs.saveAsTiff(k_folder.getName()+"/_"+(i-1)+".tiff");
        }
        
        // Save Gabor filtering results
        File g_folder = new File("gabor");
        if(!g_folder.exists())
            g_folder.mkdir();
        
        // Fourrier Filtering
        if(opt_gabor == FOURIER){
            ImageProcessor fft1, fft2, result; //original, kernel, and result
            float[] pixels;
            
            //Window size fft
            int fft_size = 1;
            int max = Math.max(width, height); 
            while(fft_size < max) //get power 2 to resize images
                fft_size*=2;
            
            //copy original image
            fft1 = doForwardFFT(original.duplicate().getProcessor().convertToFloatProcessor(), fft_size);
            
            ImageProcessor k_temp;
            for (int i=0; i<kernels.size(); i++){
                
                //Flipping required for convolution (transpose)
                k_temp = kernels.getProcessor(i+1);
                k_temp.flipHorizontal();
                k_temp.flipVertical();
                
                fft2 = doForwardFFT(k_temp, fft_size);
                result = multiply(fft1, fft2, fft_size, width, height);
                
                //normalize signal to sqrt 2*PI
                //TODO: Update normalization method as a function of kernel parameters
                pixels = (float[])result.getPixels();
                for(int k=0; k<pixels.length; k++)
                    pixels[k] = pixels[k]/(float)(Math.sqrt(2*Math.PI));
                result.setPixels(pixels);
                
                fs = new FileSaver(new ImagePlus("_"+i, result));
                fs.saveAsTiff(g_folder.getName()+"/_"+i+".tiff");
                is.addSlice("angle="+i+ " gamma="+A_RATIO+ " psi="+PSHIFT+ " freq="+FREQ_CEN, result);
                
            }
        }
        
        // Spacial Domain Filtering
        else if(opt_gabor == SPATIAL){
            float[] kernel;
            ImageProcessor ip = original.getProcessor();
            FloatProcessor fp;
            Convolver c = new Convolver();
            
            int kw, kh;
            kw = kernels.getWidth()%2 == 0 ? kernels.getWidth()-1 : kernels.getWidth();
            kh = kernels.getHeight()%2 == 0 ? kernels.getHeight()-1 : kernels.getHeight();
            
            for (int i=0; i<kernels.size(); i++){
                kernel = (float[]) kernels.getProcessor(i+1).getPixels();
                fp = ip.convertToFloatProcessor();
                c.convolveFloat(fp, kernel, kw, kh);      

                is.addSlice("angle="+i+ " gamma="+A_RATIO+ " psi="+PSHIFT+ " freq="+FREQ_CEN, fp);
                fs = new FileSaver(new ImagePlus("_"+i,fp));
                fs.saveAsTiff(g_folder.getName()+"/_"+i+".tiff");
            }
        }

        //Thresold values for segmentation - (make_resp_hist)
        float val;
        for(int k=1 ; k <= is.getSize(); k++){
            for(int i=0; i<mask.getWidth(); i++){
                for(int j=0; j<mask.getHeight(); j++){
                    val = mask.getProcessor().getPixelValue(i, j);
                    if(val == 0)
                        is.getProcessor(k).putPixelValue(i, j, Float.MAX_VALUE); 
                }
            }
        }
        
        return is;
    }
    
    //make_resp_hist.pro
    public static double[][] make_resp_hist(ImageProcessor proc){
        
        float[] pixels = (float[])proc.getPixels();
        
        List<Float> list = new ArrayList<>();
        for(int i=0; i<pixels.length; i++){ //segment valid histogram values
            if(pixels[i] != Float.MAX_VALUE)
                list.add(pixels[i]);
        }

        float[] RESP = new float[list.size()];
        float[] RESPS_SRT = new float[list.size()];
        for(int i=0; i<RESP.length; i++){
            RESP[i] = list.get(i);
            RESPS_SRT[i] = list.get(i);
        }
        
        //Sort Array
        Arrays.sort(RESPS_SRT);
        
        int NRESP = RESP.length;
        long IH_LO = Math.round(0.5*EXC_FRAC*NRESP);
        long IH_HI = Math.round(Math.round((1.0-0.5*EXC_FRAC)*NRESP));
        
        float RESP_LO = RESPS_SRT[(int)IH_LO];
        float RESP_HI = RESPS_SRT[(int)IH_HI];
        
        double DH = (double)(RESP_HI-RESP_LO)/(double)NH;
        
        double[] IRESPS_SRT = new double[NRESP];
        for(int i=0; i<IRESPS_SRT.length; i++)
            IRESPS_SRT[i] = Math.floor((double)(RESPS_SRT[i] - RESP_LO)/DH);
        
        double[][] H = new double[NH+1][2];
        
        for(int i=0; i<NH; i++)
            H[i][0] = RESP_LO + DH*(0.5 + i);
        
        H[NH][0] = DH; //last element
        
        long NHCNT = 0;
        for(int i=0; i<NH; i++){
            int NH_i = 0;
            for(int j=0; j<IRESPS_SRT.length; j++){
                if(IRESPS_SRT[j] == i)
                    NH_i++;
            }
                
            H[i][1] = NH_i;
            NHCNT+= NH_i;
        }
        
        H[NH][1] = NRESP - NHCNT;
        
        return H;
        
    }
    
    //lfe_fit_histo.pro - MAIN
    public static void lfe_fit_histo(double[][] H, int k){
        
        double[] INFORM_GL = new double[2]; //gausian, laplacian
        double DHX = H[NH][0];
        
        H_LO = new double[NH];
        H_HI = new double[NH];
        
        for(int i=0; i<NH; i++){
            H_LO[i] = H[i][0] - 0.5*DHX;
            H_HI[i] = H[i][0] + 0.5*DHX;
        }
        
        double[] P_HISTO = new double[NH+1];
        
        double sum=0;
        for(int i=0; i<H.length; i++)
            sum+=H[i][1];
        
        for(int i=0; i<H.length; i++){
            P_HISTO[i] = H[i][1]/sum;
            //System.out.println(P_HISTO[i]);
        }
        
        P_H = P_HISTO;
        
        double H_MN = 0;
        for(int i=0; i<NH; i++)
            H_MN += H[i][0]*P_H[i];
        
        double H_SD = 0;
        for(int i=0; i<NH; i++)
            H_SD += Math.pow(H[i][0]-H_MN,2)*P_H[i];
        
        H_SD = Math.sqrt(H_SD);
        
        double[] X = {H_MN, H_SD};
        
        double[] P_GAUSS = new double[NH+1];
        
        int NPRM = 2;
        double[][] XBND = new double[NPRM][2];
        XBND[0][0] = H[0][0];
        XBND[0][1] = H[NH-1][0];
        XBND[1][0] = DHX;
        XBND[1][1] = 0.5*NH*DHX;
        
        MultivariateFunction lfe_func_gauss = ((double[] X1) -> {
            double EPS = 1.0E-09;
            
            //P_GAUSS[0:NH-1] = (GAUSS_PDF((H_HI-X[0])/X[1])-GAUSS_PDF((H_LO-X[0])/X[1])) > 0.0
            double sum1 = 0;
            for (int i = 0; i<NH; i++) {
                P_GAUSS[i] = cdf((H_HI[i] - X1[0]) / X1[1]) - cdf((H_LO[i] - X1[0]) / X1[1]);
                sum1 += P_GAUSS[i];
            }
            P_GAUSS[NH] = 1.0 - sum1;
            //RELENT = TOTAL(P_H*ALOG((P_H+EPS)/(P_GAUSS+EPS)))/ALOG(2.0)
            double RELENT=0;
            for(int i=0; i<P_H.length; i++)
                RELENT+= (P_H[i]*Math.log((P_H[i]+EPS)/(P_GAUSS[i]+EPS)))/Math.log(2.0);
            
            return RELENT;
            
        });
        
        PowellOptimizer opt;
        BOBYQAOptimizer boby;
        SimplexOptimizer optimizer;
        PointValuePair minGauss=null, minLap=null;
        
        switch(opt_gauss){
            case SIMPLEX:
                optimizer = new SimplexOptimizer(1e-10, 1e-30);
                AbstractSimplex simplex = new MultiDirectionalSimplex(2);
                minGauss = optimizer.optimize(new MaxEval(1000),new ObjectiveFunction(lfe_func_gauss),simplex,GoalType.MINIMIZE,new InitialGuess(X));
                break;
            
            case BOBYQA:
                boby = new BOBYQAOptimizer(6); //not sure what number of interpolation means
                minGauss = boby.optimize(new MaxEval(1000),new MaxIter(1000),new ObjectiveFunction(lfe_func_gauss),GoalType.MINIMIZE,new InitialGuess(X),new SimpleBounds(new double[]{XBND[0][0], XBND[1][0]}, new double[]{XBND[0][1], XBND[1][1]}));
                break;
        
            case POWELL:
                opt = new PowellOptimizer(1e-15, 1e-17, 1e-15, 1e-15);
                minGauss = opt.optimize(new MaxEval(1000),new ObjectiveFunction(lfe_func_gauss),GoalType.MINIMIZE,new InitialGuess(X));
                break;
                
            default:
                System.out.println("Warning: GAUSS_ID UNKNOWN");
                boby = new BOBYQAOptimizer(6); //not sure what number of interpolation means
                minGauss = boby.optimize(new MaxEval(1000),new MaxIter(1000),new ObjectiveFunction(lfe_func_gauss),GoalType.MINIMIZE,new InitialGuess(X),new SimpleBounds(new double[]{XBND[0][0], XBND[1][0]}, new double[]{XBND[0][1], XBND[1][1]}));
                break;
        }
        

        double RE_GAUSS = lfe_func_gauss(minGauss.getPoint(), P_GAUSS);
        
        P_H = P_GAUSS;
        
        double[] X_LAP = minGauss.getPoint();
        
        X_LAP[1] = 0.8*X_LAP[1]/Math.sqrt(2.0);
        
        double[] P_LAP = new double[NH+1];
        
        MultivariateFunction lfe_func_lap = ((double[] X1) -> {
            double EPS = 1.0E-09;
            double[] Q_LO = new double[NH];
            for(int i=0; i<Q_LO.length; i++) Q_LO[i] = 1.0; //init Q_LO

            ArrayList<Integer> INDS_LO = new ArrayList<>(); //get lo indices
            for(int i=0; i< H_LO.length; i++){
                if(H_LO[i] < X1[0])
                    INDS_LO.add(i);
            }

            for(int i=0; i<INDS_LO.size(); i++)
                Q_LO[INDS_LO.get(i)] = 0.0;     //zero lo indices

            double[] Q_HI = new double[NH];
            for(int i=0; i<Q_HI.length; i++) Q_HI[i] = 1.0; //init Q_HI

            ArrayList<Integer> INDS_HI = new ArrayList<>();
            for(int i=0; i< H_HI.length; i++){
                if(H_HI[i] < X1[0])
                    INDS_HI.add(i);
            }

            for(int i=0; i<INDS_HI.size(); i++)
                Q_HI[INDS_HI.get(i)] = 0.0;     //zero hi indices

            double[] CDF_LO = new double[NH];
            double[] CDF_HI = new double[NH];

            double sum1=0;
            for(int i=0; i< CDF_LO.length; i++){
                CDF_LO[i] = 0.5*Math.exp((H_LO[i] - X1[0])/X1[1])*(1.0-Q_LO[i]) + (1.0-0.5*Math.exp(-(H_LO[i] - X1[0])/X1[1]))*Q_LO[i];
                CDF_HI[i] = 0.5*Math.exp((H_HI[i] - X1[0])/X1[1])*(1.0-Q_HI[i]) + (1.0-0.5*Math.exp(-(H_HI[i] - X1[0])/X1[1]))*Q_HI[i];
                P_LAP[i] = CDF_HI[i] - CDF_LO[i];
                sum1 += P_LAP[i];
            }

            P_LAP[NH] = 1.0 - sum1;

            double RELENT=0;
            for(int i=0; i<P_LAP.length; i++)
                RELENT+= (P_LAP[i]*Math.log((P_LAP[i]+EPS)/(P_H[i]+EPS)))/Math.log(2.0);

            return RELENT;
            
        });
        
        switch(opt_lap){
            case SIMPLEX:
                optimizer = new SimplexOptimizer(1e-10, 1e-30);
                AbstractSimplex simplex = new MultiDirectionalSimplex(2);
                minLap = optimizer.optimize(new MaxEval(1000),new ObjectiveFunction(lfe_func_lap),simplex,GoalType.MINIMIZE,new InitialGuess(X_LAP));
                break;
            
            case BOBYQA:
                boby = new BOBYQAOptimizer(6); //not sure what number of interpolation means
                minLap = boby.optimize(new MaxEval(1000),new MaxIter(1000),new ObjectiveFunction(lfe_func_lap),GoalType.MINIMIZE,new InitialGuess(X_LAP),new SimpleBounds(new double[]{XBND[0][0], XBND[1][0]}, new double[]{XBND[0][1], XBND[1][1]}));
                break;
        
            case POWELL:
                opt = new PowellOptimizer(1e-15, 1e-17, 1e-15, 1e-15);
                minLap = opt.optimize(new MaxEval(1000),new ObjectiveFunction(lfe_func_lap),GoalType.MINIMIZE,new InitialGuess(X_LAP));
                break;
                
            default:
                System.out.println("Warning: LAP_ID UNKNOWN");
                boby = new BOBYQAOptimizer(6); //not sure what number of interpolation means
                minLap = boby.optimize(new MaxEval(1000),new MaxIter(1000),new ObjectiveFunction(lfe_func_lap),GoalType.MINIMIZE,new InitialGuess(X_LAP),new SimpleBounds(new double[]{XBND[0][0], XBND[1][0]}, new double[]{XBND[0][1], XBND[1][1]}));
                break;
        }
        
        
        double RE_LAP = lfe_func_lap(minLap.getPoint(), P_LAP);
        
        double[] PRMS = new double[]{RE_GAUSS, RE_LAP, 
                                    minGauss.getPoint()[0], minGauss.getPoint()[1],
                                    minLap.getPoint()[0], minLap.getPoint()[1]};
           
        //report
        String name = imageName.replace(".raw", "_result_"+String.format("%.2f",angle[k])+"_" + FREQ_CEN);
        try {
            
            PrintWriter writer = new PrintWriter(new File(name + ".csv"));
            
            writer.println("H[0],H[1],P_HISTO,P_GAUSS,P_LAP");
            for(int i=0; i<P_HISTO.length; i++)
                writer.println(H[i][0] + "," + H[i][1] + "," + P_HISTO[i] + "," + P_GAUSS[i] + "," + P_LAP[i]);
            
            //Save summary on last line
            writer.println();
//            //DEBUG
//            System.out.println(imageName.replace(".raw", "") + "\t" + maskName.replace(".raw", "") + "\t" + String.format("%.2f",angle[k]) + "\t" + String.format("%.2f",(100.0*PRMS[0]/PRMS[1])) + "\t" +
//                            PRMS[0] + "\t" + PRMS[1] + "\t" + PRMS[2] + "\t" + PRMS[3] + "\t" + PRMS[4] + "\t" + PRMS[5]);
            writer.println(imageName.replace(".raw", "") + "," + maskName.replace(".raw", "") + "," + String.format("%.2f",angle[k]) + "," + String.format("%.2f",(100.0*PRMS[0]/PRMS[1])) + "," +
                           PRMS[0] + "," + PRMS[1] + "," + PRMS[2] + "," + PRMS[3] + "," + PRMS[4] + "," + PRMS[5]);
            
            writer.close();
            
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        //Plot Graphics: Optional
        if(flag_plot){
            LFEChart chart = new LFEChart("LFE", P_HISTO, P_GAUSS, P_LAP, PRMS, H);
            savePlot(chart, name);
            chart.dispose();
        }
        
    }
    
    //implemented for initial guess
    public static double lfe_func_gauss(double[] X, double[] P_GAUSS){
        
        double EPS = 1.0E-09;
        //P_GAUSS[0:NH-1] = (GAUSS_PDF((H_HI-X[0])/X[1])-GAUSS_PDF((H_LO-X[0])/X[1])) > 0.0
        double sum=0;
        for(int i=0; i<NH; i++){
            P_GAUSS[i] = cdf((H_HI[i]-X[0])/X[1])-cdf((H_LO[i]-X[0])/X[1]);
            sum += P_GAUSS[i];
        }
        
        P_GAUSS[NH] = 1.0 - sum;
        
        //RELENT = TOTAL(P_H*ALOG((P_H+EPS)/(P_GAUSS+EPS)))/ALOG(2.0)
        double RELENT=0;
        for(int i=0; i<P_H.length; i++)
            RELENT+= (P_H[i]*Math.log((P_H[i]+EPS)/(P_GAUSS[i]+EPS)))/Math.log(2.0);
        
        return RELENT;
        
    }
    
    //implemented for initial guess
    public static double lfe_func_lap(double[] X, double[] P_LAP){
        
        double EPS = 1.0E-09;
        double[] Q_LO = new double[NH];
        for(int i=0; i<Q_LO.length; i++) Q_LO[i] = 1.0; //init Q_LO
        
        ArrayList<Integer> INDS_LO = new ArrayList<>(); //get lo indices
        for(int i=0; i< H_LO.length; i++){
            if(H_LO[i] < X[0])
                INDS_LO.add(i);
        }
        
        for(int i=0; i<INDS_LO.size(); i++)
            Q_LO[INDS_LO.get(i)] = 0.0;     //zero lo indices
        
        double[] Q_HI = new double[NH];
        for(int i=0; i<Q_HI.length; i++) Q_HI[i] = 1.0; //init Q_HI
        
        ArrayList<Integer> INDS_HI = new ArrayList<>();
        for(int i=0; i< H_HI.length; i++){
            if(H_HI[i] < X[0])
                INDS_HI.add(i);
        }
        
        for(int i=0; i<INDS_HI.size(); i++)
            Q_HI[INDS_HI.get(i)] = 0.0;     //zero hi indices
        
        double[] CDF_LO = new double[NH];
        double[] CDF_HI = new double[NH];
        
        double sum=0;
        for(int i=0; i< CDF_LO.length; i++){
            CDF_LO[i] = 0.5*Math.exp((H_LO[i] - X[0])/X[1])*(1.0-Q_LO[i]) + (1.0-0.5*Math.exp(-(H_LO[i] - X[0])/X[1]))*Q_LO[i];
            CDF_HI[i] = 0.5*Math.exp((H_HI[i] - X[0])/X[1])*(1.0-Q_HI[i]) + (1.0-0.5*Math.exp(-(H_HI[i] - X[0])/X[1]))*Q_HI[i];
            P_LAP[i] = CDF_HI[i] - CDF_LO[i];
            sum += P_LAP[i];
        }
        
        P_LAP[NH] = 1.0 - sum;
        
        double RELENT=0;
        for(int i=0; i<P_LAP.length; i++)
            RELENT+= (P_LAP[i]*Math.log((P_LAP[i]+EPS)/(P_H[i]+EPS)))/Math.log(2.0);
        
        return RELENT;
        
    }
    
    // return pdf(x) = standard Gaussian pdf
    public static double pdf(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }

    // return pdf(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
    public static double pdf(double x, double mu, double sigma) {
        return pdf((x - mu) / sigma) / sigma;
    }
    
    // return cdf(z, mu, sigma) = Gaussian cdf with mean mu and stddev sigma
    public static double cdf(double z, double mu, double sigma) {
        return cdf((z - mu) / sigma);
    }
    
    // return cdf(z) = standard Gaussian cdf using Taylor approximation
    public static double cdf(double z) {
        if (z < -8.0) return 0.0;
        if (z >  8.0) return 1.0;
        double sum = 0.0, term = z;
        for (int i = 3; sum + term != sum; i += 2) {
            sum  = sum + term;
            term = term * z * z / i;
        }
        return 0.5 + sum * pdf(z);
    }
    
    public static void savePlot(LFEChart chart, String name) {
        chart.pack();
        Container c = chart.getContentPane();
        BufferedImage im = new BufferedImage(c.getWidth(), c.getHeight(), BufferedImage.TYPE_INT_ARGB);
        c.paint(im.getGraphics());
        try {
            ImageIO.write(im, "PNG", new File(name+".png"));
        } catch (IOException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private static ImageProcessor doForwardFFT(ImageProcessor ip, int size){
        
        int shift_x = (int)(size/2 - ip.getWidth()/2);
        int shift_y = (int)(size/2 - ip.getHeight()/2);
        
        //padding with zeros to fft size
        FloatProcessor fp = new FloatProcessor(size, size);
        for(int i=0; i<ip.getWidth(); i++){
            for(int j=0; j<ip.getHeight(); j++)
                fp.putPixelValue(i+shift_x, j+shift_y, ip.getf(i, j));
        }
        
        //init complex values
        float [] complex = new float[size*size*2];
        float [] padded = (float[]) fp.getPixels();
        int count =0;
        for(int i=0; i<padded.length; i++){
            complex[count] = padded[i];
            complex[count+1] = 0;
            count+=2;
        }
        
        FloatFFT_2D fft2d = new FloatFFT_2D(size, size);
        fft2d.complexForward(complex);
        
        fp = new FloatProcessor(complex.length/2, complex.length/2);
        fp.setPixels(complex);
        
        return fp;
        
    }
    
    private static ImageProcessor multiply(ImageProcessor ip1, ImageProcessor ip2, int size, int width, int height){
        
        ImageProcessor ip3=null;
        
        float[] pix1 = (float[]) ip1.getPixels();
        float[] pix2 = (float[]) ip2.getPixels();
        float[] pix3 = new float[pix1.length]; //result
        
        //(a+bi)*(c+di) = (ac-bd) + (bc+ad)i
        for(int i=0; i<pix1.length; i+=2){
            pix3[i] = (pix1[i]*pix2[i] - pix1[i+1]*pix2[i+1]);
            pix3[i+1] = (pix1[i+1]*pix2[i] + pix1[i]*pix2[i+1]);
        }
        
        FloatFFT_2D fft2d = new FloatFFT_2D(size, size);
        fft2d.complexInverse(pix3, true);
        
        int count =0;
        float[] padded = new float[size*size];
        for(int i=0; i<padded.length; i++){
            padded[i] = pix3[count] + pix3[count+1];
            count +=2;
        }
        
        FloatProcessor fp = new FloatProcessor(size, size, padded);
        swapQuadrants(fp);
        
        //croping
        Roi roi = new Roi((int)(size/2 - width/2), (int)(size/2 - height/2), width, height);
        fp.setRoi(roi);
        
        ImagePlus spec = new ImagePlus("FFT Convolution", fp.crop());
        //spec.show();
            
        return fp.crop();
    }
    
    public static void swapQuadrants (ImageProcessor ip) {
        
            ImageProcessor t1, t2;
            int size = ip.getWidth()/2;
            ip.setRoi(size,0,size,size);
            t1 = ip.crop();
            ip.setRoi(0,size,size,size);
            t2 = ip.crop();
            ip.insert(t1,0,size);
            ip.insert(t2,size,0);
            ip.setRoi(0,0,size,size);
            t1 = ip.crop();
            ip.setRoi(size,size,size,size);
            t2 = ip.crop();
            ip.insert(t1,size,size);
            ip.insert(t2,0,0);
            
    }
}