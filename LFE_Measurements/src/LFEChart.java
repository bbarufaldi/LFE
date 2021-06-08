
import java.awt.Color;
import java.util.Arrays;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Administrator
 */
public class LFEChart extends ApplicationFrame {
        
        public LFEChart(String applicationTitle, double[] P_HISTO, double[] P_GAUSS, double[] P_LAP, double[] PRMS, double[][] H) {
        super(applicationTitle);
        
            //Plot serie
            XYSeries histo = new XYSeries("P_HISTO");
            XYSeries gauss = new XYSeries("P_GAUSS");
            XYSeries lap = new XYSeries("P_LAP");
            
            double[] x_range = new double[P_HISTO.length];
            for(int i=0; i<P_HISTO.length; i++)
                x_range[i] = H[i][0];
            
            Arrays.sort(x_range);

            for(int i=0; i<P_HISTO.length; i++){
                histo.add(P_HISTO[i], x_range[i]);
                gauss.add(P_GAUSS[i], x_range[i]);
                lap.add(P_LAP[i], x_range[i]);
                
                //System.out.println(x_range[i] + "\t" + P_HISTO[i] + "\t" + P_GAUSS[i] + "\t" + P_LAP[i]);
            }

            XYSeriesCollection dataset = new XYSeriesCollection();
            dataset.addSeries(histo);
            dataset.addSeries(gauss);
            dataset.addSeries(lap);

            //get max
            double P_MAX=0;
            for(int i=0; i<P_HISTO.length; i++){
                double max_histo = P_HISTO[i];
                double max_gauss = P_GAUSS[i];
                double max_lap = P_LAP[i];

                if(P_MAX < max_histo) P_MAX = max_histo;
                if(P_MAX < max_gauss) P_MAX = max_gauss;
                if(P_MAX < max_lap) P_MAX = max_lap;
            }

            JFreeChart xylineChart = ChartFactory.createXYLineChart("Distribution Plots, LFE = "+String.format("%.2f",(100.0*PRMS[0]/PRMS[1]))+"%",
                                                                "Bin Probability", "Bin Position",
                                                                dataset,
                                                                PlotOrientation.HORIZONTAL,
                                                                true,true,false);
            ChartPanel chartPanel = new ChartPanel(xylineChart);
            
            chartPanel.setPreferredSize( new java.awt.Dimension(750 , 500));
            XYPlot plot = xylineChart.getXYPlot();
            ValueAxis domain = plot.getDomainAxis();
            //ValueAxis range = plot.getRangeAxis();
            domain.setRange(0.00, Math.round(P_MAX*100.0)/100.0); //set Y range
            //range.setRange(0, P_HISTO.length-2); //set X range

            XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
            renderer.setSeriesPaint(0 , Color.RED);
            renderer.setSeriesPaint(1 , Color.BLUE);
            renderer.setSeriesPaint(2 , Color.MAGENTA);
            
            //renderer.setSeriesStroke(0 , new BasicStroke( 4.0f ));
            //renderer.setSeriesStroke(1 , new BasicStroke( 3.0f ));
            //renderer.setSeriesStroke(2 , new BasicStroke( 2.0f ));
            renderer.setSeriesLinesVisible(0, false);
            renderer.setSeriesLinesVisible(1, false);
            renderer.setSeriesLinesVisible(2, false);
            
            plot.setRenderer(renderer); 
            setContentPane(chartPanel);

            //chartPanel.setVisible(true);
            
        }
}
