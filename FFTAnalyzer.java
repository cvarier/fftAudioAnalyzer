/**
* A simple audio analyzer which uses the Fast Fourier Transform algorithm

* Author: Chaitanya Varier
* Description: Samples an input signal, computes its FT and its discretized frequencies and power spectrum.
			   The program then filters it by a specified frequency, freq1, and its first two harmonics, 
			   freq2 and freq3, in the frequency domain. Next, the inverse FT is computed from the cleaned
			   signal and used to produce an audio file. Finally, the plots for the input signal and 
			   computed power spectrum are outputted.
*/

import java.io.*;
import java.io.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import org.math.plot.*;
import org.math.plot.plotObjects.*;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.transform.DftNormalization;
import WavFile.*;


public class Fourier {
	
	public static double PI = Math.PI;

	public static void main(String[] Args) {
			
		double duration = 4.0; 
		long sampleRate = 100000; // Samples/s
		int N = (int) (duration*sampleRate); // Number of samples
		double freq1 = 523.0, freq2 = freq1*3, freq3 = freq1*7; // Frequency in s^-1
		double T1 = 1.0 / freq1, T2 = 1.0 / freq2, T3 = 1.0 / freq3; // Period in s
		double dt = 1.0 / sampleRate; // Time step
		double fNyq = 1.0 / (2*dt); // Nyquist frequency
		
		//////////////////////////////////////
		// Sample Input Signal
		//////////////////////////////////////
		
		// Create array to hold the signal values
		double [] time = new double[N];
		double[] signal = new double[N];

		time[0] = 0;
		signal[0] = 0; // sin(0) = 0

		// Sample the signal and fill the arrays
		for (int i = 1; i < N; i++) {
			time[i] = time[i - 1] + dt;
			signal[i] = 10*Math.sin(2*PI / T1*time[i]) + 2*Math.sin(2*PI / T2*time[i]) + 2*Math.sin(2*PI / T3*time[i]);
		}
		
		//////////////////////////////////////
		// Pad Signal
		//////////////////////////////////////
		
		// Calcuate the next highest power of 2 after the signal length
		int n = 0;
		
		while ((int)Math.pow(2,n) < signal.length) {
			n++;
		}
		
		int paddedSize = (int)Math.pow(2,n);
		
		double[] paddedSignal = new double[paddedSize];
		
		// Copy and pad the signal
		for (int i = 0; i < paddedSignal.length; i++) {
			if (i < signal.length){
				paddedSignal[i] = signal[i];
			} else {
				paddedSignal[i] = 0;
			}
		}
		
		//////////////////////////////////////
		// Compute Fourier Transform
		//////////////////////////////////////
		
		// Create a FFT object and obtain the FT of the signal
		FastFourierTransformer FFT = new FastFourierTransformer(DftNormalization.STANDARD);
		Complex[] signalFT = FFT.transform(paddedSignal,TransformType.FORWARD);
		
		//////////////////////////////////////
		// Compute Power spectrum
		//////////////////////////////////////
		
		// Create arrays for signal frequency and power
		double[] frequency = new double[signalFT.length / 2];
		double[] power = new double[signalFT.length / 2];
		
		// Calculate the signal's frequencies and powers
		for (int m = 0; m < frequency.length; m++) {
			frequency[m] = m / dt / signalFT.length;
			power[m] = Math.pow(signalFT[m].getReal(),2) + Math.pow(signalFT[m].getImaginary(),2);
		}
		
		//////////////////////////////////////
		// Clean Signal
		//////////////////////////////////////
		
		// Find the indeces with the frequency closest to freq1 and its next two harmonics
		int freq1I = 0, freq2I = 0, freq3I = 0;
		double minDiffF1 = Math.abs(freq1-frequency[0]);
		double minDiffF2 = Math.abs(freq2-frequency[0]);
		double minDiffF3 = Math.abs(freq3-frequency[0]);
		
		for (int m = 0; m < frequency.length; m++) {
			if (Math.abs(freq1-frequency[m]) < minDiffF1) {
				minDiffF1 = Math.abs(freq1 - frequency[m]);
				freq1I = m;
			}
			if (Math.abs(freq2-frequency[m]) < minDiffF2) {
				minDiffF2 = Math.abs(freq2 - frequency[m]);
				freq2I = m;
			}
			if (Math.abs(freq3-frequency[m]) < minDiffF3) {
				minDiffF3 = Math.abs(freq3 - frequency[m]);
				freq3I = m;
			}
		}
		
		// Remove all the components from signalFT except those that correspond to the above frequencies
		double[] signalFTReClean = new double[signalFT.length];
		double[] signalFTImClean = new double[signalFT.length];
		
		for (int m = 0; m < signalFT.length; m++) {
			if (m != freq1I && m != freq2I && m != freq3I) {
				signalFTReClean[m] = 0;
				signalFTImClean[m] = 0;
				if (m < frequency.length) {
					frequency[m] = 0;
				}
			} else {
				signalFTReClean[m] = signalFT[m].getReal();
				signalFTImClean[m] = signalFT[m].getImaginary();
			}
		}
		
		//////////////////////////////////////
		// Obtain Inverse Fourier Transform
		//////////////////////////////////////
		
		// Create a FFT object and obtain the inverse FT of the clean signal
		FastFourierTransformer FFT2 = new FastFourierTransformer(DftNormalization.STANDARD);
		Complex[] cleanSignalFT = FFT.transform(signalFTReClean,TransformType.INVERSE);
		
		// Obtain the clean signal
		double[] cleanSignalFTWrite = new double[cleanSignalFT.length];
		
		for (int m = 0; m < cleanSignalFT.length; m++) {
			cleanSignalFTWrite[m] = cleanSignalFT[m].getReal();
		}
		
		//////////////////////////////////////
		// Synthesize Sound File
		//////////////////////////////////////
	
		try {
			long sampleRateClean = sampleRate;
			// Number of frames required for the specified duration
			long numFrames = N;
			// Create a wav file with the name specified as the first argument
			WavFile wavFileOut = WavFile.newWavFile(new File("C5FundamentalOvertones.wav"), 1, numFrames, 16, sampleRateClean);
			// Write the buffer
			wavFileOut.writeFrames(cleanSignalFTWrite, cleanSignalFTWrite.length);
			wavFileOut.close();
		}
			
		catch (Exception e) {
			System.err.println(e);
		}
		  
		//////////////////////////////////////		
		// Plotting Signal Power Spectrum
		//////////////////////////////////////

		Font plotFont = new Font(Font.MONOSPACED,Font.PLAIN,12);
		
		Plot2DPanel plot1 = new Plot2DPanel();
		
		// Plot the signal
		plot1.addLinePlot("Signal",Color.magenta,time,signal);
		plot1.setFixedBounds(0,0,signal.length*dt); //Manually set bounds on x-axis	
		plot1.setAxisLabels("Time", "Amplitude");
		plot1.getAxis(0).setLabelPosition(0.5, -0.1);
		plot1.getAxis(0).setLabelFont(plotFont);
		plot1.getAxis(1).setLabelPosition(-0.15,0.5);
		plot1.getAxis(1).setLabelFont(plotFont);
		BaseLabel title1 = new BaseLabel("Signal from " + wavFileName,Color.BLACK,0.5, 1.1);
		title1.setFont(plotFont);
			plot1.addPlotable(title1);

		JFrame frame1 = new JFrame("Output 1");
		frame1.setSize(1024,576);
		frame1.setContentPane(plot1);
		frame1.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		frame1.setVisible(true);

		// Find maximum power
		double PMax = 0.0;
		for (int j = 0; j < power.length; j++) {
			if (power[j] > PMax) {
				PMax = power[j];
			}
		} 
		
		Plot2DPanel plot2 = new Plot2DPanel();
				
		// Plot the power spectrum
		plot2.addLinePlot("Power", frequency, power);
		plot2.setFixedBounds(0,0,fNyq);	//Manually set bounds on x-axis
		plot2.setFixedBounds(1,0,1.01*PMax); // Manually set bounds on y-axis
		plot2.setAxisLabels("Frequency", "Power");
		plot2.getAxis(0).setLabelPosition(0.5,-0.1);
		plot2.getAxis(0).setLabelFont(plotFont);
		plot2.getAxis(1).setLabelPosition(-0.15,0.5);
		plot2.getAxis(1).setLabelFont(plotFont);
		BaseLabel title2 = new BaseLabel("Power spectrum of " + wavFileName, Color.BLACK, 0.5, 1.1);
		title2.setFont(plotFont);
			plot2.addPlotable(title2);
		
		JFrame frame2 = new JFrame("Output 2");
		frame2.setSize(1024,576);
		frame2.setLocationRelativeTo(null);
		frame2.setContentPane(plot2);
		frame2.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		frame2.setVisible(true);
	}
}
