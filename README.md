## Description

This is a simple audio analyzer written in Java which uses the **Fast Fourier Transform algorithm** to filter the harmonics of a specified frequency from an input audio signal and compute its power spectrum. Here's a step-by-step description of what the program does:

1. The program first samples an input signal, computes its FT, its discretized frequencies and power spectrum.
2. It then filters it by a specified frequency, `freq1`, and its first two harmonics, `freq2` and `freq3`, in the frequency domain. 
3. Next, the inverse FT is computed from the cleaned signal and used to produce an audio file. 
4. Finally, the plots for the input signal and computed power spectrum are outputted.

## Screenshots

Here's an example of a signal the program would read:

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/42006370-6571f06e-7a47-11e8-982e-0d7b4608a7f8.png" width = "600" height = "300"/>

The program will then output the corresponding power spectrum of that signal from its Fourier Transform:

<p align="center"><img src ="https://user-images.githubusercontent.com/16710726/42006332-35564efc-7a47-11e8-9a6b-00bfeff5e059.png" width = "600" height = "300"/>
