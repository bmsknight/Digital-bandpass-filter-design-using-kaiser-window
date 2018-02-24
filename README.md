# Digital-bandpass-filter-design-using-kaiser-window

The code implements an FIR filter using MATLAB for prescribed specifications using the windowing method in conjunction with the Kaiser window.

* Maximum passband ripple 0.13 dB

* Minimum stopband attenuation  56 dB

* Lower passband edge 300  rad/s

* Upper passband edge 600  rad/s

* Lower stopband edge 200  rad/s

* Upper stopband edge 750  rad/s

* Sampling frequency  2000  rad/s

1. Using the windowing method in conjunction with the Kaiser window, an
FIR bandpass digital filter is designed that will satisfy the specifications given above.

2. The magnitude response of the digital filter obtained for the frequency range
0 to ωs/2 rad/s is plotted.

3.The amplitude response for frequencies in the passband is plotted.

4. The impulse response is plotted.

5. The operation of the filter is checked by plotting the time-domain response of the digital filter to an excitation
![equation](http://i64.tinypic.com/2tnyf.png)

where
ω1 is the middle frequency of the lower stopband,
ω2 is middle frequency of the passband, and
ω3 is middle frequency of the upper stopband. 
The response  achieved is compared with that expected if an ideal bandpass filter is used, i.e.,
one that has a gain of 1 in the passband and 0 in the stopbands. 

6. By using the DFTs of the input and the output signals (estimated using an FFT algorithm), it is demonstrated that the output signal is a filtered version of the
input signal and that the correct frequencies have been passed.
