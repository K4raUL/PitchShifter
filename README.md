# PitchShifter
Multiply sound frequencies in .wav files by a custom coefficient (using https://github.com/jdupuy/dj_fft)

## Usage (gcc)  
`g++ pshifter.cpp -O2 -o shift -w`  
`./shift "sounds/example.wav"` - minimal syntax, coefficient = 1., FFT size = 2048  
`./shift "sounds/example.wav" 1.28` - coefficient = 1.28, FFT size = 2048  
`./shift "sounds/example.wav" 0.714 4096` - coefficient = 0.714, FFT size = 4096

##  
You can try different FFT sizes to increase result quality, but it **must** be power of 2.  
Default FFT size is 2048, because it is close to 50ms chunk on 44100 Hz sample rate (2205).  
FFT size = 4096 often shows better results, but dj_fft algorithms do not always work correctly on audio data  
vectors greater than 2048 elements (result differs from original, doing FFT forward and FFT backward after),  
causing volume loss.  
Program supports different .wav file types (sample rate, sample size, channel count, byte order, signed/unsigned).

