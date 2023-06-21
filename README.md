# Filtered Backpropagation for Coherent Plane-Wave Compounding

<!-- shields -->
[![GitHub][license-shield]][license-url]
![GitHub repository size][size-shield]
![GitHub language count][languages-shield]
![GitHub stargazers][stars-shield]
![GitHub forks][forks-shield]
[![ko-fi][ko-fi-shield]][ko-fi-url]

[license-shield]: https://img.shields.io/badge/license-citationware-blue
[license-url]: https://github.com/mschiffn/fbp/blob/main/LICENSE
[size-shield]: https://img.shields.io/github/repo-size/mschiffn/fbp
[languages-shield]: https://img.shields.io/github/languages/count/mschiffn/fbp
[stars-shield]: https://img.shields.io/github/stars/mschiffn/fbp.svg
[forks-shield]: https://img.shields.io/github/forks/mschiffn/fbp.svg
[ko-fi-shield]: https://img.shields.io/badge/ko--fi-Donate%20a%20coffee-yellowgreen
[ko-fi-url]: https://ko-fi.com/L4L7CCWYS

<!-- content -->
Simple
[MATLAB](https://mathworks.com/products/matlab.html) implementation of
the filtered backpropagation
[[1]](#SchiffnerAI2012),
[[2]](#DevaneyUI1982) for
coherent plane-wave compounding.

## Background

The filtered backpropagation (FBP) derives from
the Fourier diffraction theorem for
steered plane waves.
This theorem resembles the Fourier slice theorem.
The algorithm operates in
the two-dimensional space under
the Born approximation.

## :gear: Getting Started

1. Clone the repository or download the release to your local hard drive.

```
git clone https://github.com/mschiffn/fbp
```

2. Switch to the directory containing the repository

3. Run the script example.m

## Folder Structure

The repository has the following structure:

    .
    ├── data_RF.mat     # measurement data from tissue phantom
    ├── example.m       # main script
    ├── fbp_pw.m        # filtered backpropagation (FBP) algorithm
    ├── LICENSE         # license file
    └── README.md       # this readme

## :notebook: References

1. <a name="SchiffnerAI2012"></a>
M. F. Schiffner and G. Schmitz,
"Plane wave pulse-echo ultrasound diffraction tomography with a fixed linear transducer array,"
Acoust. Imaging, ser. Acoust. Imaging, A. Nowicki, J. Litniewski, and T. Kujawska, Eds., vol. 31, Springer Netherlands, 2012, pp. 19–30.
[![DOI:10.1007/978-94-007-2619-2_3](https://img.shields.io/badge/DOI-10.1007%2F978--94--007--2619--2__3-blue)](https://doi.org/10.1007/978-94-007-2619-2_3)

2. <a name="DevaneyUI1982"></a>
A. J. Devaney,
"A filtered backpropagation algorithm for diffraction tomography,"
Ultrasonic Imaging, vol. 4, no. 4, pp. 336–350, Oct. 1982.
[![DOI:10.1016/0161-7346(82)90017-7](https://img.shields.io/badge/DOI-10.1016%2F0161--7346(82)90017--7-blue)](https://doi.org/10.1016/0161-7346(82)90017-7)

