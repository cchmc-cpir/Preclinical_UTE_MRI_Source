# Instructions for using pjn_MUTE

This code is based on the native sequence, UTE3D. In addition to the standard functionality of that sequence, we have added a variety of useful tools, most of which are outlined below.
The primary differences between UTE3D and this sequence are:
* Golden Means Acquisition Ordering
* Collection of points prior to gradients in order to mitigate digital filter corruption (These points are discarded in reconstruction).

A comprehensive reconstruction package for this sequence is contained in this repository, in the folder Preclinical_UTE_MRI_Source/pjn_MUTE/Recon

## Echo Time Considerations
* The TE that is set in the routine card accounts for the extra points acquired at the beginning of readout before gradients are turned on. Digital filter corruption endures for ~70 us, so the minimum TE should be set no shorter than that value.
* Multiple TEs can be imaged by increasing the Number of TEs parameter.
* Multiple TE images are acquired in an interleaved fashion, with projection 1 acquired for TE1, TE2, ..., TEN prior to moving on to projection 2

## Considerations for setting Imaging Parameters
* Bruker, unfortunately, has issues with digital filter corruption, which causes issues for quantitative imaging. The way we work around that is by delaying when we turn on gradients and acquiring some points in the absence of gradients. Then, those extra points at the beginning can be thrown out. Importantly, the digital filter corruption lasts for a certain amount of time, so the number of points you use is going to depend on the bandwidth. The corruption seems to last ~70 µs, so I usually try to have at least 80 µs worth of Acq shift points (i.e. minimum echo time of ~80 µs).
* Extra points: You can specify a number of “extra points” to acquire at the end of your radial acquisition. This protects from accidental undersampling due to gradient delays.
*	Spoiling is done by specifying how many multiples of the readout gradient you want to spoil with. This is found in the sequence tab, and is specified by the parameters “Spoiling”, and “Spoiler Amp”. For example, if your readout gradient has an area A, and you specify that you want 400% spoiling, your spoiler gradient will have an area 4×A. “Spoiler Amp” determines the amplitude of this gradient – Higher amplitude leads to faster spoiler.

## Triggering
* This sequence does not include the proper code to use continuous RF pulsing in between gating periods. As such, it is not possible to maintain steady state magnetization during a gated acquisition, which makes quantitative imaging impossible.
* It is recommended to acquire images without gating and reconstruct using retrospective gating techniques - A basic method is included in the associated reconstruction package.
