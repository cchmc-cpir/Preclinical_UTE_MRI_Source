# Instructions for using cpir_UTE_V5

This code is based on the native sequence, UTE3D. In addition to the standard functionality of that sequence, we have added a variety of useful tools, most of which are outlined below.
The primary differences between UTE3D and this sequence are:
* Golden Means Acquisition Ordering
* Collection of points prior to gradients in order to mitigate digital filter corruption (These points are discarded in reconstruction).
* Continuous RF pulsing during respiratory gating in order to maintain steady state magnetization

A comprehensive reconstruction package for this sequence can be found at: https://github.com/cchmc-cpir/preclinical-reconstruction

## Echo Time Considerations
* The TE that is set in the routine card accounts for the extra points acquired at the beginning of readout before gradients are turned on. Digital filter corruption endures for ~70 us, so the minimum TE should be set no shorter than that value.
* Multiple TEs can be imaged by increasing the Number of TEs parameter.
* By default, multiple TE images are acquired in an interleaved fashion, with projection 1 acquired for TE1, TE2, ..., TEN prior to moving on to projection 2
* If "Flyback" is checked, multiple TE images are instead acquired by reversing gradients at the end of the radial readout. In this fashion, all desired echo times are acquired with each RF pulse.
* When "Flyback" is checked, TEs are calculated based on image size, bandwidth, and gradient limitations. If "Flyback" is not checked, the user can select their desired TE values.

## Imaging Diffusion
* Using the parameters in the Contrast->Diffusion Tab, you can acquire multiple diffusion weighted images in an interleaved fashion.
* Note that b-value units are in s/cm^2 to be more applicable to xenon imaging. b-values for 1H imaging are more commonly reported in s/mm^2

## Considerations for setting Imaging Parameters
* Bruker, unfortunately, has issues with digital filter corruption, which causes issues for quantitative imaging. The way we work around that is by delaying when we turn on gradients and acquiring some points in the absence of gradients. Then, those extra points at the beginning can be thrown out. Importantly, the digital filter corruption lasts for a certain amount of time, so the number of points you use is going to depend on the bandwidth. The corruption seems to last ~70 µs, so I usually try to have at least 80 µs worth of Acq shift points (i.e. minimum echo time of ~80 µs).
* Extra points: You can specify a number of “extra points” to acquire at the end of your radial acquisition. This protects from accidental undersampling due to gradient delays.
* Slab Selection: You can turn slab selection on or off. If on, you must set the slab thickness (in mm). When you do this, it is also good to change your pulse shape and duration. The default for radial is to use a block pulse that is really short (~0.005 ms). If doing slab selection, probably change to “calculated”, “gauss”, or “sinc3”, and a longer duration (usually I use a minimum of 0.25 ms).
* Rewinder: You can select whether you want your sequence to use a rewinder (return to k0 at the end of readout) or not. This is under the contrast tab. It defaults to “yes”.
* RF Spoiling is optional (Default "Yes"). This option is found in the contrast tab
*	Spoiling is done by specifying how many multiples of the readout gradient you want to spoil with. This is found in the sequence tab, and is specified by the parameters “Spoiling”, and “Spoiler Amp”. For example, if your readout gradient has an area A, and you specify that you want 400% spoiling, your spoiler gradient will have an area 4×A. “Spoiler Amp” determines the amplitude of this gradient – Higher amplitude leads to faster spoiler.

## Triggering
* If you are imaging 1H, the sequence will automatically trigger off of trigger input 1 (usually an SAI monitoring system), and it will automatically trigger off of trigger input 4 (a homebuilt ventilator, in our case) for 129Xe. 
* When the sequence is triggering off of input 4, you need to specify the number of projections you want to acquire per trigger. That is, the ventilator will give off a short trigger signal, then, the scanner will collect the number of projections that you specify. 
* When using radial for quantitative imaging, you need to maintain a steady-state magnetization, which means that you need to be consistently pulsing even when triggering. That is, you should pulse and acquire data when the trigger is high, and should pulse but not acquire when the trigger is low. However, Bruker doesn’t have a good away to do this, and you need to do some kludge-y stuff in the ppg to make it work. This, along with some of the extra gradients involved in slab selection and diffusion weighting require arbitrary lengthening of the TR. Empirically, it seems that a basic triggered radial needs at least 7 ms for TR. That value is longer for slab selective or diffusion weighted sequences. I think that I’ve coded it such that the system forces you to use a TR that will work. However, if you get an error when running a triggered acquisition, just try lengthening the TR a little bit, and it should run. 
