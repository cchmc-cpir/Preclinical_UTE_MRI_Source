# Preclinical_UTE_MRI_Source
Source code for radial UTE sequences developed for the CCHMC Bruker 7T animal MRI system.

Author: Peter J. Niedbalski

Contact: pniedbalski@kumc.edu

## cpir_UTE_V5
This sequence is a 3D radial MRI sequence that uses Golden Means Projection Ordering to encode images. Several layers of functionality have been added, including
* Shifted Acquisition Window to mitigate points corrupted by digital filter
* Multiple, interleaved TE acquisition
* Multiple TE acquisition using Radial Flyback acquisition
* Diffusion Weighted Imaging
* Steady-state RF pulsing during gated acquisition
* Hyperpolarized 129Xe MRI, including adpated gating methods for imaging while animals are on a ventilator

## pjn_MUTE
This sequence is a 3D radial MRI sequence that uses Golden Means Projection Ordering to encode images. This is a more basic sequence that has less added functionality. This sequence is designed to be used for structural imaging of 1H in the lungs using a non-gated acquisition. Images can be subsequently retrospectively gated using any of a variety of post-acquisiton techniques. Specific functionality includes:
* Shifted Acquisition Window to mitigate points corrupted by digital filter
* Multiple, interleaved TE acquisition
