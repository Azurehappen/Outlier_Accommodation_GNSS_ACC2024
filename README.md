# Outlier Accommodation for Multi-GNSS Precise Point Positioning using Risk-Averse Performance-Specified (RAPS) Approach

An extension of [Multi GNSS Repo](https://github.com/Azurehappen/Multi_GNSS_PPP_Tool)

This repo is the implementation of the paper that has been accepted by the 2024 American Control Conference.

If you are interested in this work, please cite our [paper](https://arxiv.org/abs/2402.01860)
```
@misc{hu2024rapsppp,
	title={Outlier Accommodation for GNSS Precise Point Positioning using Risk-Averse State Estimation}, 
	author={Wang Hu and Jean-Bernard Uwineza and Jay A. Farrell},
	year={2024},
	eprint={2402.01860},
	archivePrefix={arXiv}
}
```

To run this implementation, please run 'MultiGNSS_Single_Main.m'. You can also change the estimator via 'p.est_mode' where the parameter for outlier accommodation using RAPS is 'p.raps_ned_est'.
(More details about running this implementation will be filled in soon.....)

Targeted GNSS signals:
GPS: L1
GAL: E1
BDS: B1I
