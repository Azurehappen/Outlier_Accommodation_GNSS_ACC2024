# Outlier Accommodation for Multi-GNSS Precise Point Positioning using Risk-Averse Performance-Specified (RAPS) Approach

An extension of [Multi GNSS Repo](https://github.com/Azurehappen/Multi_GNSS_PPP_Tool)

This repo is the implementation of the paper that has been accepted by the 2024 American Control Conference.

To run this implementation, please run 'MultiGNSS_Single_Main.m'. You can also change the estimator via 'p.est_mode' where the parameter for outlier accommodation using RAPS is 'p.raps_ned_est'.

Targeted GNSS signals:
GPS: L1
GAL: E1
BDS: B1I
