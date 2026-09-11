# Sophisticated-Biomedical-signal-Processing-of-Brain-EEG
Signal Processing GUI and MATLAB code
# EEG Signal Processing and Analysis GUI

## MSc Project

**Project Title:** Sophisticated Biomedical Signal Processing of Brain Electroencephalogram (EEG) Data

**Author:** Oghenekome Allison Uvo-Akpos  
**Degree:** MSc Embedded Systems and IC Design  
**Institution:** Liverpool John Moores University  
**Year:** 2026

---

## Overview

This repository contains the MATLAB application and validation scripts developed for the MSc project **“Sophisticated Biomedical Signal Processing of Brain Electroencephalogram (EEG) Data.”**

Electroencephalogram (EEG) signals are non-stationary biomedical signals whose frequency characteristics can vary over time. The project investigates and compares three complementary signal-processing techniques for analysing and visualising EEG data:

- Filter–Hilbert analysis
- Short-Time Fourier Transform (STFT)
- Multitaper spectral analysis using Discrete Prolate Spheroidal Sequences (DPSS)

The techniques are integrated into a MATLAB App Designer graphical user interface (GUI), allowing EEG data to be loaded, processing parameters to be configured, results to be visualised, and the outputs of the different methods to be compared.

---

## Repository Structure

```text
EEG-Signal-Processing-GUI/
│
├── EEGAnalysis.mlapp
├── README.md
│
└── validation/
    ├── validation_hz.m
    └── figure11_timevarying_validation.m
