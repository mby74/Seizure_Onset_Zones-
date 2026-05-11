# Identifying Seizure Onset Zones via Nonlinear Dynamical Geometry of iEEG Signals

**UMKC – Saint Luke’s Hospital**  
**Date:** April 4, 2026  

---

## Abstract

Seizure onset zone (SOZ) localization remains a central challenge in epilepsy research and clinical decision-making. Conventional machine learning approaches often rely on large labeled datasets and lack interpretability, limiting their clinical applicability.  

This repository presents a mechanism-driven nonlinear dynamical systems framework for analyzing intracranial EEG (iEEG) data. By reconstructing latent brain dynamics via delay embedding and quantifying the geometry of the resulting attractor, we identify early transitions from complex (chaotic) to ordered (ictal-like) states.  

The results suggest that specific grid electrodes, particularly G7, exhibit the earliest transition, followed by rapid recruitment of mesial temporal and orbitofrontal regions. This supports a hypothesis of localized neocortical onset with rapid network propagation.

---

## Background and Significance

Traditional AI/ML approaches in epilepsy often behave as black boxes, requiring large annotated datasets and offering limited interpretability. These methods frequently fail to generalize across patients due to variability in brain dynamics and recording conditions.  

In contrast, the brain can be modeled as a nonlinear dynamical system, where seizures correspond to phase transitions in state space. This framework allows us to move beyond pattern recognition and instead investigate the underlying mechanism of seizure generation.  

The central idea is that normal brain activity is high-dimensional and chaotic, whereas seizure activity corresponds to a more ordered, lower-dimensional dynamical regime. Detecting the earliest transition toward this ordered state provides both temporal and spatial localization of seizure onset.

---

## Mathematical Framework

### Delay Embedding (Takens’ Theorem)

The hidden neural state is reconstructed using delay coordinates:

$$
\mathbf{X}(t) = [x(t), x(t-\tau), x(t-2\tau), \dots, x(t-(m-1)\tau)]
$$

This produces a phase-space representation that preserves the topology of the underlying system.

---

### Geometric Features

The reconstructed attractor is quantified using the following measures.

The participation dimension measures how many effective directions the data occupy in state space:

$$
D_p = \frac{\left(\sum_i \lambda_i\right)^2}{\sum_i \lambda_i^2}
$$

The correlation dimension characterizes scaling behavior of pairwise distances:

$$
C(r) \sim r^{D_2}
$$

The Lyapunov proxy measures local divergence of nearby trajectories:

$$
d(t) \approx d(0)e^{\lambda t}
$$

Determinism quantifies the predictability of recurrence structure:

$$
DET = \frac{\text{diagonal recurrence points}}{\text{total recurrence points}}
$$

The radius measures geometric spread of the attractor:

$$
R = \sqrt{\frac{1}{N}\sum_{i=1}^{N}\left\|\mathbf{X}_i - \bar{\mathbf{X}}\right\|^2}
$$

---

### Order Score

To quantify seizure-like ordering, we define:

$$
S = -z(R) - z(D_p) - z(D_2) - z(\lambda) + z(DET)
$$

Higher values indicate more ordered dynamics, corresponding to seizure-like behavior.

---

### Onset Detection

For each channel \(c\), the candidate onset time is defined as:

$$
t_c = \min \left( t \;:\; \Delta S_c(t) > \theta_c \right)
$$

Channels are ranked based on the earliest transition toward this ordered state.

---

## Methods

The dataset consists of multichannel intracranial EEG (iEEG) recordings obtained at Saint Luke’s Hospital and processed at UMKC. The data include recordings from grid electrodes, anterior hippocampal (AH), cingulate (Cing), posterior hippocampal (PH), and orbitofrontal (BO) regions.  

The sampling rate is approximately 500 Hz, and the recording duration is approximately 2 seconds, yielding a short but high-resolution window suitable for exploratory dynamical analysis. A total of 70 channels were analyzed.  

Each channel was treated as a one-dimensional observable of an underlying high-dimensional neural system. Using delay embedding, we reconstructed the attractor for each channel. Geometric features were computed over sliding time windows, allowing us to track the evolution of dynamical structure over time.  

The order score was computed for each channel and window, smoothed, and used to detect transitions toward more ordered dynamics. Candidate onset times were defined as the earliest significant increase in the order score.

---

## Results

### Figure 1: Time Series of Top Candidate Channels

![Figure 1](figure1_top_channel_timeseries.png)

**Figure 1.** Normalized raw signals from the top-ranked channels (G7, G24, G30, G31, G36, AH6).

These time series show structured, non-random dynamics. While not strictly periodic, the signals exhibit organized fluctuations consistent with deterministic behavior rather than noise. This suggests that the underlying system is governed by nonlinear dynamics that can be reconstructed in phase space.

---

### Figure 2: Delay-Embedded Attractor (G7)

![Figure 2](figure2_top_channel_attractor.png)

**Figure 2.** Delay-embedded attractor for channel G7 using a 3D reconstruction.

The attractor exhibits a structured but irregular geometry, indicating a system that is neither random nor purely periodic. This intermediate structure is characteristic of nonlinear dynamical systems approaching a transition. The geometry suggests evolving constraints in the system's degrees of freedom.

---

### Figure 3: Order Score Heatmap

![Figure 3](figure3_order_score_heatmap.png)

**Figure 3.** Heatmap of order score across channels and time (higher values indicate more ordered dynamics).

This figure shows the spatiotemporal evolution of dynamical ordering. Early increases in order score are concentrated in specific grid channels, followed by propagation to other regions. This pattern supports the hypothesis of localized onset followed by network recruitment.

---

### Figure 4: Feature Trajectories for Top Channel (G7)

![Figure 4](figure4_top_channel_feature_trajectories.png)

**Figure 4.** Temporal evolution of geometric features and order score for channel G7.

The participation and correlation dimensions decrease near the transition, indicating reduced complexity. Determinism increases prior to the transition, reflecting more predictable dynamics. The Lyapunov proxy shows fluctuations consistent with changing stability. The order score rises sharply around 1.5–1.6 seconds, marking the transition toward a more ordered state.

---

### Figure 5: Candidate Onset Ranking by Channel

![Figure 5](figure5_onset_ranking.png)

**Figure 5.** Sorted candidate onset times across all channels.

The earliest transitions occur in channels G7, G24, G30, G31, and G36, all at approximately 0.1991 seconds. These channels represent the strongest candidates for seizure onset based on dynamical criteria.

---

### Figure 6: Group-Level Onset Summary

![Figure 6](figure6_group_onset_summary.png)

**Figure 6.** Mean onset time by anatomical group.

At the group level, anterior hippocampal (AH) channels show the earliest average onset, followed by cingulate and posterior hippocampal regions. This indicates early involvement of mesial temporal structures at the network level.

---

## Discussion

The results demonstrate that seizure onset can be identified as a transition toward a more ordered dynamical state, characterized by reduced dimensionality, reduced divergence, and increased determinism.  

At the channel level, the earliest transitions occur in specific grid electrodes, particularly G7. This suggests a focal neocortical origin of seizure activity. However, group-level analysis shows that mesial temporal structures, especially the anterior hippocampus, are involved very early on average. This indicates rapid propagation from a focal onset region to a broader network.  

These findings highlight the importance of distinguishing between **earliest individual channel activation** and **average group-level involvement**. The former suggests the initial source, while the latter reflects network recruitment dynamics.  

The study is limited by the short duration of the recording and the absence of clinically annotated seizure onset times. As a result, the findings should be considered exploratory. Future work should include longer recordings, multiple seizures, and integration with clinical EEG interpretation and anatomical mapping.  

Despite these limitations, the approach demonstrates that nonlinear dynamical geometry provides an interpretable and data-efficient framework for seizure analysis. Unlike black-box machine learning models, this method directly reveals where, when, and how seizures emerge in the brain.

---

## Citation

If you use this work, please cite:

Bani Yaghoub, M. (2026)  
Identifying Seizure Onset Zones via Nonlinear Dynamical Geometry of iEEG Signals  
UMKC – Saint Luke’s Hospital
