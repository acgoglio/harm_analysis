# harm_analysis 
# repository:
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
# 
# The aim of this code is to perform a validation of EAS system outputs in terms of tidal harmonic analysis. 
#
# The harmonic analysis implemented in this package is based on the Foreman fit method [Foreman et al.; 1991] 
# and employs the ttide package [Pawlowicz et al; 2002].
#
# The validation is composed by several analysis and comparisons including comparisons with respect to TPXO model [Egbert et al; 2002] 
# and with respect to results from the following papers [Palma et al. 2020], [Tsimplis et al; 1995], [Agresti; 2018] and [Arabelos et al,2011] 
#
# The validation procedure, along with the code structure, is diveded in two sections:
# 1) The PUNCTUAL harmonic analysis, concerning tide-gauge mod/obs/literature comparisons in tide-gauge locations
# 2) The AREAL harmonic analysis dealing with the comparison on the whole model domain
#
# For both this sections is provided a specific description of tasks, usage and outputs.
# The tasks of each section are summarized in the following lines:
#
# PUNCTUAL TIDAL HARMONIC ANALYSIS 
# ==============================
# TASKS:
#
# - Pre-processing of tide-gauges sea level data from EMODnet (counts and fills the GAPs in the time series; removes time-series with a num of gaps higher then a user-fixed threshold )
# - Diagnostic analysis on tide-gauge data (plot comparisons between modeled and observed time-series, spectra and empirical distribution; quality flag checks;)
# - Evaluation of harmonic analysis fit skill through signal-to-noise analysis and comparison between two different fit methods
# - Maps and Lists the resulting tide-gauges datasets 
# - Computation of Amplitude and Phase in tide-gauges location by means of Foreman fit method 
# - mod/obs/TPXO model comparison in terms of harmonic analysis results 
# - Comparisons of the results with respect to literature and to TPXO model in terms of 
#    - amplitude differences statistics
#    - phase differences statistic 
#    - linear regressions
#    - vectorial distances
#    - RMSmisfits 
#
#
# AREAL TIDAL HARMONIC ANALYSIS 
# ==============================
# TASKS:
#
# - Extraction of sea level field time-serie on the whole domain from EAS system hourly outputs
# - Areal harmonic analysis: Foreman fit applied to the whole domain
# - Visualization of the areal harmonic analysis results (Amplitude and phase maps )
#   (Includes the possibility to apply the same palette as literature [Palma et al., Agresti et al., Arabelos et al.])
# - Comparisons with respect to TPXO model in terms of:
#    - Maps of TPXO model amplitude and phase on the TPXO grid
#    - Maps of TPXO model amplitude and phase on the EAS system grid (TPXO is interpolated on the med 24 grid)
#    - Amplitude difference (Maps per tidal component) 
#    - Root Mean square differences between amplitudes computation (for points with bathymetry higher than a user-fixed threshold)
#    - Vectoria differences (Maps per tidal component)
#    - Root Mean Square misfits computation
#    - Do-Seong F, E and Ea Envelope Form Factors (Do-Seong; 2020) (Maps for both datasets)
#    - Comparison between eas system and TPXO model bathymetries 
#
