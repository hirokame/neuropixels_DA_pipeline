# neuropixels_DA_pipeline
The pipeline for the Neuropixels (spikeGLX) data and DA photometry (TDT) data.


## Analysis Pipeline Overview

Primary workflow (spike sorting → waveform metrics → classification):
1) Spike sorting: run `SpikeSorting.ipynb` to generate Kilosort outputs (`kilosort4` folder).
2) Analyzer creation: run `pipeline/spikeinterface_waveform_extraction/extract_si2.ipynb` to build `analyzer_beh` / `analyzer_tag` and export spikeinterface metrics/waveforms.
3) Waveform metrics export: run `savewaveform_info.ipynb` to compute and save `template_metrics.csv`, ACC/ISI files, and waveform npys into each session’s `kilosort4` directory.
4) Cell classification: run `Cellclassify.ipynb` to ingest the saved metrics/waveforms and perform unit tagging/quality analysis.

Photometry (dopamine) workflow:
- Run `Matlab4TDT/TDT_test.m` (calls `TDT_demod` then `TDT_dFF_stage2`) to generate `_dFF.mat` outputs for each session.

Other analyses (LFP, behavior, etc.) are handled by their respective notebooks/scripts, but the base pipelines above are the standard sequence for spikes and dopamine photometry.
