# lagrangian-precision-wola

MATLAB and VHDL source for the paper:
> J. C. Maravilla, "A Perceptually Guided Lagrangian Framework for Fixed-Point Precision in Hearing-Aid Filterbanks," *IEEE Trans. Circuits Syst. II*, submitted 2026.

The paper proposes a psychoacoustic Lagrangian framework that allocates fixed-point arithmetic precision heterogeneously across the 64 subbands of a WOLA filterbank, concentrating bits where they contribute most to speech intelligibility. The framework is validated on a Xilinx Artix-7 FPGA (XC7A200T) using post-implementation SAIF-annotated power estimation.

---

## Structure

```
.
├── matlab/
│   ├── core/           shared functions used by all experiments
│   ├── experiments/    numbered experiment scripts (E1-E11)
│   └── utils/          listening and diagnostic tools
│
└── hdl/
    ├── src/            synthesizable VHDL source files
    ├── tb/             VHDL testbenches for SAIF power measurement
    └── ip/
        └── xfft_0/     Xilinx FFT IP core configuration (xfft_0.xci)
```

---

## MATLAB

### Core (`matlab/core/`)

**allocators.m** - the three bit allocation strategies: uniform (equal bits per band), linear (water-filling), and proposed (psychoacoustic Lagrangian with hardware cost model). Also contains the Newton solver for the per-band KKT condition and a budget verification utility.

**local_functions.m** - shared WOLA signal processing: analysis FFT, synthesis IFFT and overlap-add, subband quantization, and envelope distortion.

**iso_budget_config.m** - shared parameters for all experiments: M=64, fs=16kHz, alpha=1.37, beta=0, B_min=4, B_max=24, budget fractions, and the psychoacoustic importance function S(k).

**discover_librispeech.m** - recursively finds all `.flac` files under one or more LibriSpeech root directories. Returns a struct array with path, speaker, chapter, and corpus tag.

### Experiments (`matlab/experiments/`)

Run in order. **E1 must run first** as it generates `E1_sigma2_raw.csv` and `E1_sigma2_eff.csv` which all downstream experiments load.

| Script | What it does | Paper |
|---|---|---|
| E1_iso_budget.m | Main distortion results: D_w and D_mse for all strategies at all four budget levels. Also generates allocation vectors and sigma2 files for downstream use. | Table II, Figure 4 |
| E2_iso_budget.m | Confirms D_w advantage is algorithmic by isolating mid-band and high-band envelope distortion per strategy. | Section VI-A |
| E6_iso_budget.m | STOI evaluation across all 5,567 LibriSpeech utterances at all four budget levels. | Section VI-B |
| E5_iso_budget.m | Wilcoxon signed-rank tests comparing uniform and linear against proposed on distortion and STOI. Requires E1 and E6 outputs. | Section VI-B |
| E7_iso_budget.m | Robustness across phonetically defined subsets: whispered, fricative-heavy, unvoiced-burst, and full corpus. | Section VI-C |
| E9_iso_budget.m | HASPI v2 and HASQI v2 evaluation on a 200-utterance stratified subset under a sensorineural audiogram. Requires the HASPI/HASQI toolbox. | Section VI-B |
| E10_iso_distortion.m | Iso-distortion analysis: binary search for the minimum bit budget each strategy needs to match uniform's D_w. | Table III |
| E11_bitdepth_sweep.m | Per-band bit-depth sensitivity sweep from B=1 to B=24. Derives the B_min=4 floor used throughout. | Section III-F |

### Utils (`matlab/utils/`)

**E_listen.m** - picks random clips, processes them through all three strategies at all four budgets, and writes labelled wav files for informal listening comparison. Not part of the paper results pipeline.

---

## VHDL

### Source files (`hdl/src/`)

**WOLA_Filterbank.vhd** - top-level entity. Structural architecture connecting all five components: analysis, precision controller, operand isolation, MAC, and synthesis.

**WOLA_Analysis.vhd** - analysis path: Hann windowing, 128-point FFT via the xfft_0 IP core, outputs 64 positive-frequency subband coefficients.

**Wola_Synthesis.vhd** - synthesis path: conjugate-symmetric IFFT, synthesis windowing, overlap-add reconstruction.

**Precision_Controller.vhd** - four-state FSM (IDLE, FETCH, EXECUTE, WAIT) that sequences through all 64 subbands per frame, reads B_k from the allocation ROM, and drives the AND mask and MAC enable signals.

**Operand_Isolation.vhd** - AND-mask datapath. Zeros the (B_MAX - B_k) LSBs of both MAC operands, eliminating partial-product toggle activity in those columns.

**Precision_MAC.vhd** - 24x24-bit registered multiplier-accumulator. Receives masked operands from Operand_Isolation.

**Precision_MAC_lut.vhd** - LUT-only array multiplier used for hardware cost characterization. Synthesized with DSP48 inference suppressed to measure the structural O(B^2) LUT count directly.

**MAC_Path_Wrapper.vhd** - testbench-facing wrapper that exposes the Precision Controller, Operand Isolation, and MAC as a single DUT with a clean frame/band handshake interface. Both testbenches instantiate this.

### Testbenches (`hdl/tb/`)

Both testbenches instantiate `MAC_Path_Wrapper` and drive it with LFSR-generated operand data across 256 frames. They exercise the precision control and MAC switching activity path. The SAIF power numbers reported in the paper come from these simulations.

**tb_filterbank_saif.vhd** - iso-budget testbench. 12 runs covering all combinations of three strategies and four budget levels. Allocation vectors sourced from `E1_iso_bit_allocation_<N>pct.csv`. Run index is a Vivado generic (G_RUN=0..11).

**tb_filterbank_saif_isodist.vhd** - iso-distortion testbench. 12 runs using the minimum bit budgets from `E10_iso_distortion_detail.csv`, the fewest bits each strategy needs to match uniform's D_w at each quality target.

### IP core (`hdl/ip/xfft_0/`)

**xfft_0.xci** - Vivado IP configuration for the 128-point FFT used in WOLA_Analysis. Targets Artix-7, configured for the transform length, data width, and scaling used in the paper. Vivado regenerates all output products from this file.

---

## Reproducing the results

### MATLAB

1. Download LibriSpeech `dev-clean` and `dev-other` from [openslr.org/12](https://www.openslr.org/12) and extract locally.
2. Update `LIBRISPEECH_ROOTS` at the top of each experiment script to point to your local copies.
3. Add `matlab/core/` to your MATLAB path.
4. Run experiments in order starting with E1:
   ```matlab
   run('matlab/experiments/E1_iso_budget.m')
   ```
5. E9 requires the HASPI v2 / HASQI v2 toolbox. Update `HASPI_DIR` in the script to your local copy.

### VHDL

1. Open Vivado 2025.2 and create a project targeting XC7A200T (FBG484, speed grade -1).
2. Add all files from `hdl/src/` as design sources.
3. Add `hdl/ip/xfft_0/xfft_0.xci` and regenerate IP output products.
4. Add `hdl/tb/tb_filterbank_saif.vhd` as a simulation source.
5. Run post-implementation simulation for each run by overriding the generic:
   ```tcl
   set_property generic {G_RUN=0} [get_filesets sim_1]
   launch_simulation -mode post-implementation
   ```
6. After each simulation, read the SAIF file and report power:
   ```tcl
   read_saif ./saif_run_0.saif
   report_power
   ```
7. Repeat for G_RUN=0 through 11 for iso-budget results, then repeat with `tb_filterbank_saif_isodist.vhd` for iso-distortion results.

---

## Dependencies

| Tool | Version |
|---|---|
| MATLAB | R2023b or later |
| MATLAB Audio Toolbox | required for stoi() and pesq() |
| MATLAB Statistics Toolbox | required for signrank() in E5 |
| HASPI v2 / HASQI v2 | required for E9 only |
| Vivado | 2025.2 WebPACK |
| Target device | Xilinx XC7A200T (FBG484, speed grade -1) |

---

## Notes

- All MATLAB experiments use alpha=1.37, beta=0, B_min=4, B_max=24. These are defined in `iso_budget_config.m` and restated at the top of each self-contained script.
- The VHDL testbenches measure FPGA dynamic power on the MAC and precision datapath only, not the full audio pipeline. These figures are not directly comparable to hearing-aid ASIC power budgets.
- LibriSpeech audio files are not included. The corpus is freely available at [openslr.org/12](https://www.openslr.org/12).

---

## License

MIT
