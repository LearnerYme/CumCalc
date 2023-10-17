# Cumulant Calculation

Version: 3.1

Author: Yige Huang

## Guide

1. Use `make cumulant` to get `runCumulant`, which do raw cumulant calculation for each RefMult3 bin and CBWC procedure.

2. Use `make cbwc` to get `cbwc`, which get CBWC results from `raw.root`.

3. Use `make duoCBWC` to get `duoCBWC`, which from U and L run raw root files get CBWC results.

## Change log

17.10.2023 by yghuang (3.1):

> duoCBWC now support up to 10 files, and will read file list / centrality bin list instead of 2 file names.

18.08.2023 by yghuang (3.0):

> Full constructed cumulant calculation system.
