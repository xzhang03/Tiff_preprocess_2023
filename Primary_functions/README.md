## Functions

#### Bleach correct the entire tiff movie (experimental):
```Matlab
movout = TiffBleachCorrectMovie(mov, varargin)
```

#### Manually matching dendrite and soma ROIs per FOV:
```Matlab
TiffDenMatchingManual(mouse, date, runnum, optotune, varargin)
```

#### Generate an opto-triggered movie:
```Matlab
movout = TiffGenerateOptoTiffMovie(mouse,date,run, pmt,sbxtype,varargin)
```

#### Make a classifier-based movie:
```Matlab
TiffMakeClassifierMovie(mousecell, datecell, runvec, optotunecell, varargin)
```

#### Pull signals core (private):
```Matlab
 cellsort = TiffPullSignalsCore(mov, cellsort, xybin)
```

#### Pull signals from cellpose segments:
```Matlab
grinface = TiffSignalsCellpose(mouse, date, runs, varargin)
```

#### Manually match ROIs across days and experiments:
```Matlab
TiffROIMatchingManual(mousecell, datecell, runs, optotunes, varargin)
```

#### Manually match ROIs across experiment types:
```Matlab
TiffROIMatchingManualXExpt(mousecell1, datecell1, runs1, optotunes1, mousecell2, datecell2, runs2, optotunes2, varargin)
```

#### Add soma ROIs to signal files:
```Matlab
TiffSignalsAddSomaROIs(mouse, date, runs, varargin)
```


#### Apply windowed-DFF calculation to signal files:
```Matlab
cellsort = TiffSignalsWindowedDFF(cellsort, fps, time_window, percentile)
```

#### Quick register a stack from center to the top/bottom:
```Matlab
TiffStrackReg(fpath,varargin)
```

#### XY rigid registration:
```Matlab
Tiffxyreg(mouse, date, run, varargin)
```

#### XY rigid registration core (private):
```Matlab
shifts = Tiffxyreg_core(im, ref, gausssize)
```

## Obesolete functions

#### Make split movies (Obsolete):
```Matlab
TiffCombineOptoSplits(varargin)
```

#### Extract traces based on morphological segmentation:
```Matlab
TiffMorphologicalFilterExtractROIs(mouse, date, run, varargin)
```

#### Regress signals to xy shifts:
```Matlab
TiffSignalsRegress(mouse, date, runs, varargin)
```
