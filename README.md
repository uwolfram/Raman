# Module to evaluate Raman spectra
## General
Two scripts are provided ramanavrgfits.py and ramanavrgspectra.py. Both deal with repeated measures per sample but the first averages nonlinear least square fits of portions of the spectra, while the second averages the spectra for each repeated measures and then determines the nonlinear least square fits.

Code evaluates n individual Raman spectra taken from bone samples and averages the resulting fits 

Code was originally developed 2015 for: 
    [1] Mirzaali, M., Schwiedrzik, J., Thaiwichai, S., Best, J., Michler, J., Zysset, P. & Wolfram, U. 
        Mechanical properties of cortical bone and their relationships with age, gender, composition and 
        microindentation properties in the elderly. Bone 93, 196–211 (2016).

Spectra are evaluated for:
    v1PO4
    v2PO4
    Amide1
    Amide3
    PYD
    Lipids
    
Boundaries for these spectral regions are hard-encoded at the moment. Check [1] for baclground info.


##Input
spectra : 
    Directory holding a set of baseline corrected spectra consisting of only the counts.
    
nbrspec:
    Number of specimens tested.
    
nbrmeas:
    Number of repeated measurements per specimen
    
shifts:
    File containing counts and shifts. Do not ask why this was setup like this.
    

## Returns

ramanavrgfits-mean.dat && ramanavrgfits-std.dat: 
    2 ascii data files providing the mean and std (where possible) over nbrmeas for
    filename (donor)
    age
    gender
    v1PO4fwhm
    v1PO4int
    v1PO4peak
    v2PO4int
    amid1int
    amid1peak
    amid3int
    pyd
    CH3int
    lipid
        
## Example run
```bash
python3 ramanavrgfits.py -d Median-subtracted -n 3 -m 10 -s S18_F70L_1_AX_01.txt
```
OR
```bash
python3 ramanavrgspectra.py -d Median-subtracted -n 3 -m 10 -s S18_F70L_1_AX_01.txt
```
