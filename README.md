# Analysis code for the R(D*) measurement with parked data in CMSSW

This code has to be run after sourceing CMSSW_10_2_3. It does not require compilation but needs the CMSSW environment to be active.

## Suggestions for installation

```
cd ~
mkdir RDstAnalysis
cd RDstAnalysis

cmsrel CMSSW_10_2_3
cd CMSSW_10_2_3/src
cmsenv
cd ../..

git clone https://github.com/ocerri/BPH_RD_Analysis
```

# Utilities

### Notebook extensions

To install notebook extensions:
``pip install jupyter_contrib_nbextensions``
``jupyter contrib nbextension install --user``
In case config parser have erros try to install the specific release
``pip install --user configparser==3.5.1``
