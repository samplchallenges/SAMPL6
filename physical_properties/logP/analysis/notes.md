## 2019/03/31

### EVALUATION OF LOGP PREDICTIONS

#### Create Conda Environment

$ conda create -n sampl6_logP python=3.6  
$ source activate sampl6_logP  
$ conda install matplotlib  
$ conda install seaborn  
$ conda install scipy


#### Correct problematic characters in submission files

- `3vqbi-1559-logP-cosmologic-2.csv`: `<96>` replaced with `-`
- `6fyg5-1559-logP-method-2.csv `: `^M` at the end of each line were deleted.
- `ahmtf-1559-logP-b3pw91_tz-2.csv`: `<91>` and `<92>` replaced with `'`.  
- `dyxbt-1559-logP-b3pw91_tz-1.csv`: `<91>` and `<92>` replaced with `'`.  
- `hf4wj-1559-logP-vohringerlab-4.csv`: Remove empty lines after `name` and `category` keywords.
- `hmz0n-1559-logP-cosmologic-1.csv`: `<95>` and `<96>` replaced with `-`.
- `o7djk-1559-logP-b3pw91_tz-3.csv`: `<91>` and `<92>` replaced with `'`.
- `pcv32-1559-logP-method-4.csv`: `^M` at the end of each line were deleted.
- `zdj0j-1559-logP-method-3.csv`: `^M` at the end of each line were deleted.
- `2ggir-1559-logP-procaccicgen-1.csv`: Add a new line after `Methods` keyword.
- `hmz0n-1559-logP-cosmologic-1.csv`: `ü` in `Bürger` was replaced with `u`.
- `odex0-1559-logP-InterX-4.csv`: Add a new line after `Name:` keyword.
- `3vqbi-1559-logP-cosmologic-2.csv`: `ü` in `Bürger` was replaced with `u`.
- `bzeez-1559-logP-procaccigaff2-1.csv`: Add a new line after `Methods` keyword.
- `5svjv-1559-logP-procaccigaff2-2.csv`: Add a new line after `Methods` keyword.
- `padym-1559-logP-InterX-5.csv`: Add a new line after `Name:` keyword.
- `6cm6a-1559-logP-InterX-2.csv`: Add a new line after `Name:` keyword.
- `fcspk-1559-logP-InterX-3.csv`: Add a new line after `Name:` keyword.
- `ggm6n-1559-logP-procacciopls-2.csv`: Add a new line after `Methods` keyword.
- `y0xxd-1559-logP-procacci_cgen-2.csv`: Add a new line after `Methods` keyword.
- `v2q0t-1559-logP-InterX-6.csv`: Add a new line after `Name:` keyword.
- `25s67-1559-logP-procacciopls-1.csv`: Add a new line after `Methods` keyword.
- `eg52i-1559-logP-InterX-1.csv`: Add a new line after `Name:` keyword.
- `vzgyt-1559-logp-shin-1.csv`: File renamed as `vzgyt-1559-logP-shin-1.csv` to fix the `logP` section in file name.  

## 2019/04/05

#### Create new Conda environment based on `sampl6_pKa` environment

Delete old environment:
$ conda env remove -n sampl6_logP

Copy requirements file from /sampl6-physicochemical-properties/analysis_of_pKa_predictions/analysis_of_typeIII_predictions

$ conda create --name sampl6_logP --file requirements.txt

UnsatisfiableError: The following specifications were found to be in conflict:
  - libgfortran 3.0.0 1
  - numpy 1.15.1 py36_blas_openblashd3ea46f_0
Use "conda info <package>" to see the dependencies for each package.

I commented out libgfortran 3.0.0 1 line.


