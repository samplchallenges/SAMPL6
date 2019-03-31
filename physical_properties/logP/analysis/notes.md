## 03/31/2019

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
