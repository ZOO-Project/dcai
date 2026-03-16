# DCAI: GeoCroissant Recipes

## What is GeoCroissant?

GeoCroissant extends the [MLCommons Croissant](http://mlcommons.org/croissant/1.1) metadata standard with geospatial concepts for GeoAI workflows. It introduces support for spatial and temporal coverage, coordinate reference systems, spatial resolution, band configuration, and time-series metadata.

## Requirements

Install the required dependency:

```bash
pip install mlcroissant
```

## Usage

Load a GeoCroissant dataset using the `mlcroissant` Python library:

```python
import mlcroissant as mlc

dataset = mlc.Dataset("geocroissant.json")
```

## Conformance

GeoCroissant datasets must declare conformance to both the Croissant and GeoCroissant specifications:

```json
"dct:conformsTo": [
  "http://mlcommons.org/croissant/1.1",
  "http://mlcommons.org/croissant/geo/1.0"
]
```

## Authors

```
Rajat Shinde, Manil Maskey, AG Stephens, Harsh Shinde, Joseph Edgerton, Dr. Tejasri Nampally.,  
Douglas Fils, Edenna Chen, Claus Weiland, Pedram Ghamisi, Gerald Fenoy,  
Yuhan Douglas Rao, Omar Benjelloun, and Elena Simperl
```

GeoCroissant Working Group · [croissant-geo@mlcommons.org](mailto:croissant-geo@mlcommons.org)

## Acknowledgements

- MLCommons GeoCroissant Working Group
- MLCommons Croissant Working Group
- Open Geospatial Consortium (OGC) GeoAI Domain Working Group

