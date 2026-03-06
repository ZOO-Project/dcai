# GeoCroissant Metadata Standard

GeoCroissant is an extension of the MLCommons Croissant metadata standard explicitly designed for Earth Observation and Geospatial Machine Learning datasets. 

It aims to make GeoAI datasets more findable, accessible, interoperable, and reusable (FAIR) by providing structured, machine-actionable metadata.

## Core Geospatial Extensions

GeoCroissant adds several key geospatial properties to the core Croissant standard:

1. **`spatialCoverage`** and **`temporalCoverage`**: Define the specific geographic extent (bounding boxes) and temporal range of the dataset. 
2. **`geocr:coordinateReferenceSystem`**: The CRS identifier (e.g., "EPSG:4326") necessary to properly align data on Earth.
3. **`geocr:spatialResolution`** & **`geocr:temporalResolution`**: Defines the physical spacing between pixels/points (e.g., `<value> 30, <unitText> 'm'`) and the time interval between observations.
4. **`geocr:bandConfiguration`**: Describes multidimensional data structures typical in remote sensing, mapping the number of bands and their logical order.
5. **`geocr:spectralBandMetadata`**: Per-band physical descriptors such as `centerWavelength` and `bandwidth`.

## Using GeoCroissant

A compliant dataset must declare conformance to both Croissant and GeoCroissant vocabularies via its context and `conformsTo` fields:

```json
"conformsTo": [
    "http://mlcommons.org/croissant/1.1",
    "http://mlcommons.org/croissant/geo/1.0"
]
```

### Example Implementation

Here is a simplified example of how optical satellite characteristics are structured using the `geocr` namespace:

```json
{
  "geocr:coordinateReferenceSystem": "EPSG:4326",
  "geocr:spatialResolution": {
    "@type": "QuantitativeValue",
    "value": 30.0,
    "unitText": "m"
  },
  "geocr:bandConfiguration": {
    "@type": "geocr:BandConfiguration",
    "geocr:totalBands": 6,
    "geocr:bandNameList": ["Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2"]
  },
  "geocr:spectralBandMetadata": [
    {
      "@type": "geocr:SpectralBand",
      "name": "Blue",
      "geocr:centerWavelength": {
        "@type": "QuantitativeValue",
        "value": 490,
        "unitText": "nm"
      },
      "geocr:bandwidth": {
        "@type": "QuantitativeValue",
        "value": 65,
        "unitText": "nm"
      }
    }
  ]
}
```

## Advanced & Responsible GeoAI Domains

Beyond standard image parameters, GeoCroissant provides hooks for specialized domains and robust AI practices:
- **`geocr:samplingStrategy`** & **`geocr:spatialBias`**: Documents biases relative to geography (e.g., heavily sampled North America vs under-sampled global south) making downstream AI more responsible.
- **Space Weather & Heliophysics**: `geocr:MultiWavelengthConfiguration` and `geocr:SolarInstrumentCharacteristics` cater to multi-channel solar observations.
- **Data Discovery**: `geocr:spatialIndex` accommodates rapid catalog indexing based on spatial tiling schemes (like DGGS) to process massive web-hosted datasets efficiently.
