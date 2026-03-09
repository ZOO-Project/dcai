import json
from rdflib import Graph, Namespace, URIRef, Literal, BNode
from rdflib.namespace import DCTERMS, DCAT, FOAF, XSD, RDF
from urllib.parse import quote


def geocroissant_to_geodcat_jsonld(geocroissant_json, output_file="geodcat.jsonld"):
    """Convert GeoCroissant JSON-LD to GeoDCAT-AP compliant format"""
    g = Graph()

    # Namespaces
    GEO = Namespace("http://www.opengis.net/ont/geosparql#")
    SCHEMA = Namespace("https://schema.org/")
    SPDX = Namespace("http://spdx.org/rdf/terms#")
    ADMS = Namespace("http://www.w3.org/ns/adms#")
    PROV = Namespace("http://www.w3.org/ns/prov#")
    GEOCR = Namespace("http://mlcommons.org/croissant/geo/")

    g.bind("dct", DCTERMS)
    g.bind("dcat", DCAT)
    g.bind("foaf", FOAF)
    g.bind("geo", GEO)
    g.bind("schema", SCHEMA)
    g.bind("spdx", SPDX)
    g.bind("adms", ADMS)
    g.bind("prov", PROV)
    g.bind("geocr", GEOCR)

    # Create dataset URI
    dataset_name = geocroissant_json.get("name", "dataset")
    # URL-encode the dataset name to handle spaces and special characters
    safe_name = quote(dataset_name, safe='')
    dataset_uri = URIRef(f"https://example.org/{safe_name}")
    
    # Basic dataset properties
    g.add((dataset_uri, RDF.type, DCAT.Dataset))
    g.add((dataset_uri, RDF.type, SCHEMA.Dataset))
    g.add((dataset_uri, DCTERMS.title, Literal(geocroissant_json["name"])))
    g.add((dataset_uri, DCTERMS.description, Literal(geocroissant_json["description"])))
    
    # License
    if "license" in geocroissant_json:
        g.add((dataset_uri, DCTERMS.license, URIRef(geocroissant_json["license"])))
    
    # Version
    if "version" in geocroissant_json:
        g.add((dataset_uri, ADMS.version, Literal(geocroissant_json["version"])))
    
    # Date published
    if "datePublished" in geocroissant_json:
        g.add((dataset_uri, DCTERMS.issued, Literal(geocroissant_json["datePublished"], datatype=XSD.date)))
    
    # ConformsTo
    for conformance in geocroissant_json.get("conformsTo", []):
        g.add((dataset_uri, DCTERMS.conformsTo, URIRef(conformance)))
    
    # Keywords
    for keyword in geocroissant_json.get("keywords", []):
        g.add((dataset_uri, DCAT.keyword, Literal(keyword)))
    
    # Spatial coverage
    spatial_coverage = geocroissant_json.get("spatialCoverage", {})
    if spatial_coverage and "geo" in spatial_coverage:
        geo_shape = spatial_coverage["geo"]
        if "box" in geo_shape:
            # Parse the bounding box (south west north east format)
            bbox = geo_shape["box"].split()
            if len(bbox) == 4:
                spatial_uri = URIRef(f"{dataset_uri}/spatial")
                g.add((dataset_uri, DCTERMS.spatial, spatial_uri))
                g.add((spatial_uri, RDF.type, DCTERMS.Location))
                
                # Create WKT polygon from bounding box
                south, west, north, east = bbox
                wkt_bbox = f"POLYGON(({west} {south}, {east} {south}, {east} {north}, {west} {north}, {west} {south}))"
                g.add((spatial_uri, GEO.asWKT, Literal(wkt_bbox, datatype=GEO.wktLiteral)))
    
    # Temporal coverage
    if "temporalCoverage" in geocroissant_json:
        temporal_coverage = geocroissant_json["temporalCoverage"]
        if "/" in temporal_coverage:
            start_date, end_date = temporal_coverage.split("/")
            temporal_uri = URIRef(f"{dataset_uri}/temporal")
            g.add((dataset_uri, DCTERMS.temporal, temporal_uri))
            g.add((temporal_uri, RDF.type, DCTERMS.PeriodOfTime))
            g.add((temporal_uri, DCAT.startDate, Literal(start_date, datatype=XSD.date)))
            g.add((temporal_uri, DCAT.endDate, Literal(end_date, datatype=XSD.date)))
    
    # GeoCroissant specific properties
    if "geocr:coordinateReferenceSystem" in geocroissant_json:
        crs_uri = URIRef(f"http://www.opengis.net/def/crs/{geocroissant_json['geocr:coordinateReferenceSystem']}")
        g.add((dataset_uri, GEOCR.coordinateReferenceSystem, crs_uri))
    
    # Spatial resolution
    if "geocr:spatialResolution" in geocroissant_json:
        spatial_res = geocroissant_json["geocr:spatialResolution"]
        if isinstance(spatial_res, dict) and "@type" in spatial_res:
            res_node = BNode()
            g.add((dataset_uri, GEOCR.spatialResolution, res_node))
            g.add((res_node, RDF.type, SCHEMA.QuantitativeValue))
            if "value" in spatial_res:
                g.add((res_node, SCHEMA.value, Literal(spatial_res["value"])))
            if "unitText" in spatial_res:
                g.add((res_node, SCHEMA.unitText, Literal(spatial_res["unitText"])))
    
    # Temporal resolution
    if "geocr:temporalResolution" in geocroissant_json:
        temporal_res = geocroissant_json["geocr:temporalResolution"]
        if isinstance(temporal_res, dict) and "@type" in temporal_res:
            res_node = BNode()
            g.add((dataset_uri, GEOCR.temporalResolution, res_node))
            g.add((res_node, RDF.type, SCHEMA.QuantitativeValue))
            if "value" in temporal_res:
                g.add((res_node, SCHEMA.value, Literal(temporal_res["value"])))
            if "unitText" in temporal_res:
                g.add((res_node, SCHEMA.unitText, Literal(temporal_res["unitText"])))
    
    # Distributions
    for dist in geocroissant_json.get("distribution", []):
        if dist.get("@type") == "cr:FileObject":
            dist_id = dist.get("@id", "distribution")
            dist_uri = URIRef(f"{dataset_uri}/distribution/{dist_id}")
            g.add((dataset_uri, DCAT.distribution, dist_uri))
            g.add((dist_uri, RDF.type, DCAT.Distribution))
            
            if "name" in dist:
                g.add((dist_uri, DCTERMS.title, Literal(dist["name"])))
            if "description" in dist:
                g.add((dist_uri, DCTERMS.description, Literal(dist["description"])))
            if "contentUrl" in dist:
                g.add((dist_uri, DCAT.accessURL, URIRef(dist["contentUrl"])))
            if "encodingFormat" in dist:
                g.add((dist_uri, DCAT.mediaType, Literal(dist["encodingFormat"])))
            if "md5" in dist:
                checksum_node = BNode()
                g.add((dist_uri, SPDX.checksum, checksum_node))
                g.add((checksum_node, RDF.type, SPDX.Checksum))
                g.add((checksum_node, SPDX.algorithm, SPDX.checksumAlgorithm_md5))
                g.add((checksum_node, SPDX.checksumValue, Literal(dist["md5"])))
        
        elif dist.get("@type") == "cr:FileSet":
            # Handle FileSet as a special type of distribution
            dist_id = dist.get("@id", "fileset")
            dist_uri = URIRef(f"{dataset_uri}/distribution/{dist_id}")
            g.add((dataset_uri, DCAT.distribution, dist_uri))
            g.add((dist_uri, RDF.type, DCAT.Distribution))
            g.add((dist_uri, RDF.type, GEOCR.FileSet))
            
            if "name" in dist:
                g.add((dist_uri, DCTERMS.title, Literal(dist["name"])))
            if "description" in dist:
                g.add((dist_uri, DCTERMS.description, Literal(dist["description"])))
            if "encodingFormat" in dist:
                g.add((dist_uri, DCAT.mediaType, Literal(dist["encodingFormat"])))
            if "includes" in dist:
                g.add((dist_uri, GEOCR.includes, Literal(dist["includes"])))
    
    # Record sets and fields (as additional metadata)
    for record_set in geocroissant_json.get("recordSet", []):
        if record_set.get("@type") == "cr:RecordSet":
            rs_id = record_set.get("@id", record_set.get("name", "recordset"))
            rs_uri = URIRef(f"{dataset_uri}/recordset/{rs_id}")
            g.add((dataset_uri, GEOCR.recordSet, rs_uri))
            g.add((rs_uri, RDF.type, GEOCR.RecordSet))
            
            if "name" in record_set:
                g.add((rs_uri, DCTERMS.title, Literal(record_set["name"])))
            if "description" in record_set:
                g.add((rs_uri, DCTERMS.description, Literal(record_set["description"])))
            
            # Handle fields
            for field in record_set.get("field", []):
                if field.get("@type") == "cr:Field":
                    field_id = field.get("@id", field.get("name", "field"))
                    field_uri = URIRef(f"{rs_uri}/field/{field_id}")
                    g.add((rs_uri, GEOCR.field, field_uri))
                    g.add((field_uri, RDF.type, GEOCR.Field))
                    
                    if "name" in field:
                        g.add((field_uri, DCTERMS.title, Literal(field["name"])))
                    if "description" in field:
                        g.add((field_uri, DCTERMS.description, Literal(field["description"])))
                    if "dataType" in field:
                        g.add((field_uri, GEOCR.dataType, Literal(field["dataType"])))

    # Serialize outputs
    g.serialize(destination=output_file, format="json-ld", indent=2)
    print(f"GeoDCAT JSON-LD metadata written to {output_file}")

    g.serialize(destination="geodcat.ttl", format="turtle")
    print("GeoDCAT Turtle metadata written to geodcat.ttl")


if __name__ == "__main__":
    with open("croissant.json", "r") as f:
        geocroissant = json.load(f)

    geocroissant_to_geodcat_jsonld(geocroissant, output_file="geodcat.jsonld")
