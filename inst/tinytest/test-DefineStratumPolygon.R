stratumFile_wkt <- system.file("testresources", "strata_sandeel_2020_firstCoverage.wkt", package = "RstoxBase")
stratumPolygon_wkt <- RstoxBase::DefineStratumPolygon(
    DefinitionMethod = "ResourceFile", 
    FileName = stratumFile_wkt
)

stratumFile_GeoJSON <- system.file("testresources", "strata_sandeel_2020_firstCoverage.geojson", package = "RstoxBase")
stratumPolygon_GeoJSON <- RstoxBase::DefineStratumPolygon(
    DefinitionMethod = "ResourceFile", 
    FileName = stratumFile_GeoJSON, 
    StratumNameLabel = "StratumName"
)

stratumFile_shape <- system.file("testresources", "strata_sandeel_2020_firstCoverage/StratumPolygon.shp", package = "RstoxBase")
stratumPolygon_shape <- RstoxBase::DefineStratumPolygon(
    DefinitionMethod = "ResourceFile", 
    FileName = stratumFile_shape, 
    StratumNameLabel = "StratumNam"
)

expect_equal(stratumPolygon_wkt, stratumPolygon_GeoJSON, check.attributes = F)
stratumPolygon_shape$StratumNam  <-   NULL
expect_equal(stratumPolygon_wkt, stratumPolygon_shape, check.attributes = F)


stratumNames <- RstoxBase::getStratumNames(stratumPolygon_wkt)

expect_equal(
    stratumNames, 
    c("AlbjoernLing", "Engelsk_Klondyke_2020", "Inner_Shoal_East_2016", "Inner_Shoal_North_2020", "Inner_Shoal_West_2018", "Nordgyden", "Ostbanken_2020", "Outer_Shoal_2020_1", "Vestbanken_North_2020", "VestbankenSouthEast", "VestbankenSouthWest", "Vikingbanken")
)
expect_equal(
    class(stratumPolygon_wkt)[1], 
    "SpatialPolygonsDataFrame"
)
